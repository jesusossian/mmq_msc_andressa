using JuMP, Gurobi, DelimitedFiles

# Solving the lot size by Lagrangian Relaxation
cp = 0.0
cr = 0.0
hp = 200
hr = 200
fp = 0.2
fr = 1.0

# Reading demand data
d, header = readdlm("demand.csv", ',', header=true)

# Reading return data
r, header = readdlm("return.csv", ',', header=true)

n = length(d)

lambdap = zeros(Float64,n)
lambdar = zeros(Float64,n)

# the length of 'd' and  the length return must match
@assert length(d) == length(r)

# making these data global so that any function can access data
global n, d, r, cp, cr, hp, hr, fp, fr, lambdap, lambdar

function optimal(n)

    SD = zeros(Int,n,n)
    SR = zeros(Int,n,n)
 
    for i=1:n
        SD[i,i] = d[i]
        SR[i,i] = r[i]
        for j=(i+1):n
            SD[i,j] = SD[i,j-1] + d[j]
            SR[i,j] = SR[i,j-1] + d[j]
        end
    end

    model = Model(Gurobi.Optimizer)
    set_silent(model)
  
    @variable(model, 0 <= xp[t=1:n] <= Inf)
    @variable(model, 0 <= xr[t=1:n] <= Inf)
    @variable(model, yp[t=1:n], Bin)
    @variable(model, yr[t=1:n], Bin)
    @variable(model, 0 <= sp[t=1:n] <= Inf)
    @variable(model, 0 <= sr[t=1:n] <= Inf)
  
    @objective(model, Min, 
        sum(cp*xp[t] + hp*sp[t] + fp*yp[t] for t=1:n) 
        + sum(cr*xr[t] + hr*sr[t] + fr*yr[t] for t=1:n)
    )

    @constraint(model, balanceP0, xp[1] + xr[1] - sp[1] == d[1])
    @constraint(model, balanceP[t=2:n], sp[t-1] + xp[t] + xr[t] - sp[t] == d[t])

    @constraint(model, balanceR0, -xr[1] - sr[1] == - r[1])
    @constraint(model, balanceR[t=2:n], sr[t-1] - xr[t] - sr[t] == - r[t])

    @constraint(model, setupP[t=1:n], xp[t] <= SD[t,n]*yp[t])

    @constraint(model, setupR[t=1:n], xr[t] <= min(SR[1,t],SD[t,n])*yr[t])
    
    @constraint(model, vin03[t=2:n,l=t:n,t<=l], 
        sp[t-1] + sum(SD[k,l]*(yp[k]+yr[k]) for k=t:l) >= SD[t,l]
    )
    
    @constraint(model, vin04[t=1:n,l=t:n,t<=l], 
        sr[l] + sum(SR[t,k]*yr[k] for k=t:l) >= SR[t,l]
    )

    @constraint(model, sp[n] == 0)

    optimize!(model)

    Z_opt = objective_value(model)

    xp_opt = value.(xp)
    yp_opt = value.(yp)
    sp_opt = value.(sp)
    xr_opt = value.(xr)
    yr_opt = value.(yr)
    sr_opt = value.(sr)

    return Z_opt, xp_opt, yp_opt, sp_opt, xr_opt, yr_opt, sr_opt

end


function lower_bound()
  
    SD = zeros(Int,n,n)
    SR = zeros(Int,n,n)
 
    for i=1:n
        SD[i,i] = d[i]
        SR[i,i] = r[i]
        for j=(i+1):n
            SD[i,j] = SD[i,j-1] + d[j]
            SR[i,j] = SR[i,j-1] + d[j]
        end
    end

    model = Model(Gurobi.Optimizer)
    set_silent(model)
  
    @variable(model, 0 <= xp[t=1:n] <= Inf)
    @variable(model, 0 <= xr[t=1:n] <= Inf)
    @variable(model, 0 <= yp[t=1:n] <= 1)
    @variable(model, 0 <= yr[t=1:n] <= 1)
    @variable(model, 0 <= sp[t=1:n] <= Inf)
    @variable(model, 0 <= sr[t=1:n] <= Inf)
  
    @objective(model, Min, 
        sum((cp + lambdap[t])*xp[t] + hp*sp[t] + (fp - lambdap[t]*SD[t,n] )*yp[t] for t=1:n) 
        + sum(cr*xr[t] + hr*sr[t] + fr*yr[t] for t=1:n)
    )

    @constraint(model, balanceP0, xp[1] + xr[1] - sp[1] == d[1])
    @constraint(model, balanceP[t=2:n], sp[t-1] + xp[t] + xr[t] - sp[t] == d[t])

    @constraint(model, balanceR0, -xr[1] - sr[1] == - r[1])
    @constraint(model, balanceR[t=2:n], sr[t-1] - xr[t] - sr[t] == - r[t])

    @constraint(model, setupP[t=1:n], xp[t] <= SD[t,n]*yp[t])

    @constraint(model, setupR[t=1:n], xr[t] <= min(SR[1,t],SD[t,n])*yr[t])
    
    @constraint(model, vin03[t=2:n,l=t:n,t<=l], 
        sp[t-1] + sum(SD[k,l]*(yp[k]+yr[k]) for k=t:l) >= SD[t,l]
    )
    
    @constraint(model, vin04[t=1:n,l=t:n,t<=l], 
        sr[l] + sum(SR[t,k]*yr[k] for k=t:l) >= SR[t,l]
    )

    @constraint(model, sp[n] == 0)

    optimize!(model)

    Z_D = objective_value(model)
    
    xp_d = value.(xp)
    yp_d = value.(yp)
    sp_d = value.(sp)
    xr_d = value.(xr)
    yr_d = value.(yr)
    sr_d = value.(sr)

    return Z_D, xp_d, yp_d, sp_d, xr_d, yr_d, sr_d

end


function upper_bound()
     
    SD = zeros(Int,n,n)
    SR = zeros(Int,n,n)
 
    for i=1:n
        SD[i,i] = d[i]
        SR[i,i] = r[i]
        for j=(i+1):n
            SD[i,j] = SD[i,j-1] + d[j]
            SR[i,j] = SR[i,j-1] + d[j]
        end
    end
    
    xp_sol = zeros(Float64,n)
    yp_sol = zeros(Float64,n)
    sp_sol = zeros(Float64,n)
    xr_sol = zeros(Float64,n)
    yr_sol = zeros(Float64,n)
    sr_sol = zeros(Float64,n)

    yp_val = zeros(Int,n)
    yr_val = zeros(Int,n)

    model = Model(Gurobi.Optimizer)
    set_silent(model)
  
    @variable(model, 0 <= xp[t=1:n] <= Inf)
    @variable(model, 0 <= xr[t=1:n] <= Inf)
    @variable(model, yp[t=1:n], Bin)
    @variable(model, yr[t=1:n], Bin)
    @variable(model, 0 <= sp[t=1:n] <= Inf)
    @variable(model, 0 <= sr[t=1:n] <= Inf)
    
    @objective(model, Min, 
        sum(cp*xp[t] + hp*sp[t] + fp*yp[t] for t=1:n) 
        + sum(cr*xr[t] + hr*sr[t] + fr*yr[t] for t=1:n)
    )

    @constraint(model, balanceP0, xp[1] + xr[1] - sp[1] == d[1])
    @constraint(model, balanceP[t=2:n], sp[t-1] + xp[t] + xr[t] - sp[t] == d[t])

    @constraint(model, balanceR0, -xr[1] - sr[1] == - r[1])
    @constraint(model, balanceR[t=2:n], sr[t-1] - xr[t] - sr[t] == - r[t])

    @constraint(model, setupP[t=1:n], xp[t] <= SD[t,n]*yp[t])

    @constraint(model, setupR[t=1:n], xr[t] <= min(SR[1,t],SD[t,n])*yr[t])
    
    @constraint(model, vin03[t=2:n,l=t:n,t<=l], 
        sp[t-1] + sum(SD[k,l]*(yp[k]+yr[k]) for k=t:l) >= SD[t,l]
    )
    
    @constraint(model, vin04[t=1:n,l=t:n,t<=l], 
        sr[l] + sum(SR[t,k]*yr[k] for k=t:l) >= SR[t,l]
    )

    @constraint(model, sp[n] == 0)
  
    # set hsrf and fsrf according to the initial parameters
    hsrf = 5 #params.horsizerf
    fsrf = 2 #params.fixsizerf

    time = 0
    alpha = 1
    beta = min(alpha+hsrf-1, n)
  
    #relax variable
    for t in 1:n
        unset_binary(yp[t])
        unset_binary(yr[t])
    end
  
    while(beta <= n)
        t1 = time_ns()
  
        if(alpha > 1)
            #fix variables yp
            for t in 1:alpha-1
                if yp_val[t] > 0.99
                    set_lower_bound(yp[t],1.0)
                    #set_start_value(yp[t],1.0)
                else
                    set_upper_bound(yp[t],0.0)
                    #set_start_value(yp[t],0.0)
                end
            end

            #fix variables yr
            for t in 1:alpha-1
                if yr_val[t] > 0.99
                    set_lower_bound(yr[t],1.0)
                else
                    set_upper_bound(yr[t],0.0)
                end
            end
        end

        #set binary variables
        for t in alpha:beta
            set_binary(yp[t])
            set_binary(yr[t])
        end
    
        #relax binary variables    
        for t in beta+1:n
            if is_binary(yp[t]) == true
                unset_binary(yp[t])
                set_lower_bound(yp[t],0.0)
                set_upper_bound(yp[t],1.0)
            end
        end
    
        for t in beta+1:n
            if is_binary(yr[t]) == true
                unset_binary(yr[t])
                set_lower_bound(yr[t],0.0)
                set_upper_bound(yr[t],1.0)
             end
        end

        # solve problem
        optimize!(model)

        yp_val = value.(yp)
        yr_val = value.(yr)

        alpha = alpha+fsrf
    
        if beta == n
            beta = n+1
        else
            beta = min(alpha+hsrf-1,n)
        end

        t2 = time_ns()
        time += (t2-t1)/1.0e9
  
    end
  
    Z = objective_value(model)

    xp_sol = value.(yp)
    yp_sol = value.(yp)
    sp_sol = value.(sp)
    xr_sol = value.(xr)
    yr_sol = value.(yr)
    sr_sol = value.(sr)

    return Z, xp_sol, yp_sol, sp_sol, xr_sol, yr_sol, sr_sol

end


function lagrangian_relaxation()
    # The maximum number of iterations allowed
    MAX_ITER = 10

    # To track the upper and lower bounds
    UB = Array{Float64}(undef, 0)
    LB = Array{Float64}(undef, 0)

    # The best-known upper and lower bounds
    Z_UB = Inf
    Z_LB = -Inf

    # The best-known feasible solutions
    xp_best = zeros(Float64,n)
    xr_best = zeros(Float64,n)
    yp_best = zeros(Float64,n)
    yr_best = zeros(Float64,n)
    sp_best = zeros(Float64,n)
    sr_best = zeros(Float64,n)

    cp_m = zeros(Float64,n)
    cr_m = zeros(Float64,n)
    hp_m = zeros(Float64,n)
    hr_m = zeros(Float64,n)

    # Initial multiplier
    lambdap = zeros(Float64,n)
    lambdar = zeros(Float64,n)
    #print(lambdap)

    gradp = zeros(Float64,n)

    SD = zeros(Int,n,n)
    SR = zeros(Int,n,n)
 
    for i=1:n
        SD[i,i] = d[i]
        SR[i,i] = r[i]
        for j=(i+1):n
            SD[i,j] = SD[i,j-1] + d[j]
            SR[i,j] = SR[i,j-1] + d[j]
        end
    end

    for k=1:MAX_ITER
      
      println("++++++++++++++++++ iteracao $k ++++++++++++++++++")
      
      # Obtaining the lower and upper bounds
      Z_D, xp_d, yp_d, sp_d, xr_d, yr_d, sr_d = lower_bound()
      Z, xp_s, yp_s, sp_s, xr_s, yr_s, sr_s = upper_bound()
      println("upper bound= $Z")
      println("lower bound = $Z_D")

      # Updating the upper bound
      if Z < Z_UB
          Z_UB = Z
          xp_best = xp_s
          yp_best = yp_s 
          sp_best = sp_s 
          xr_best = xr_s 
          yr_best = yr_s
          sr_best = sr_s 
      end

      # Updating the lower bound
      if Z_D > Z_LB
          Z_LB = Z_D
      end

      # Adding the bounds from the current iteration to the record
      push!(UB, Z)
      push!(LB, Z_D)

      # subgradiente
      for t=1:n
          gradp[t] = xp_best[t] - SD[t,n]*yp_best[t]
      end

      # Qgrad
      Qgrad = 0
      for t=1:n
          Qgrad += gradp[t] * gradp[t]
      end

      # Determining the step size and updating the multiplier
      theta = 1.0
      t = theta * (Z_UB - Z_D) / Qgrad
      println("Qgrad = $Qgrad, Z_UB = $Z_UB, Z_D = $Z_D, t = $t")

      for t=1:n
          lambdap = max(0,t*gradp[t])
      end

      # Computing the optimality gap
      opt_gap = (Z_UB-Z_LB) / Z_UB
      
      println("LB = $(round(Z_LB,digits=2)) || UB = $(Z_UB) || gap = $(round(opt_gap,digits=4))) || t = $(round(t,digits=2)))")
      
      if opt_gap < 0.000001
        println("opt gap small")
        break
      end
      
      if t < 0.000001
        println("step t small")
        break
      end

    end

    return Z_UB, xp_best, yp_best, sp_best, xr_best, yr_best, sr_best, UB, LB
    
end


# Finding the exact optimal solution
Z_opt, xp_opt, yp_opt, sp_opt, xr_opt, yr_opt, yr_opt = optimal(n)

# Finding a solution by Lagrangian relaxation
Z_UB, xp_best, yp_best, sp_best, xr_best, yr_best, sr_best, UB, LB = lagrangian_relaxation()

#println("Z_opt = $(round(Z_opt,digits=2))")
#println("Z_UB = $(Z_UB)")
