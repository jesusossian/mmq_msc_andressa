import gurobipy as gp
from gurobipy import GRB
import numpy as np

MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

def fix_and_optimize(partp, partr, yp_sol, yr_sol, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):

    try:
        yp_val = np.zeros(N)
        yr_val = np.zeros(N)

        # create model
        model = gp.Model("CLSR")

        # create variables
        xp = model.addVars(list(range(N)), lb =0.0, ub = float('inf'),vtype=GRB.CONTINUOUS, name="xp")
        yp = model.addVars(list(range(N)), lb =0.0, ub = 1.0,vtype=GRB.BINARY, name="yp")
        sp = model.addVars(list(range(N)), lb =0.0, ub = float('inf'),vtype=GRB.CONTINUOUS, name="sp")
        xr = model.addVars(list(range(N)), lb =0.0, ub = float('inf'),vtype=GRB.CONTINUOUS, name="xr")
        yr = model.addVars(list(range(N)), lb =0.0, ub = 1.0,vtype=GRB.BINARY, name="yr")
        sr = model.addVars(list(range(N)), lb =0.0, ub = float('inf'),vtype=GRB.CONTINUOUS, name="sr")

        #if (partp is not None):
        #    print(f"partp_fo: {partp, max(partp)}")
			
        #if (partr is not None):
        #    print(f"partr_fo: {partr, max(partr)}")
			
        #input("\ntecle <enter> para continuar")

        for i in range(N):
            if (partp is not None) and (i not in partp):
                yp[i].lb = yp_sol[i]
                yp[i].ub = yp_sol[i]
            else:
                yp[i].start = yp_sol[i]

        for i in range(N):
            if (partr is not None) and (i not in partr):
                yr[i].lb = yr_sol[i]
                yr[i].ub = yr_sol[i]
            else:
                yr[i].start = yr_sol[i]
                
        #for i in range(N):
        #    xp[i].start = xp_sol[i]
        #    xr[i].start = xr_sol[i]
        #    sp[i].start = sp_sol[i]
        #    sr[i].start = sr_sol[i]
        
        model.update()
        
        # set objective
        model.setObjective(gp.quicksum(
            PP[i]*xp[i] + HP[i]*sp[i] + FP[i]*yp[i] +
            PR[i]*xr[i] + HR[i]*sr[i] + FR[i]*yr[i] for i in range(N)) , sense = GRB.MINIMIZE)

        # add constraints
        model.addConstr(xp[0] + xr[0] - sp[0] == D[0])
        
        model.addConstrs(sp[i-1] + xp[i] + xr[i] - sp[i] == D[i] for i in range(N) if i > 0 )
        
        model.addConstr(R[0] - xr[0] - sr[0] == 0)
        
        model.addConstrs(sr[i-1] + R[i] - xr[i] - sr[i] == 0 for i in range(N) if i > 0)
        
        model.addConstrs(xp[i] - yp[i]*min(C,SD[i][N-1]) <= 0 for i in range(N))
        
        model.addConstrs(xr[i] - yr[i]*min(SR[0][i], SD[i][N-1],C) <= 0 for i in range(N))
        
        model.addConstrs(xp[i] + xr[i] <= C for i in range(N))
       
        #model.write(file_name+"_model.lp")

        # parameters 
        model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
        model.setParam(GRB.Param.MIPGap, EPSILON)
        model.setParam(GRB.Param.Threads,1)
        #model.setParam(GRB.Param.Cuts, -1)
        #model.setParam(GRB.Param.Presolve,-1)

        # optimize model
        model.optimize()

        yp_val = [yp[i].X for i in range(N)]
        yr_val = [yr[i].X for i in range(N)]

        objval = model.ObjVal

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    return objval, yp_val, yr_val
