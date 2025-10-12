import gurobipy as gp
from gurobipy import GRB
import numpy as np

MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

def clsr_std_math(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C, yp_sol, yr_sol):

	yp_val = np.zeros(N)
	yr_val = np.zeros(N)

	try:

		# create model
		model = gp.Model("clsr")

		# create variables
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		sp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="sp")
		sr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="sr")
		yp = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yr")

		# fix variables yp, yr
		for i in range(N):
			yp[i].start = yp_sol[i]
			yr[i].start = yr_sol[i]
		
		model.update()
		
		# set objective
		model.setObjective(gp.quicksum(
			PP[i]*xp[i] + sp[i]*HP[i] + yp[i]*FP[i] + 
			xr[i]*PR[i] + sr[i]*HR[i] + yr[i]*FR[i] for i in range(N)), sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(xp[0] + xr[0] - sp[0] == D[0])

		model.addConstrs(sp[i-1] + xp[i] + xr[i] - sp[i] == D[i] for i in range(N) if i > 0 )
		
		model.addConstr(R[0] - xr[0] - sr[0] == 0)
		
		model.addConstrs(sr[i-1] + R[i] - xr[i] - sr[i] == 0 for i in range(N) if i > 0)
		
		model.addConstrs(xp[i] - yp[i]*min(C,SD[i][N-1]) <= 0 for i in range(N))
		
		model.addConstrs(xr[i] - yr[i]*min(SR[0][i], SD[i][N-1],C) <= 0 for i in range(N))
		
		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))
	   
	    #model.write(file_name+"_model.lp")

		# set parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads,1)
		#model.setParam(GRB.Param.Cuts, -1)
		#model.setParam(GRB.Param.Presolve,-1)

		# relax model
		#for v in model.getVars():
		#	v.setAttr('vtype', 'C')

		# optimize model
		model.optimize()

		tmp = 0
		if model.status == GRB.OPTIMAL:
			tmp = 1

		objval = model.ObjVal
		objbound = model.ObjBound
		mipgap = model.MIPGap
		runtime = model.Runtime
		nodecount = model.NodeCount
	
	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, objbound, mipgap, runtime, nodecount, tmp

	
def clsr_std_mip(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):
	try:

		# create model
		model = gp.Model("clsr_std_mip")

		# create variables
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		sp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="sp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")
		sr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="sr")
		
		model.update()

		# set objective
		model.setObjective(gp.quicksum(
			PP[i]*xp[i] + HP[i]*sp[i] + FP[i]*yp[i] +  
			PR[i]*xr[i] + HR[i]*sr[i] + yr[i]*FR[i] for i in range(N)), sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(xp[0] + xr[0] - sp[0] == D[0])

		model.addConstrs(sp[i-1] + xp[i] + xr[i] - sp[i] == D[i] for i in range(N) if i > 0 )
		
		model.addConstr(R[0] - xr[0] - sr[0] == 0)
		
		model.addConstrs(sr[i-1] + R[i] - xr[i] - sr[i] == 0 for i in range(N) if i > 0)
		
		model.addConstrs(xp[i] - yp[i]*min(SD[i][N-1],C) <= 0 for i in range(N))
		
		model.addConstrs(xr[i] - yr[i]*min(SR[0][i],SD[i][N-1],C) <= 0 for i in range(N))
		
		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))
	    
		# export .lp
		#model.write(file_name+"_model.lp")

		# parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads, 1)
		#model.setParam(GRB.Param.Cuts, -1)
		#model.setParam(GRB.Param.Presolve, -1)

		# relax model
		#for v in model.getVars():
		#	v.setAttr('vtype', 'C')

		# optimize model
		model.optimize()
		
		tmp = 0
		if model.status == GRB.OPTIMAL:
			tmp = 1

		#xp_sol = [xp[i].X for i in range(N)]
		#xr_sol = [xr[i].X for i in range(N)]
		#sp_sol = [sp[i].X for i in range(N)]
		#sr_sol = [sr[i].X for i in range(N)]
		#yp_sol = [yp[i].X for i in range(N)]
		#yr_sol = [yr[i].X for i in range(N)]
	
		objval = model.ObjVal
		objbound = model.ObjBound 
		mipgap = model.MIPGap
		runtime = model.Runtime
		nodecount = model.NodeCount 

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, objbound, mipgap, runtime, nodecount, tmp


def clsr_std_lp(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):
	try:

		#desvio = [0]*N
		#for i in range(N):
		#	desvio[i] = C - D[i]

		#sumDesvio = sum(desvio)
		#sumD = sum(D)

		# create model
		model = gp.Model("clsr_std_lp")

		# create variables
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		sp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="sp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")
		sr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="sr")
		
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

		model.addConstrs(xp[i] - yp[i]*min(SD[i][N-1], C) <= 0 for i in range(N))

		model.addConstrs(xr[i] - yr[i]*min(SR[0][i],SD[i][N-1], C) <= 0 for i in range(N))

		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))

		# relax model
		for v in model.getVars():
			v.setAttr('vtype', 'C')
	    
		# export .lp
		#model.write(file_name+"_model.lp")

		# parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads, 1)
		#model.setParam(GRB.Param.Cuts, -1)
		#model.setParam(GRB.Param.Presolve, -1)

		# optimize model
		model.optimize()
				
		objval = model.ObjVal
		runtime = model.Runtime

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, runtime
