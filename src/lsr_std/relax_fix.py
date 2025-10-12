import gurobipy as gp
from gurobipy import GRB
import numpy as np

MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

def relax_fix(partp, partr, yp_sol , yr_sol, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):

	yp_val = np.zeros(N)
	yr_val = np.zeros(N)

	try:
		# create model
		model = gp.Model("clsr")

		# create variables
		xp = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="xp")
		yp = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yp")
		sp = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="sp")
		xr = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="xr")
		yr = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yr")
		sr = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="sr")
		
		#if (partp is not None):
		#	print(f"partp_rf: {partp, max(partp)}")
			
		#if (partr is not None):
		#	print(f"partr_rf: {partr, max(partr)}")
			
		#input("\ntecle <enter> para continuar")
		
		for i in range(N):
			if (partp is not None) and (i > max(partp)):
				yp[i].VType = gp.GRB.CONTINUOUS
				yp[i].lb = 0.0
				yp[i].ub = 1.0
			elif (partp is not None) and (i in partp):
				yp[i].lb = 0.0
				yp[i].ub = 1.0
			else:
				yp[i].lb = yp_sol[i]
				yp[i].ub = yp_sol[i]

		for i in range(N):
			if (partr is not None) and (i > max(partr)):
				yr[i].VType = gp.GRB.CONTINUOUS
				yr[i].lb = 0.0
				yr[i].ub = 1.0
			elif (partr is not None) and (i in partr):
				yr[i].lb = 0.0
				yr[i].ub = 1.0
			else:
				yr[i].lb = yr_sol[i]
				yr[i].ub = yr_sol[i]
		
		model.update()

		# set objective
		model.setObjective(gp.quicksum(
			PP[i]*xp[i] + HP[i]*sp[i] + FP[i]*yp[i] +
			PR[i]*xr[i] + HR[i]*sr[i] + FR[i]*yr[i] for i in range(N)), sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(xp[0] + xr[0] - sp[0] == D[0])
		
		model.addConstrs(sp[i-1] + xp[i] + xr[i] - sp[i] == D[i] for i in range(N) if i > 0 )
		
		model.addConstr(R[0] - xr[0] - sr[0] == 0)
		
		model.addConstrs(sr[i-1] + R[i] - xr[i] - sr[i] == 0 for i in range(N) if i > 0)
		
		model.addConstrs(xp[i] - yp[i]*min(C,SD[i][N-1]) <= 0 for i in range(N))
		
		model.addConstrs(xr[i] - yr[i]*min(SR[0][i], SD[i][N-1]) <= 0 for i in range(N))
		
		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))
		
		#model.write(file_name+".lp")

		# Parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads,1)
		#model.setParam(GRB.Param.Cuts,-1)
		#model.setParam(GRB.Param.Presolve,-1)

		# Optimize model
		model.optimize()

		objval = model.ObjVal

		yp_val = [yp[i].X for i in range(N)]
		yr_val = [yr[i].X for i in range(N)]

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, yp_val, yr_val
