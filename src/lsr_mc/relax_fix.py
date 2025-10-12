import gurobipy as gp
from gurobipy import GRB
import numpy as np
MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

def relax_fix(partp, partr, yp_sol , yr_sol, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):

	#xp_sol1 = [0]*N
	#xr_sol1 = [0]*N
	#wp_sol1 = (np.zeros((N,N))).tolist()
	#wr_sol1 = (np.zeros((N,N))).tolist()
	#vor_sol1 = (np.zeros((N,N))).tolist()
	yp_val = np.zeros(N)
	yr_val = np.zeros(N)

	#objval = 0.0
	#bestbound = 0.0
	#numnode = 0.0
	#elapsed = 0.0
	#gap = 0.0

	CP = np.zeros(N)
	CR = np.zeros(N)
	KP = 0
	KR = 0 

	for i in range(N):
		CP[i] = PP[i]

		for j in range(i,N):
			CP[i] = CP[i] + HP[j]

	for i in range(N):
		CR[i] = PR[i]
		for j in range(i,N):
			CR[i] += HP[j]
		for j in range(i,N):
			CR[i] -= HR[j]

	for i in range(N):
		aux = 0 
		for j in range(i+1):
			aux += D[j]

		KP += HP[i]*aux

	for i in range(N):
		aux = 0

		for j in range(i+1):
			aux += R[j]

		KR += HR[i]*aux

	K = KR - KP

	try:
		# create model
		model = gp.Model("clsr")

		# create variables
		wp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wp")
		wr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wr")
		vor = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="vor")
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")
		
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
		#obj = gp.quicksum(xp[i]*CP[i] for i in range(N)) 
		#obj += gp.quicksum(yp[i]*FP[i] for i in range(N))
		#obj += gp.quicksum(xr[i]*CR[i] for i in range(N)) 
		#obj += gp.quicksum(yr[i]*FR[i] for i in range(N)) 
		#obj += K
		
		obj = 0.0
		for i in range(N):
			obj += xp[i]*CP[i]
			obj += yp[i]*FP[i]
			obj += xr[i]*CR[i]
			obj += yr[i]*FR[i]
		obj += K

		# add constraints
		model.setObjective(obj, sense = GRB.MINIMIZE)
		
		# add constraints
		for i in range(N):
			ctr = 0.0
			for j in range(i+1):
				ctr += wp[j,i]
				ctr += wr[j,i]
			model.addConstr(ctr >= D[i])

		for i in range(N):
			ctr = 0.0
			for j in range(i+1):
				ctr+= vor[j,i]
			for j in range(i,N):
				ctr+=(-wr[i,j])
			model.addConstr(ctr == 0)

		for i in range(N):
			ctr = 0.0
			for j in range(i,N):
				ctr += vor[i,j]
			model.addConstr(ctr <= R[i])

		for i in range(N):
			for j in range(i,N):
				model.addConstr(wp[i,j] + yp[i]*(-D[j]) <= 0)

		for i in range(N):
			for j in range(i,N):
				model.addConstr(wr[i,j] + yr[i]*(-min(SR[0][i],D[j])) <= 0)

		for i in range(N):
			for j in range(i,N):
				model.addConstr(vor[i,j] + yr[j]*(-R[i]) <= 0)

		for i in range(N):
			ctr = 0.0 
			ctr += xp[i]
			for j in range(i,N):
				ctr += (-wp[i,j])
			model.addConstr(ctr == 0)

		for i in range(N):
			ctr = 0.0
			ctr += xr[i]
			for j in range(i,N):
				ctr += (-wr[i,j])
			model.addConstr(ctr == 0)

		#model.addConstrs(xp[i] - yp[i]*min(C,SD[i][N-1]) <= 0 for i in range(N))

		#model.addConstrs(xr[i] - yr[i]*min(SR[0][i],SD[i][N-1],C) <= 0 for i in range(N))

		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))
		
		#model.write(file_name+".lp")

		# Parameters 
		model.setParam(GRB.Param.TimeLimit,MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap,EPSILON)
		model.setParam(GRB.Param.Threads,1)
		#model.setParam(GRB.Param.Cuts,-1)
		#model.setParam(GRB.Param.Presolve,-1)

		# optimize model
		model.optimize()

		objval = model.ObjVal
		#objbound = model.ObjBound
		#nodecount = model.NodeCount
		#runtime = model.Runtime
		#mipgap = model.MIPGap
	
		yp_val = [yp[i].X for i in range(N)]
		yr_val = [yr[i].X for i in range(N)]

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, yp_val, yr_val
