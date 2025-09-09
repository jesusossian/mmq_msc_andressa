import gurobipy as gp
from gurobipy import GRB
import numpy as np

MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

lsdbar = 1

def clsr_mc(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C, yp_sol, yr_sol):

	CP = np.zeros(N)
	CR = np.zeros(N)
	KP = 0
	KR = 0 

	for i in range(N):
		CP[i] = PP[i]
		for j in range(i,N):
			CP[i] = CP[i] + HP[j]

#	for l in range(0,N):
#		if l < 1:
#			rrl[l] = R[l]
#		else:
#			rrl[l] = R[l] + max(0.0, rrl[l-1]-D[l-1])
#
#		dl[l] = max(0.0,D[l]-rrl[l])

#	for k in range(0,N):
#		sdl[k][k] = dl[k]
#		for j in range(k+1,N):
#			sdl[k][j] = sdl[k][j-1] + dl[j]

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
		yp = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yr")

		for i in range(N):
			yp[i].start = yp_sol[i]
			yr[i].start = yr_sol[i]

		model.update()

		# set objective
		#obj = 0.0
		#obj += gp.quicksum(xp[i]*CP[i] for i in range(N)) 
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
		
		model.setObjective(obj, sense=GRB.MINIMIZE)

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
			ctr =0.0
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

		# add constraints
		#model.addConstrs(gp.quicksum(wp[j,i]+wr[j,i] for j in range(i+1)) >= D[i] for i in range(N))

		#model.addConstrs(gp.quicksum(vor[j,i] for j in range(i+1))
		#    + gp.quicksum((-wr[i,j]) for j in range(i,N)) >= 0 for i in range(N) )

		#model.addConstrs(gp.quicksum(vor[i,j] for j in range(i,N)) <= R[i] for i in range(N))

		#model.addConstrs(wp[i,j] <= yp[i]*D[j]  for i in range(N) for j in range(i,N) )

		#model.addConstrs(wr[i,j] <= min(SR[0][i],D[i])  for j in range(i,N) for i in range(N) )

		#model.addConstrs(vor[i,j] + yr[j]*(-R[i])  <= 0 for j in range(i,N) for i in range(N))

		#model.addConstrs(xp[i]+ gp.quicksum((-wp[i,j]) for j in range(i,N)) == 0 for i in range(N))

		#model.addConstrs(xr[i]+ gp.quicksum((-wr[i,j]) for j in range(i,N)) == 0 for i in range(N))

		#model.addConstrs(xp[i] - yp[i]*min(C,SD[i][N-1]) <= 0 for i in range(N))
		
		#model.addConstrs(xr[i] - yr[i]*min(SR[0][i], SD[i][N-1]) <= 0 for i in range(N))

		#model.addConstrs(xp[i]+xr[i] <= C for i in range(N))

		#if lsdbar == 1:
		#	model.addConstrs((gp.quicksum(xp[l] for l in range(i)) 
		# 		+ gp.quicksum(yp[l]*sdl[l][j] for l in range(i,j+1))) \
		#		>= sdl[0][j] for i in range(N) for j in range(i,N))

		# parameters 
		model.setParam(GRB.Param.TimeLimit,MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap,EPSILON)
		model.setParam(GRB.Param.Threads,1)
		#model.setParam(GRB.Param.Cuts,-1)
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
	
def clsr_mc_mip(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):

	CP = np.zeros(N)
	CR = np.zeros(N)
	#rrl = np.zeros(N)
	#dl = np.zeros(N)
	#sdl = (np.zeros((N,N))).tolist()
	KP = 0
	KR = 0 

	for i in range(N):
		CP[i] = PP[i]
		for j in range(i,N):
			CP[i] = CP[i] + HP[j]

#	for l in range(0,N):
#		if l < 1:
#			rrl[l] = R[l]
#		else:
#			rrl[l] = R[l] + max(0.0, rrl[l-1]-D[l-1])
#
#		dl[l] = max(0.0,D[l]-rrl[l])
#
#	for k in range(0,N):
#		sdl[k][k] = dl[k]
#		for j in range(k+1,N):
#			sdl[k][j] = sdl[k][j-1] + dl[j]

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
		model = gp.Model("clsr_mc_mip")

		# create variables		
		wp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wp")
		wr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wr")
		vor = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="vor")
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")

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
		
		model.setObjective(obj, sense = GRB.MINIMIZE)

		for i in range(N): # 37
			ctr = 0.0
			for j in range(i+1):
				ctr += wp[j,i]
				ctr += wr[j,i]
			model.addConstr(ctr >= D[i])

		for i in range(N): # 38
			ctr = 0.0
			for j in range(i+1):
				ctr+= vor[j,i]
			for j in range(i,N):
				ctr+=(-wr[i,j])
			model.addConstr(ctr == 0)

		for i in range(N): #39
			ctr =0.0
			for j in range(i,N):
				ctr += vor[i,j]
			model.addConstr(ctr <= R[i])

		for i in range(N): # 40
			for j in range(i,N):
				model.addConstr(wp[i,j] + yp[i]*(-D[j]) <= 0)

		for i in range(N): # 41
			for j in range(i,N):
				model.addConstr(wr[i,j] + yr[i]*(-min(SR[0][i],D[j])) <= 0)

		for i in range(N): # 42
			for j in range(i,N):
				model.addConstr(vor[i,j] + yr[j]*(-R[i]) <= 0)

		for i in range(N):
			ctr = xp[i]
			for j in range(i,N):
				ctr += (-wp[i,j])
			model.addConstr(ctr == 0)

		for i in range(N): # 44
			#ctr = 0.0
			ctr = xr[i]
			for j in range(i,N):
				ctr += (-wr[i,j])
			model.addConstr(ctr == 0)

		#model.addConstrs(xp[i] <= yp[i]*SD[i][N-1] for i in range(N))

		#model.addConstrs(xr[i] <= yr[i]*min(SR[0][i],SD[i][N-1]) for i in range(N))

		model.addConstrs(xp[i] <= yp[i]*min(SD[i][N-1],C) for i in range(N))

		model.addConstrs(xr[i] <= yr[i]*min(SR[0][i],SD[i][N-1],C) for i in range(N))

		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))

		#if lsdbar == 1:
		#	model.addConstrs(
		# 		(gp.quicksum(xp[l] for l in range(i)) 
		# 		+ gp.quicksum(yp[l]*sdl[l][j] for l in range(i,j+1))) >= 
		# 		sdl[0][j] for i in range(N) for j in range(i,N))

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

		# get status
		tmp = 0
		if model.status == GRB.OPTIMAL:
			tmp = 1

		# get value
		objval = model.ObjVal
		objbound = model.ObjBound
		mipgap = model.MIPGap
		runtime = model.Runtime
		nodecount = model.NodeCount

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, objbound, mipgap, runtime, nodecount, tmp


def clsr_mc_lp(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):

	CP = np.zeros(N)
	CR = np.zeros(N)
	#rrl = np.zeros(N)
	#dl = np.zeros(N)
	#sdl = (np.zeros((N,N))).tolist()
	KP = 0
	KR = 0 

	for i in range(N):
		CP[i] = PP[i]
		for j in range(i,N):
			CP[i] = CP[i] + HP[j]

	#for l in range(0,N):
	#	if l < 1:
	#		rrl[l] = R[l]
	#	else:
	#		rrl[l] = R[l] + max(0.0, rrl[l-1]-D[l-1])
	#	dl[l] = max(0.0,D[l]-rrl[l])

	#for k in range(0,N):
	#	sdl[k][k] = dl[k]
	#	for j in range(k+1,N):
	#		sdl[k][j] = sdl[k][j-1] + dl[j]

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
		model = gp.Model("clsr_mc_lp")

		# create variables		
		wp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wp")
		wr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wr")
		vor = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="vor")
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")

		# set objective
		obj = 0.0
		for i in range(N):
			obj += xp[i]*CP[i]
			obj += yp[i]*FP[i]
			obj += xr[i]*CR[i]
			obj += yr[i]*FR[i]
		obj += K
		
		model.setObjective(obj, sense = GRB.MINIMIZE)

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
			ctr =0.0
			for j in range(i,N):
				ctr += vor[i,j]
			model.addConstr(ctr <= R[i])

		for i in range(N):
			for j in range(i,N):
				model.addConstr(wp[i,j] + yp[i]*(-D[j]) <=0)

		for i in range(N):
			for j in range(i,N):
				model.addConstr(wr[i,j] + yr[i]*(-min(SR[0][i],D[j])) <=0)

		for i in range(N):
			for j in range(i,N):
				model.addConstr(vor[i,j] + yr[j]*(-R[i]) <= 0)

		for i in range(N):
			ctr = xp[i]
			for j in range(i,N):
				ctr += (-wp[i,j])
			model.addConstr(ctr == 0)

		for i in range(N):
			ctr = xr[i]
			for j in range(i,N):
				ctr += (-wr[i,j])
			model.addConstr(ctr == 0)

		#model.addConstrs(xp[i] - yp[i]*min(C,SD[i][N-1]) <= 0 for i in range(N))
		
		#model.addConstrs(xr[i] - yr[i]*min(SR[0][i],SD[i][N-1],C) <= 0 for i in range(N))

		model.addConstrs(xp[i] + xr[i] <= C for i in range(N))

		#if lsdbar == 1:
		#	model.addConstrs(
		#		gp.quicksum(xp[l] for l in range(i)) + 
		#		gp.quicksum(yp[l]*sdl[l][j] for l in range(i,j+1))) >= 
		#		sdl[0][j] for i in range(N) for j in range(i,N)
		#)

		# parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads, 1)
		#model.setParam(GRB.Param.Cuts, -1)
		#model.setParam(GRB.Param.Presolve,-1)

		# relax model
		for v in model.getVars():
			v.setAttr('vtype', 'C')

		# optimize model
		model.optimize()

		objval = model.ObjVal
		runtime = model.Runtime

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, runtime
