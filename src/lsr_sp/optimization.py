import gurobipy as gp
from gurobipy import GRB
import numpy as np

MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

def clsr_sp_math(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C, yp_sol, yr_sol):

	yp_val = np.zeros(N)
	yr_val = np.zeros(N)

	CSP = (np.zeros((N,N))).tolist()
	CSR = (np.zeros((N,N))).tolist()
	CR = (np.zeros((N,N))).tolist()
	CL = np.zeros(N)

	for i in range(N):
		for j in range(i,N):
			CR[i][j] = sum(HR[t]*SR[i][t] for t in range(i,j))
			CSP[i][j] = PP[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))
			CSR[i][j] = PR[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))

	for i in range(N):
		CL[i] = sum(HR[j]*SR[i][j] for j in range(i,N))

	try:

		# create model
		model = gp.Model("clsr")

		# create variables
		zsp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_sp")
		zsr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_sr")
		zr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_r")
		l = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="l")
		yp = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.BINARY, name="yr")

		# fix variables yp, yr
		for i in range(N):
			yp[i].start = yp_sol[i]
			yr[i].start = yr_sol[i]

		# set objective
		#obj = 0.0
		#for i in range(N):
		#	obj += yp[i]*FP[i] + yr[i]*FR[i] + l[i]*CL[i] 
		#   obj += gp.quicksum(zsp[i,j]*CSP[i][j] + zsr[i,j]*CSR[i][j] + zr[i,j]*CR[i][j] for j in range(i,N))
		
		obj = 0.0
		for i in range(N):
			obj += yp[i]*FP[i] + yr[i]*FR[i] + l[i]*CL[i] 
			for j in range(i,N):
				obj += zsp[i,j]*CSP[i][j] + zsr[i,j]*CSR[i][j] + zr[i,j]*CR[i][j] 

		model.setObjective(obj, sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(gp.quicksum(zsp[0,j] + zsr[0,j] for j in range(N)) == 1)
		
		model.addConstrs(
			gp.quicksum(zsp[i,t-1] + zsr[i,t-1] for i in range(t)) - 
			gp.quicksum(zsp[t,j] + zsr[t,j] for j in range(t, N)) == 0  for t in range(1,N) 
			)
				
		model.addConstrs(gp.quicksum(zsp[t,j] for j in range(t,N)) <= yp[t] for t in range(N))
			
		model.addConstrs(gp.quicksum(zsr[t,j] for j in range(t,N)) <= yr[t] for t in range(N))
			
		model.addConstr(gp.quicksum(zr[0,j] for j in range(N)) + l[0] == 1)
					
		model.addConstrs(
			gp.quicksum(zr[i,t-1] for i in range(0,t)) == 
			gp.quicksum(zr[t,j]  for j in  range(t,N)) + l[t] for t in range(1,N)
			)       
				
		model.addConstrs(gp.quicksum(zr[i,t] for i in range(0,t+1)) <= yr[t] for t in range(N))
			
		model.addConstrs(
			gp.quicksum(SR[i][t]*zr[i,t] for i in range(t+1) ) ==
			gp.quicksum(SD[t][j]*zsr[t,j] for j in range(t,N)) for t in range(N)
			)

       #model.addConstrs(gp.quicksum(SD[i][t]*zsp[i,t] for t in range(i,N)) <= yp[i]*min(C,SD[i][N-1]) for i in range(N))

	   #model.addConstrs(
	   #     gp.quicksum(SD[i][t]*zsp[i,t] for t in range(i,N)) <= yp[i]*SD[i][N-1] for i in range(N)
	   #    )
		
		#model.addConstrs(
		#    gp.quicksum(SD[i][t]*zsr[i,t] for t in range(i,N)) <= yr[i]*SD[i][N-1] for i in range(N)
		#   )
		
		model.addConstrs(
			gp.quicksum(SD[t][k]*zsp[t,k] for k in range(t,N)) + 
			gp.quicksum(SD[t][k]*zsr[t,k] for k in range(t,N)) <= C for t in range(N)
			)

		# Parameters 
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
	

def clsr_sp_mip(N, PP, PR, FP, FR, HP, HR, D, R, SD,SR,C):

	CSP = (np.zeros((N,N))).tolist()
	CSR = (np.zeros((N,N))).tolist()
	CR = (np.zeros((N,N))).tolist()
	CL = np.zeros(N)

	for i in range(N):
		for j in range(i,N):
			CR[i][j]  = sum(HR[t]*SR[i][t] for t in range(i,j))
			CSP[i][j] = PP[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))
			CSR[i][j] = PR[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))

	for i in range(N):
		CL[i] = sum(HR[j]*SR[i][j] for j in range(i,N))	

	try:
		# create model
		model = gp.Model("clsr_sp_mip")

		# add variables		
		zsp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="zsp")
		zsr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="zsr")
		zr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="zr")
		l = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="l")
		yp = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.BINARY, name="yr")

		# set objective
		#obj = 0.0
		#for i in range(N):
		#	obj += yp[i]*FP[i] + yr[i]*FR[i] + l[i]*CL[i]
		#	obj += gp.quicksum(zsp[i,j]*CSP[i][j] + zsr[i,j]*CSR[i][j] + zr[i,j]*CR[i][j] for j in range(i,N))

		obj = 0.0
		for i in range(N):
			obj += yp[i]*FP[i] + yr[i]*FR[i] + l[i]*CL[i] 
			for j in range(i,N):
				obj += zsp[i,j]*CSP[i][j] + zsr[i,j]*CSR[i][j] + zr[i,j]*CR[i][j] # 16

		model.setObjective(obj, sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(gp.quicksum(zsp[0,j] + zsr[0,j] for j in range(N)) == 1)
		
		model.addConstrs(
			gp.quicksum(zsp[i,t-1] + zsr[i,t-1] for i in range(t)) - 
			gp.quicksum(zsp[t,j] + zsr[t,j] for j in range(t, N)) == 
			0  for t in range(1,N) 
		)
				
		model.addConstrs(
			gp.quicksum(zsp[t,j] for j in range(t,N)) <= 
			yp[t] for t in range(N)
		)
			
		model.addConstrs(
			gp.quicksum(zsr[t,j] for j in range(t,N)) <= 
			yr[t] for t in range(N)
		)
			
		model.addConstr(gp.quicksum(zr[0,j] for j in range(N)) + l[0] == 1)
					
		model.addConstrs(
			gp.quicksum(zr[i,t-1] for i in range(0,t)) == 
			gp.quicksum(zr[t,j]  for j in  range(t,N)) + l[t] for t in range(1,N)
		)       
				
		model.addConstrs(
			gp.quicksum(zr[i,t] for i in range(0,t+1)) <= 
			yr[t] for t in range(N)
		)    
			
		model.addConstrs(
			gp.quicksum(SR[i][t]*zr[i,t] for i in range(t+1) ) ==
			gp.quicksum(SD[t][j]*zsr[t,j] for j in range(t,N)) for t in range(N)
		)

		model.addConstrs(
		    gp.quicksum(SD[i][t]*zsp[i,t] for t in range(i,N)) <= 
		    yp[i]*min(SD[i][N-1],C) for i in range(N)
		)

		#model.addConstrs(
		#	gp.quicksum(SD[i][t]*zsp[i,t] for t in range(i,N)) <= 
		#	yp[i]*SD[i][N-1] for i in range(N)
		#)
		
		#model.addConstrs(
		#	gp.quicksum(SD[i][t]*zsr[i,t] for t in range(i,N)) <= 
		#	yr[i]*min(SR[0][i],SD[i][N-1],C) for i in range(N)
		#)

		#model.addConstrs(
		#	gp.quicksum(SD[i][t]*zsr[i,t] for t in range(i,N)) <= 
		#	yr[i]*SD[i][N-1] for i in range(N)
		#)
		
		#model.addConstrs(
		#	gp.quicksum(SD[t][k]*zsp[t,k] for k in range(t,N)) + 
		#	gp.quicksum(SD[t][k]*zsr[t,k] for k in range(t,N)) <= 
		#	C for t in range(N)
		#)
		
		# parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads, 1)
		#model.setParam(GRB.Param.Cuts, -1)
		#model.setParam(GRB.Param.Presolve, -1)

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


def clsr_sp_lp(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C):

	CSP = (np.zeros((N,N))).tolist()
	CSR = (np.zeros((N,N))).tolist()
	CR = (np.zeros((N,N))).tolist()
	CL = np.zeros(N)

	for i in range(N):
		for j in range(i,N):
			CR[i][j]  = sum(HR[t]*SR[i][t] for t in range(i,j))
			CSP[i][j] = PP[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))
			CSR[i][j] = PR[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))

	for i in range(N):
		CL[i] = sum(HR[j]*SR[i][j] for j in range(i,N))

	try:

		# create model
		model = gp.Model("clsrp_sp_lp")

		# create variables
		zsp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_sp")
		zsr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_sr")
		zr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_r")
		l = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="l")
		yp = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="yp")
		yr = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="yr")

		# set objective
		obj = 0.0
		for i in range(N):
			obj += yp[i]*FP[i] + yr[i]*FR[i] + l[i]*CL[i] 
			for j in range(i,N):
				obj += zsp[i,j]*CSP[i][j] + zsr[i,j]*CSR[i][j] + zr[i,j]*CR[i][j] 

		model.setObjective(obj, sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(gp.quicksum(zsp[0,j] + zsr[0,j] for j in range(N)) == 1)
		
		model.addConstrs(
			gp.quicksum(zsp[i,t-1] + zsr[i,t-1] for i in range(t)) - 
			gp.quicksum(zsp[t,j] + zsr[t,j] for j in range(t, N)) == 0  for t in range(1,N)
		)
				
		model.addConstrs(gp.quicksum(zsp[t,j] for j in range(t,N)) <= yp[t] for t in range(N))
			
		model.addConstrs(gp.quicksum(zsr[t,j] for j in range(t,N)) <= yr[t] for t in range(N))
			
		model.addConstr(gp.quicksum(zr[0,j] for j in range(N)) + l[0] == 1)
					
		model.addConstrs(
			gp.quicksum(zr[i,t-1] for i in range(0,t)) == 
			gp.quicksum(zr[t,j]  for j in  range(t,N)) + l[t] for t in range(1,N)
		)       
				
		model.addConstrs(gp.quicksum(zr[i,t] for i in range(0,t+1)) <= yr[t] for t in range(N))    
			
		model.addConstrs(
			gp.quicksum(SR[i][t]*zr[i,t] for i in range(t+1) ) ==
			gp.quicksum(SD[t][j]*zsr[t,j] for j in range(t,N)) for t in range(N)
		)

		#model.addConstrs(
		#	gp.quicksum(SD[i][t]*zsp[i,t] for t in range(i,N)) <= 
		#	yp[i]*SD[i][N-1] for i in range(N)
		#)

		#model.addConstrs(
		#	gp.quicksum(SD[i][t]*zsr[i,t] for t in range(i,N)) <= 
		#	yr[i]*min(SR[0][i],SD[i][N-1]) for i in range(N)
	    #)

		#model.addConstrs(
		#    gp.quicksum(SD[i][t]*zsp[i,t] for t in range(i,N)) <= 
		#    yp[i]*min(SD[i][N-1],C) for i in range(N)
		#)
		
		#model.addConstrs(
		#	gp.quicksum(SD[i][t]*zsr[i,t] for t in range(i,N)) <= 
		#	yr[i]*min(SR[0][i],SD[i][N-1],C) for i in range(N)
		#)


		model.addConstrs(
			gp.quicksum(SD[t][k]*zsp[t,k] for k in range(t,N)) +
			gp.quicksum(SD[t][k]*zsr[t,k] for k in range(t,N)) <= C for t in range(N)
		)
		
		# parameters 
		model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap, EPSILON)
		model.setParam(GRB.Param.Threads,1)
		#model.setParam(GRB.Param.Cuts, -1)
		#model.setParam(GRB.Param.Presolve,-1)

		# optimize model
		model.optimize()
		
		objval = model.ObjVal
		runtime = model.Runtime

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, runtime
