import gurobipy as gp
from gurobipy import GRB

MAX_CPU_TIME = 3600.0
EPSILON = 1e-6

# Define a callback function
def mycallback(model, where):
    if where == GRB.Callback.POLLING:
        # ignore POLLING callback
        pass
    elif where == GRB.Callback.PRESOLVE:
        # Presolve callback
        pass # Ignore presolve callback
        cdels = model.cbGet(GRB.Callback.PRE_COLDEL)
        rdels = model.cbGet(GRB.Callback.PRE_ROWDEL)
        #if cdels or rdels:
        #    print(f'{cdels} columns and {rdels} rows are removed')
    elif where == GRB.Callback.SIMPLEX:
        # Simplex callback
        pass # ignore SIMPLEX callback
        itcnt = model.cbGet(GRB.Callback.SPX_ITRCNT)
        if itcnt - model._lastiter >= 100:
            model._lastiter = itcnt
            obj = model.cbGet(GRB.Callback.SPX_OBJVAL)
            ispert = model.cbGet(GRB.Callback.SPX_ISPERT)
            pinf = model.cbGet(GRB.Callback.SPX_PRIMINF)
            dinf = model.cbGet(GRB.Callback.SPX_DUALINF)
            if ispert == 0:
                ch = ' '
            elif ispert == 1:
                ch = 'S'
            else:
                ch = 'P'
            print(f'Simplex callback: {itcnt}, {obj}, {ch}, {pinf}, {dinf}')
    elif where == GRB.Callback.MIP:
        # General MIP callback
        pass # Ignore MIP callback
        nodecnt = model.cbGet(GRB.Callback.MIP_NODCNT)
        objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
        objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
        solcnt = model.cbGet(GRB.Callback.MIP_SOLCNT)
        #if nodecnt - model._lastnode >= 100:
        #    model._lastnode = nodecnt
        #    actnodes = model.cbGet(GRB.Callback.MIP_NODLFT)
        #    itcnt = model.cbGet(GRB.Callback.MIP_ITRCNT)
        #    cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)
        #    print(f'{nodecnt}, {actnodes}, {itcnt}, {objbst}, {objbnd}, {solcnt}, {cutcnt}')
        #if abs(objbst - objbnd) < 0.1 * (1.0 + abs(objbst)):
        #    print('Stop early - 10% gap achieved')
        #    model.terminate()
        #if nodecnt >= 10000 and solcnt:
        #    print('Stop early - 10000 nodes explored')
        #    model.terminate()
        if solcnt == 1:
            runtime = model.cbGet(GRB.Callback.RUNTIME)
            print(f'Solution: node {nodecnt}, obj {round(objbst,2)}, bound {round(objbnd,2)}, numsol {solcnt}, time {round(runtime,4)}')
            # store solution
            model._first_incumbent = objbst
            model._time_incumbent = runtime
            model.terminate()
    elif where == GRB.Callback.MIPSOL:
        # MIP solution callback
        #pass # ignore MIPSOL callback
        nodecnt = model.cbGet(GRB.Callback.MIPSOL_NODCNT)
        obj = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        solcnt = model.cbGet(GRB.Callback.MIPSOL_SOLCNT)
        x = model.cbGetSolution(model._vars)
        runtime = model.cbGet(GRB.Callback.RUNTIME)
        if solcnt == 1:
            print(f'Solution: node {nodecnt}, obj {obj}, numsol {solcnt}, time {round(runtime,2)}')
            print('Stop early - first incumbent solution')
            model.terminate()
    elif where == GRB.Callback.MIPNODE:
        # MIP node callback
        pass # ignore MIPNODE callback
        #print('**** New node ****')
        if model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
            x = model.cbGetNodeRel(model._vars)
            model.cbSetSolution(model.getVars(), x)
    elif where == GRB.Callback.BARRIER:
        # Barrier callback
        pass # ignore BARRIER callback
        itcnt = model.cbGet(GRB.Callback.BARRIER_ITRCNT)
        primobj = model.cbGet(GRB.Callback.BARRIER_PRIMOBJ)
        dualobj = model.cbGet(GRB.Callback.BARRIER_DUALOBJ)
        priminf = model.cbGet(GRB.Callback.BARRIER_PRIMINF)
        dualinf = model.cbGet(GRB.Callback.BARRIER_DUALINF)
        cmpl = model.cbGet(GRB.Callback.BARRIER_COMPL)
        print(f'Barrier callback: {itcnt}, {primobj}, {dualobj}, {priminf}, {dualinf}, {cmpl}')
    elif where == GRB.Callback.MESSAGE:
        # Message callback
        pass #ignore MESSAGE callback
        msg = model.cbGet(GRB.Callback.MSG_STRING)
        model._logfile.write(msg)


def clsr_std(N, PP, PR, FP, FR, HR, HP, D, R, SD, SR, C):
    
    try:

        # create model
        model = gp.Model("clsr_std")

        # define variables
        xp = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="xp")
        yp = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yp")
        sp = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="sp")
        xr = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="xr")
        yr = model.addVars(list(range(N)), lb=0.0, ub=1.0, vtype=GRB.BINARY, name="yr")
        sr = model.addVars(list(range(N)), lb=0.0, ub=float('inf'), vtype=GRB.CONTINUOUS, name="sr")

        model.update()

        # set objective
        model.setObjective(gp.quicksum(
            PP[i]*xp[i] + PR[i]*xr[i] + 
            HP[i]*sp[i] + HR[i]*sr[i] + 
            FP[i]*yp[i] + FR[i]*yr[i] for i in range(N)), sense = GRB.MINIMIZE)

        # add constraints
        model.addConstr(xp[0] + xr[0] - sp[0] == D[0])

        model.addConstrs(sp[i-1] + xp[i] + xr[i] - sp[i] == D[i] for i in range(N) if i > 0 )

        model.addConstr(R[0] - xr[0] - sr[0] == 0)

        model.addConstrs(sr[i-1] + R[i] - xr[i] - sr[i] == 0 for i in range(N) if i > 0)

        model.addConstrs(xp[i] - (min(SD[i][N-1],C))*yp[i] <= 0 for i in range(N))

        model.addConstrs(xr[i] - (min(SR[0][i],SD[i][N-1],C))*yr[i] <= 0 for i in range(N))

        model.addConstrs(xp[i] + xr[i] <= C for i in range(N))
        
        # export .lp
        #model.write(file_name+"_model.lp")

        # parameters 
        model.setParam(GRB.Param.TimeLimit, MAX_CPU_TIME)
        model.setParam(GRB.Param.MIPGap, EPSILON)
        model.setParam(GRB.Param.Threads, 1)
        #model.setParam(GRB.Param.Cuts, -1)
        #model.setParam(GRB.Param.Presolve, -1)
        model.setParam(GRB.Param.SolutionLimit, 1)

        model.setParam('OutputFlag', 0)
        #model.setParam('Heuristics', 0)
        model.setParam("LogToConsole", 0)
        
        logfile = open('callback.log', 'w')

        # Pass data into my callback function
        model._lastiter = -GRB.INFINITY
        model._lastnode = -GRB.INFINITY
        model._logfile = logfile
        
        # Initialize a custom attribute to store variables for the callback
        model._vars = model.getVars()

        # optimize model
        #model.optimize(mycallback)
        model.optimize()
        
        status = 0
        if model.status == GRB.OPTIMAL:
            status = 1

        objval = model.ObjVal
        objbound = model.ObjBound
        mipgap = model.MIPGap
        runtime = model.Runtime
        nodecount = model.NodeCount

        #print(f"obj: {model.ObjVal}")

        # After optimization, you can access the stored first incumbent solution        
        if hasattr(model, '_first_incumbent'):
            print(f"objval: {model._first_incumbent}")
            print(f"time: {model._time_incumbent}")
        else:
            print("\nNo incumbent solution found during optimization.")
                
        logfile.close()
        
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    return objval, objbound, mipgap, runtime, nodecount, status
