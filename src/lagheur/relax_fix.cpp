double relax_fix (double *HP, double *HR, double *PP, double *PR,double *FP, double *FR, double *D, double *R, int tam_sub,int iteration, int tam_sub_inter) {


  double C;
  double soma;
  double fator;
  double objval;
  double elapsed_;
  int optimstatus;
  bool cap;
  //int tam_sub_inter;
  int status;
  
  cap = true;
  status = 0;
  //tam_sub_inter = 1;

  /* Capacity */
  soma = 0;
  for (int i=0; i<N; i++) {
    soma += D[i];
  }
  fator = 1.5; //1.5, 1.75, 2.0 
  C = (soma * fator)/N;
  
  try {
    // Create an environment
      GRBEnv env = GRBEnv(true);
      env.set("LogFile", "mip.log");
      env.start();
    
      //Create an empty model
      GRBModel model = GRBModel(env);
    
      //Create variables
      GRBVar sp[N];
      GRBVar xp[N];
      GRBVar yp[N];
      GRBVar sr[N];
      GRBVar xr[N];
      GRBVar yr[N];
    
      for(int i=0; i<N; i++) {
        sp[i] = model.addVar(0.0,GRB_INFINITY,0.0, GRB_CONTINUOUS,"sp("+IntToString(i)+")");
        xp[i] = model.addVar(0.0,GRB_INFINITY,0.0, GRB_CONTINUOUS,"xp("+IntToString(i)+")");
        sr[i] = model.addVar(0.0,GRB_INFINITY,0.0, GRB_CONTINUOUS,"sr("+IntToString(i)+")");
        xr[i] = model.addVar(0.0,GRB_INFINITY,0.0, GRB_CONTINUOUS,"xr("+IntToString(i)+")");

          if (tam_sub_inter == N){
             yp[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yp("+IntToString(i)+")");
             yr[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yr("+IntToString(i)+")");

          }
          else{
             if (i < tam_sub_inter){
                yp[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yp("+IntToString(i)+")");
                yr[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yr("+IntToString(i)+")");
              }
              else{
                yp[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS,"yp("+IntToString(i)+")");
               yr[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS,"yr("+IntToString(i)+")");
         
              }
          }
        
      }
    
      /* objective function */
      GRBLinExpr FO = 0;
      GRBLinExpr expr = 0;
      for(int i=0; i<N; i++) {
        FO += xp[i]*PP[i] + sp[i]*HP[i] + xr[i]*PR[i] + sr[i]*HR[i] + yp[i]*FP[i] + yr[i]*FR[i];
      }
      model.setObjective(FO, GRB_MINIMIZE);
    
      /* constraints */      
      for (int i=0; i<N; i++) {
        expr = 0;
        if (i>0) expr += sp[i-1];
        expr +=  xp[i] + xr[i];
        if (i<N-1) expr -= sp[i];
        model.addConstr(expr == D[i]);
      }
      
      for (int j=0; j <(iteration - 1)*tam_sub; j++ ){
          model.addConstr(yp[j] == yp_int[j]);
          model.addConstr(yr[j] == yr_int[j]);
      }

      for (int i=0; i<N; i++) {
        expr=0;
        if (i>0) expr += sr[i-1];
        expr += R[i] - xr[i] - sr[i] ;
        model.addConstr(expr == 0);
      }
    
      if (cap == true) {
        for (int i=0; i<N; i++) {
          expr = 0;
          expr = xp[i]  - yp[i]*SD[i][N-1];
          model.addConstr(expr <=0);
        }
          
        for (int i=0; i<N; i++) {
          expr = 0;
          expr = xr[i] - yr[i]*min(SR[0][i], SD[i][N-1]);
          model.addConstr(expr <= 0);
        }

        for (int i=0; i<N; i++) {
          expr = 0;
          expr = xp[i] + xr[i] - C;
          model.addConstr(expr <= 0);
        }      
      }
    
      /* Parameters */
      model.set(GRB_DoubleParam_TimeLimit, MAX_CPU_TIME);
      //model.set(GRB_DoubleParam_MIPGap, EPSILON);
      model.set(GRB_IntParam_Threads, 1);
      model.set(GRB_IntParam_Cuts, -1);
      model.set(GRB_IntParam_Presolve, -1);
    

      /* solve problem */
      model.optimize();
    
      /* status problem */
      optimstatus = model.get(GRB_IntAttr_Status);
    
      objval = model.get(GRB_DoubleAttr_ObjVal);
      bestbound = model.get(GRB_DoubleAttr_ObjBound);
      numnode = model.get(GRB_DoubleAttr_NodeCount);
      elapsed = model.get(GRB_DoubleAttr_Runtime);
      //gap = model.get(GRB_DoubleAttr_MIPGap);

      for (int i = 0; i < N; i++){
        cout<<yp[i].get(GRB_DoubleAttr_X)<<"  ";
      } 
      cout<<endl;

      cout<<"yr"<<endl;

      for (int i = 0; i < N; i++) {
        cout<<yr[i].get(GRB_DoubleAttr_X)<<"  ";
      } 
      cout<<endl;

      for (int i = 0; i < N; i++) {
        yp_int[i] = yp[i].get(GRB_DoubleAttr_X);
        yr_int[i] = yr[i].get(GRB_DoubleAttr_X);
        
      }

      model.write("model.lp");
    
      if (optimstatus == GRB_OPTIMAL) {
        status = 1;
        cout << "Model is optimal" << endl;
      } else if (optimstatus == GRB_INF_OR_UNBD) {
        cout << "Model is infeasible or unbounded" << endl;
      } else if (optimstatus == GRB_INFEASIBLE) {
        cout << "Model is infeasible" << endl;
      } else if (optimstatus == GRB_UNBOUNDED) {
        cout << "Model is unbounded" << endl;
      } else {
        cout << "Optimization was stopped with status = " << optimstatus << endl;
      }
    
    }
    catch (GRBException e) {
      cout << "Error code = " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    }
    catch(...) {
      cout << "Exception during optimization" << endl;
    }
  
  return objval;
  
}

string IntToString (int a)
{
  ostringstream temp;
  temp<<a;
  return temp.str();
}
