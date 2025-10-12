double calc_soma(double *v){
  double soma = 0;
  for(int i = 0 ; i < N; i++){
    soma = soma + pow(v[i],2); 
  }

  return soma;
}

double maximo(double a, double b){
  if(a > b){
    return a;
  }
  else{
    return b;
  }
}

double subgradiente(double *HP, double *HR, double *PP, double *PR,double *FP, double *FR, double *D, double *R){
  double C;
  double soma;
  double fator;
  double objval_relax_fix;
  double objval_lagrange;
  double elapsed_;
  int optimstatus;
  bool cap;
  // int tam_sub_inter;
  int status;
  int contador = 0;
  
  double z_max = -100000000000;

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
  
  //Cota superior

  int tam_sub_conj = 20;
  int iteration = 1;
  int tam_sub_iter = 20;
  
  bool sair = false;
  while (sair == false){

    if ((iteration*tam_sub_conj)+  tam_sub_conj > N){
        tam_sub_iter = tam_sub_iter + (N - (iteration*tam_sub_conj));
        objval_relax_fix = relax_fix (HP,HR,PP,PR,FP,FR, D, R, tam_sub_conj, iteration,tam_sub_iter);
        sair = true;
    }
    else{
        objval_relax_fix = relax_fix (HP,HR,PP,PR,FP,FR, D, R, tam_sub_conj, iteration,tam_sub_iter);
       

        tam_sub_iter = tam_sub_iter + tam_sub_conj;
    }

    iteration = iteration + 1;
    
  }

  //relaxação lagrangeana

  u = new double[N];
  sgv = new double[N];

  for(int i =0; i< N ; i++){
    u[i] = 0;
    sgv[i] = 0;
  }
   
  double alfa = 1.95;
  double step=0;

  bool sair_busca = false;
  int iteracao_sub = 1;
  while(sair_busca == false){
    cout<<"iteracao :"<<iteracao_sub <<endl;

    objval_lagrange = lagrangian(HP,HR,PP,PR,FP,FR, D, R);
    cout<<"calcula passo"<<endl;
    step = alfa*((objval_relax_fix - objval_lagrange)/calc_soma(sgv));
    cout<<"passo: "<<step<<endl;
    for(int i=0;i < N;i++){
      u[i] = maximo(0,u[i]+step*sgv[i]);
    }
    for(int i=0;i < N;i++){
      cout<<u[i]<< " ";
    }
    cout<<endl;
    int maior_zero = 0;
    for(int i=0; i < N; i++){
      if (sgv[i] > 0){
        maior_zero =1;
        break;
      }
    }

    if(maior_zero==1){
      if(objval_lagrange > z_max){
        z_max = objval_lagrange;
      }
    }
    else{

      contador = contador + 1;

    }
    if(iteracao_sub == 4){
      sair_busca=true;
    }
    if (contador == 30){
      alfa = alfa/2;
    }

    if((abs(objval_lagrange - objval_relax_fix) < 0.001) || (alfa < 0.005)){

      cout<<"saiu"<<endl;
      sair_busca = true;
    }
    iteracao_sub= iteracao_sub+1;

  }

  return z_max;
  
}
