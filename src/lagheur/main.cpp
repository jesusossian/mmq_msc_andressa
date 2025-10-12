#include "cabecalho.h"

int main (int argc, char** argv){

  double *HP, *HR, *PP, *FP, *PR, *FR, *D, *R;
  int i, j, t;
  FILE* fp;
  double objval;
  FILE* fpout;
  time_t start, end;
  int elapsed;
  double setup;

  /*mip or lp*/
  opc.solve_lp = 0;
  opc.solve_mip = 0;

  /*instances*/
  opc.instance_sifaleras = 0;
  opc.instance_mathijn = 0;

  N = 0;
  objval = 0.0;
  bestbound = 0.0;
  gap = 0.0;
  elapsed = 0.0;
  
  if (argc < 1 || argc > 10)
    {
      printf("numero de parametros de entrada errado \n");
      return 0;
    }
  
  fp = fopen(argv[1], "r");
  if (fp == NULL)
    {
      fprintf(stderr, "nao pode abrir %s \n", argv[1]);
      return 0;
    }
  
  for (i=2; i<argc; i++)
    {
      if (strcmp(argv[i], "-lp") == 0)
	{
	  opc.solve_lp = 1;
	}
      else if (strcmp(argv[i], "-mip") == 0)
	{
	  opc.solve_mip = 1;
	}
      else if (strcmp(argv[i], "-isifa") == 0)
	{
	  opc.instance_sifaleras = 1;
	}
      else if (strcmp(argv[i], "-imath") == 0)
	{
	  opc.instance_mathijn = 1;
	}
      else if (opc.instance_mathijn == 1 && N==0)
	{
	  N = atoi(argv[i]);
	}
      else if (opc.instance_mathijn == 1 && N!=0)
	{
	  setup = atoi(argv[i]);
	}
      else
	{
	  printf("parametros de entrada errado \n");
	}
    }
  
  if (opc.instance_sifaleras == 1)
    {
      if (fscanf(fp, "%d", &N) != 1)
	{
	  fprintf(stderr, "error : falta entrada  \n");
	  exit(1);
	}
    }


  yp_int = new double[N];
  yr_int = new double[N];

  for (int i =0; i < N; i++){
    yp_int[i] = -1;
    yr_int[i] = -1;
  }
  
  HP = calloc_vetor_double(N);
  PP = calloc_vetor_double(N);
  FP = calloc_vetor_double(N);
  HR = calloc_vetor_double(N);
  PR = calloc_vetor_double(N);
  FR = calloc_vetor_double(N);
  D = calloc_vetor_double(N);
  R = calloc_vetor_double(N);
  SD = calloc_matriz_double(N,N);
  SR = calloc_matriz_double(N,N);
  ssp = calloc_vetor_double(N);
  sxp = calloc_vetor_double(N);
  syp = calloc_vetor_double(N);
  ssr = calloc_vetor_double(N);
  sxr = calloc_vetor_double(N);
  syr = calloc_vetor_double(N);

  lerdados(fp, HP, PP, FP, HR, PR, FR, D, R, setup);

  fclose(fp);
  fp = NULL;
  /*
    printf("N = %d \n", N);
    printf("HP = [ "); for (i=0; i<N; i++) printf("%.2f ", HP[i]); printf("] \n");
    printf("PP = [ "); for (i=0; i<N; i++) printf("%.2f ", PP[i]); printf("] \n");
    printf("FP = [ "); for (i=0; i<N; i++) printf("%.2f ", FP[i]); printf("] \n");
    printf("HR = [ "); for (i=0; i<N; i++) printf("%.2f ", HR[i]); printf("] \n");
    printf("PR = [ "); for (i=0; i<N; i++) printf("%.2f ", PR[i]); printf("] \n");
    printf("D = [ ");  for (i=0; i<N; i++) printf("%.2f ", D[i]); printf("] \n");
    printf("R = [ "); for (i=0; i<N; i++) printf("%.2f ", R[i]); printf("] \n");
  */

  for (i=0; i<N; i++)
    {
      SD[i][i] = D[i];
      SR[i][i] = R[i];
      for (j=i+1; j<N; j++)
        {
	  SD[i][j] = SD[i][j-1] + D[j];
	  SR[i][j] = SR[i][j-1] + R[j];
        }
    }
 
  int tam_sub_conj = 10;

  int iteration = 1;
  int tam_sub_iter = 10;
  
  bool sair = false;
  while (sair == false){



    if ((iteration*tam_sub_conj)+  tam_sub_conj > N){
        tam_sub_iter = tam_sub_iter + (N - (iteration*tam_sub_conj));
        objval = relax_fix (HP,HR,PP,PR,FP,FR, D, R, tam_sub_conj, iteration,tam_sub_iter);
        sair = true;
    }
    else{
        objval = relax_fix(HP,HR,PP,PR,FP,FR, D, R, tam_sub_conj, iteration,tam_sub_iter);
       

        tam_sub_iter = tam_sub_iter + tam_sub_conj;
    }

    iteration = iteration + 1;
    
}

  objval = subgradiente(HP, HR, PP, PR,FP, FR, D, R);

cout<<"FO: "<<objval<<endl;

  /*
  printf("SOLUCAO : \n");
  printf("objval = %.2f \n", objval);
  printf("sp = [ "); for (i=0; i<N; i++) printf("%.2f ", ssp[i]); printf(" ] \n");
  printf("xp = [ "); for (i=0; i<N; i++) printf("%.2f ", sxp[i]); printf(" ] \n");
  printf("yp = [ "); for (i=0; i<N; i++) printf("%.2f ", syp[i]); printf(" ] \n");
  printf("sr = [ "); for (i=0; i<N; i++) printf("%.2f ", ssr[i]); printf(" ] \n");
  printf("xr = [ "); for (i=0; i<N; i++) printf("%.2f ", sxr[i]); printf(" ] \n");
  printf("sr = [ "); for (i=0; i<N; i++) printf("%.2f ", syr[i]); printf(" ] \n");
  */
  
  //fpout = fopen ("table.txt", "a");
  //fprintf (fpout, "%s;%.1f;%.1f;%1.f;%d;%d;%d \n", argv[1], objval,bestbound, gap,numnode, elapsed,status);


  free_vetor_double(HP);
  free_vetor_double(PP);
  free_vetor_double(FP);
  free_vetor_double(HR);
  free_vetor_double(PR);
  free_vetor_double(FR);
  free_vetor_double(D);
  free_vetor_double(R);
  free_matriz_double(N, SD);
  free_matriz_double(N, SR);
  free_vetor_double(ssp);
  free_vetor_double(sxp);
  free_vetor_double(syp);
  free_vetor_double(ssr);
  free_vetor_double(sxr);
  free_vetor_double(syr);

  delete [] yp_int;
  delete [] yr_int;
  
  return 0;

}
