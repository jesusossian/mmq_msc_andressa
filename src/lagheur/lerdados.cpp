void lerdados (FILE *fp, double *HP, double *PP, double *FP, double *HR, double *PR, double *FR, double *D, double *R, double setup) {
  int i, t;
  char mystring[111][200];
  
  if (opc.instance_sifaleras == 1)
    {
      if (fscanf(fp, "%lf", &(FR[0])) != 1)
	{
	  fprintf(stderr, "error : falta entrada \n");
	  exit(1);
	}
      for (i=1; i<N; i++) FR[i] = FR[0];
      
      if (fscanf(fp, "%lf", &(FP[0])) != 1)
	{
	  fprintf(stderr, "error : falta entrada \n");
	  exit(1);
	}
      for (i=1; i<N; i++) FP[i] = FP[0];
      
      if (fscanf(fp, "%lf", &(HR[0])) != 1)
	{
	  fprintf(stderr, "error : falta entrada  \n");
	  exit(1);
	}
      for (i=1; i<N; i++) HR[i] = HR[0];
      
      if (fscanf(fp, "%lf", &(HP[0])) != 1)
	{
	  fprintf(stderr, "error : falta entrada \n");
	  exit(1);
	}
      for (i=1; i<N; i++) HP[i] = HP[0];
      
      for (i=0; i<N; i++)
	{
	  if (fscanf(fp, "%lf", &(D[i])) != 1)
	    {
	      fprintf(stderr, "error : falta entrada \n");
	      exit(1);
	    }
	}
      
      for (i=0; i<N; i++)
	{
	  if (fscanf(fp, "%lf", &(R[i])) != 1)
	    {
	      fprintf(stderr, "error : falta entrada \n");
	      exit(1);
	    }
	}
      
    }
   
  if (opc.instance_mathijn == 1)
    {
      
      for (i=0; i<11; i++)
	{
	  if(fgets(mystring[i], 200, fp)==NULL)
	    {
	      fprintf(stderr, "error : entrada 0 \n");
	      exit(1);
	    }
	}
      
      for (i=0; i<N; i++)
	{
	  if (fscanf(fp, "%d", &(t)) != 1)
	    {
	      fprintf(stderr, "error : entrada 1 \n");
	      exit(1);
	    }
	  t -= 1;
	  
	  if (fscanf(fp, "%lf", &(D[t])) != 1)
	    {
	      fprintf(stderr, "error : entrada 2 \n");
	      exit(1);
	    }
	  
	  if (fscanf(fp, "%lf", &(R[t])) != 1)
	    {
	      fprintf(stderr, "error : entrada 3 \n");
	      exit(1);
	    }
	  
	  if(fgets(mystring[i+10], 200, fp)==NULL)
	    {
	      fprintf(stderr, "error : entrada 4\n");
	      exit(1);
	    }
	}
      
      for (i=0; i<N; i++) FP[i]=setup;
      for (i=0; i<N; i++) FR[i]=setup;
      for (i=0; i<N; i++) HP[i]=1.0;
      for (i=0; i<N; i++) HR[i]=1.0;
      for (i=0; i<N; i++) PP[i]=0.0;
      for (i=0; i<N; i++) PR[i]=0.0;
      
    }  

}
