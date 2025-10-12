double** malloc_matriz_double (const int m, const int n)
{
  double **v;
  int i;
  
  v = (double **) malloc (m * sizeof(double));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  for (i = 0; i < m; i++)
    {
      v[i] = (double *) malloc ((n) * sizeof(double));
      if (v[i] == NULL)
	{
	  printf ("error: memoria insuficiente \n");
	  exit (EXIT_FAILURE);
	}
    }
  
  return v;
  
}

int** malloc_matriz_int (const int m, const int n)
{
  int **v;
  int i;
  
  v = (int **) malloc (m * sizeof(int *));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  for (i = 0; i < m; i++)
    {
      v[i] = (int*) malloc ((n) * sizeof(int));
      if (v[i] == NULL)
	{
	  printf ("error: memoria insuficiente \n");
	  exit (EXIT_FAILURE);
	}
    }
  
  return v;
  
}

bool* malloc_vetor_bool (const int n)
{
  bool *v;
  
  v = (bool *) malloc(n * sizeof(bool));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  return v;
  
}

double* malloc_vetor_double (const int n)
{
  double *v;
  
  v = (double *) malloc(n * sizeof(double));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  return v;
  
}

int* malloc_vetor_int (const int n)
{
  int *v;
  
  v = (int *) malloc(n * sizeof(int));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  return v;
  
}

int* calloc_vetor_int (const int n)
{
  int *v;
  
  v = (int *) calloc(n, sizeof(int));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  return v;
  
}

double* calloc_vetor_double (const int n)
{
  double *v;
  
  v = (double *) calloc(n, sizeof(double));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  return v;
  
}

int** calloc_matriz_int (const int m, const int n)
{
  int **v;
  int i;
  
  v = (int **) calloc (m, sizeof(int *));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  for (i = 0; i < m; i++)
    {
      v[i] = (int*) calloc (n, sizeof(int));
      if (v[i] == NULL)
	{
	  printf ("error: memoria insuficiente \n");
	  exit (EXIT_FAILURE);
	}
    }
  
  return v;
  
}

double** calloc_matriz_double (const int m, const int n) {
  double **v;
  int i;
  
  v = (double **) calloc (m, sizeof(double *));
  if (v == NULL)
    {
      printf ("erro: memoria insuficiente \n");
      exit (EXIT_FAILURE);
    }
  
  for (i = 0; i < m; i++)
    {
      v[i] = (double*) calloc (n, sizeof(double));
      if (v[i] == NULL)
	{
	  printf ("error: memoria insuficiente \n");
	  exit (EXIT_FAILURE);
	}
    }
  
  return v;
  
}

void free_vetor_int (int *vetor)
{
  if (vetor)
    {
      free(vetor);
      vetor = NULL;
    }
}

void free_vetor_double (double *vetor)
{
  if (vetor)
    {
      free(vetor);
      vetor = NULL;
    }
}

void free_matriz_int (int n, int **matriz)
{
  int i;
  
  for (i = 0; i < n; i++)
    {
      if (matriz[i])
	{
	  free(matriz[i]);
	  matriz[i] = NULL;
	}
    }
  
  if (matriz)
    {
      free(matriz);
      matriz = NULL;
    }
  
}

void free_matriz_double (int n, double **matriz)
{
  int i;
  
  for (i = 0; i < n; i++)
    {
      if (matriz[i])
	{
	  free(matriz[i]);
	  matriz[i] = NULL;
	}
    }
  
  if (matriz)
    {
      free(matriz);
      matriz = NULL;
    }
  
}
