/* Name: Shahzeb Siddiqui
 * Class: CS 311 High Performance Computing 1
 * Project: PLU Decomposition with PAPI
 * KAUST ID: 124211
 * Description: Apply the LU decomposition method using Doolittle algorithm on a square input matrix
 * 		which decomposes into lower triangle and upper triangle. These two matrices should decompose
 * 		in manner such that LT * UT = A 
 * 		where: 
 * 			LT: Lower Triangle Matrix 
 * 			UT: Upper Triangle Matrix
 *			 A: Input Matrix 
 */
// library definition
#include <stdio.h>
#include <stdlib.h>
#include <papi.h> 
#define tolerance 0.0001

// function definition
void printMatrix(double **array);
void filewriteMatrix(double **P,double **L,double **U);
void ReadMatrix();
double **allocate_array(int row, int col);
int validateLU(double **A, double **Aorig, int N);

int N;													// size of square matrix
double **A;												// input matrix

int main()
{
  int i,j,k,l;												// iteration counters
  double **Aorig;
  double **Anew1,**Anew2;										// array to check results after PLU decomposition
  double **P,**LT,**UT;											// LU decomposition arrays Lower Triangle, Upper Triangle
 
  const PAPI_hw_info_t *hwinfo; 
 
  int retval;
  
  // read matrix from file and store in matrix A
  ReadMatrix();
  
  // PAPI Events
  //int Events[3] = {PAPI_TOT_CYC, PAPI_TOT_INS,PAPI_FP_OPS};			// Total Cycles, Total Instructions, Floating Point Operations
  //int Events[3] = {PAPI_L1_DCA,PAPI_L1_DCH,PAPI_L1_DCM};				// L1 Data Cache Access, L1 Data Cache Hits, L1 Data Cache Miss
  //int Events[3] = {PAPI_L1_ICA,PAPI_L1_ICH,PAPI_L1_ICM};				// L1 Instruction Cache Access, L1 Instruction Cache Hit, L1 Instruction Cache Miss
  //int Events[3] = {PAPI_L2_TCM,PAPI_L2_TCH,PAPI_L2_TCA};				// L2 Total Cache Miss, L2 Total Cache Hit, L2 Total Cache Access
  //int Events[3] = {PAPI_L2_ICM,PAPI_L2_ICH,PAPI_L2_ICA};				// L2 Instruction Cache Miss, L2 Instruction Cache Hit, L2 Instruction Cache Access
  //int Events[2] = {PAPI_BR_NTK,PAPI_BR_TKN};  						// Branch Not-Taken, Branch Taken
  int Events[1] = {PAPI_RES_STL};										// Stall on any Resource
  long_long values[3];


  // allocating matrices NxN for P,LT,UT,Anew for computing below
  LT = allocate_array(N,N);
  UT = allocate_array(N,N);
  P = allocate_array(N,N);
  Anew1 = allocate_array(N,N); 
  Anew2 = allocate_array(N,N);
  Aorig = allocate_array(N,N);

  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
		Aorig[i][j] = A[i][j];
  }
  if (PAPI_start_counters(Events,1) != PAPI_OK)
  {
	printf("unable to start PAPI counters\n");
	exit(1);
  }
  for (i = 0; i < N; i++)
    P[i][i] = 1.0;

  // PLU Algorithm
  for (i = 0; i < N; i++)
  {
	for (j = i; j < i+1; j++)
   	{
		
		double maxNumInCol = A[j][j];			
		if (A[j][j] == 0.f)
		{
			int rowindex;
			double tempA, tempP;
			// loop column where 0 exist in diagonal and swap with max value in that column by searching in lower column of diagonal
			for (k = j+1; k < N; k++)
			{
				if (A[k][i] > maxNumInCol)
				{
					maxNumInCol = A[k][i]; 
					rowindex = k;
				}
			}
			if (maxNumInCol != 0.f)
			{
			// swap rows i and rowindex for matrix A and permutation matrix P
			for (l = 0; l < N; l++)
			{
				tempA = A[i][l];
				A[i][l] = A[rowindex][l];
				A[rowindex][l] = tempA;
				
				tempP = P[i][l];
				P[i][l] = P[rowindex][l];
				P[rowindex][l] = tempP;
			}
			}
			
		}
		if (maxNumInCol != 0.f)
		{
		// loop u
		for (k = i+1; k < N; k++)
		{
			A[k][i] = A[k][i] / A[j][j];
		}
		//calculate A(N-1) - uv'/a
		for (k = i+1; k < N; k++)
		{
			for (l = i+1; l < N; l++)
			{
				A[k][l] = A[k][l] - A[k][i]*A[i][l];
			}
		}
		}

	}
  }
  // assign arrays UT and LT values in there region based on A. LT contains the diagonal 1s
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
	{
		if (i > j)
		{
			LT[i][j] = A[i][j];
		}
		else if(i == j)
		{
			LT[i][j] = 1.0;
			UT[i][j] = A[i][j];
		}
		else
		{
			UT[i][j] = A[i][j];
		}
	}
  }

     if (PAPI_read_counters(values,1) != PAPI_OK)
     {
	printf("unable to read PAPI counters");
	exit(1);
     }


  FILE *fp = fopen ("test1","aw");
  
  fprintf(fp,"Task 2: Hardware Performance Counters\n");              
  fprintf(fp,"--------------------------------------\n");
  fprintf(fp,"PLU Statistics\n");
/*  
  fprintf(fp,"Total Cycles 			: %ld\n",(long)values[0]);
  fprintf(fp,"Total Instructions		: %ld\n",(long)values[1]);
  fprintf(fp,"Total Floating point operations 	: %ld\n",(long)values[2]);
*/

//  fprintf(fp,"\n\nL1 Cache Statistics\n");
/*
  fprintf(fp,"L1 data cache accesses		: %ld\n",(long)values[0]);
  fprintf(fp,"L1 data cache hits		: %ld\n",(long)values[1]);
  fprintf(fp,"L1 data cache misses		: %ld\n",(long)values[2]);
*/
/*
  fprintf(fp,"L1 instruction cache accesses             : %ld\n",(long)values[0]);
  fprintf(fp,"L1 instruction cache hits                 : %ld\n",(long)values[1]);
  fprintf(fp,"L1 instruction cache misses               : %ld\n",(long)values[2]);
*/

 
//  fprintf(fp,"\n\nL2 Cache Statistics\n");
/*  
  fprintf(fp,"L2 total cache miss             : %ld\n",(long)values[0]);
  fprintf(fp,"L2 total cache hits             : %ld\n",(long)values[1]);
  fprintf(fp,"L2 total cache access		: %ld\n",(long)values[2]);
*/
/*  
  fprintf(fp,"Level 2 instruction cache misses 	: %ld\n",(long)values[0]);
  fprintf(fp,"Level 2 instruction cache hits		: %ld\n",(long)values[1]);
  fprintf(fp,"Level 2 instruction cache accesses	: %ld\n",(long)values[2]);
*/  
/*
  fprintf(fp,"\n\nBranch Instruction Statistics\n");
  fprintf(fp,"Branch Instructions Not Taken		: %ld\n",(long)values[0]);
  fprintf(fp,"Branch Instructions Taken			: %ld\n",(long)values[1]);
*/


//  fprintf(fp,"Cycles stalled on any resource			: %ld\n",(long)values[0]);

/*  if ((retval=PAPI_flops(&irealtime2,&iproctime2,&iflpops2,&imflops2)) < PAPI_OK)
  {
        printf("Error acquring floating point events\n");
        exit(1);
  }
 */
  
  // Verification: Anew2 = A by doing Anew2 = P * LT * UT 
  //  Anew = P *LT
  for (i = 0; i < N; i++)
  {
  	for (j = 0; j < N; j++)
        {
        	for (k = 0; k < N; k++)
                	Anew1[i][j] += P[i][k] * LT[k][j];
        }
  }
  // Anew2 = P * LT * UT 
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
	{
		for (k = 0; k < N; k++)
			Anew2[i][j] += Anew1[i][k] * UT[k][j];
	}
  }
  if (PAPI_stop_counters(values,1) != PAPI_OK)
  {
	exit(1);
  }
   fprintf(fp,"Matrix Multiplication Statistics\n");
  
/* fprintf(fp,"Total Cycles                        : %ld\n",(long)values[0]);
   fprintf(fp,"Total Instructions                  : %ld\n",(long)values[1]);
   fprintf(fp,"Total Floating point operations     : %ld\n",(long)values[2]);
*/
//  fprintf(fp,"\n\nL1 Cache Statistics\n");
/*

  fprintf(fp,"L1 data cache accesses            : %ld\n",(long)values[0]);
  fprintf(fp,"L1 data cache hits                : %ld\n",(long)values[1]);
  fprintf(fp,"L1 data cache misses              : %ld\n",(long)values[2]);      
*/
/*
  fprintf(fp,"L1 instruction cache accesses             : %ld\n",(long)values[0]);
  fprintf(fp,"L1 instruction cache hits 	        : %ld\n",(long)values[1]);
  fprintf(fp,"L1 instruction cache misses               : %ld\n",(long)values[2]);
*/


//  fprintf(fp,"\n\nL2 Cache Statistics\n");
/*  
  fprintf(fp,"L2 total cache miss           : %ld\n",(long)values[0]);
  fprintf(fp,"L2 total cache hits             : %ld\n",(long)values[1]);
  fprintf(fp,"L2 total cache access             : %ld\n",(long)values[2]);
*/
/*
  fprintf(fp,"Level 2 instruction cache misses          : %ld\n",(long)values[0]);
  fprintf(fp,"Level 2 instruction cache hits            : %ld\n",(long)values[1]);
  fprintf(fp,"Level 2 instruction cache accesses        : %ld\n",(long)values[2]);
*/
/*
  fprintf(fp,"\n\nBranch Instruction Statistics\n");
  fprintf(fp,"Branch Instructions Not Taken		  : %ld\n",(long)values[0]);
  fprintf(fp,"Branch Instructions Taken 	   	  : %ld\n",(long)values[1]);
*/

//  fprintf(fp,"Cycles stalled on any resource                     : %ld\n",(long)values[0]);


  fclose(fp);
//  printf("LU decomposition Results. Anew2 = P * LT * UT\n");
//  printf("\n  Anew2:\n"); printMatrix(Anew2); 

  retval = validateLU(Anew2,Aorig,N);
  if (retval)
    printf("\nProblem with matrix, LU decomposition not within tolerance level of %f\n",tolerance);
  
  filewriteMatrix(P,LT,UT);
  
  FILE *fp1 = fopen("verificationMatrix","w");
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
		fprintf(fp1,"%.2lf\t",Anew2[i][j]);
	fprintf(fp1,"\n");
  }
  fclose(fp1);
  
  // PAPI System Info
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
	exit(1);
  if ((hwinfo = PAPI_get_hardware_info())==NULL)
	exit(1);

  printf("\nTASK 3: System Information via PAPI:\n");
  printf("-----------------------------------\n");
  printf("SMP CPU count: %d\n",hwinfo->ncpu);
  printf("Node count: %d\n",hwinfo->nnodes);
  printf("Total-CPUs: %d\n",hwinfo->totalcpus);
  printf("Vendor ID: %d\nVendor: %s\n",hwinfo->vendor,hwinfo->vendor_string);
  printf("Model Number: %d\nModelString: %s\n",hwinfo->model,hwinfo->model_string);
  printf("CPU-frequency: %f Mhz.\n",hwinfo->mhz);
return 0;
}
double **allocate_array(int row_dim, int col_dim)
{
      double **result;											// matrix
      int i,j;												// iteration counters

      // Allocate an array of pointers to hold pointers to the rows of the array 
     result=(double **)malloc(row_dim*sizeof(double *));

     // The first pointer is in fact a pointer to the entire array 
     result[0]=(double *)malloc(row_dim*col_dim*sizeof(double));

     // The remaining pointers are just pointers into this array, offset in units of col_dim
     for(i=1; i<row_dim; i++)
     {
           result[i]=result[i-1]+col_dim;
     }

     // initialize array to 0s
     for(i=0;i<row_dim;i++)
	for(j=0;j<col_dim;j++)
		result[i][j] = 0.f;

     return result;
 }

void printMatrix(double ** array)
{
 
      int i,j;												// iteration counters
      for (i = 0; i < N; i++)
      {
            for (j = 0; j < N; j++)
                    printf("%0.4f ",array[i][j]);
            printf("\n");
      }
}
void filewriteMatrix (double **P, double **L, double **U)
{
  char *filename = "outputfile.txt";
  FILE *fp = fopen(filename,"w");
  int i,j;
  // writing matrix P
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
	fprintf(fp,"%.2lf ",P[i][j]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  // writing matrix L
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
		fprintf(fp,"%.2lf ",L[i][j]);
  	fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  // writing matrix U
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
		fprintf(fp,"%.2lf ",U[i][j]);
	fprintf(fp,"\n");
  }
  fclose(fp);
}
void ReadMatrix()
{
	int size;											// variable to store 1st line in file that is size of matrix
	int row = 0, col = 0;										// variable used for indexing 2D array while reading matrix 
	FILE * fp;											// file pointer used for reading file
	const char* filename = "sample_LU_512.in";								// filename of text file
	fp = fopen(filename,"r");
	if (fp == NULL)
		perror("Error opening file\n");
   
	fscanf(fp,"%d",&size);
 	A = allocate_array(size,size);
	N = size;
	// read until end of file
	while (feof(fp) == 0)
	{
		fscanf(fp,"%lf",&A[row][col]);
		// next row set column to 0
		if (col == size-1)
		{
			col = 0; 
			row++;
		}
		else
			col++;
	}
 	// close file
  	fclose(fp);
}
int validateLU (double **Anew, double **Aorig, int N)
{
  int i,j;
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
	{
		if (fabsf(Anew[i][j] - Aorig[i][j]) >= tolerance)
		{
			printf("tolerance exceed for Anew[%d][%d] = %f, \t Aorig[%d][%d] = %f\n",i,j,Anew[i][j],i,j,Aorig[i][j]);
 			return 1;
		}
	}
  }
  return 0;
}		
 

