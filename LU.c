/* Name: Shahzeb Siddiqui
 * Class: CS 311 High Performance Computing 1
 * Project: LU Decomposition
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
 
// function definition
void printMatrix(float **array);
void ReadMatrix();
float **allocate_array(int row, int col);

int N;													// size of square matrix
float **A;												// input matrix

int main()
{
  int i,j,k,l;												// iteration counters
  float **Anew;												// array to check results after LU decomposition
  float **LT,**UT;											// LU decomposition arrays Lower Triangle, Upper Triangle
  
  // read matrix from file and store in matrix A
  ReadMatrix();
  
  printf("   InputMatrix\n");
  printMatrix(A);

  // allocating matrices NxN for LT,UT,Anew for computing below
  LT = allocate_array(N,N);
  UT = allocate_array(N,N);
  Anew = allocate_array(N,N); 

  // Doolittle Algorithm
  for (i = 0; i < N; i++)
  {
	for (j = i; j < i+1; j++)
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
  // Verification: Anew = A by doing Anew = LT * UT 
  for (i = 0; i < N; i++)
  {
	for (j = 0; j < N; j++)
	{
		for (k = 0; k < N; k++)
			Anew[i][j] += LT[i][k] * UT[k][j];
	}
  }
  printf("LU decomposition Results. Anew = LT * UT\n");
  printf("\n  LT:\n"); printMatrix(LT);
  printf("\n  UT:\n"); printMatrix(UT);
  printf("\n  Anew:\n"); printMatrix(Anew);  

return 0;
}
float **allocate_array(int row_dim, int col_dim)
{
      float **result;											// matrix
      int i,j;												// iteration counters

      // Allocate an array of pointers to hold pointers to the rows of the array 
     result=(float **)malloc(row_dim*sizeof(float *));

     // The first pointer is in fact a pointer to the entire array 
     result[0]=(float *)malloc(row_dim*col_dim*sizeof(float));

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

void printMatrix(float ** array)
{
 
      int i,j;												// iteration counters
      for (i = 0; i < N; i++)
      {
            for (j = 0; j < N; j++)
                    printf("%0.4f ",array[i][j]);
            printf("\n");
      }
}
void ReadMatrix()
{
	int size;											// variable to store 1st line in file that is size of matrix
	int row = 0, col = 0;										// variable used for indexing 2D array while reading matrix 
	FILE * fp;											// file pointer used for reading file
	const char* filename = "matrixfile";								// filename of text file
	fp = fopen(filename,"r");
	if (fp == NULL)
		perror("Error opening file\n");
   
	fscanf(fp,"%d",&size);
 	A = allocate_array(size,size);
	N = size;
	// read until end of file
	while (feof(fp) == 0)
	{
		fscanf(fp,"%f",&A[row][col]);
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

