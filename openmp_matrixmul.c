#include <stdio.h>
#include "timer.h"
#include <omp.h>
#define SIZE 2048 

float a[SIZE][SIZE];
float b[SIZE][SIZE];
float c[SIZE][SIZE];

int main()
{
  StartTimer();
  
  int i,j,k;
  omp_set_num_threads(2);
  // Initialize matrices.
   #pragma omp parallel for
   for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      a[i][j] = (float)i + j;
      b[i][j] = (float)i - j;
      c[i][j] = 0.0f;
    }
  }
  // Compute matrix multiplication.
  #pragma omp parallel
  {
  #pragma omp for
  for (i = 0; i < SIZE; ++i) {
    for (j = 0; j < SIZE; ++j) {
      for (k = 0; k < SIZE; ++k) {
	c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  }

  double runtime = GetTimer();
  #pragma omp parallel
  {
    if(omp_get_thread_num() == 0)
	printf("execution time: %f s with %d threads\n", runtime / 1000.f,omp_get_num_threads());
  }
  /*
  for (i = 0; i < SIZE; i++)
  {
    for (j = 0; j < SIZE; j++)
	printf("%0.6f\t",c[i][j]);
	printf("\n");
  }*/
  return 0;
}
