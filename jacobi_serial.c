#include <stdio.h>

#define epsilon 0.0001

main(int argc, char* argv[]) 
{
     int n;
     float *Xi_old,*Xi_new;//Xi in each iteration
     short convergence_flag;
     
     int i,j;
     FILE *inputFile;
     
     inputFile = fopen ("inputFile.txt","r");
     fscanf(inputFile, "%d", &n);
     
	       
     short **a_b = (short**)malloc(n*sizeof(short*));//Coefficients a[i,j] followed by b[i]
     for(i=0; i<n; i++)
     {
	  a_b[i] = (short*)malloc((n+1)*sizeof(short));
	  for(j=0; j<n; j++)
	  {
	       fscanf(inputFile, "%d", &a_b[i][j]);//a
	  }
	  fscanf(inputFile, "%d", &a_b[i][n]);//b
     }
     fclose(inputFile);

     
     Xi_old = (float*)malloc(n*sizeof(float));
     Xi_new = (float*)malloc(n*sizeof(float));
     for(j=0; j<n; j++)
     {
	  Xi_old[j] = 0;//Initialization 
     }
     
     //Computations
     do
     {
	  convergence_flag = 0; 
	  for(i=0; i<n; i++)
	  {
	       float sum = 0;
	       for(j=0; j<n; j++)
	       {
		    sum += a_b[i][j]*Xi_old[j]; 
	       }
	       sum -= a_b[i][i]*Xi_old[i];
	       
	       Xi_new[i] = (1.0/a_b[i][i]) * (a_b[i][n]/* b[i] */-sum);
	       
	       if(fabs(Xi_old[i]-Xi_new[i]) <= epsilon)
		    convergence_flag++;
	  }
	  
	  for(i=0; i<n; i++)
	  {
	       Xi_old[i] = Xi_new[i];
	  }
     }while(convergence_flag<n);
     
     
     
     printf("\n");
     for (i = 1; i <n+1; i++) 
     {
	  printf(" X%d = %f\t", i, Xi_new[i-1]);    
     }
     printf("\n\n");

     free(a_b);
     free(Xi_old);
     free(Xi_new);
}
