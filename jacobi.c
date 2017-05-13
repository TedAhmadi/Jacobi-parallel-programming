#include <stdio.h>
#include "mpi.h"

#define epsilon 0.0001

main(int argc, char* argv[]) 
{
     int n, my_rank, p;
     float *Xi_old,*Xi_new;//Xi in each iteration
     short *a_b;
     int tag;
     short convergence_flag;
     MPI_Comm slave_comm;
     MPI_Group world_group, slave_group;
     MPI_Status status;
     
     MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
     MPI_Comm_size(MPI_COMM_WORLD, &p);
     
	  
	  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	  int master_ranks[] = {0};
	  MPI_Group_excl(world_group, 1, master_ranks, &slave_group);
	  MPI_Comm_create(MPI_COMM_WORLD, slave_group, &slave_comm);

	  //Checking whether there are correct number of processes needed.
	  FILE * inputFile = fopen ("inputFile.txt","r");
	  fscanf(inputFile, "%d", &n);
	  if((n+1)!=p)
	  {
	       printf("There must be %d processes (= problemDimention + 1). Terminating process %d!\n",n+1,my_rank);
	       fclose(inputFile);
	       MPI_Finalize();
	       exit(0);
	  }
	  fclose(inputFile);
	  
	  
	  if(my_rank == 0)
	  {
	       int i,j;
	       FILE *inputFile;
	       
	       inputFile = fopen ("inputFile.txt","r");
	       fscanf(inputFile, "%d", &n);
	       
	              
	       short **a_b_complete = (short**)malloc(n*sizeof(short*));//Coefficients a[i,j] followed by b[i]
	       for(i=0; i<n; i++)
	       {
		    a_b_complete[i] = (short*)malloc((n+1)*sizeof(short));
		    for(j=0; j<n; j++)
		    {
			 fscanf(inputFile, "%d", &a_b_complete[i][j]);//a
		    }
		    fscanf(inputFile, "%d", &a_b_complete[i][n]);//b
	       }
	       fclose(inputFile);
	                     
	       //Distributing coefficients
	       //MPI_Scatter (a_b_complete[0][0], (n+1), coefficients, a_b, (n+1), coefficients, 0, MPI_COMM_WORLD);
	       for(i=1; i<p; i++)
	       {
		    tag = 1;
		    MPI_Send(&n, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
		    
		    tag = 2;
		    MPI_Send(&a_b_complete[i-1][0], n+1, MPI_SHORT, i, tag, MPI_COMM_WORLD);
	       }
	       
	       //Collecting final Xi for all i
	       float *Xi;
	       Xi = (float*)malloc(n*sizeof(float));
	       tag = 4;
	       for (i = 1; i <p; i++) 
	       {
		    MPI_Recv(&Xi[i-1], 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &status);
	       }printf("\n");
	       for (i = 1; i <p; i++) 
	       {
		    printf(" X%d = %f\t", i, Xi[i-1]);		    
	       }
	       printf("\n\n");
	       
	       
	       
	       for(i=0; i<n; i++)
	       {
		    free(a_b_complete[i]);
	       }
	       free(a_b_complete);
	       free(Xi);
	  }
	  
	  else
	  {
	       int i;
	       
	       tag = 1;
	       MPI_Recv(&n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
	       a_b = (short*)malloc((n+1)*sizeof(short));
	              
	       tag = 2;
	       
	       MPI_Recv(a_b, n+1, MPI_SHORT, 0, tag, MPI_COMM_WORLD, &status);
	       
	       Xi_old = (float*)malloc(n*sizeof(float));
	       Xi_new = (float*)malloc(n*sizeof(float));
	       for(i=0; i<n; i++)
	       {
		    Xi_old[i] = 0;//Initialization 
	       }
	       
	       //Computations
	       i = my_rank-1;
	       do
	       {
		    float sum = 0;
		    int j;
		    for(j=0; j<n; j++)
		    {
			 sum += a_b[j]*Xi_old[j]; 
		    }
		    sum -= a_b[i]*Xi_old[i];
		    
		    Xi_new[i] = (1.0/a_b[i]) * (a_b[n]/* b */-sum);
		    
		    if(fabs(Xi_old[i]-Xi_new[i]) <= epsilon)
			 convergence_flag = 1;
		    else
			 convergence_flag = 0; 
		    
		    Xi_old[i] = Xi_new[i];
		    
		    
		    MPI_Barrier(slave_comm);
		    
		    
		    MPI_Allreduce (&convergence_flag, &convergence_flag, 1, MPI_SHORT, MPI_MIN, slave_comm);//Convergence check
		
		    tag = 3;
		    for(j=1; j<p; j++)
		    {
			 if(j != my_rank)
			 {
			      MPI_Send(&Xi_old[i], 1, MPI_FLOAT, j, tag, MPI_COMM_WORLD);
			      MPI_Recv(&Xi_old[j-1], 1, MPI_FLOAT, j, tag, MPI_COMM_WORLD, &status);
			 }
		    }      
		    //printf("%d ---> %f ,   %f %f %f\n", my_rank, Xi_new[i],Xi_old[0],Xi_old[1],Xi_old[2]);
	       }while(convergence_flag==0);
	       tag = 4;
	       MPI_Send(&Xi_new[i], 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
	       
	       free(a_b);
	       free(Xi_old);
	       free(Xi_new);
	  }
     
     MPI_Finalize();
}
