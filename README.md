# Jacobi-parallel-programming
"Jacobi iteration method for solving a system of linear equation"

 
compile:  mpicc jacobi.c -o jacobi

run:      mpirun -np P jacobi   /*P is the number of processes*/


"inputFile.txt":
n
a11 a12 . . . a1n  b1
a21 a22 . . . a2n  b2
.		.   .
.		.   .
.		.   .
an1 an2 . . . ann  bn


The number of processes "P" must be equal the problem dimension "n"
