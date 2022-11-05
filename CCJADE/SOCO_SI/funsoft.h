#ifndef _FUNSOFTCOMPUTING_H

#define _FUNSOFTCOMPUTING_H 1

#ifndef _TFITNESS 
#define _TFITNESS 1
typedef long double tFitness;
#endif 

#define MAX_D 1000

extern tFitness Shifted_Sphere( int dim ,  double* x );
extern tFitness Schwefel_Problem(int dim, double* x);
extern tFitness Shifted_Rosenbrock(int dim, double* x);
extern tFitness Shifted_Rastrigin(int dim, double * x);
extern tFitness Shifted_Griewank(int dim, double * x);
extern tFitness Shifted_Ackley(int dim, double* x);
extern tFitness f_Schwefel2_22(int dim, double *s);
extern tFitness f_Schwefel1_2(int dim, double *s);
extern tFitness Extended_f_10(int dim, double *x);
extern tFitness f_Bohachevsky(int dim, double *s);
extern tFitness f_Schaffer(int dim, double *s);
extern tFitness f_Hybrid_12(int dim, double *s);
extern tFitness f_Hybrid_13(int dim, double *s);
extern tFitness f_Hybrid_14(int dim, double *s);
extern tFitness f_Hybrid_15(int dim, double *s);
extern tFitness f_Hybrid_16new(int dim, double *s);
extern tFitness f_Hybrid_17new(int dim, double *s);
extern tFitness f_Hybrid_18new(int dim, double *s);
extern tFitness f_Hybrid_19new(int dim, double *s);

#endif
