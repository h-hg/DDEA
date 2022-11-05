#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include "funsoft.h"
#include "cec08data.h"
#include "f7f11data.h"
#include <assert.h>
#define abss(a)     (a<0 ? (-a) : a)
#
static double pi = M_PI;
static double e = M_E;


tFitness Shifted_Sphere( int dim , double* x ){
    int i;
    tFitness z;
    tFitness F = 0;
    for(i=0;i<dim;i++){
        z = x[i] - sphere[i];
        F += z*z;
    }
    return F + f_bias[0];
}

tFitness Schwefel_Problem( int dim , double * x ){
    int i;
    tFitness z;
    tFitness F = abss(x[0]);
    for(i=1;i<dim;i++){
    	  z = x[i] - schwefel[i];
        F = (F > abss(z) ? F : abss(z));
    }
    return F + f_bias[1]; 
}

tFitness Shifted_Rosenbrock( int dim , double* x ){
    int i;
	double z[MAX_D];
    tFitness F = 0;

    for(i=0;i<dim;i++) z[i] = x[i] - rosenbrock[i] + 1;   

    for(i=0;i<dim-1;i++){    
        F = F + 100*( pow((pow(z[i],2)-z[i+1]) , 2) ) + pow((z[i]-1) , 2);
    }
    return F + f_bias[2]; 
}

tFitness Shifted_Rastrigin( int dim , double * x )
{
    int i;
    double z;
    tFitness F = 0;
    for(i=0;i<dim;i++){  
        z = x[i] - rastrigin[i];
        F = F + ( pow(z,2) - 10*cos(2*pi*z) + 10);
    }
    return F + f_bias[3]; 
}

tFitness Shifted_Griewank( int dim , double * x ){
    int i;
    tFitness z;
    tFitness F1 = 0;
    tFitness F2 = 1;
    for(i=0;i<dim;i++){       
        z = x[i] - griewank[i];
        F1 = F1 + ( pow(z,2) / 4000 );
        F2 = F2 * ( cos(z/sqrt(i+1)));

    }
    return (F1 - F2 + 1 + f_bias[4]); 
}

tFitness Shifted_Ackley( int dim , double* x ){
    int i;
    double z;
    double Sum1 = 0;
    double Sum2 = 0;
    double F = 0;
    for(i=0;i<dim;i++){   
        z = x[i] - ackley[i];
        Sum1 = Sum1 + pow(z , 2 );
        Sum2 = Sum2 + cos(2*pi*z);
    }
    F = -20*exp(-0.2*sqrt(Sum1/dim)) -exp(Sum2/dim) + 20 + e + f_bias[5];

    return F; 
}

tFitness f_Schwefel2_22(int dim, double *s) 
{
    tFitness sum, currentGen, prod;

    sum = 0.0;
    prod = 1.0;

    for (int i = 0; i < dim; i++) 
    {
        currentGen = fabs(s[i]-f7[i]);
        sum += currentGen;
        prod *= currentGen;
    }

    return sum + prod;
}

tFitness f_Schwefel2_22NoDesplazamiento(int dim, double *s) 
{
    tFitness sum, currentGen, prod;

    sum = 0.0;
    prod = 1.0;

    for (int i = 0; i < dim; i++) 
    {
        currentGen = fabs(s[i]);
        sum += currentGen;
        prod *= currentGen;
    }

    return sum + prod;
}



tFitness f_Schwefel1_2(int dim, double *s)
{
   	tFitness Sum=0.0, Val=0.0;

	for (int i = 0; i < dim; i++)
    {  
        Val += s[i]-f8[i];
	    Sum += Val * Val;
    }

	return Sum;
}

tFitness f_10(double x, double y)
{
     double p, z, t;

     p=(x*x+y*y);
     
     z=pow(p, 0.25);
     t=sin(50.0*pow(p, 0.1));
     t=t*t+1.0;
     
     return z*t;
}

tFitness Extended_f_10(int dim, double *x)
{
     double suma=0.0;

     for(int i=0; i<dim-1; i++)
           suma+=f_10(x[i]-f9[i], x[i+1]-f9[i+1]);

     suma+=f_10(x[dim-1]-f9[dim-1], x[0]-f9[0]);

     return suma;
}

tFitness Extended_f_10NoDesplazamiento(int dim, double *x)
{
     double suma=0.0;

     for(int i=0; i<dim-1; i++)
           suma+=f_10(x[i], x[i+1]);

     suma+=f_10(x[dim-1], x[0]);

     return suma;
}

tFitness f_Bohachevsky(int dim, double *s) 
{   
    const double PI = 3.141592653589793;
    tFitness sum = 0.0;
    int i;
    double currentGen;
    double nextGen;

    currentGen = s[0]-f10[0];

    for (i = 1; i < dim; i++) 
    {
        nextGen = s[i]-f10[i];
        sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
	double c1=cos(3.0 * PI * currentGen);
	double c2=cos(4.0 * PI * nextGen);
        sum += 0.7 -(0.3*c1+0.4*c2);
        currentGen = nextGen;
    }

    return sum;
}

tFitness f_BohachevskyNoDesplazamiento(int dim, double *s) 
{   
    const double PI = 3.141592653589793;
    tFitness sum = 0.0;
    int i;
    double currentGen;
    double nextGen;

    currentGen = s[0];

    for (i = 1; i < dim; i++) 
    {
        nextGen = s[i];
        sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
	double c1=cos(3.0 * PI * currentGen);
	double c2=cos(4.0 * PI * nextGen);
        sum += 0.7 -(0.3*c1+0.4*c2);
        currentGen = nextGen;
    }

    return sum;
}


tFitness f_Schaffer(int dim, double *s) 
{
    int i;
    tFitness sum;
    double aux, aux2;
    double currentGen, nextGen;

    sum = 0.0;
    currentGen = s[0]-f11[0];
    currentGen = currentGen * currentGen;

    for (i = 1; i < dim; i++) 
    {
        nextGen = s[i]-f11[i];
        nextGen = nextGen * nextGen;
        aux = currentGen + nextGen;
        currentGen = nextGen;
        aux2 = sin(50. * pow(aux, 0.1));
        sum += pow(aux, 0.25) * (aux2 * aux2 + 1.0);
    }

    return sum;
}

static void divideFunctions(int dim, double *s, double *part1, double *part2, double m, int *psize1, int *psize2) {
    int shared;
    int rest, i, total;
    double *partrest;

    if (m <= 0.5) {
       partrest = part2;
    }
    else {
       partrest = part1;
       m = 1-m;
    }

    shared = (int) floor(dim*m);
    rest = 2*shared;

    for (i = 0; i < shared; i++) {
	part1[i] = s[i*2];
	part2[i] = s[i*2+1];
    }
    total = dim-shared;
    
    for (i = 0; i < total-shared; i++) {
	partrest[i+shared] = s[i+rest];
    }

    *psize1 = shared;
    *psize2 = dim-shared;

    if (partrest == part1) {
        int temp = *psize1;
	*psize1 = *psize2;
	*psize2 = temp;
    }
}

tFitness f_Hybrid_12(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
    tFitness f1, f2;
    divideFunctions(dim,s, part1, part2, 0.25, &size1, &size2);

    f1=Extended_f_10NoDesplazamiento(size1, part1);
    f2=Shifted_Sphere(size2,part2)-f_bias[0];
    assert(f1 >= 0);
    assert(f2 >= 0);

    return f1+f2;
}

tFitness f_Hybrid_13(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
    tFitness f1, f2;

    divideFunctions(dim,s, part1, part2, 0.25, &size1, &size2);
    f1=Extended_f_10NoDesplazamiento(size1, part1);
    f2=Shifted_Rosenbrock(size2,part2)-f_bias[2];

    return f1+f2;
}

tFitness f_Hybrid_14(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
    tFitness f1, f2;
    divideFunctions(dim,s, part1, part2, 0.25, &size1, &size2);

    f1=Extended_f_10NoDesplazamiento(size1, part1);
    f2=Shifted_Rastrigin(size2, part2)-f_bias[3];

    return f1+f2; 
}

tFitness f_Hybrid_15(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
	double desp[MAX_D];
    int size1, size2;
    tFitness f1, f2;

    for (int i = 0; i < dim; ++i) {
	desp[i] = s[i] - f15[i];
    }

    divideFunctions(dim, desp, part1, part2, 0.25, &size1, &size2);

    f1=f_BohachevskyNoDesplazamiento(size1, part1);
    f2=f_Schwefel2_22NoDesplazamiento(size2, part2);
    return f1+f2; 
}

tFitness f_Hybrid_16new(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
    tFitness f1, f2;
    divideFunctions(dim,s, part1, part2, 0.5, &size1, &size2);

    f1=Extended_f_10NoDesplazamiento(size1, part1);
    assert(f1 >= 0);
    f2=Shifted_Sphere(size2, part2)-f_bias[0];
    assert(f2 >= 0);

    return f1+f2;
}

tFitness f_Hybrid_17new(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
    tFitness f1, f2;

    divideFunctions(dim,s, part1, part2, 0.75, &size1, &size2);
    f1=Extended_f_10NoDesplazamiento(size1, part1);
    f2=Shifted_Rosenbrock(size2,part2)-f_bias[2];

    return f1+f2;
}

tFitness f_Hybrid_18new(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
    tFitness f1=0;
    tFitness f2=0;

    divideFunctions(dim,s, part1, part2, 0.75, &size1, &size2);

    f1=Extended_f_10NoDesplazamiento(size1, part1);
    f2=Shifted_Rastrigin(size2, part2)-f_bias[3];
    assert(isfinite(f1));
    assert(isfinite(f2));

    return f1+f2;
}

tFitness f_Hybrid_19new(int dim, double *s) 
{
	double part1[MAX_D], part2[MAX_D];
    int size1, size2;
	double desp[MAX_D];
    tFitness f1, f2;

    for (int i = 0; i < dim; ++i) {
	desp[i] = s[i] - f19[i];
    }

    divideFunctions(dim, desp, part1, part2, 0.75, &size1, &size2);

    f1=f_BohachevskyNoDesplazamiento(size1, part1);
    f2=f_Schwefel2_22NoDesplazamiento(size2, part2);

    return f1+f2;
}
