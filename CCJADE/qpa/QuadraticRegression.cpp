//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : QuadraticRegression.cpp
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
// More details on the following paper:
//
// De Falco, I., Della Cioppa, A., Trunfio, G.A. 'Surrogate-assisted Cooperative Coevolution
// for Large-Scale Optimization of Computationally Expensive Objective Functions, submitted'
//
//=============================================================================================
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <conio.h>
#include <cmath>
#include <ctime>
#include <chrono>
#include "QuadraticRegression.h"

using namespace std;


int QuadraticRegression::maxModInColumnindex(double **A, int rows, int column, int starting_column)
{
    int maxvindex = 0;
    double maxv = -1.0e+10;

    for(int i=starting_column; i<rows; i++)
    {
        if(maxv < fabs(A[i][column]))
        {
            maxv = fabs(A[i][column]);
            maxvindex = i;
        }
    }
    return maxvindex;
}

void QuadraticRegression::solve(double **A, double *b, int N)
{
    int k,pivot;
    int pivotflag = 1;

    // Steps (1) and (2) (decomposition and solution of Ly = b)
    switch(pivotflag)
    {
    case 1: // Case in which pivoting is employed

        for(k=0; k<N-1; k++)
        {
            pivot = maxModInColumnindex(A, N, k, k);

            double *tmp = A[pivot];
            A[pivot] = A[k];
            A[k] = tmp;

            double t = b[pivot];
            b[pivot] = b[k];
            b[k] = t;

            for(int i=k+1; i<N; i++)
            {
                double l_ik = A[i][k]/A[k][k];
                for(int j=k; j<N; j++)
                    A[i][j] = A[i][j] - l_ik*A[k][j];
                b[i] = b[i] - l_ik*b[k];
            }
        }
        break;

    case 0:  // Case 0/default in which no pivoting is used

    default:
        for(k=0; k<N-1; k++)
            for(int i=k+1; i<N; i++)
            {
                double l_ik = A[i][k]/A[k][k];
                for(int j=k; j<N; j++)
                    A[i][j] = A[i][j] - l_ik*A[k][j];
                b[i] = b[i] - l_ik*b[k];
            }

    }

    // Step (3) (backsolveing to solve Ux=y)
    b[N-1] = b[N-1]/A[N-1][N-1];
    for(k=N-2; k>=0; k--)
    {
        for(int j=k+1; j<N; j++)
            b[k] -= A[k][j]*b[j];
        b[k] = b[k]/A[k][k];
    }
}


QuadraticRegression::QuadraticRegression(unsigned _dim, unsigned _numPoints)
{
    dim = _dim;
    numPoints = _numPoints;

    //numero di coefficienti incogniti
    nt = (dim + 1)*(dim + 2) / 2;

    if (numPoints < nt)
    {
        cerr << "Insufficient number of points for QPA" << endl;
        exit(1);
    }

    //mean = new double[dim];
    //loBound = new double[dim];
    //hiBound = new double[dim];

    X = new double*[numPoints];    // matrice numPoints x nt
    for (int h = 0; h<numPoints; h++)
        X[h] = new double[nt];

    B = new double*[nt];  // X^T X -> matrice nt x nt
    for (int h = 0; h<nt; h++)
        B[h] = new double[nt];

    c = new double[nt];

    ready = false;
}

QuadraticRegression::QuadraticRegression(std::vector< std::vector<double> > &points, std::vector< tFitness > &values)
{
    numPoints = points.size();
    dim = points[0].size();
    QuadraticRegression(dim, numPoints);
    setData(points, values);
}


void QuadraticRegression::setData(std::vector< std::vector<double> > &points, std::vector< tFitness > &values)
{
    if (numPoints != points.size())
    {
        cerr << "Error in QPA" << endl;
        exit(1);
    }

    if (dim != points[0].size())
    {
        cerr << "Error in QPA" << endl;
        exit(1);
    }

    nt = (dim + 1)*(dim + 2) / 2;

    if (numPoints != nt)
    {
        cerr << "Error in QPA" << endl;
        exit(1);
    }

    /*
    for(int k=0; k<dim; ++k)
    {
    mean[k] = 0.0;
    loBound[k] = 1.0E08;
    hiBound[k] = -1.0E08;
    }
    */

    for (int i = 0; i<numPoints; ++i)
    {
        X[i][0] = 1.0;

        for (int k = 0; k<dim; ++k)
        {
            X[i][1 + k] = points[i][k]; //linear terms

            /*
            mean[k] += points[i][k]/numPoints;

            if( points[i][k]>hiBound[k] )
            hiBound[k] = points[i][k];
            if( points[i][k]<loBound[k] )
            loBound[k] = points[i][k];
            */
        }

        int r = 0;
        for (int k = 0; k<dim; ++k)
            for (int q = k; q<dim; ++q, ++r)
                X[i][1 + dim + r] = points[i][k] * points[i][q]; //quadratic terms
    }

    for (int i = 0; i<nt; ++i)
        for (int j = 0; j<nt; ++j)
        {
            B[i][j] = 0;
            for (int k = 0; k<numPoints; ++k)
                B[i][j] += X[k][i] * X[k][j];
        }


    for (int i = 0; i<nt; ++i)
    {
        c[i] = 0;
        for (int k = 0; k<numPoints; ++k)
            c[i] += X[k][i] * values[k];
    }

    solve(B, c, nt);

    ready = true;
}


QuadraticRegression::~QuadraticRegression()
{
    for(int h=0; h<numPoints; h++)
        delete [] X[h];
    delete [] X;

    for(int h=0; h<nt; h++)
        delete [] B[h];
    delete [] B;

    delete [] c;

    // delete [] mean;
    //delete [] hiBound;
    //delete [] loBound;
}


tFitness QuadraticRegression::evaluate(std::vector<double> &point)
{
    tFitness v = c[0];
    int r = 1;



    for(int k=0; k<dim; ++k, ++r)
        v += c[r]*point[k]; //termini lineari

    for(int k=0; k<dim; ++k)
        for(int q=k; q<dim; ++q, ++r)
            v += c[r]*point[k]*point[q]; //termini quadratici

    return v;
}

std::vector<double> QuadraticRegression::evaluateGradient(std::vector<double> &point)
{
    std::vector<double> g;

    for(int i=0; i<dim; ++i)
    {
        double gi = c[1+i];

        int r = dim + 1;

        for(int k=0; k<dim; ++k)
            for(int q=k; q<dim; ++q, ++r)
            {
                //c[r]*[ D point[k]/D xi * point[q] + point[k] * D point[q]/D xi ]; //termini quadratici
                if( k==i )
                    gi += c[r]*point[q];
                if( q==i )
                    gi += c[r]*point[k];
            }

        g.push_back(gi);
    }

    return g;
}

/*
void readData(std::string sFilename, std::vector<std::vector<double>> &input, std::vector<double> &output, int nInputCols)
{
  ifstream is;
  is.open(sFilename.c_str(), ios::binary);

  // get length of file:
  is.seekg (0, std::ios::end);
  long length = is.tellg();
  is.seekg (0, std::ios::beg);

  // allocate memory:
  char *buffer = new char [length];

  // read data as a block:
  is.read (buffer,length);

  // create string stream of memory contents
  // NOTE: this ends up copying the buffer!!!
  string bs(buffer);
  std::istringstream iss(bs);

  // delete temporary buffer
  delete [] buffer;

  // close filestream
  is.close();

  std::string s;
  input.clear();
  output.clear();
  int r = 0;
  while ( std::getline(iss, s) )
  {
	 std::istringstream ss(s);
	 std::string field;
	 std::vector<double> rv;
	 input.push_back(rv);
	 int col = 0;
	 while ( std::getline(ss, field, ',' ) )
	 {
		 std::istringstream i(field);
         double val;
		 i >> val;
		 if( col<nInputCols )
		   input[r].push_back(val);
		 else
		   output.push_back(val);
		 ++col;
	 }
	 ++r;
  }
}

int _tmain(int argc, _TCHAR* argv[])
{
	std::vector< std::vector<double> > points;
	std::vector< double > values;
	readData("C:\\VCProjects\\SCIARA-SHARKSurrogateCalibration\\test.dat", points, values, 2);
	printf("Creating quadratic regression on %d points\n", points.size());

	QuadraticRegression qr(points, values);

	FILE *f=fopen("testqr.dat", "wt");
	for(int i=0; i<(int)points.size(); ++i)
	{
       for(int k=0; k<points[0].size(); ++k)
		 fprintf(f, "%lf ", points[i][k]);
	   fprintf(f, "%lf\n", qr.evaluate(points[i]));
	}
	fclose(f);

	system("pause");

	return 0;
}
*/
