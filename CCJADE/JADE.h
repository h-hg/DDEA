//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Jade.h
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
// More details on the following paper:
//
// De Falco, I., Della Cioppa, A., Trunfio, G.A. 
// 'Investigating Surrogate-assisted Cooperative Coevolution for Large-Scale Global Optimization', 
// submitted'
//=============================================================================================

#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <vector>
#include <deque>
#include <map>
#include <cstddef>
#include "funsoft.h"
#include "dlib/optimization.h"
#include "dlib/global_optimization.h"
#include "dlib/svm.h"
#include "QuadraticRegression.h"
#include "RBFNetwork.h"
#include "gp.h"
#include "cg.h"
#include "gp_utils.h"
#include "rprop.h"

class Decomposer;

//using namespace dlib;
using namespace std;

typedef dlib::matrix<double, 0, 1> column_vector;

typedef enum { sNone = 0, sGP, sQPA, sRBFN, sSVR} typeOfSurrogate;

struct Pattern {
	vector<double> point;
	double fitness;
	double dist;
	Pattern(vector<double> & p, tFitness f) : point(p), fitness(f) {};
};


class JADE
{
	Decomposer &decomposer;

	typedef dlib::matrix<double, 0, 1> sample_type;
	typedef dlib::radial_basis_kernel<sample_type> kernel_type;
	


	struct doCompareIndividuals
	{
		doCompareIndividuals(const tFitness *_f) : f(_f) { }
		const tFitness *f;

		bool operator()(const int & i1, const int & i2)
		{
			return f[i1] < f[i2];
		}
	};

	struct doCompareIndividualsByVariance
	{
		doCompareIndividualsByVariance(const tFitness *_v) : v(_v) { }
		const tFitness *v;

		bool operator()(const int & i1, const int & i2)
		{
			return v[i1] > v[i2];
		}
	};

	struct doComparePatterns
	{
		doComparePatterns(){};

		bool operator()(const vector<Pattern>::iterator & i1, const vector<Pattern>::iterator & i2)
		{
			return i1->dist < i2->dist;
		}
	};

public:
	JADE(unsigned _dimension, unsigned _numberOfIndividuals, Decomposer &_group, typeOfSurrogate _sType);
	~JADE();
	void setCoordinates(vector<unsigned> &_coordinates);
	void setCoordinates(unsigned *coordinates, unsigned numOfCoordinates);
	void SAUpdate(typeOfSurrogate sType);	
	void update();
	void updateContextVector();
	void updateContextVector(vector<double> &cv, vector<unsigned> &coords, unsigned &vi);
	void sortPopulation(vector<tFitness> &fitness, vector<int> &sortIndex);
	void evaluateOffsprings(vector<bool> &toEvaluate, vector<bool> &hasTrueFitness);	
	int evaluateParents();	
	tFitness calculateFitnessValue(vector<double> &p, bool updateArchive=true);
	tFitness calculateSurrogateFitnessValue(vector<double> &p, typeOfSurrogate sType, bool &offspringHasTrueFitness);	
	tFitness calculateGPSurrogatePredictionVariance(vector<double> &p);	
	void optimize(int iterations);
	void updateIndexOfBest();
	void loadIndividuals(vector< vector<double> > &population);
	void updateIndividuals(vector< vector<double> > &population);
	void storeIndividuals(vector< vector<double> > &population);
	void setParentFitness(vector<tFitness> &fitnessValues);
	void addElementToArchive(vector<double> &individual, tFitness trueFitness);
	void createGP();	
	void optimizeGPParameters();
	void trainGlobalSurrogate();	
	void emptyArchive();
	vector<double> &getCollaborator();	
	double GPLogLikelihood(const column_vector &p);
	double GPLogLikelihoodD(double *p);
	const column_vector GPLogLikelihoodGradient(const column_vector &p);
	void findSVRparameters(std::vector<sample_type> &samples, std::vector<double> &targets, double &gamma, double &c);

	unsigned nfe;
	vector<unsigned> coordinates;	
	unsigned int dimension;
	unsigned int numGPPars;
	vector<double> gpPars;
	vector<double> gpPars_l, gpPars_u;

	void setQuadraticRegression(QuadraticRegression *q);
	void setRBFN(RBFNetwork *r);
	void setArchive(vector<Pattern> *a);

	unsigned int ite;
	double	JADE_mu_cr;
	double	JADE_mu_ff;

	///array containing the positions of all individuals
	vector< vector<double> > parents;
	vector< vector<double> > offsprings;
	vector<int> sortIndex;
	vector<double> FF;
	vector<double> CR;
	vector<double> SSFF;  // successful F values
	vector<double> SSCR;  // successful CR values 
	vector<int> binaryVector;
	vector< tFitness > parentsFitness;///array containing the current fitness of all particles 
	vector< bool > parentHasTrueFitness;
	vector< tFitness > offspringsFitness;
	vector< tFitness > offspringsVariance;

	tFitness bestFitness;

	///array containing the index of the best position attained so far
	unsigned indexOfBest;

	unsigned numberOfIndividuals;

	///buffer
	vector< double > xp;

	uniform_real_distribution<double> unifRandom;

	bool internalArchive;
	vector< Pattern > *archive;
		
	libgp::GaussianProcess *gp;
	libgp::RProp rprop;
	libgp::CG cg;
	string GPCovType;	
	
	tFitness currentMax;
	tFitness currentMin;
	vector<double> minCoordInArchive;
	vector<double> maxCoordInArchive;
	bool surrogateIsValid;		
	unsigned minNumberOfPatterns;
	unsigned maxNumberOfPatterns;	
	typeOfSurrogate sType;
	bool internalQR;
	QuadraticRegression *qr;
	bool internalRBFN;
	RBFNetwork *rbfn;

	dlib::svr_trainer<kernel_type> svrTrainer;
	dlib::decision_function<kernel_type> svr;
};

