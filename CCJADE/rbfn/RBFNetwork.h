#pragma once
#include<vector>
#include<math.h>
#include <iostream>     
#include <algorithm>    
#include "KmeansPP.h"
#include "UtilRBFN.h"



typedef std::vector<double> datapoint;
typedef enum { rbfntOneSigmaPerCentre = 0, rbfntOneSigmaPerDir } rbfnType;
class RBFNetwork
{
public:
	RBFNetwork();
	RBFNetwork(const std::vector<datapoint> &training_data, const std::vector<tFitness> &training_values);
	RBFNetwork(const std::vector<datapoint> &training_data, const std::vector<tFitness> &training_values,
		       const std::vector<datapoint> &testing_data, const std::vector<tFitness> &testing_values);
	~RBFNetwork(void);

	void reset();

	/* Start Training the Radial Basis Function network
		Takes the number of RBF centroids, the learning rate, the number of iteration and a print flag as input
		Saves the output model to be used in testing and single predictions 
		return accuracy and mse (by reference)
		*/
	tFitness training(int num_rbf_units, double learning_rate, int num_iterations, tFitness toll, bool print_flag = false);
	//tFitness training(int num_rbf_units, double learning_rate, int num_iterations, tFitness toll, bool print_flag = false);

	/* Start Testing the RBF Network to make sure it's not overfitting
		(Should be done after training) */
	void testing(const std::vector<datapoint> &testing_data, const std::vector<tFitness> &testing_values);

	
	/* Predict a single data point*/
	tFitness predictValue(const datapoint &data_point);

	//TODO
	void saveModel();
	void loadModel();

	std::vector<datapoint> training_data;
	std::vector<tFitness> training_values;
	std::vector<datapoint> testing_data;
	std::vector<tFitness> testing_values;

private:	
	rbfnType type;
	std::vector<double> sigma;	
	std::vector< std::vector<tFitness> > hxc;
	std::vector<tFitness> layer2_weights;
	std::vector<datapoint> rbf_centroids;
	//std::vector<double>total_centroids_dist;

	// Random Number seed devices/engines/distributions
    std::random_device rd;
	std::default_random_engine random_engine;
	std::uniform_real_distribution<double> random_real_gen;

	void build_hxc();
	void initSigma();
	double sqrdistance(const datapoint &a ,const datapoint &b);
	double sqrdistance_sigmadir(const datapoint &a, const datapoint &b);
	double sqrdistance_sigmadir_center(const datapoint &a, const datapoint &b, const int centroidIndex);
	tFitness basisFunction(const datapoint &data_point, const int centroidIndex);
	
};

