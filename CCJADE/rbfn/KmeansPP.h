#pragma once

#include <vector>
#include <random>
#include <assert.h>

typedef std::vector<double> datapoint;
using namespace std;

class KmeansPP
{	

public:

	//typedef std::vector<double> datapoint;

	/*	Takes the input data where each data point is in a vector<double> format */
	KmeansPP(const std::vector<datapoint> &input_data);
	
	~KmeansPP(void);

	/*	Run the K-Means++ (Plus Plus) Algorithm
		Takes the number of desired clusters K and passed vector of datapoints (vector of doubles) as input
		Returns a vector for each cluster which includes the indices of their corresponding data points	*/
	void RunKMeansPP(int K, std::vector<datapoint> &centroids);

	
private:
	
	//Store the input data for many future runs without reinitializing
	std::vector<datapoint> input_data;

	std::vector<int> nearest_cluster_idx;
	std::vector<double> nearest_cluster_dist;
	std::vector<datapoint>initial_centroids_;
	std::vector<datapoint> cur_centroids_ ;
	std::vector<datapoint> prev_centroids_;

	// Random Number seed devices/engines/distributions
	//mt19937 mt;
    std::random_device rd;
	std::default_random_engine random_engine;
	std::uniform_int_distribution<int> random_index_gen;
	std::uniform_real_distribution<double> random_real_gen;


	void init();
	void updateNearestClusters(const std::vector<datapoint> &centroids_);
	void updateCentroids(const std::vector<datapoint> &centroids_);
	int getNextInitialCentroidIndex();
	int getClosestCentroidIndex(int data_point_idx, const std::vector<datapoint> &centroids_);
	double distance(const datapoint &a , const datapoint &b);
	bool equalCentroids(const std::vector<datapoint> &a , const std::vector<datapoint> &b);
	datapoint getMeanCentroid(const std::vector<datapoint> &centroids_);

};

