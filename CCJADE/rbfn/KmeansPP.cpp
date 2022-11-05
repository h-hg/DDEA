#include "KmeansPP.h"

using namespace std;

KmeansPP::KmeansPP(const std::vector<datapoint> &input_data) 
	: input_data(input_data), random_index_gen(0, input_data.size()-1), random_real_gen(0, 1), random_engine(rd())
{}

KmeansPP::~KmeansPP(void)
{}

void KmeansPP::RunKMeansPP(int K, vector<datapoint> &centroids)
{
	assert(K<=input_data.size());
	vector<vector<int> > clusters_vec(K);
	// First: initalize the initial centroids according to K-Means Plus Plus Algorithm
	init();
	int first_centroid = random_index_gen(random_engine);
	initial_centroids_.push_back(input_data[first_centroid]);
	for(int i=1 ; i<K ; i++)
	{
		updateNearestClusters(initial_centroids_);
		initial_centroids_.push_back(input_data[getNextInitialCentroidIndex()]);
	}
	// Second: Continue as in the regular K-means clustering algoithm
	cur_centroids_ = initial_centroids_;
	do
	{
		prev_centroids_ = cur_centroids_;
		updateNearestClusters(cur_centroids_);
		updateCentroids(cur_centroids_);
	}
	while(!equalCentroids(cur_centroids_,prev_centroids_));

	// Push resuts into the clusters vector
	for(int i = 0 ; i<input_data.size() ; i++)
		clusters_vec[nearest_cluster_idx[i]].push_back(i);
	centroids = cur_centroids_;
}

int KmeansPP::getNextInitialCentroidIndex()
{
	/*	Adding the rest of the points according to the probability D(x)/SIGMA(D(x))
	where D(x) is the distance between a datapoint x and it's nearest cluster */

	// Total Error i.e SIGMA(D(x)
	double total_distance=0;
	for (int i = 0; i < input_data.size() ; ++i)
		total_distance += nearest_cluster_dist[i];

	// The probability D(x)/SIGMA(D(x))
	vector<double> cumm_prob(input_data.size(),0);
	for (int i = 0; i < input_data.size() ; ++i)
		cumm_prob[i] = (nearest_cluster_dist[i] / total_distance);

	// Cummulating
	for (int i = 1; i < input_data.size() ; ++i)
		cumm_prob[i] += cumm_prob[i-1];

	// Choosing the next point with a probabilty D(x)/SIGMA(D(x))
	int rand_num = random_real_gen(random_engine);
	for (int i = 0; i < cumm_prob.size() ; ++i)
		if (rand_num < cumm_prob[i])
			return i;

	return cumm_prob.size()- 1;
}

void KmeansPP::updateNearestClusters(const vector<datapoint> &centroids_)
{
	for (int i = 0; i < input_data.size() ; ++i)
	{
		int idx = getClosestCentroidIndex(i, centroids_);		
		nearest_cluster_idx[i] = idx;
		nearest_cluster_dist[i] = distance(centroids_[idx],input_data[i]);
	}
}

void KmeansPP::updateCentroids(const vector<datapoint> &centroids_)
{
	vector<int>freq(centroids_.size(),0);
	vector<datapoint>new_centroids(centroids_.size(),vector<double>(centroids_[0].size(),0));
	for (int i = 0; i < input_data.size() ; ++i) // training set count
	{
		++freq[nearest_cluster_idx[i]];
		for(int j=0 ; j< input_data[i].size() ; ++j) // feature count
			new_centroids[nearest_cluster_idx[i]][j] += input_data[i][j];
	}

	for (int i = 0; i < centroids_.size() ; ++i)
		if(freq[i])
			for (int j = 0; j < centroids_[i].size() ; ++j)
				new_centroids[i][j] *= (1.0/(double)freq[i]);
	cur_centroids_ = new_centroids;
}

int KmeansPP::getClosestCentroidIndex(int data_point_idx, const std::vector<datapoint> &centroids_)
{
	int closest_cluster=-1;
	double minDistance=1e18;
	for(int i = 0 ; i<centroids_.size() ; ++i)
	{
		double dist_dp_centroid = distance(centroids_[i],input_data[data_point_idx]);
		if(minDistance > dist_dp_centroid)
		{
			minDistance = dist_dp_centroid;
			closest_cluster = i;
		}
	}
	assert(closest_cluster!=-1);
	return closest_cluster;	
}

double KmeansPP::distance(const datapoint &a , const datapoint &b)
{
	double dist=0;
	for(int i = 0 ; i< a.size() ; i++)
		dist += (a[i]-b[i])*(a[i]-b[i]);
	return dist;
}

bool KmeansPP::equalCentroids(const vector<datapoint> &centroid_vec_a , const vector<datapoint> &centroid_vec_b)
{
	double total_dist=0;
	for(int i = 0 ; i< centroid_vec_a.size() ; i++)
		total_dist += distance(centroid_vec_a[i],centroid_vec_b[i]);
	return (total_dist<1e-9);
}


void KmeansPP::init()
{
	nearest_cluster_idx.assign(input_data.size(),-1);
	nearest_cluster_dist.assign(input_data.size(),0);
}

