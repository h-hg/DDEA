#include "RBFNetwork.h"

using namespace std;

RBFNetwork::RBFNetwork(const vector<datapoint> &training_data, const vector<tFitness> &training_values)
    : training_data(training_data), training_values(training_values), random_real_gen(-1, 1), random_engine(rd())
{
    type = rbfntOneSigmaPerCentre;
}

RBFNetwork::RBFNetwork()
    :  random_real_gen(-1, 1), random_engine(rd())
{
    type = rbfntOneSigmaPerDir;
}

RBFNetwork::RBFNetwork(const vector<datapoint> &training_data, const vector<tFitness> &training_values,
                       const vector<datapoint> &testing_data, const vector<tFitness> &testing_values)
    : training_data(training_data), training_values(training_values),
      testing_data(testing_data), testing_values(testing_values),
      random_real_gen(-1, 1), random_engine(rd())
{
    type = rbfntOneSigmaPerCentre;
}

RBFNetwork::~RBFNetwork(void)
{
    type = rbfntOneSigmaPerCentre;
}


void RBFNetwork::reset()
{
    training_data.clear();
    testing_data.clear();
    training_values.clear();
    testing_values.clear();
}

tFitness RBFNetwork::training(int num_rbf_units, double learning_rate, int num_iterations, tFitness toll, bool print_flag)
{


    if (training_data.size() == 0)
    {
        cout << "Not enough patterns in rbfn training" << endl;
        exit(1);
    }

    int dimension = training_data[0].size();

    if (print_flag)
    {
        printf("Starting RBF Network Training with %d units and learning rate=%f...\n", num_rbf_units, learning_rate);
        printf("Getting RBF Centroids using K-means++...\n");
    }

    // Calculate RBF Centroids
    KmeansPP KMPP(training_data);
    rbf_centroids.clear();
    KMPP.RunKMeansPP(num_rbf_units, rbf_centroids);

    initSigma();

    if (print_flag)
        printf("Building the basis function values...\n");

    // Init the Second Layer with values [-1,1]
    layer2_weights.assign(num_rbf_units, 0);
for (auto &w : layer2_weights)
        w = random_real_gen(random_engine);

    vector<int> indexes(training_data.size());
    for (int i = 0; i < training_data.size(); ++i)
        indexes[i] = i;

    build_hxc();

    // Train the second layer weights
    tFitness mse = 0, prev_testing_mse = 1.0E100;
    for (int iter = 0; iter < num_iterations; iter++)
    {
        random_shuffle(indexes.begin(), indexes.end());

        for (int jj = 0; jj < training_data.size(); jj++)
        {
            int j = indexes[jj];

            for (int i = 0; i < rbf_centroids.size(); i++)
                hxc[j][i] = basisFunction(training_data[j], i);

            //value of RBFN for the j-th training point
            tFitness Fj = 0.0;
            for (int i = 0; i < rbf_centroids.size(); i++)
                Fj += hxc[j][i]*layer2_weights[i];

            //error for the i-th training pattern
            tFitness error_j = training_values[j] - Fj;

            for (int q = 0; q < dimension; ++q)
                for (int i = 0; i < num_rbf_units; ++i)
                {
                    double deltaSigmai = 0.0;
                    double sigmai = sigma[i*dimension + q];
                    double sigmai2 = sigmai * sigmai;
                    double sigmai3 = sigmai2 * sigmai;
                    double d = (training_data[j][q] - rbf_centroids[i][q]);
                    double h = error_j * layer2_weights[i] * hxc[j][i] * d;
                    sigma[dimension*i + q] += learning_rate * h*d / sigmai3;
                    rbf_centroids[i][q] += learning_rate * h / sigmai2;
                }

            for (int i = 0; i < num_rbf_units; i++)
                layer2_weights[i] += learning_rate * hxc[j][i] * error_j;

            mse += error_j*error_j;
        }


        if (testing_data.size()>0)
        {
            tFitness testing_mse = 0.0;
            for (int i = 0; i<testing_data.size(); i++)
            {
                tFitness error = testing_values[i] - predictValue(testing_data[i]);
                testing_mse += error * error;
            }
            testing_mse *= (double)(1.0 / (double)testing_data.size());

            if (testing_mse > prev_testing_mse)
                break;
            prev_testing_mse = testing_mse;
        }


        if (print_flag)
        {
            // Gathering Statistics
            mse = 0.0;
            for (int i = 0; i<training_data.size(); i++)
            {
                tFitness error = training_values[i] - predictValue(training_data[i]);
                mse += error * error;
            }
            mse *= (double)(1.0 / (double)training_data.size());
            printf("Training (%*d/%d), MSE=[%.3f], Progress [%.2f]\n", 2, (iter + 1), num_iterations, mse, ((double)((double)(iter + 1) / (double)num_iterations) * 100.0));
        }

        //if ( mse < toll )
        //break;
    }
    if (print_flag)
        printf("\n----------------------------\n");

    return mse;
}



void RBFNetwork::build_hxc()
{
    //total_centroids_dist.clear();
    //total_centroids_dist.assign(rbf_centroids.size(),0);
    //// Calculate total distance for each centroid
    //for(int i = 0 ; i<rbf_centroids.size() ; i++)
    //	for(auto &data_point : training_data)
    //		total_centroids_dist[i] += distance(data_point,rbf_centroids[i]);

    // Build RBF Units
    hxc.assign(training_data.size(),vector<tFitness>());
    for(int i = 0 ; i<training_data.size() ; i++)
    {
        for(int j = 0 ; j<rbf_centroids.size() ; j++)
            hxc[i].push_back(basisFunction(training_data[i], j));
    }
}

tFitness RBFNetwork::basisFunction(const datapoint &data_point, const int centroidIndex)
{
    if ( type==rbfntOneSigmaPerCentre )
        return exp(-sqrdistance(data_point, rbf_centroids[centroidIndex]) / (2 * sigma[centroidIndex] * sigma[centroidIndex]));
    else if (type == rbfntOneSigmaPerDir)
    {
        //return exp(-sqrdistance_sigmadir(data_point, rbf_centroids[centroidIndex]));
        return exp(-sqrdistance_sigmadir_center(data_point, rbf_centroids[centroidIndex], centroidIndex));
    }
    else
    {
        cout << "Unknown rbfn type" << endl;
        exit(1);
    }

}

double RBFNetwork::sqrdistance(const datapoint &a, const datapoint &b)
{
    double dist=0;
    for(int i = 0 ; i< a.size() ; i++)
        dist += (a[i]-b[i])*(a[i]-b[i]);
    return dist;
}

double RBFNetwork::sqrdistance_sigmadir(const datapoint &a, const datapoint &b)
{
    double dist = 0;
    for (int i = 0; i< a.size(); i++)
        dist += (a[i] - b[i])*(a[i] - b[i])/(2*sigma[i]*sigma[i]);
    return dist;
}

double RBFNetwork::sqrdistance_sigmadir_center(const datapoint &x, const datapoint &c, const int centroidIndex)
{
    double dist = 0;
    int dimension = x.size();
    for (int i = 0; i < dimension; i++)
    {
        double s = sigma[centroidIndex*dimension + i];
        dist += (x[i] - c[i])*(x[i] - c[i]) / (2*s*s);
    }
    return dist;
}

tFitness RBFNetwork::predictValue(const datapoint &data_point)
{
    vector<tFitness> cur_rbf_unit;
    for(int j = 0 ; j<rbf_centroids.size() ; j++)
        cur_rbf_unit.push_back(basisFunction(data_point, j));
    tFitness v = Utility::multiplyVectors(cur_rbf_unit, layer2_weights);
    return v;
}

void RBFNetwork::initSigma()
{
    double variance = -1;
    for (int i = 0; i < rbf_centroids.size() ; ++i)
        for (int j = i + 1; j < rbf_centroids.size() ; ++j)
            variance = std::max(variance, sqrdistance(rbf_centroids[i], rbf_centroids[j]));
    variance *= (1.0 / ((double)rbf_centroids.size()));

    if ( type == rbfntOneSigmaPerCentre )
        sigma.assign(rbf_centroids.size(), sqrt(variance));
    else if (type == rbfntOneSigmaPerDir)
        //sigma.assign(rbf_centroids[0].size(), sqrt(variance));
        sigma.assign(rbf_centroids.size()*rbf_centroids[0].size(), sqrt(variance));
    else
    {
        cout << "Unknown rbfn type" << endl;
        exit(1);
    }
}


void RBFNetwork::testing(const std::vector<datapoint> &testing_data, const std::vector<tFitness> &testing_values)
{
    printf("Testing...\n");
    tFitness mse = 0, err = 0;
    for(int i = 0 ; i<testing_data.size() ; i++)
    {
        err = predictValue(testing_data[i]) - testing_values[i];
        mse += err*err;
    }
    mse *= (1.0/(double)testing_data.size());
    printf("Testing Results MSE = %.6f\n", mse);
    printf("------------------------------\n");
}
