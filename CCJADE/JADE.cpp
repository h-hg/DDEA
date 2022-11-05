//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Jade.cpp
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

#include "JADE.h"
#include "CCDE.h"
#include <chrono>

using namespace std;

//******************************************************************************************/
//
//
//
//******************************************************************************************/
JADE::JADE(unsigned _dimension, unsigned _numberOfIndividuals, Decomposer &_group, typeOfSurrogate _sType) :
    decomposer(_group), dimension(_dimension), numberOfIndividuals(_numberOfIndividuals), sType(_sType)
{
    JADE_mu_cr = 0.5;
    JADE_mu_ff = 0.5;
    nfe = 0;
    parentsFitness.resize(numberOfIndividuals, 0);
    parentHasTrueFitness.resize(numberOfIndividuals, false);
    offspringsFitness.resize(numberOfIndividuals, 0);
    FF.resize(numberOfIndividuals, 0);
    CR.resize(numberOfIndividuals, 0);
    SSFF.resize(numberOfIndividuals, 0);
    SSCR.resize(numberOfIndividuals, 0);
    indexOfBest = (unsigned)(numberOfIndividuals * unifRandom(decomposer.eng));
    binaryVector.resize(dimension, 0);
    xp.resize(decomposer.CCOptimizer.problemDimension);

    minCoordInArchive.resize(dimension, 0);
    maxCoordInArchive.resize(dimension, 0);

    numGPPars = dimension + 3;
    gpPars.resize(numGPPars);
    gpPars_l, gpPars_u;
    for (int i = 0; i < dimension; ++i)
    {
        gpPars_l.push_back(log(1.0E-02));
        gpPars_u.push_back(log(10.0));
    }
    gpPars_l.push_back(log(1.0E-03));
    gpPars_u.push_back(log(1.0));

    gpPars_l.push_back(log(1.0E-03));
    gpPars_u.push_back(log(1.0));

    gpPars_l.push_back(log(1.0E-9));
    gpPars_u.push_back(log(1.0E-2));

    for (int i = 0; i < numGPPars; ++i)
        gpPars[i] = (gpPars_l[i] + gpPars_u[i]) / 2.0;

    gp = NULL;
    rbfn = NULL;
    qr = NULL;
    surrogateIsValid = false;

    if (sType == sQPA || sType==sNone )
    {
        minNumberOfPatterns = ((dimension + 1)*(dimension + 2)) / 2;
        maxNumberOfPatterns = 3 * minNumberOfPatterns;
    }
    else if (sType == sRBFN || sType == sSVR )
    {
        minNumberOfPatterns = numberOfIndividuals; //minimum size of archive to be used for generating RBFN surrogates
        maxNumberOfPatterns = 10000; //use all patterns
    }
    else
    {
        minNumberOfPatterns = numberOfIndividuals; //minimum size of archive to be used for generating GP surrogates
        maxNumberOfPatterns = 80; //older patterns are ignored
    }


    internalQR = internalRBFN = internalArchive = true;
    if ( sType == sQPA ) qr = new QuadraticRegression(dimension, minNumberOfPatterns);
    else if (sType == sRBFN) rbfn = new RBFNetwork();
    archive = new vector< Pattern >;

}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
JADE::~JADE()
{
    if ( internalQR )
        delete qr;
    if ( internalRBFN )
        delete rbfn;
    if ( internalArchive )
        delete archive;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
vector<double> &JADE::getCollaborator()
{
    return parents[indexOfBest];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setQuadraticRegression(QuadraticRegression *q)
{
    if (internalQR)
        delete qr;
    qr = q;
    internalQR = false;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setRBFN(RBFNetwork *r)
{
    if ( internalRBFN )
        delete rbfn;
    rbfn = r;
    internalRBFN = false;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setArchive(vector<Pattern> *a)
{
    if (internalArchive)
        delete archive;
    archive = a;
    internalArchive = false;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateContextVector()
{
    for (unsigned ld = 0; ld < coordinates.size(); ld++)
        decomposer.contextVector[coordinates[ld]] = parents[indexOfBest][ld];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateContextVector(vector<double> &cv, vector<unsigned> &coords, unsigned &vi)
{
    for (unsigned ld = 0; ld < coordinates.size(); ld++)
    {
        cv[vi + ld] = parents[indexOfBest][ld];
        coords[vi + ld] = coordinates[ld];
    }
    vi += coordinates.size();
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::createGP()
{
    delete gp;
    GPCovType = "CovSEard";
    gp = new libgp::GaussianProcess(dimension, GPCovType);
    int np = gp->covf().get_param_dim();
    Eigen::VectorXd params(np);

    for (int i = 0; i < np; ++i)
        params(i) = (params(i) + params(i)) / 2.0;

    gp->covf().set_loghyper(params);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setCoordinates(vector<unsigned> &_coordinates)
{
    coordinates = _coordinates;
};


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setCoordinates(unsigned *_coordinates, unsigned numOfCoordinates)
{
    if (coordinates.size() == numOfCoordinates)
    {
        for (unsigned i = 0; i < numOfCoordinates; ++i)
            coordinates[i] = _coordinates[i];
        return;
    }

    coordinates.clear();
    for (unsigned i = 0; i < numOfCoordinates; ++i)
        coordinates.push_back(_coordinates[i]);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::loadIndividuals(vector< vector<double> > &population)
{
    parents.clear();

    for (unsigned i = 0; i < numberOfIndividuals && i < population.size(); ++i)
    {
        vector< double > position;

        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            position.push_back(population[i][coordinates[ld]]);

        parents.push_back(position);
    }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateIndividuals(vector< vector<double> > &population)
{
    for (unsigned i = 0; i < numberOfIndividuals && i < population.size(); ++i)
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            parents[i][ld] = population[i][coordinates[ld]];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::storeIndividuals(vector< vector<double> > &population)
{
    for (unsigned i = 0; i < numberOfIndividuals && i < population.size(); ++i)
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            population[i][coordinates[ld]] = parents[i][ld];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::sortPopulation(vector<tFitness> &fitness, vector<int> &sortIndex)
{
    sortIndex.resize(fitness.size());

    for (unsigned int j = 0; j < fitness.size(); j++)
        sortIndex[j] = j;

    sort(&sortIndex[0], &sortIndex[0] + sortIndex.size(), doCompareIndividuals(&fitness[0]));
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::SAUpdate(typeOfSurrogate sType)
{
    //Sort the population from best to worst
    sortPopulation(parentsFitness, sortIndex);
    //Generate the CR and F values based on Gaussian and Cauchy distribution, respectively
    cauchy_distribution<double> cauchy(JADE_mu_ff, 0.1);

    normal_distribution<double> gaussian(JADE_mu_cr, 0.1);

    for (unsigned int i = 0; i < parents.size(); i++)
    {
        do
        {
            FF[i] = cauchy(decomposer.eng);
        } while (FF[i] <= 0.0);

        if (FF[i] > 1.0) FF[i] = 1.0;

        CR[i] = gaussian(decomposer.eng);

        if (CR[i] < 0.0) CR[i] = 0.0;

        if (CR[i] > 1.0) CR[i] = 1.0;
    }

    offsprings.clear();
    vector<bool> toEvaluate;
    for (unsigned int i = 0; i < parents.size(); i++)
    {
        unsigned r1, r2, r3;
        //Generate the mutant vector
        //Randomly choose the p_best individual
        unsigned p_index = unifRandom(decomposer.eng) * parents.size() * decomposer.CCOptimizer.JADE_p;
        p_index = sortIndex[p_index];

        //Select three parents randomly
        do
        {
            r1 = unifRandom(decomposer.eng) * parents.size();
        } while (r1 == i);

        do
        {
            r2 = unifRandom(decomposer.eng) * parents.size();
        } while (r2 == i || r2 == r1);

        do
        {
            r3 = unifRandom(decomposer.eng) * parents.size();
        } while (r3 == i || r3 == r2 || r3 == r1);

        vector<double> child(dimension, 0);

        for (unsigned int j = 0; j < dimension; j++)
        {
            if (decomposer.CCOptimizer.JADE_mutationStrategy == 1)
            {
                child[j] = parents[i][j] +
                           FF[i] * (parents[p_index][j] - parents[i][j]) +
                           FF[i] * (parents[r1][j] - parents[r2][j]);
            }
            else if (decomposer.CCOptimizer.JADE_mutationStrategy == 2)
            {
                child[j] = parents[r1][j] +
                           FF[i] * (parents[p_index][j] - parents[r1][j]) +
                           FF[i] * (parents[r2][j] - parents[r3][j]);
            }

            if (child[j] < decomposer.CCOptimizer.lowerLimit || child[j] > decomposer.CCOptimizer.upperLimit)
                child[j] = decomposer.CCOptimizer.lowerLimit + unifRandom(decomposer.eng) * (decomposer.CCOptimizer.upperLimit - decomposer.CCOptimizer.lowerLimit);
        }

        offsprings.push_back(child);
        toEvaluate.push_back(false);

        //Generate the binary vector based on the binomial crossover
        unsigned j_rnd = dimension * unifRandom(decomposer.eng);

        for (unsigned j = 0; j < dimension; j++)
        {
            if (unifRandom(decomposer.eng) < CR[i] || j == j_rnd)
                binaryVector[j] = 1;
            else
                binaryVector[j] = 0;
        }

        //Repair the crossover rate with its binary vector generated above
        unsigned int tt = 0;

        for (unsigned int j = 0; j < dimension; j++)
            tt += binaryVector[j];

        CR[i] = (double)tt / ((double)dimension);

        //Generate the trial vector based on the binary vector and the mutant vector
        for (unsigned int j = 0; j < dimension; j++)
            if (binaryVector[j] == 0)
                offsprings[i][j] = parents[i][j];
            else
                toEvaluate[i] = true;
    }

    //Evaluate the offspring population
    vector<bool> offspringHasTrueFitness(offsprings.size(), true);
    vector<bool> alreadyEvaluated(offsprings.size(), false);
    bool failureInSurrogate = false;

    offspringsVariance.resize(offsprings.size());

    if (archive->size() <  minNumberOfPatterns)
    {
        evaluateOffsprings(toEvaluate, offspringHasTrueFitness);

        for (int j = 0; j < offsprings.size(); ++j)
            alreadyEvaluated[j] = true;
    }
    else
    {
        if (sType == sGP || sType == sRBFN || sType == sSVR )
            trainGlobalSurrogate();

        unsigned nSurrogateEvals = 0;
        for (unsigned id = 0; id < offsprings.size(); ++id)
        {
            if (toEvaluate[id])
            {
                if (!alreadyEvaluated[id])
                {
                    bool isTrueFitness = false;
                    offspringsFitness[id] = calculateSurrogateFitnessValue(offsprings[id], sType, isTrueFitness);
					//double tf = calculateFitnessValue(offsprings[id], true);
					//cout << tf << " " << fabs(tf - offspringsFitness[id]) << endl;
                    offspringHasTrueFitness[id] = isTrueFitness;
                    alreadyEvaluated[id] = true;

                    if (isinf(offspringsFitness[id]) || isnan(offspringsFitness[id]))
                    {
                        failureInSurrogate = true;
						cout << "failure in surrogate" << endl;
                        offspringsFitness[id] = calculateFitnessValue(offsprings[id], true);
                        offspringHasTrueFitness[id] = true;
                        offspringsVariance[id] = 0.0;
                    }
                    else
                    {
                        tFitness var = 0.0;
                        if ( sType==sGP )
                            var = calculateGPSurrogatePredictionVariance(offsprings[id]);
                        offspringsVariance[id] = var;
                    }
                }
            }
            else
            {
                offspringsFitness[id] = parentsFitness[id];
                offspringsVariance[id] = 0.0;
                offspringHasTrueFitness[id] = parentHasTrueFitness[id];
            }

        }

    }

    vector<int> sortIndex;
    sortIndex.resize(offspringsVariance.size());
    for (unsigned int j = 0; j < offspringsVariance.size(); j++)
        sortIndex[j] = j;


    //Evaluate with the true fitness the offspring individual with the highest variance
    if ( sType == sGP )
    {
        int imv = 0;
        tFitness maxV = offspringsVariance[0];
        for (int q=1; q<sortIndex.size(); ++q)
            if (offspringsVariance[q]>maxV)
            {
                maxV = offspringsVariance[q];
                imv = q;
            }

        if (!offspringHasTrueFitness[imv])
        {
            offspringsFitness[imv] = calculateFitnessValue(offsprings[imv]);
            offspringHasTrueFitness[imv] = true;
        }
    }


    //Evaluate with the exact fitness the best individual
	
	sort(&sortIndex[0], &sortIndex[0] + sortIndex.size(), doCompareIndividuals(&offspringsFitness[0]));
    while (!offspringHasTrueFitness[sortIndex[0]] )
    {
        offspringsFitness[sortIndex[0]] = calculateFitnessValue(offsprings[sortIndex[0]]);
        offspringHasTrueFitness[sortIndex[0]] = true;
        sort(&sortIndex[0], &sortIndex[0] + sortIndex.size(), doCompareIndividuals(&offspringsFitness[0]));
    }

    //Selection and save the successful parameters
    SSFF.clear();
    SSCR.clear();

    for (unsigned int i = 0; i < parents.size(); i++)
    {
        if (offspringsFitness[i] <= parentsFitness[i])
        {

            for (unsigned int j = 0; j < dimension; j++)
                parents[i][j] = offsprings[i][j];

            parentsFitness[i] = offspringsFitness[i];
            parentHasTrueFitness[i] = offspringHasTrueFitness[i];

            //Save the successful CR and F values
            SSFF.push_back(FF[i]);
            SSCR.push_back(CR[i]);
        }
    }

    //Update mu_CR and mu_F based on the successful CRs and Fs
    if (SSCR.size())
    {
        //Update the mu_CR values
        double mean_cr = 0.0;
        for (unsigned int i = 0; i < SSCR.size(); i++)
            mean_cr += SSCR[i];
        mean_cr = mean_cr / ((double)SSCR.size());

        JADE_mu_cr = (1 - decomposer.CCOptimizer.JADE_c) * JADE_mu_cr + decomposer.CCOptimizer.JADE_c * mean_cr;

        //Update the mu_F value
        double mean_ff = 0.0;
        double t1 = 0.0;
        double t2 = 0.0;

        for (unsigned int i = 0; i < SSFF.size(); i++)
        {
            t1 += SSFF[i] * SSFF[i];
            t2 += SSFF[i];
        }

        mean_ff = t1 / t2;

        //Lehmer mean
        JADE_mu_ff = (1 - decomposer.CCOptimizer.JADE_c) * JADE_mu_ff + decomposer.CCOptimizer.JADE_c * mean_ff;
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::update()
{
    //Sort the population from best to worst
    sortPopulation(parentsFitness, sortIndex);

    //Generate the CR and F values based on Gaussian and Cauchy distribution, respectively
    cauchy_distribution<double> cauchy(JADE_mu_ff, 0.1);

    normal_distribution<double> gaussian(JADE_mu_cr, 0.1);

    for (unsigned int i = 0; i < parents.size(); i++)
    {
        do
        {
            FF[i] = cauchy(decomposer.eng);
        } while (FF[i] <= 0.0);

        if (FF[i] > 1.0) FF[i] = 1.0;

        CR[i] = gaussian(decomposer.eng);

        if (CR[i] < 0.0) CR[i] = 0.0;

        if (CR[i] > 1.0) CR[i] = 1.0;
    }

    offsprings.clear();
    vector<bool> toEvaluate;
    for (unsigned int i = 0; i < parents.size(); i++)
    {
        unsigned r1, r2, r3;
        //Generate the mutant vector
        //Randomly choose the p_best individual
        unsigned p_index = unifRandom(decomposer.eng) * parents.size() * decomposer.CCOptimizer.JADE_p;
        p_index = sortIndex[p_index];

        //Select three parents randomly
        do
        {
            r1 = unifRandom(decomposer.eng) * parents.size();
        } while (r1 == i);

        do
        {
            r2 = unifRandom(decomposer.eng) * parents.size();
        } while (r2 == i || r2 == r1);

        do
        {
            r3 = unifRandom(decomposer.eng) * parents.size();
        } while (r3 == i || r3 == r2 || r3 == r1);

        vector<double> child(dimension, 0);

        for (unsigned int j = 0; j < dimension; j++)
        {
            if (decomposer.CCOptimizer.JADE_mutationStrategy == 1)
            {
                child[j] = parents[i][j] +
                           FF[i] * (parents[p_index][j] - parents[i][j]) +
                           FF[i] * (parents[r1][j] - parents[r2][j]);
            }
            else if (decomposer.CCOptimizer.JADE_mutationStrategy == 2)
            {
                child[j] = parents[r1][j] +
                           FF[i] * (parents[p_index][j] - parents[r1][j]) +
                           FF[i] * (parents[r2][j] - parents[r3][j]);
            }

            if (child[j] < decomposer.CCOptimizer.lowerLimit || child[j] > decomposer.CCOptimizer.upperLimit)
                child[j] = decomposer.CCOptimizer.lowerLimit + unifRandom(decomposer.eng) * (decomposer.CCOptimizer.upperLimit - decomposer.CCOptimizer.lowerLimit);
        }

        offsprings.push_back(child);
        toEvaluate.push_back(false);

        //Generate the binary vector based on the binomial crossover
        unsigned j_rnd = dimension * unifRandom(decomposer.eng);

        for (unsigned j = 0; j < dimension; j++)
        {
            if (unifRandom(decomposer.eng) < CR[i] || j == j_rnd)
                binaryVector[j] = 1;
            else
                binaryVector[j] = 0;
        }

        //Repair the crossover rate with its binary vector generated above
        unsigned int tt = 0;

        for (unsigned int j = 0; j < dimension; j++)
            tt += binaryVector[j];

        CR[i] = (double)tt / ((double)dimension);

        //Generate the trial vector based on the binary vector and the mutant vector
        for (unsigned int j = 0; j < dimension; j++)
            if (binaryVector[j] == 0)
                offsprings[i][j] = parents[i][j];
            else
                toEvaluate[i] = true;
    }

    //Evaluate the child population
    vector<bool> offspringHasTrueFitness(offsprings.size(), true);
    evaluateOffsprings(toEvaluate, offspringHasTrueFitness);


    //Selection and save the successful parameters
    SSFF.clear();
    SSCR.clear();

    for (unsigned int i = 0; i < parents.size(); i++)
    {
        if (offspringsFitness[i] <= parentsFitness[i])
        {
            for (unsigned int j = 0; j < dimension; j++)
                parents[i][j] = offsprings[i][j];

            parentsFitness[i] = offspringsFitness[i];

            //Save the successful CR and F values
            SSFF.push_back(FF[i]);
            SSCR.push_back(CR[i]);
        }
    }

    //Update mu_CR and mu_F based on the successful CRs and Fs
    if (SSCR.size())
    {
        //Update the mu_CR values
        double mean_cr = 0.0;

        for (unsigned int i = 0; i < SSCR.size(); i++)
            mean_cr += SSCR[i];

        mean_cr = mean_cr / ((double)SSCR.size());
        JADE_mu_cr = (1 - decomposer.CCOptimizer.JADE_c) * JADE_mu_cr + decomposer.CCOptimizer.JADE_c * mean_cr;

        //Update the mu_F value
        double mean_ff = 0.0;
        double t1 = 0.0;
        double t2 = 0.0;

        for (unsigned int i = 0; i < SSFF.size(); i++)
        {
            t1 += SSFF[i] * SSFF[i];
            t2 += SSFF[i];
        }

        mean_ff = t1 / t2;

        //Lehmer mean
        JADE_mu_ff = (1 - decomposer.CCOptimizer.JADE_c) * JADE_mu_ff + decomposer.CCOptimizer.JADE_c * mean_ff;
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setParentFitness(vector<tFitness> &fitnessValues)
{
    parentsFitness = fitnessValues;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::evaluateOffsprings(vector<bool> &toEvaluate, vector<bool> &hasTrueFitness)
{
    //use real fitness function
    for (unsigned d = 0; d < decomposer.CCOptimizer.problemDimension; ++d)
        xp[d] = decomposer.contextVector[d];

    offspringsFitness.resize(offsprings.size());

    for (unsigned i = 0; i < offsprings.size(); i++)
    {
        if ( !toEvaluate[i] )
        {
            offspringsFitness[i] = parentsFitness[i];
            hasTrueFitness[i] = parentHasTrueFitness[i];
            continue;
        }

        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = offsprings[i][ld];

        offspringsFitness[i] = decomposer.CCOptimizer.computeFitnessValue(xp);
        nfe++;
        hasTrueFitness[i] = true;
        addElementToArchive(offsprings[i], offspringsFitness[i]);

        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = decomposer.contextVector[coordinates[ld]];
    }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
int JADE::evaluateParents()
{
    //use real fitness function
    for (unsigned d = 0; d < decomposer.CCOptimizer.problemDimension; ++d)
        xp[d] = decomposer.contextVector[d];

    parentsFitness.resize(parents.size());

    for (unsigned i = 0; i < parents.size(); i++)
    {
        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = parents[i][ld];

        parentsFitness[i] = decomposer.CCOptimizer.computeFitnessValue(xp);
        nfe++;
        parentHasTrueFitness[i] = true;
        addElementToArchive(parents[i], parentsFitness[i]);

        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = decomposer.contextVector[coordinates[ld]];
    }

    updateIndexOfBest();

    return parents.size();
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness JADE::calculateFitnessValue(vector<double> &p, bool updateArchive)
{
    for (unsigned d = 0; d < decomposer.CCOptimizer.problemDimension; ++d)
        xp[d] = decomposer.contextVector[d];

    for (unsigned ld = 0; ld < coordinates.size(); ld++)
        xp[coordinates[ld]] = p[ld];

    tFitness f = decomposer.CCOptimizer.computeFitnessValue(xp);

    nfe++;

    if ( updateArchive )
        addElementToArchive(p, f);
    return f;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
double JADE::GPLogLikelihood(const column_vector &p)
{
    int np = gp->covf().get_param_dim();
    Eigen::VectorXd params(np);
    for (int i = 0; i < np; ++i)
        params(i) = exp(p(i, 0));
    gp->covf().set_loghyper(params);
    return gp->log_likelihood();
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
double JADE::GPLogLikelihoodD(double *p)
{
    int np = gp->covf().get_param_dim();
    Eigen::VectorXd params(np);
    for (int i = 0; i < np; ++i)
        params(i) = exp(p[i]);
    gp->covf().set_loghyper(params);
    return -gp->log_likelihood();
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
const column_vector JADE::GPLogLikelihoodGradient(const column_vector &p)
{
    int np = gp->covf().get_param_dim();
    Eigen::VectorXd params(np);
    for (int i = 0; i < np; ++i)
        params(i) = exp(p(i, 0));
    gp->covf().set_loghyper(params);
    Eigen::VectorXd grad = gp->log_likelihood_gradient();

    column_vector g(np, 1);
    for (int i = 0; i < np; ++i)
        g(i, 0) = grad(i);

    return g;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::optimizeGPParameters()
{
    using namespace dlib;
    using namespace std::placeholders;

    int np = gp->covf().get_param_dim();

    column_vector xc(np, 1);
    column_vector lowerGPPars(np, 1), upperGPPars(np, 1);

    Eigen::VectorXd params = gp->covf().get_loghyper();

    for (int i = 0; i < np; ++i)
    {
        lowerGPPars(i, 0) = gpPars_l[i];
        upperGPPars(i, 0) = gpPars_u[i];
    }

    for (int i = 0; i < np; ++i)
        gpPars[i] = (gpPars_l[i] + gpPars_u[i])*unifRandom(decomposer.eng);

    for (int i = 0; i < np; ++i)
        xc(i, 0) = gpPars[i];

    std::function<double(const column_vector &)> func = std::bind(&JADE::GPLogLikelihood, this, _1);
    std::function<const column_vector(const column_vector &)> der = std::bind(&JADE::GPLogLikelihoodGradient, this, _1);
    surrogateIsValid = true;
    try
    {
        find_max_box_constrained(lbfgs_search_strategy(5), objective_delta_stop_strategy(1e-08, 100), func, der, xc, lowerGPPars, upperGPPars);
    }
    catch (...)
    {
        surrogateIsValid = false;
    }

    for (int i = 0; i < np; ++i)
        params(i) = exp(xc(i, 0));

    gp->covf().set_loghyper(params);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateIndexOfBest()
{
    bestFitness = std::numeric_limits<tFitness>::infinity();
    for (unsigned i = 0; i < parents.size(); ++i)
        if (parentsFitness[i] < bestFitness)
        {
            indexOfBest = i;
            bestFitness = parentsFitness[i];
        }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::emptyArchive()
{
    archive->clear();
    surrogateIsValid = false;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::addElementToArchive(vector<double> &individual, tFitness trueFitness)
{
    //check for duplicates
    for (auto it = archive->begin(); it != archive->end(); ++it)
    {
        double d = 0;
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            d += (individual[ld] - it->point[ld])*(individual[ld] - it->point[ld]);
        if ( d < 1.0E-18 )
        {
            return;
        }
    }

    archive->push_back(Pattern(individual, trueFitness));
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness JADE::calculateSurrogateFitnessValue(vector<double> &p, typeOfSurrogate sType, bool &isTrueFitness)
{
    isTrueFitness = false;

    //Avoid using the surrogate if already converged
    if (fabs(decomposer.CCOptimizer.globalBestFitness - decomposer.CCOptimizer.optimum) <1.0E-16)
    {
        isTrueFitness = true;
		cout << "Already converged" << endl;
        return calculateFitnessValue(p);
    }

    clock_t startTime = clock();

    if (sType == sGP || sType == sRBFN || sType== sSVR )
    {
        if (surrogateIsValid && sType == sGP)
        {
            double *x;
            x = new double[p.size()];
            for (unsigned ld = 0; ld < coordinates.size(); ld++)
                x[ld] = -1.0 + 2.0*(p[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);

            double f = currentMin + gp->f(x)*(currentMax - currentMin);

            delete[] x;

            return f;

        }
        else if (surrogateIsValid && sType == sRBFN)
        {
            vector<double> x(p.size());
            for (unsigned ld = 0; ld < coordinates.size(); ld++)
                x[ld] = -1.0 + 2.0*(p[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);

            tFitness f = currentMin + rbfn->predictValue(x)*(currentMax - currentMin);

            return f;
        }
		else if (surrogateIsValid && sType == sSVR)
		{
			sample_type x;
			x.set_size(p.size());
			for (unsigned ld = 0; ld < coordinates.size(); ld++)
				x(ld, 0) = -1.0 + 2.0*(p[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);

			tFitness f = currentMin + svr(x)*(currentMax - currentMin);

			return f;
		}
        else
        {
            isTrueFitness = true;
			cout << "Surrogate not valid: using true fitness" << endl;
            return calculateFitnessValue(p);
        }
    }
    else if (sType == sQPA)
    {
        if (archive->size() < minNumberOfPatterns)
        {
            isTrueFitness = true;

            return calculateFitnessValue(p);
        }

        vector< vector<Pattern>::iterator > sortArchive;
        for (auto it = archive->begin(); it != archive->end(); ++it)
        {
            double d = 0;
            for (unsigned ld = 0; ld < coordinates.size(); ld++)
                d += (p[ld] - it->point[ld])*(p[ld] - it->point[ld]);
            it->dist = d;
            sortArchive.push_back(it);
        }

        sort(sortArchive.begin(), sortArchive.end(), doComparePatterns());

        vector< vector<double> > points;
        vector<tFitness> values;

        currentMin = currentMax = sortArchive[0]->fitness;
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
        {
            minCoordInArchive[ld] = sortArchive[0]->point[ld];
            maxCoordInArchive[ld] = sortArchive[0]->point[ld];
        }

        for (int i = 0; i < minNumberOfPatterns; ++i)
        {
            if (currentMax < sortArchive[i]->fitness)
                currentMax = sortArchive[i]->fitness;
            if (currentMin > sortArchive[i]->fitness)
                currentMin = sortArchive[i]->fitness;

            for (unsigned ld = 0; ld < coordinates.size(); ld++)
            {
                if (minCoordInArchive[ld] > sortArchive[i]->point[ld])
                    minCoordInArchive[ld] = sortArchive[i]->point[ld];
                if (maxCoordInArchive[ld] < sortArchive[i]->point[ld])
                    maxCoordInArchive[ld] = sortArchive[i]->point[ld];
            }
        }

        for (int i = 0; i < minNumberOfPatterns; ++i)
        {
            vector<double> x(sortArchive[i]->point.size());
            for (unsigned ld = 0; ld < coordinates.size(); ld++)
            {
                x[ld] = -1.0 + 2.0*(sortArchive[i]->point[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);
            }
            points.push_back(x);

            tFitness f = (sortArchive[i]->fitness - currentMin) / (currentMax - currentMin);

            values.push_back(f);
        }

        qr->setData(points, values);

        vector<double> x(p.size());
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            x[ld] = -1.0 + 2.0*(p[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);

        startTime = clock();

        tFitness f = qr->evaluate(x);

        return currentMin + f*(currentMax - currentMin);
    }
    else
    {
        cerr << "Unknown surrogate type in calculateSurrogateFitnessValue" << endl;
        exit(1);
    }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness JADE::calculateGPSurrogatePredictionVariance(vector<double> &p)
{
    if (surrogateIsValid)
    {
        double *x;
        clock_t startTime = clock();
        x = new double[p.size()];
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            x[ld] = -1.0 + 2.0*(p[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);

        double v = currentMin + gp->var(x)*(currentMax - currentMin);

        delete[] x;

        return v;
    }
    else return 0;
}





void JADE::findSVRparameters(std::vector<sample_type> &samples, std::vector<double> &targets, double &bestGamma, double &bestC)
{
	using namespace dlib;

	std::vector<int> ind;
	for (int i = 0; i < samples.size(); ++i) ind.push_back(i);
	std::random_shuffle(ind.begin(), ind.end());

	std::vector<sample_type> cvSamples[2];
	std::vector<double> cvTargets[2];
	int i = 0;
	for (; i < 0.8*samples.size(); ++i)
	{
		cvSamples[0].push_back(samples[ind[i]]);
		cvTargets[0].push_back(targets[ind[i]]);
	}
	for (; i < samples.size(); ++i)
	{
		cvSamples[1].push_back(samples[ind[i]]);
		cvTargets[1].push_back(targets[ind[i]]);
	}


	dlib::svr_trainer<kernel_type> svrTrainer;

	std::vector<double> gammas;
	std::vector<double> cs;
	double minErr = std::numeric_limits<double>::infinity();

	gammas.push_back(0.01);
	gammas.push_back(0.1);
	//gammas.push_back(0.25);
	gammas.push_back(0.5);
	gammas.push_back(1.0);

	cs.push_back(1);
	cs.push_back(2);
	cs.push_back(5);
	cs.push_back(10);
	//cs.push_back(100);


	for (int i = 0; i < gammas.size(); ++i)
	{
		double gamma = gammas[i];
		for (int j = 0; j < cs.size(); ++j)
		{
			double c = cs[j];

			svrTrainer.set_c(c);
			svrTrainer.set_kernel(kernel_type(gamma));
			svrTrainer.set_epsilon_insensitivity(1.0E-8);

			svr = svrTrainer.train(cvSamples[0], cvTargets[0]);

			double err = 0;
			for (int q = 0; q < cvSamples[1].size(); ++q)
				err += fabs(svr(cvSamples[1][q]) - cvTargets[1][q]);

			if (err < minErr)
			{
				bestGamma = gamma;
				bestC = c;
				minErr = err;
			}
		}
	}
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::trainGlobalSurrogate()
{
	

	if (archive->size() < minNumberOfPatterns)
    {
        cerr << "archive->size() < minArchiveSize in trainGlobalSurrogate" << endl;
        exit(1);
    }

    //Avoid using the surrogate if already converged
    if (fabs(decomposer.CCOptimizer.globalBestFitness - decomposer.CCOptimizer.optimum)<1.0E-16)
    {
        surrogateIsValid = false;
        return;
    }

    clock_t startTime = clock();

    //forget the latest archive->size() - maxNumberOfPatterns patterns
    unsigned startIndex = max(0, (int)(archive->size() - maxNumberOfPatterns));

    currentMin = currentMax = (*archive)[startIndex].fitness;
    for (unsigned ld = 0; ld < coordinates.size(); ld++)
    {
        minCoordInArchive[ld] = (*archive)[startIndex].point[ld];
        maxCoordInArchive[ld] = (*archive)[startIndex].point[ld];
    }

    for (int i = startIndex; i < archive->size(); ++i)
    {
        if ( currentMax < (*archive)[i].fitness )
            currentMax = (*archive)[i].fitness;
        if ( currentMin > (*archive)[i].fitness )
            currentMin = (*archive)[i].fitness;

        for (unsigned ld = 0; ld < coordinates.size(); ld++)
        {
            if (minCoordInArchive[ld] > (*archive)[i].point[ld])
                minCoordInArchive[ld] = (*archive)[i].point[ld];
            if (maxCoordInArchive[ld] < (*archive)[i].point[ld])
                maxCoordInArchive[ld] = (*archive)[i].point[ld];
        }
    }

	currentMax += 1.0E-06;

    if ( sType == sGP )
    {
        createGP();

        gp->clear_sampleset();
        for (int i = startIndex; i < archive->size(); ++i)
        {
            double *x = new double[coordinates.size()];
            for (unsigned ld = 0; ld < coordinates.size(); ld++)
            {
                x[ld] = -1.0 + 2.0*((*archive)[i].point[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);
            }

            gp->add_pattern(x, ((*archive)[i].fitness - currentMin) / (currentMax - currentMin));

            delete[] x;
        }

        optimizeGPParameters();
    }
    else if (sType == sRBFN)
    {
        rbfn->reset();

        unsigned numPatterns = archive->size() - startIndex;
        vector<int> indexes(numPatterns);
        for (int i = 0; i < numPatterns; ++i)
            indexes[i] = i;

        random_shuffle(indexes.begin(), indexes.end());
        bool hold_out = false;
        for (int ii = 0; ii<numPatterns; ++ii)
        {
            int i = startIndex + indexes[ii];
            vector<double> x(coordinates.size());
            for (unsigned ld = 0; ld < coordinates.size(); ld++)
            {
                x[ld] = -1.0 + 2.0*((*archive)[i].point[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);
            }

            //tFitness f = -1.0 + 2.0*((*archive)[i].fitness - currentMin) / (currentMax - currentMin);
            tFitness f = ((*archive)[i].fitness - currentMin) / (currentMax - currentMin);

            if (ii < 2 * numPatterns / 3 || !hold_out)
            {
                rbfn->training_data.push_back(x);
                rbfn->training_values.push_back(f);
            }
            else
            {
                rbfn->testing_data.push_back(x);
                rbfn->testing_values.push_back(f);
            }
        }
        int num_rbf_units = rbfn->training_data.size() / 5;
        double learning_rate = 0.2;
        int num_iterations = 30;
        tFitness toll = 1.0E-09;
        tFitness mse = rbfn->training(num_rbf_units, learning_rate, num_iterations, toll);
        surrogateIsValid = true;
    }
	else if (sType == sSVR)
	{
		using namespace dlib;

		std::vector<sample_type> samples;
		std::vector<double> targets;
		for (int i = startIndex; i < archive->size(); ++i)
		{
			sample_type x;
			x.set_size(coordinates.size());
			for (unsigned ld = 0; ld < coordinates.size(); ld++)
			{
				x(ld, 0) = -1.0 + 2.0*((*archive)[i].point[ld] - minCoordInArchive[ld]) / (maxCoordInArchive[ld] - minCoordInArchive[ld]);
			}

			samples.push_back(x);
			targets.push_back(((*archive)[i].fitness - currentMin) / (currentMax - currentMin));
			//cout << (currentMax - currentMin) << endl;
		}


		// Now we are making a typedef for the kind of kernel we want to use.  I picked the
		// radial basis kernel because it only has one parameter and generally gives good
		// results without much fiddling.


		// Now setup a SVR trainer object.  It has three parameters, the kernel and
		// two parameters specific to SVR.  
		double gamma = 0.1;
		double c = 10;

		findSVRparameters(samples, targets, gamma, c);
	
		svrTrainer.set_kernel(kernel_type(gamma));

		// This parameter is the usual regularization parameter.  It determines the trade-off 
		// between trying to reduce the training error or allowing more errors but hopefully 
		// improving the generalization of the resulting function.  Larger values encourage exact 
		// fitting while smaller values of C may encourage better generalization.
		svrTrainer.set_c(c);

		// Epsilon-insensitive regression means we do regression but stop trying to fit a data 
		// point once it is "close enough" to its target value.  This parameter is the value that 
		// controls what we mean by "close enough".  In this case, I'm saying I'm happy if the
		// resulting regression function gets within 0.001 of the target value.
		svrTrainer.set_epsilon_insensitivity(1.0E-8);

		// Now do the training and save the results
		svr = svrTrainer.train(samples, targets);
		surrogateIsValid = true;
	}
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::optimize(int iterations)
{
    if ( qr )
        qr->reset();

    for (ite = 0; ite < iterations; ++ite)
    {
        if ( decomposer.CCOptimizer.numberOfEvaluations + this->nfe >= decomposer.CCOptimizer.maxNumberOfEvaluations )
            break;

        if (sType == sNone)
            update();
        else if (sType == sGP || sType == sQPA || sType==sRBFN || sType==sSVR )
        {
            SAUpdate(sType);
        }
        else
        {
            cerr << "Unknown type of surrogate" << endl;
            exit(1);
        }
    }
    updateIndexOfBest();
}
