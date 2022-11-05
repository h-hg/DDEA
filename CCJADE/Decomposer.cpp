//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Decomposer.cpp
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


#include "Decomposer.h"
#include "CCDE.h"
#include "sobol.hpp"
#include <random>


//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer::Decomposer(CCDE &_CCOptimizer,
                       unsigned seed,
                       vector<unsigned> &_coordinates,
                       unsigned _sizeOfSubcomponents,
                       unsigned _individualsPerSubcomponent,
                       vector< vector<double> > &_population,
                       vector<double>  &_contextVector,
                       tFitness _expectedOptimum,
                       bool RG,
                       typeOfSurrogate sType,
                       bool allocateOptimizers) : CCOptimizer(_CCOptimizer), sizeOfSubcomponents(_sizeOfSubcomponents),
    individualsPerSubcomponent(_individualsPerSubcomponent), applyRandomGrouping(RG), sType(sType)
{
    eng.seed(seed);
    expectedOptimum = _expectedOptimum;
    bestAchievedFitness = std::numeric_limits<tFitness>::infinity();
    coordinates = _coordinates;

    for (unsigned i = 0; i < individualsPerSubcomponent; ++i)
        population.push_back(_population[i]);

    unsigned d = 0, size = sizeOfSubcomponents;
    while ( d<coordinates.size() )
    {
        if ( d + size > coordinates.size() )
            size = coordinates.size() - d;

        unsigned *optCoord = &(coordinates[d]);

        baseCoordIndex.push_back(d);
        sizes.push_back(size);

        int s = tau_sobol(size);
        double *x = i8_sobol_generate(size, individualsPerSubcomponent, s);
        for (int j = 0; j < individualsPerSubcomponent; j++)
        {
            for (int k = 0; k < size; ++k)
            {
                population[j][optCoord[k]] = CCOptimizer.lowerLimit + x[size*j + k] * (CCOptimizer.upperLimit - CCOptimizer.lowerLimit);
                //population[j][optCoord[k]] = CCOptimizer.lowerLimit + unifRandom(eng) * (CCOptimizer.upperLimit - CCOptimizer.lowerLimit);
            }
        }
        delete[] x;

        if (allocateOptimizers)
        {
            JADE *optimizer = new JADE(size, individualsPerSubcomponent, *this, sType);
            optimizer->setCoordinates(optCoord, size);
            optimizers.push_back(optimizer);
        }

        d += size;
    }
    numberOfSubcomponents = sizes.size();

    optimizers.resize(numberOfSubcomponents);

    contextVector.resize(CCOptimizer.problemDimension);
    _contextVector.resize(CCOptimizer.problemDimension);
    for (unsigned i = 0; i < CCOptimizer.problemDimension; ++i)
        contextVector[i] = _contextVector[i] = population[unifRandom(eng)*population.size()][i];
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
vector< JADE* >  Decomposer::allocateOptimizers(vector<unsigned> &indexes)
{
    vector< JADE* > bunchOfOptimizers;

    for (int i = 0; i < indexes.size(); ++i)
    {
        unsigned j = indexes[i];
        JADE *optimizer = new JADE(sizes[j], individualsPerSubcomponent, *this, sType);
        optimizer->setCoordinates(&coordinates[baseCoordIndex[j]], sizes[j]);
        optimizers[j] = optimizer;
        bunchOfOptimizers.push_back(optimizer);
        optimizer->loadIndividuals(population);
    }
    return bunchOfOptimizers;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
JADE* Decomposer::allocateOptimizer()
{
    JADE *optimizer = new JADE(sizes[0], individualsPerSubcomponent, *this, sType);
    optimizer->setCoordinates(&coordinates[baseCoordIndex[0]], sizes[0]);
    optimizer->loadIndividuals(population);
    return optimizer;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void  Decomposer::setOptimizersCoordinatesAndEvaluatePopulation(vector<unsigned> &indexes)
{
    for (int i = 0; i < indexes.size(); ++i)
    {
        unsigned j = indexes[i];
        JADE *optimizer = optimizers[j];
        optimizer->setCoordinates(&coordinates[baseCoordIndex[j]], sizes[j]);
        optimizer->updateIndividuals(population);
        optimizer->evaluateParents();
        optimizer->updateIndexOfBest();
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void  Decomposer::setOptimizersCoordinates(vector<unsigned> &indexes)
{
    for (int i = 0; i < indexes.size(); ++i)
    {
        unsigned j = indexes[i];
        JADE *optimizer = optimizers[j];
        optimizer->setCoordinates(&coordinates[baseCoordIndex[j]], sizes[j]);
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setPopulation(vector< vector<double> > &_population)
{
    population.clear();
    for (unsigned i = 0; i < individualsPerSubcomponent; ++i)
        population.push_back(_population[i]);
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setOptimizersCoordinatesAndEvaluatePopulation()
{
    unsigned d = 0, k = 0, size = sizeOfSubcomponents;
    while ( d<coordinates.size() )
    {
        if (d + size > coordinates.size())
            size = coordinates.size() - d;

        JADE *optimizer = optimizers[k];

        optimizer->setCoordinates(&(coordinates[d]), size);

        optimizer->loadIndividuals(population);
        optimizer->evaluateParents();
        optimizer->updateIndexOfBest();

        d += size;
        k++;
    }

}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer::~Decomposer()
{
    for (unsigned i = 0; i < optimizers.size(); ++i)
        delete optimizers[i];
};



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setSeed(unsigned seed)
{
    eng.seed(seed);
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setSubcomponentsOfEqualSize(unsigned newSizeOfSubcomponents, vector<tFitness> &fitnessValues, typeOfSurrogate sType)
{
    for (unsigned i = 0; i < optimizers.size(); ++i)
        delete optimizers[i];
    optimizers.clear();

    sizeOfSubcomponents = newSizeOfSubcomponents;
    unsigned numberOfSubcomponents = coordinates.size() / newSizeOfSubcomponents;

    individualsPerSubcomponent = fitnessValues.size();

    for (unsigned i = 0; i<numberOfSubcomponents; ++i)
    {
        optimizers.push_back(new JADE(sizeOfSubcomponents, individualsPerSubcomponent, *this, sType));

        optimizers[i]->setCoordinates(&(coordinates[i*sizeOfSubcomponents]), sizeOfSubcomponents);

        optimizers[i]->loadIndividuals(population);
        optimizers[i]->setParentFitness(fitnessValues);

        //optimizers[i]->evaluateParents();

        // optimizers[i]->updateIndexOfBest();
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setCoordinates(vector<unsigned> &_coordinates)
{
    coordinates = coordinates;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::updateContextVector(JADE *optimizer)
{
    vector<double> v = optimizer->getCollaborator();
    tFitness newBestCandidate = optimizer->calculateFitnessValue(v);
    CCOptimizer.numberOfEvaluations++;
    if ( newBestCandidate < bestAchievedFitness )
    {
        for (unsigned ld = 0; ld<v.size(); ld++)
            contextVector[optimizer->coordinates[ld]] = v[ld];

        bestAchievedFitness = newBestCandidate;
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::buildContextVector()
{
    for (unsigned j = 0; j<optimizers.size(); ++j)
    {
        optimizers[j]->nfe = 0;
        vector<double> v = optimizers[j]->getCollaborator();
        tFitness newBestCandidate = optimizers[j]->calculateFitnessValue(v);
        if ( newBestCandidate<bestAchievedFitness )
        {
            for (unsigned ld = 0; ld<v.size(); ld++)
                contextVector[optimizers[j]->coordinates[ld]] = v[ld];
            bestAchievedFitness = newBestCandidate;
        }
    }
    bestAchievedFitness = CCOptimizer.computeFitnessValue(contextVector);
    CCOptimizer.numberOfEvaluations++;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::randomGrouping()
{
    if ( optimizers.size() && this->applyRandomGrouping )
    {
        shuffle(coordinates.begin(), coordinates.end(), eng);
        //setOptimizersCoordinatesAndEvaluatePopulation();

        unsigned numOfCoordinatesPerSubgroup = coordinates.size() / optimizers.size();
        for (unsigned i = 0; i < optimizers.size(); ++i)
        {
            optimizers[i]->setCoordinates(&(coordinates[i*numOfCoordinatesPerSubgroup]), numOfCoordinatesPerSubgroup);
            optimizers[i]->loadIndividuals(population);
            optimizers[i]->evaluateParents();
            CCOptimizer.numberOfEvaluations += optimizers[i]->nfe;
            optimizers[i]->nfe = 0;
        }
    }
}
