//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Decomposer.h
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
#include "JADE.h"
#include <vector>
#include <queue>
#include <map>
#include <random>
#include <functional>
#include <numeric>

class CCDE;

using namespace std;

typedef mt19937 RandomEngine;

class Decomposer
{
public:
    vector<unsigned> coordinates;
    CCDE &CCOptimizer;
    vector< JADE* > optimizers;
    int sizeOfSubcomponents;
    bool applyRandomGrouping;
    tFitness bestAchievedFitness;
    RandomEngine eng;
    uniform_real_distribution<double> unifRandom;
    unsigned individualsPerSubcomponent;
    unsigned numberOfSubcomponents;

    //Current population
    vector< vector<double> > population;

    //Final global best position and context vector
    vector<double> contextVector;

    tFitness expectedOptimum;

    //Fitnesses of population
    vector< tFitness > fitnessValues;

    typeOfSurrogate sType;

    vector<unsigned> sizes;
    vector<unsigned> baseCoordIndex;

    Decomposer(CCDE &_CCOptimizer, unsigned seed, vector<unsigned> &_coordinates,
               unsigned _sizeOfSubcomponents,
               unsigned _individualsPerSubcomponent,
               vector< vector<double> > &_population,
               vector<double>  &_contextVector,
               tFitness _expectedOptimum,
               bool RG,
               typeOfSurrogate sType, bool allocateOptimizers=true);
    ~Decomposer();
    vector< JADE* >  allocateOptimizers(vector<unsigned> &indexes);
    JADE*  allocateOptimizer();
    void setPopulation(vector< vector<double> > &_population);
    void setSubcomponentsOfEqualSize(unsigned newSizeOfSubcomponents, vector<tFitness> &fitnessValues, typeOfSurrogate sType);
    void setCoordinates(vector<unsigned> &_coordinates);
    void updateContextVector(JADE *optimizer);
    void buildContextVector();
    void randomGrouping();
    void setSeed(unsigned seed);
    void setOptimizersCoordinatesAndEvaluatePopulation();
    void setOptimizersCoordinatesAndEvaluatePopulation(vector<unsigned> &indexes);
    void setOptimizersCoordinates(vector<unsigned> &indexes);
};
