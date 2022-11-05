//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : main.cpp
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

#define _USE_MATH_DEFINES
#include "CCDE.h"
#include "tclap/CmdLine.h"

using namespace TCLAP;

using namespace std;

FunctionCallback functions[] = { Shifted_Sphere, Schwefel_Problem, Shifted_Rosenbrock, Shifted_Rastrigin,
                                 Shifted_Griewank, Shifted_Ackley, f_Schwefel2_22, f_Schwefel1_2, Extended_f_10, f_Bohachevsky, f_Schaffer, f_Hybrid_12,
                                 f_Hybrid_13, f_Hybrid_14, f_Hybrid_15, f_Hybrid_16new, f_Hybrid_17new, f_Hybrid_18new, f_Hybrid_19new
                               };

char* functionsNames[] = { "Shifted_Sphere", "Schwefel_Problem", "Shifted_Rosenbrock", "Shifted_Rastrigin",
                           "Shifted_Griewank", "Shifted_Ackley", "Schwefel2_22", "Schwefel1_2", "Extended_f_10", "Bohachevsky", "Schaffer", "Hybrid_12",
                           "Hybrid_13", "Hybrid_14", "Hybrid_15", "Hybrid_16", "Hybrid_17", "Hybrid_18", "Hybrid_19"
                         };

double functionsDomainBounds[] = {  100.0,  100.0, 100.0,    5.0,  600.0,   32.0, 10.0, 65.536, 100.0, 15.0, 100.0, 100.0, 100.0, 5.0, 10.0, 100.0, 100.0, 5.0, 10.0};
tFitness  functionsOptima[] = { -450.0, -450.0, 390.0, -330.0, -180.0, -140.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void optimization(int argc, char* argv[])
{   
	unsigned int functionIndex;
	typeOfSurrogate sType; //allowed: {sNone, sGP, sQPA, sRBFN, sSVR}
	unsigned int numRep;
	unsigned int numItePerCycle;
	unsigned int problemDimension;
	unsigned int sizeOfSubcomponents;
	unsigned int numOfIndividuals;
	unsigned int numberOfEvaluations;
    
	vector<int> seeds;
    unsigned maxNumRep = 100;
    for (unsigned i = 0; i < maxNumRep; ++i)
        seeds.push_back(i);

	//parsing SACCJADE command line
	try
	{
		CmdLine cmd("SACCJADE - Surrogate-assisted Cooperative Coevolution for Large-Scale Optimization of Computationally Expensive Objective Functions", ' ', "1.0");

		ValueArg<unsigned int> functionArg("f", "function", "function to optimize [1-19]", false, 1, "int");
		cmd.add(functionArg);

		ValueArg<unsigned int> surrogateArg("m", "metamodel", "type of fitness metamodel [0->none; 1->GP; 2->QPA; 3->RBFN; 4->SVR]", false, 4, "int");
		cmd.add(surrogateArg);		

		ValueArg<unsigned int> repArg("r", "repetitions", "number of independent repetitions [1-100]", false, 1, "int");
		cmd.add(repArg);

		ValueArg<unsigned int> iteArg("i", "iterations", "number of JADE iterations per cycle", false, 6, "int");
		cmd.add(iteArg);		

		ValueArg<unsigned int> dimArg("d", "dimension", "problem dimension", false, 1000, "int");
		cmd.add(dimArg);		

		ValueArg<unsigned int> sdimArg("s", "subdim", "size of subcomponents", false, 4, "int");
		cmd.add(sdimArg);		

		ValueArg<unsigned int> npopArg("p", "popsize", "number of individuals in each subcomponent", false, 25, "int");
		cmd.add(npopArg);		

		ValueArg<unsigned int> feArg("e", "fevals", "allowed number of exact fitness evaluations", false, 500000, "int");
		cmd.add(feArg);		

		cmd.parse(argc, argv);
		
		functionIndex = functionArg.getValue();
		sType = (typeOfSurrogate)surrogateArg.getValue();
		numRep = repArg.getValue();
		numItePerCycle = iteArg.getValue();
		problemDimension = dimArg.getValue();
		sizeOfSubcomponents = sdimArg.getValue();
		numOfIndividuals = npopArg.getValue();
		numberOfEvaluations = feArg.getValue();

	}
    catch (ArgException& e)
    {
	   cout << "ERROR: " << e.error() << " " << e.argId() << endl;
    }

		
	if (sType == sNone)
        cout << "Using CCJADE" << endl;
    else if (sType == sGP)
        cout << "SACCJADE with Gausiann Process" << endl;
    else if (sType == sQPA)
        cout << "SACCJADE with Quadratic Polynomial Local Approsimation" << endl;
    else if (sType == sRBFN)
        cout << "SACCJADE with Radial Basis Function Network" << endl;
	else if (sType == sSVR)
		cout << "SACCJADE with Support Vector Regression" << endl;
    else
    {
        cerr << "unknown surrogate" << endl;
        exit(1);
    }

	if ( functionIndex < 1 || functionIndex>19 )
	{
		cerr << "function index out of allowed bounds [1..19]" << endl;
		exit(1);
	}
	functionIndex--;

	if (problemDimension<1 || problemDimension>1000)
	{
		cerr << "problem dimension must be in [1..1000]" << endl;
		exit(1);
	}

	if (sizeOfSubcomponents<1 || sizeOfSubcomponents>problemDimension)
	{
		cerr << "problem dimension must be in [1..problem dimension]" << endl;
		exit(1);
	}

    cout << "Optimizing f" << functionIndex+1 << " (" << functionsNames[functionIndex] << ")" << endl;
	cout << "Problem dimension = " << problemDimension << endl;
    cout << "Number of iterations per cycle = " << numItePerCycle << endl;    	
	cout << "The problem is decomposed in " << problemDimension / sizeOfSubcomponents << " equal subcomponents of size " << sizeOfSubcomponents << endl;
    cout << "Number of individuals per subcomponent = " << numOfIndividuals << endl;
	cout << "Number of repetitions = " << numRep << endl;
	cout << "Allowed number of exact function evaluations = " << numberOfEvaluations << endl;

    double time = 0;
    vector< vector<ConvPlotPoint> > convergences;
    for (unsigned k = 0; k < numRep; ++k)
    {
        vector<ConvPlotPoint> convergence;
        CCDE ccde;
        int seed = seeds[k];
        ccde.optimize(functions[functionIndex], problemDimension, functionsDomainBounds[functionIndex],
                      functionsOptima[functionIndex], numberOfEvaluations, sizeOfSubcomponents, numOfIndividuals,
                      convergence, seed, sType, numItePerCycle);
        time = ccde.elapsedTime;
        convergences.push_back(convergence);
    }

    char fName[256];
    FILE *file;

    if (convergences.size() == 0) return;
    if (convergences[0].size() == 0) return;

    if (sType == sNone)
		sprintf_s(fName, "convplot_f%d_dec%d_popsize%d_CCJADE.csv", functionIndex + 1, sizeOfSubcomponents, numOfIndividuals);
    else if (sType == sGP)
		sprintf_s(fName, "convplot_f%d_dec%d_popsize%d_GP_SACCJADE.csv", functionIndex + 1, sizeOfSubcomponents, numOfIndividuals);
    else if (sType == sQPA)
		sprintf_s(fName, "convplot_f%d_dec%d_popsize%d_QPA_SACCJADE.csv", functionIndex + 1, sizeOfSubcomponents, numOfIndividuals);
    else if (sType == sRBFN)
		sprintf_s(fName, "convplot_f%d_dec%d_popsize%d_RBFN_SACCJADE.csv", functionIndex + 1, sizeOfSubcomponents, numOfIndividuals);
	else if (sType == sSVR)
		sprintf_s(fName, "convplot_f%d_dec%d_popsize%d_RBFN_SACCJADE.csv", functionIndex + 1, sizeOfSubcomponents, numOfIndividuals);

    fopen_s(&file, fName, "wt");
    vector<ConvPlotPoint> averageConvergence;

    int maxSize = 0;
    int idOfMaxSize = 0;
    for (unsigned q = 0; q < convergences.size(); ++q)
        if (convergences[q].size() > maxSize)
        {
            maxSize = convergences[q].size();
            idOfMaxSize = q;
        }

    for (unsigned q = 0; q < maxSize; ++q)
    {
        for (unsigned k = 0; k < numRep; ++k)
        {
            if (q < convergences[k].size())
            {
                fprintf(file, "%d; %.8Le; ", convergences[k][q].nfe, convergences[k][q].f);
            }
            else
            {
                fprintf(file, "; ;");
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
int main(int argc, char* argv[])
{
    optimization(argc, argv);
}

