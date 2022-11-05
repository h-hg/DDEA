#ifndef _UTIL_H_
#define _UTIL_H_


#include <vector>
#include <assert.h>
#include "funsoft.h"

class Utility
{

public:

	/*
	Multiply equal sized vectors and return the result in double type
	*/
	static tFitness multiplyVectors(const std::vector<tFitness> &a, const std::vector<tFitness> &b)
	{
		assert(a.size() == b.size());
		tFitness res = 0;
		for(int i = 0 ; i<a.size() ; i++)
			res += a[i] * b[i];
		return res;
	}

	/*
		Multiplies a vector with a given constant and return the result vector
	*/
	static std::vector<tFitness> multiplyVecConst(const std::vector<tFitness>&input, const tFitness c)
	{
		std::vector<tFitness>res(input.size());
		for(int i = 0 ; i<input.size() ; i++)
				res[i] = input[i] * c;
		return res;
	}

	/*
		adds a vector to a given input (passed by reference)
	*/
	static void AddVectors(std::vector<tFitness>&input, const std::vector<tFitness>&addend)
	{
		assert(input.size() == addend.size());
		for(int i = 0 ; i<input.size() ; i++)
				input[i] += addend[i];
	}

	/*
		calculates the covariance between two vectors
	*/
	static tFitness coVariance(datapoint x, datapoint y)
{
	assert(x.size() == y.size());
	tFitness a = 0, b = 0, c = 0;
	for (int i = 0; i < x.size(); ++i)
	{
		a += x[i];
		b += y[i];
		c += x[i] * y[i];
	}

	return (tFitness)(c / x.size()) - (tFitness)((a / x.size()) * (b / x.size()));
}

};

#endif