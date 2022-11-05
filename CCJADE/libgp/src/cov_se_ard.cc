// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include "cov_se_ard.h"
#include <cmath>

namespace libgp
{
  
  CovSEard::CovSEard() {}
  
  CovSEard::~CovSEard() {}
  
  bool CovSEard::init(int n)
  {
    input_dim = n;
    param_dim = n+3;
    ell.resize(input_dim);
    loghyper.resize(param_dim);
    return true;
  }
  
  double CovSEard::get(const Eigen::VectorXd &x1, const Eigen::VectorXd &x2)
  {  
    double z = (x1-x2).cwiseQuotient(ell).squaredNorm();
	double noise = 0;
	if (&x1 == &x2) noise = loghyper(input_dim + 2);
	return sf2*exp(-0.5*z) + loghyper(input_dim + 1) + noise;
  }
  
  void CovSEard::grad(const Eigen::VectorXd &x1, const Eigen::VectorXd &x2, Eigen::VectorXd &grad)
  {
	double z = (x1 - x2).cwiseQuotient(ell).squaredNorm();
    double k = exp(-0.5*z);

	for (int i = 0; i < input_dim; ++i)
		grad(i) = sf2 * k * (x1(i) - x2(i))*(x1(i) - x2(i)) / (loghyper(i)*loghyper(i));  
	grad(input_dim) = k*sf2;
	grad(input_dim + 1) = loghyper(input_dim + 1);
	if (&x1 == &x2) 
	  grad(input_dim + 2) = loghyper(input_dim + 2);
	else 
	  grad(input_dim + 2) = 0.0;		
  }

  Eigen::VectorXd CovSEard::get_characteristic_length()
  {
	  Eigen::VectorXd l;
	  l.resize(input_dim);
	  for (int i = 0; i < input_dim; ++i)
		  l(i) = ell(i)*ell(i);
	  return l;
  }
  
  void CovSEard::set_loghyper(const Eigen::VectorXd &p)
  {	  
	CovarianceFunction::set_loghyper(p);	
	for(size_t i = 0; i < input_dim; ++i) ell(i) = loghyper(i);
	sf2 = loghyper(input_dim);
  }
  
  std::string CovSEard::to_string()
  {
    return "CovSEard";
  }
}

