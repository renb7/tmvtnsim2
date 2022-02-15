#include "RcppArmadillo.h"

#ifndef __UTILITIES__
#define __UTILITIES__

double norm_rej(const double a, const double b);
double unif_rej(const double a, const double b);
double halfnorm_rej(const double a, const double b);
double exp_rej(const double a, const double b);
arma::colvec rtnormcpp(const arma::colvec& mean, 
                       const double sd, 
                       const arma::colvec& lower, 
                       const arma::colvec& upper);

#endif // __UTILITIES__
