#include <RcppEigen.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "dptf.h"
#include "tf_dp.h"

// void tf_dp (int n, double *y, double lam, double *beta);

// [[Rcpp::export]]
Eigen::VectorXd dptf(Eigen::VectorXd y, double lam) {
  Rcpp::NumericVector ny = evec_to_nvec(y);
  int n = y.size();
  Rcpp::NumericVector beta(n);
  tf_dp(n, ny.begin(), lam, beta.begin());  // get pointers to the vectors
  Eigen::VectorXd ebeta = nvec_to_evec(beta);
  return ebeta;
}

