#ifndef __UTILS_H
#define __UTILS_H

Rcpp::NumericVector evec_to_nvec(Eigen::VectorXd evec);
Eigen::VectorXd nvec_to_evec(Rcpp::NumericVector nvec);
Eigen::SparseMatrix<double> get_Dtil(int k, Rcpp::NumericVector xd);
Eigen::SparseMatrix<double> get_D(int k, Rcpp::NumericVector xd);
Eigen::VectorXd doDv(Eigen::VectorXd v, int k, Rcpp::NumericVector xd);
Eigen::VectorXd doDtv(Eigen::VectorXd v, int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDtv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
void create_lambda(Rcpp::NumericVector& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol);

#endif
