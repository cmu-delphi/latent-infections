#include <RcppEigen.h>
#include <Eigen/Sparse>
#include "admm.h"
#include "utils.h"
#include "dptf.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qr;

using namespace Rcpp;

// [[Rcpp::export]]
List estim_path_single(Eigen::VectorXd y,
                       NumericVector x,
                       Eigen::SparseMatrix<double> Cmat,
                       int korder,
                       NumericVector lambda,
                       double lambdamax = -1,
                       double lambdamin = -1,
                       int nsol = 100,
                       double rho = -1,
                       int maxadmm_iter = 1e3,
                       double tolerance = 1e-3,
                       double lambda_min_ratio = 1e-8,
                       int verbose = 0) {
  int n = Cmat.cols();
  
  if (korder < 1) stop("korder must be at least 1.");
  
  // Placeholders for solutions
  Eigen::MatrixXd thetas(n, nsol);
  NumericVector niter(nsol);
  
  // Build D matrices as needed
  Eigen::SparseMatrix<double> D;
  Eigen::SparseMatrix<double> Dk;
  Eigen::SparseMatrix<double> DkDk;
  D = get_D(korder, x);
  qr.compute(D.transpose());
  
  Dk = get_Dtil(korder, x);
  DkDk = Dk.transpose() * Dk;
  int mm = Dk.rows();
  
  
  // Generate lambda sequence if necessary
  if (abs(lambda[nsol - 1]) < tolerance / 100 && lambdamax <= 0) {
    VectorXd b(n - korder);
    VectorXd Cy = Cmat.transpose() * y;
    b = qr.solve(Cy);
    NumericVector bp = evec_to_nvec(b);
    lambdamax = max(abs(bp));
  }
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);
  
  // ADMM parameters
  double _rho;
  
  // ADMM variables
  Eigen::VectorXd beta(n);
  Eigen::VectorXd alpha(mm);
  Eigen::VectorXd u(mm);
  beta.setZero();
  alpha.setZero();
  u.setZero();
  
  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {
    if (verbose > 0) Rcout << ".";
    Rcpp::checkUserInterrupt();
    
    // rho for uneven spacing
    double x1 = x[0];
    double xn = x[x.size() - 1];
                            
    _rho = (rho < 0) ? lambda[i] * pow((xn - x1) / xn, korder) : rho;
    double lamz = lambda[i] / _rho;
    admm_gauss(maxadmm_iter, korder, y, x, Cmat, beta, alpha, u, _rho,
               lamz, DkDk, tolerance);
    
    // Store solution
    thetas.col(i) = beta;
    
    // Verbose handlers
    if (verbose > 1) Rcout << niter[i];
    if (verbose > 2) Rcout << "(" << lambda[i] << ")";
    if (verbose > 0) Rcout << std::endl;
  }
  
  // Return
  List out = List::create(
    Named("thetas") = thetas,
    Named("lambda") = lambda,
    Named("korder") = korder
  );
  
  return out;
}
