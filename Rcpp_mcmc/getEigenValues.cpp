#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;               	// 'maps' rather than copies 
using Eigen::MatrixXcd;                  // variable size matrix, double precision
using Eigen::VectorXcd;                  // variable size vector, double precision
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::ComplexEigenSolver;    // one of the eigenvalue solvers
using Eigen::EigenSolver;    // one of the eigenvalue solvers

// [[Rcpp::export]]
VectorXcd getEigenValues(Map<MatrixXd> M) {
  ComplexEigenSolver<MatrixXcd> es(M);
  return es.eigenvalues();
}
// [[Rcpp::export]]
MatrixXcd getEigenVectors(Map<MatrixXd> M) {
  ComplexEigenSolver<MatrixXcd> es(M);
  return es.eigenvectors();
}