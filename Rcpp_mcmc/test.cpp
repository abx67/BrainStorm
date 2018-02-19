#include <Rcpp.h>
//#include <RcppEigen.h>
//#include <Eigen/Eigenvalues> 
using namespace Rcpp;
//using namespace RcppEigen;


NumericMatrix mmult(NumericMatrix m1, NumericMatrix m2){
  if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);
    }
  }
  return out;
};

int my_rmultinomF(NumericVector a){
  Rcpp::Environment global = Rcpp::Environment::global_env();
  Function rmultinomF = global["rmultinomF"];
  NumericVector x = rmultinomF(a);
  return x[0];
};

NumericMatrix transpose(NumericMatrix X){
  NumericMatrix ans (X.ncol(), X.nrow()); 
  for (int i = 0; i < X.nrow(); i++){ 
    for (int j = 0; j < X.ncol(); j++){ 
      ans(j,i) = X(i, j); 
    } 
  } 
  return(ans);
};

NumericVector first_abs_eigenvec(NumericMatrix a){
  Rcpp::Environment global = Rcpp::Environment::global_env();
  Function getEigenVectors = global["getEigenVectors"];
  ComplexMatrix EigenVec = getEigenVectors(a);
  return Mod(EigenVec(_,0));
};

// [[Rcpp::export]]
SEXP bb(NumericMatrix mat) {
  Rcpp::Environment global = Rcpp::Environment::global_env();
  Function getEigenVectors = global["getEigenVectors"];
  Function getEigenValues = global["getEigenValues"];
  int k =5;
  NumericMatrix x(k,k);
  std::fill(x.begin(),x.end(),2);
  NumericVector y(k);
  std::fill(y.begin(),y.end(),0.5);
  NumericMatrix m(y.size() , 1 , y.begin() );
  //y.attr("dim") = Dimension(k, 1);
  Range rowind(0, k-2);
  Range ran(0,k-1);
  //for (int i=0;i<k;i++){x(_,i)=exp(x(_,i));}
  NumericMatrix a(1,k-1);
  NumericVector n=wrap(m(rowind,_));
  //x(_,1)=n;
  //SelfAdjointEigenSolver<MatrixXd> es(x);
  //MatrixXd Y = RcppEigen::as<MatrixXd>(x);
  //ComplexMatrix result = getEigenValues(mat);
  //NumericMatrix res=wrap(Re(result));
  //return wrap(first_real_eigenvec(mat));
  return wrap(first_abs_eigenvec(mat));
}


