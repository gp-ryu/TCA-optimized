// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace std;
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, 
                  Eigen::MatrixXd B, 
                  int n_cores){
  
  Eigen::setNbThreads(n_cores);
  //qDebug()  << Eigen::nbThreads( );
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B, 
                      int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenT(MatrixXd A){
  MatrixXd B = A.transpose().eval();
  
  return Rcpp::wrap(B);
}

// [[Rcpp::export]]
SEXP eigenW_norms(MatrixXd W, MatrixXd sigmas_hat, double tau_hat, int n_cores){
  setNbThreads(n_cores);
  MatrixXd C = W * sigmas_hat;
  MatrixXd D = (C.array()  + pow(tau_hat,2)).matrix();
  return wrap(D.cwiseSqrt());
}

// [[Rcpp::export]]
SEXP eigenArrayProduct(ArrayXXd A, ArrayXXd B){
    return wrap(A * B);
}

// [[Rcpp::export]]
SEXP eigenProduct(MatrixXd A, MatrixXd B){
  return wrap(A*B);
}

// [[Rcpp::export]]
SEXP eigenSqrt(MatrixXd A){
  return wrap(A.cwiseSqrt());
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
