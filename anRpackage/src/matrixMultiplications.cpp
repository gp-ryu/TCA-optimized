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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
