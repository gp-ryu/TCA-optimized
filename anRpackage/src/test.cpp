#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#define THREAD_NUM 16

using namespace std;
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXi> MapMati;
typedef Map<MatrixXd> MapMatd;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
SEXP tcp(MatrixXd& A, MatrixXd& B) {
  MatrixXd res = A * B.adjoint();
  return wrap(res);
}

// [[Rcpp::export]]
MatrixXd logsum(MatrixXd& A, MatrixXd& B){
  return (A.transpose() * B).diagonal();
}

// [[Rcpp::export]]
SEXP sumlog(MatrixXd& A){
  return wrap(A.array().log().sum());
}

// [[Rcpp::export]]
SEXP sumlog_p(MatrixXd& A){
  int m = A.cols();
  VectorXd res(m);
#pragma omp parallel for 
  for(int i = 0; i < m; i++){
    res(i) = A.col(i).array().log().sum();
  }
  return wrap(res.sum());
}

// [[Rcpp::export]]
SEXP sum_p(MatrixXd& A){
  int n = A.rows();
  int m = A.cols();
  int mass = (n>m)?n:m;
  int rowcol = (n>m)?0:1;
  VectorXd res(mass);
  
#pragma omp parallel for 
  for(int i=0; i < mass; i++){
    if(rowcol == 0){
      res(i) = A.row(i).sum();
    }else{
      res(i) = A.col(i).sum();
    }
  }
  return wrap(res.sum());
}



//[[Rcpp::export]]
double test1(MatrixXd& A){
  VectorXd B = A.colwise().lpNorm<2>();
  cout << B << endl;
  cout << B.inverse() << endl;
  return 1;
  // return wrap(A.array().colwise() * B.inverse().array());
  // return wrap(A.array() * B.cwiseInverse().array()).matrix();
}

//[[Rcpp::export]]
SEXP AtV(MatrixXd& A, VectorXd& B){
  VectorXd res;
  res = A.adjoint() * B;
  return wrap(res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
