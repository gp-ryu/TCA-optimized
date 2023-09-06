// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>

using namespace std;
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> MapMatd;
typedef Map<VectorXd> MapVecd;

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

// [[Rcpp::export]]
SEXP AAt(const MapMatd& A){
  int n = A.rows();
  int m = A.cols();
  MatrixXd AAt = MatrixXd(n,n).setZero().selfadjointView<Lower>()
                              .rankUpdate(A);
  return wrap(AAt);
}

// [[Rcpp::export]]
SEXP AtA(const MapMatd& A){
  int n = A.rows();
  int m = A.cols();
  MatrixXd AtA = MatrixXd(m,m).setZero().selfadjointView<Lower>()
                              .rankUpdate(A.adjoint());
  return wrap(AtA);
}