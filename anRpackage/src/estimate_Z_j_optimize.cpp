#include <iostream>
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#define THREAD_NUM 16

using namespace std;
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
int mp(){
  
  omp_set_num_threads(THREAD_NUM);
  
  #pragma omp parallel 
  {
    usleep(5000*omp_get_thread_num());
    cout << "Number of available threads: " << omp_get_num_threads() << endl;
    cout << "Currnet threand number: " << omp_get_thread_num() << endl;
    cout << "Hello world" << endl;
  }
  return 0;
}


// [[Rcpp::export]]
NumericMatrix my_matrix(int I, int J, int nthreads) {
  NumericMatrix A(I,J);
  // create a thread safe accessor for A
  RcppParallel::RMatrix<double> a(A);
  int tid;
  omp_set_num_threads(nthreads);
#pragma omp parallel for private(tid)
  for(int j = 0; j < J; j++) {
    for(int i = 0; i < I; i++) {
      tid = omp_get_thread_num();
      a(i, j) = tid ;
    }
  }
  
  return A;
}

// // [[Rcpp::export]]
// MatrixXf estimate_Z_j_optimize(const MatrixXf& W, const VectorXf& mus_hat_j, VectorXf& sigmas_hat_j, 
//                                const float& tau_hat, const MatrixXf& C1,  VectorXf& gammas_hat_j, 
//                                const List W_prime, const VectorXf& C2_prime_j, bool scale){
//   
//   
//   
//   return 0;
// }
