#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>


#define THREAD_NUM 16
using namespace std;
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> MapMatd;
typedef Map<VectorXd> MapVecd;


// [[Rcpp::export]]
MatrixXd createDiagonalMatrix(const VectorXd& diagonalValues) {
  DiagonalMatrix<double, Dynamic> diagMatrix(diagonalValues.size());
  diagMatrix.diagonal() = diagonalValues;
  return diagMatrix;
}

// [[Rcpp::export]]
MatrixXd outerprodV(const VectorXd & A) {
  int n = A.size();
  MatrixXd AtA =  MatrixXd(n,n).setZero().selfadjointView<Lower>()
                               .rankUpdate(A.adjoint());
  return AtA;
}
// [[Rcpp::export]]
MatrixXd outerprodM(const MapMatd& A){
  int n = A.cols();
  MatrixXd AtA = MatrixXd(n,n).setZero().selfadjointView<Lower>()
                              .rankUpdate(A.adjoint());
  return AtA;
}

// [[Rcpp::export]]
MatrixXd outerprodMM(MatrixXd& A, MatrixXd& B) {
  MatrixXd res = A.adjoint() * B;
  return res;
}
// [[Rcpp::export]]
MatrixXd outerprodMM_p(MatrixXd& A, MatrixXd& B) {
  initParallel();
  MatrixXd res = A.adjoint() * B;
  return res;
}
// [[Rcpp::export]]
MatrixXd outerprodMV_p(MatrixXd& A, VectorXd& B) {
  initParallel();
  MatrixXd res = A.adjoint() * B;
  return res;
}

// [[Rcpp::export]]
MatrixXd innerprodMM(MatrixXd& A, MatrixXd& B) {
  MatrixXd res = A * B.adjoint();
  return res;
}





// [[Rcpp::export]]
SEXP estimate_Z_j_cpp(const MapMatd& X, const MapMatd& W, const MapVecd& mus_hat_j, const MapVecd& sigmas_hat_j, 
                          double& tau_hat, const MapMatd& C1, const MapVecd& gammas_hat_j, 
                          const MapMatd& W_prime, const MapVecd& C2_prime_j, bool scale ){
  Clock clock;
  int n = W.rows();
  int k = W.cols();
  int p1= C1.cols();
  MatrixXd Sig_j_orig = createDiagonalMatrix(sigmas_hat_j);
  
  // MatrixXd gammas_hat_j_reshape;
  // MatrixXd C1_prime = C1 * Map<MatrixXd>(gammas_hat_j.data(), p1, k);
  
  MatrixXd Z_j_hat = MatrixXd(n,k).setZero();
  VectorXd Z_j_hat_left = (Sig_j_orig.inverse() * mus_hat_j);
  // cout << "Sig_j_oirg\n" << Sig_j_orig << "\n" 
  //      << "W_prime_outprod\n" << outerprodV(W_prime.row(0)) << "\n" 
  //      << endl;
  
  clock.tick("forloop");
  
  // #pragma omp parallel for num_threads(2)
  for(int i = 0; i < n; i++){
    // MatrixXd BA_inv = as<MatrixXd>(W_prime[i]) * Sig_j_orig;
    MatrixXd BA_inv = outerprodV(W_prime.row(i)) * Sig_j_orig;    
    double g = BA_inv.diagonal().sum();
    MatrixXd S_ij = Sig_j_orig - ((Sig_j_orig * BA_inv).array()*(1/(1+g))).matrix() ;
    
    // VectorXd C1_prime_j = C1_prime.row(i);
    
    Z_j_hat.row(i) = (S_ij * (Z_j_hat_left  + (W.row(i) * C2_prime_j(i)).transpose() ) ).array() ;
    // cout <<  "S_ij_orig\n" << S_ij << "S_ij_transpose\n" << S_ij.transpose()  << "\n\n"
    //      << "Z_j_hat\n" << ( (Sig_j_orig.inverse() * mus_hat_j) + (W.row(i) * C2_prime_j(i)).transpose() ) << "\n\n"
    //      << "Z_j_hat_left\n" <<  (Sig_j_orig.inverse() * mus_hat_j) << "\n\n"
    //      << endl;
  }
  clock.tock("forloop");
  clock.stop("naptimes");
  return wrap(Z_j_hat);
}

// [[Rcpp::export]] 
void write_tensor_header(StringVector cells, StringVector TOE_IDs){
  ofstream file;
  for(int h = 0; h < cells.size(); h++){
    file.open(cells(h));
    file << "probe";
    for(int i = 0; i < TOE_IDs.size(); i++){
      file << ',' << TOE_IDs(i);
    }
    file << "\n"; file.close();
  }
}

// [[Rcpp::export]]
void write_tensor_list(StringVector& cells, StringVector& TOE_IDs, StringVector& cpgs, List& L){

  for(int i = 0; i < cpgs.size(); i++){
    
    MatrixXd A; A= as<MatrixXd>(L.at(i));
    // Matrix <double, TOE_IDs.size(),cells.size()> A = Map<Matrix<double, TOE_IDs.size(),cells.size()>> (A_v.array(), TOE_IDs.size(),cells.size());
    for(int h = 0; h < cells.size(); h++) { // by cell types
      
      ofstream file(cells(h), std::fstream::in | std::fstream::out | std::fstream::app);
      // file.open(cells(h));
      file << cpgs(i);
      for(int j = 0; j < TOE_IDs.size(); j++){ //one cpg line with TOE_IDs
        
        file << ',' << A(j,h);
      }
      file << "\n"; file.close();
    }
  }
}
// [[Rcpp::export]]
void write_tensor(StringVector& cells, StringVector& TOE_IDs, StringVector& cpgs, MapMatd& A, int& i){

  for(int h = 0; h < cells.size(); h++) { // by cell types
    ofstream file(cells(h), std::fstream::in | std::fstream::out | std::fstream::app);
    // file.open(cells(h));
    file << cpgs(i);
    for(int j = 0; j < TOE_IDs.size(); j++){ //one cpg line with TOE_IDs
      file << ',' << A(j,h);
    }
    file << "\n"; file.close();
  }
}


/*** R
# estimate_Z_j <- function(W, mus_hat_j, sigmas_hat_j, tau_hat, C1, gammas_hat_j, W_prime, C2_prime_j, scale){
#   n <- nrow(W)
#   k <- ncol(W)
#   p1 <- ncol(C1)
#   Z_j_hat <- matrix(0,n,k)
#   Sig_j_orig <- diag(sigmas_hat_j**2)
#   Sig_j <- matrix.inverse(Sig_j_orig)
#   C1_prime <- tcrossprod(C1, t(Reshape(gammas_hat_j,p1,k)))
#   for (i in 1:n){
#     #S_ij_inv <- W_prime[[i]]+Sig_j
#     #S_ij <- matrix.inverse(S_ij_inv)
#     ## the above two lines are the straightforward (slower) way to calculate S_ij; calculate 'matrix.inverse(W_prime[[i]]+Sig_j)' using the lemma from Miller 1981:
#     BA_inv <- W_prime[[i]] %*% Sig_j_orig
#     g <- sum(diag(BA_inv))
#     S_ij <- Sig_j_orig - ((1/(1+g))*(Sig_j_orig %*% BA_inv))
#     Z_j_hat[i,] = crossprod(S_ij, ( tcrossprod(Sig_j,t(mus_hat_j+C1_prime[i,])) + W[i,]*C2_prime_j[i] ) )
#     if (scale) Z_j_hat[i,] <- Z_j_hat[i,] / diag(S_ij)**0.5
#   }
#   return(Z_j_hat)
# }

# tmp1 <- estimate_Z_j_standalone(input_X, , mus_hat, sigmas_hat, tau_hat, scale)

*/



