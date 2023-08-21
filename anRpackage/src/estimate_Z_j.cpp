#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>

using namespace std;
using namespace Rcpp;
using namespace Eigen;
//using Eigen::Map;                       // 'maps' rather than copies
//using Eigen::MatrixXd;                  // variable size matrix, double precision
//using Eigen::VectorXd;                  // variable size vector, double precision
//using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/



// [[Rcpp::export]]
MatrixXd createDiagonalMatrix(const VectorXd& diagonalValues) {
  DiagonalMatrix<double, Dynamic> diagMatrix(diagonalValues.size());
  diagMatrix.diagonal() = diagonalValues;
  return diagMatrix;
}

// [[Rcpp::export]]
MatrixXd estimate_Z_j_cpp(MatrixXd& W, VectorXd& mus_hat_j, VectorXd& sigmas_hat_j, 
                      double& tau_hat, MatrixXd& C1, VectorXd& gammas_hat_j, 
                      List W_prime, VectorXd& C2_prime_j, bool scale ){
  int n = W.rows();
  int k = W.cols();
  int p1= C1.cols();
  MatrixXd Sig_j_orig = createDiagonalMatrix(sigmas_hat_j);
  
  MatrixXd Sig_j = Sig_j_orig.inverse();
  MatrixXd gammas_hat_j_reshape;
  gammas_hat_j_reshape = Map<MatrixXd>(gammas_hat_j.data(), p1, k);
  MatrixXd C1_prime = C1 * gammas_hat_j_reshape;
  
  MatrixXd Z_j_hat(n,k);
  for(int i = 0; i < n; i++){
    //MatrixXd W_tmp = Map<MatrixXd>(W_prime[i].begin(), W_prime[i].nrow(), W_prime[i].ncol());
    MatrixXd W_tmp = as<MatrixXd>(W_prime[i]);
    MatrixXd BA_inv = W_tmp * Sig_j_orig;
    double g = BA_inv.diagonal().sum();
    MatrixXd S_ij = Sig_j_orig - ((1/(1+g)))*(Sig_j_orig * BA_inv);
    VectorXd C1_prime_j = C1_prime.row(i);
    
    VectorXd S_ij_tmptmp = W.row(i) * C2_prime_j(i);
    MatrixXd S_ij_tmp = (Sig_j * (mus_hat_j + C1_prime_j) + S_ij_tmptmp);
    VectorXd S = S_ij.transpose().eval() * S_ij_tmp;
    Z_j_hat.row(i) = S;
    /*
    cout << "W_prime_j < Sig_j_orig" << endl;
    cout << W_tmp << "\n" << Sig_j_orig << endl;
    
    cout << "BA_inv : g : S_ij : S_ij_tmp" << endl;
    cout << BA_inv << "\n" << g << "\n" << S_ij << "\n" << S_ij_tmp << endl;
    cout << endl;
    
    cout << mus_hat_j.rows() << "\t" << C1_prime_j.rows() << endl;
    cout << "mus_hat_j + C1_prime[1,] < W[1,] * C2_prime_j" << endl;
    cout << mus_hat_j + C1_prime_j << "\n" <<  W.row(i) * C2_prime_j(i) << endl;
    
    cout << "Sig_j %*% (mus_hat_j + C1_prime[1,]) >> W[1,] * C2_prime[1,1]" << endl;
    cout << Sig_j * (mus_hat_j + C1_prime_j) << "\n" << W.row(i) * C2_prime_j(i) << '\n' << endl ;
    
    cout << S_ij_tmptmp << "\n" << S_ij_tmp << "\n" << S_ij << '\n' << endl;
    cout << S << endl;
     */
  }
  return Z_j_hat;
}



/*** R
# n <- nrow(W)
# k <- ncol(W)
# p1 <- ncol(C1)
# Z_j_hat <- matrix(0,n,k)
# Sig_j_orig <- diag(sigmas_hat_j**2)
# Sig_j <- matrix.inverse(Sig_j_orig)
# C1_prime <- tcrossprod(C1, t(Reshape(gammas_hat_j,p1,k)))
# for (i in 1:n){
#   #S_ij_inv <- W_prime[[i]]+Sig_j
#   #S_ij <- matrix.inverse(S_ij_inv)
#   ## the above two lines are the straightforward (slower) way to calculate S_ij; calculate 'matrix.inverse(W_prime[[i]]+Sig_j)' using the lemma from Miller 1981:
#   BA_inv <- W_prime[[i]] %*% Sig_j_orig
#   g <- sum(diag(BA_inv))
#   S_ij <- Sig_j_orig - ((1/(1+g))*(Sig_j_orig %*% BA_inv))
#   Z_j_hat[i,] = crossprod(S_ij, ( tcrossprod(Sig_j,t(mus_hat_j+C1_prime[i,])) + W[i,]*C2_prime_j[i] ) )
#   if (scale) Z_j_hat[i,] <- Z_j_hat[i,] / diag(S_ij)**0.5
# }
# return(Z_j_hat)


*/



