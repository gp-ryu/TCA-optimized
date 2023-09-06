#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <vector>

// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
// // // [[Rcpp::depends(RcppParallel)]]
// #include <RcppParallel.h>
//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>
#include <thread>
  
#define THREAD_NUM 16
using namespace std;
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> MapMatd;
typedef Map<VectorXd> MapVecd;
typedef vector<double> vectord;


// [[Rcpp::export]]
MatrixXd outerproduct(const VectorXd & x) {
  MatrixXd m = x * x.transpose();
  return m;
}

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

MatrixXd outerprodM(const MapMatd& A){
  int n = A.rows();
  MatrixXd AtA = MatrixXd(n,n).setZero().selfadjointView<Lower>()
                              .rankUpdate(A.adjoint());
  return AtA;
}


// [[Rcpp::export]]
SEXP estimate_Z_j_standalone(MatrixXd& X, MapMatd& W, MapMatd& mus_hat, MapMatd& sigmas_hat, double& tau_hat, bool scale) {
  Clock clock;
  // int p1 = C1.cols();
  // int p2 = C2.cols();
  int m = mus_hat.rows();
  int n = W.rows();
  int k = W.cols();
  int p1 = 0;
  int p2 = 0;
  clock.tick("res");
  unordered_map <int, vectord > cellTypes;
  
  // MatrixXd B= Eigen::MatrixXd::Zero(m,n); double* pt_B = &B(0);
  // MatrixXd NK = Eigen::MatrixXd::Zero(m,n); double* pt_NK = &NK(0);
  // MatrixXd CD4T= Eigen::MatrixXd::Zero(m,n); double* pt_CD4T = &CD4T(0);
  // MatrixXd CD8T= Eigen::MatrixXd::Zero(m,n); double* pt_CD8T = &CD8T(0);
  // MatrixXd Mono= Eigen::MatrixXd::Zero(m,n); double* pt_Mono = &Mono(0);
  // MatrixXd Neutro= Eigen::MatrixXd::Zero(m,n); double* pt_Neutro = &Neutro(0);
  // MatrixXd Eosino= Eigen::MatrixXd::Zero(m,n); double* pt_Eosino = &Eosino(0);
  
  // cellTypes[0](pt_B, m,n);
  
  // vectord B; vectord* ptr_B = &B;
  // vectord CD4T; vectord* ptr_CD4T = &CD4T;
  // vectord CD8T; vectord* ptr_CD8T = &CD8T;
  // vectord NK;   vectord* ptr_NK = &NK;
  // vectord Mono; vectord* ptr_Mono = &Mono;
  // vectord Neutro; vectord* ptr_Neutro = &Neutro;
  // vectord Eosino; vectord* ptr_Eosino = &Eosino;
  
  // cellTypes[0]  = &B;
  // cellTypes[1]  = &NK;
  // cellTypes[2]  = &CD4T;
  // cellTypes[3]  = &CD8T;
  // cellTypes[4]  = &Mono;
  // cellTypes[5]  = &Neutro;
  // cellTypes[6]  = &Eosino;
  
  // NumericMatrix B(m,n);
  // NumericMatrix CD4T(m,n);
  // NumericMatrix CD8T(m,n);
  // NumericMatrix NK(m,n);
  // NumericMatrix Mono(m,n);
  // NumericMatrix Neutro(m,n);
  // NumericMatrix Eosino(m,n);
  
  clock.tock("res");
  
  // clock.tick("W_prime");
  // MatrixXd W_prime = (W)/tau_hat;
  // clock.tock("W_prime");
  // if(p2 != 0 && deltas_hat.cols() != 0){
    // MatrixXf C2_prime(n,m); C2_prime = (X - (C2 * deltas_hat)/pow(tau_hat, 2));
    // }else{
      clock.tick("C2_prime");
      
      double* pt = &X(0);
      Map<const MatrixXd> C2_prime (pt, m,n);
      // MatrixXd C2_prime; C2_prime.noalias() = X;
      clock.tock("C2_prime");
      // }
  if(p1 != 0){
    // MatrixXf C1_prime(m,k); C1_prime = C1 * (gammas_hat_j.data(), p1, k);
  }
  clock.tick("Sig_j_orig_mat");
  MatrixXd Sig_j_orig_mat(m,k); Sig_j_orig_mat = sigmas_hat.cwiseProduct(sigmas_hat);
  clock.tock("Sig_j_orig_mat");
  
  // cout << "W_prime[[j]]\n" << outerproduct( W_prime.row(0) ) << endl;
  // cout << "BA_inv\n" << outerproduct( W_prime.row(0) ) *createDiagonalMatrix(Sig_j_orig_mat.row(0))  << endl;
  
  #ifdef _OPENMP
  REprintf("Number of threads=%i\n", omp_get_max_threads());
  #endif
  clock.tick("forloop");
  
#pragma omp parallel for num_threads(omp_get_max_threads())
  for(int j = 0; j < 200; j++){
    // MatrixXd Z_j_hat = MatrixXd::Zero(m,k);
    MatrixXd Sig_j_orig(k,k); Sig_j_orig = createDiagonalMatrix(Sig_j_orig_mat.row(j));
    
        for(int i = 0; i < n; i++){
      #pragma omp critical
      {
          // if (Progress::check_abort() )
          //   return -1.0;
          if(p1 == 0 && p2 == 0){
            MatrixXd W_tmp(k,k); W_tmp = outerproduct( W.row(i)/tau_hat );//* W_prime.row(i).transpose());
            
            clock.tick("BA_inv_outerprod");
            MatrixXd BA_inv(k,k); BA_inv = W_tmp * Sig_j_orig;
            clock.tock("BA_inv_outerprod");
            double g = BA_inv.diagonal().sum();
            
            clock.tick("S");
            MatrixXd S_ij = Sig_j_orig - ((1/(1+g)))*(Sig_j_orig * BA_inv);
            // Z_j_hat.row(i) = S_ij.transpose().eval() * (Sig_j_orig.inverse() * mus_hat.row(j) + W.row(i) * C2_prime(i,j));
            VectorXd res_tmp = S_ij.transpose().eval() * (Sig_j_orig.inverse() * mus_hat.row(j) + W.row(i) * C2_prime(i,j));
            clock.tock("S");
            for(int h = 0; h < k; h++){
              cellTypes[h].push_back(res_tmp(h));
            }
            // Rcout << "S_ij\n" << S_ij << '\n' << "Sig_j\n" << Sig_j << '\n' << "W_tmp\n" << W_tmp << '\t' << "mus_hat_1 / C2_prime_1,1\n" << mus_hat.row(1) << '\t' << C2_prime(1,1) << endl;
            // Rcout << Z_j_hat.row(j) << endl;
            // if(scale){cout << endl;}
          }
        }
        // cout << Z_j_hat.block<4,6>(0,0) << endl;
        // return wrap(Z_j_hat);
      }
      
// #pragma omp for
//         for(int h = 0; h < k; h++){
//           for(int q = 0; q < n; q++){
//             cellTypes[h].push_back(Z_j_hat(q,h));
//           }
//         }
//       
  REprintf("\r| ---- ---- %d %% ---- ---- |  %d th", j/m, j);
  }
  REprintf("\n");
  clock.tock("remap");


  clock.tock("forloop");
  clock.tick("Listrization");
  // List res = List::create(Named("B") = wrap(cellTypes[0]), Named("NK") = wrap(cellTypes[1]),
  //                         Named("CD4T") = wrap(cellTypes[2]), Named("CD8T") = wrap(cellTypes[3]),
  //                         Named("Mono") = wrap(cellTypes[4]), Named("Neutro") = wrap(cellTypes[5]), Named("Eosino") = wrap(cellTypes[6]));
  clock.tock("Listrization");
  clock.tick("matrixization");
  
  vectord tmp;tmp.reserve(cellTypes[0].size() * k); for(int i = 0; i < k; i++){tmp.insert(tmp.end(), cellTypes[i].begin(), cellTypes[i].end());}
  NumericVector v_B; v_B = wrap(tmp);
  v_B.attr("dim") = Dimension(200,n,k);
  clock.tock("matrixization");
  clock.stop("naptime");
  
  return v_B;
}
