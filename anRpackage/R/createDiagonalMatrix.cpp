#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
MatrixXd createDiagonalMatrix(const Eigen::VectorXd& diagonalValues) {
  DiagonalMatrix<double, Dynamic> diagMatrix(diagonalValues.size());
  diagMatrix.diagonal() = diagonalValues;
  return diagMatrix;
}

/*** R
# Example usage
diagonalValues <- c(1.0, 2.0, 3.0)
result <- createDiagonalMatrix(diagonalValues)
print(result)
*/