#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector fast_compute_cut_values(Rcpp::NumericVector col) {
  // Sort and get unique values
  NumericVector col2 = Rcpp::unique(col);
  std::sort(col2.begin(), col2.end());

  int n = col2.size();
  if (n < 2) return col2; // return unique values if < 2

  // Compute midpoints
  Rcpp::NumericVector cuts(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    cuts[i] = 0.5 * (col2[i] + col2[i + 1]);
  }

  return cuts;
}
