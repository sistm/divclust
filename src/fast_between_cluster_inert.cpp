#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fast_between_cluster_inert(const NumericMatrix& Z,
                                         const NumericVector& cut_vals,
                                         const IntegerVector& indices,
                                         const NumericVector& X_col,
                                         const NumericVector& w,
                                         const NumericVector& D) {

  int n_cuts = cut_vals.size();
  NumericVector results(n_cuts);
  for (int c = 0; c < n_cuts; ++c) {
    double l = cut_vals[c];

    // Partition indices based on the cutoff value l
    IntegerVector indices_1, indices_2;

    // Compute weights and means for the first cluster
    double mu1 = 0.0;
    NumericVector g1(Z.ncol());

    // Compute weights and means for the second cluster
    double mu2 = 0.0;
    NumericVector g2(Z.ncol());

    for (int i = 0; i < indices.size(); ++i) {

      int idx = indices[i] - 1;  // Adjusting for 1-based indexing in R
      if (X_col[idx] <= l) {
        indices_1.push_back(idx);
        mu1 += w[idx];
        g1 += w[idx] * Z.row(idx);

      } else {
        indices_2.push_back(idx);
        mu2 += w[idx];
        g2 += w[idx] * Z.row(idx);
      }
    }

    g1 = g1 / mu1;

    g2 = g2 / mu2;

    if (indices_1.size() == 0 || indices_2.size() == 0) {
      stop("One of the two clusters is empty!");
    }

    // Compute between-cluster inertia
    double inert = 0.0;
    for (int j = 0; j < Z.ncol(); ++j) {
      double diff = g1[j] - g2[j];
      inert += D[j] * diff * diff;
    }

    results[c] = (mu1 * mu2) / (mu1 + mu2) * inert;
  }

  return results;
}
