#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd fast_between_cluster_inert_eigen(const Eigen::MatrixXd& Z,
                                                 const Eigen::VectorXd& cut_vals,
                                                 const Eigen::VectorXi& indices,
                                                 const Eigen::VectorXd& X_col,
                                                 const Eigen::VectorXd& w,
                                                 const Eigen::VectorXd& D) {

  int n_cuts = cut_vals.size();
  Eigen::VectorXd results(n_cuts);

  for (int c = 0; c < n_cuts; ++c) {
    double l = cut_vals[c];

    // Precompute sizes for two clusters
    //std::vector<int> indices_1, indices_2;
    //indices_1.reserve(indices.size());
    //indices_2.reserve(indices.size());

    // Compute weights and means for the first cluster
    double mu1 = 0.0;
    Eigen::VectorXd g1 = Eigen::VectorXd::Zero(Z.cols());

    // Compute weights and means for the second cluster
    double mu2 = 0.0;
    Eigen::VectorXd g2 = Eigen::VectorXd::Zero(Z.cols());

    for (int i = 0; i < indices.size(); ++i) {
      int idx = indices[i] - 1;  // Adjust for 1-based indexing in R
      if (X_col[idx] <= l) {
        //indices_1.push_back(idx);
        double w_idx = w[idx];
        mu1 += w_idx;
        g1 += w_idx * Z.row(idx).transpose();
      } else {
        //indices_2.push_back(idx);
        double w_idx = w[idx];
        mu2 += w_idx;
        g2 += w_idx * Z.row(idx).transpose();
      }
    }


    g1 /= mu1;
    g2 /= mu2;

    //if (indices_1.empty() || indices_2.empty()) {
    //  Rcpp::stop("One of the two clusters is empty!");
    //}

    // Compute between-cluster inertia
    Eigen::VectorXd diff = g1 - g2;
    double inert = diff.cwiseProduct(D).dot(diff);

    results[c] = (mu1 * mu2) / (mu1 + mu2) * inert;
  }

  return results;
}
