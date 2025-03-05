#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string fastSubstring_cut1(const std::string& s) {
  // Return the substring starting from the 3rd character
  if (s.size() <= 2) return ""; // Handle edge cases
  return s.substr(2); // substr uses zero-based indexing
}
