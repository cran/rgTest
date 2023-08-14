#' get the test statistic and p-value based on permutation using robust generalized edge-count test 
#' 
#' @importFrom stats cov
#' @param R1_list list of permuted weighted within-sample edge-counts of sample 1
#' @param R2_list list of permuted weighted within-sample edge-counts of sample 2
#' @param R1_test weighted within-sample edge-counts of sample 1
#' @param R2_test weighted within-sample edge-counts of sample 2
#' @param n_per number of permutations
#' 
#' @return The p-value based on permutation distribution using robust generalized graph-based test.
#' 
#' @keywords internal
#' @export
#' 
permu_gen <- function(R1_list, R2_list, R1_test, R2_test, n_per){
  mu1 = mean(R1_list)
  mu2 = mean(R2_list)
  sigma = cov(cbind(R1_list, R2_list))
  sigma_inverse = solve(sigma)
  temp = c(R1_test - mu1, R2_test - mu2) %*% sigma_inverse
  gen_t = temp %*% c(R1_test - mu1, R2_test - mu2)
  
  null_gen = rep(0, n_per)
  for (i in 1:n_per) {
    temp1 = c(R1_list[i] - mu1, R2_list[i] - mu2) %*% sigma_inverse
    null_gen[i] = temp1 %*% c(R1_list[i] - mu1, R2_list[i] - mu2)
  }
  return(sum(null_gen > gen_t[1, 1]))
}
