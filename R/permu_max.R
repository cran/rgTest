#' get the test statistic and p-value based on permutation using robust max-type edge-count test 
#' 
#' @importFrom stats cov
#' @param R1_list list of permuted weighted within-sample edge-counts of sample 1
#' @param R2_list list of permuted weighted within-sample edge-counts of sample 2
#' @param R1_test weighted within-sample edge-counts of sample 1
#' @param R2_test weighted within-sample edge-counts of sample 2
#' @param n1 number of observations in sample 1
#' @param n2 number of observations in sample 2
#' @param n_per number of permutations
#' 
#' @return The p-value based on permutation distribution using robust max-type graph-based test.
#' 
#' @keywords internal
#' @export
#' 
permu_max <- function(R1_list, R2_list, R1_test, R2_test, n1, n2, n_per){
  mu1 = mean(R1_list)
  mu2 = mean(R2_list)
  sigma = cov(cbind(R1_list, R2_list))
  max_diff_t = ((R1_test - R2_test) - (mu1 - mu2)) / 
    (sqrt(sigma[1, 1] + sigma[2, 2] - 2*sigma[1, 2]))
  p_w = (n1-1)/(n1+n2-2)
  q_w = 1 - p_w
  max_w_t = ((q_w*R1_test + p_w*R2_test) - (q_w*mu1 + p_w*mu2)) / 
    (sqrt(q_w^2*sigma[1, 1] + p_w^2*sigma[2, 2] + 2*p_w*q_w*sigma[1, 2]))
  max_t = max(abs(max_diff_t), max_w_t)
  
  null_max = rep(0, n_per)
  deno1 = sqrt(sigma[1, 1] + sigma[2, 2] - 2*sigma[1, 2])
  deno2 = sqrt(q_w^2*sigma[1, 1] + p_w^2*sigma[2, 2] + 2*p_w*q_w*sigma[1, 2])
  for (i in 1:n_per) {
    max_diff_temp = ((R1_list[i] - R2_list[i]) - (mu1 - mu2))/deno1
    max_w_temp = ((q_w*R1_list[i] + p_w*R2_list[i]) - (q_w*mu1 + p_w*mu2))/deno2
    null_max[i] = max(abs(max_diff_temp), max_w_temp)
  }
  return(sum(null_max > max_t))
}
