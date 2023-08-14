#' get the approximate test statistic and p-value based on asymptotic theory using robust max-type edge-count test 
#' 
#' @importFrom stats pnorm
#' @param asy_res analytic expressions of expectations, variances and covariances
#' @param R1_test weighted within-sample edge-counts of sample 1
#' @param R2_test weighted within-sample edge-counts of sample 2
#' @param n1 number of observations in sample 1
#' @param n2 number of observations in sample 2
#' 
#' @return A list containing the following components:
#' \item{test_statistic}{the asymptotic test statistic using robust max-type graph-based test.}
#' \item{p_value}{the asymptotic p-value using robust max-type graph-based test.}
#' 
#' @keywords internal
#' @export
#' 
asy_max <- function(asy_res, R1_test, R2_test, n1, n2){
  max_diff_t = ((R1_test - R2_test) - (asy_res$mu1 - asy_res$mu2)) / 
    (sqrt(asy_res$sig11 + asy_res$sig22 - 2*asy_res$sig12))
  p_w = (n1-1)/(n1+n2-2)
  q_w = 1 - p_w
  max_w_t = ((q_w*R1_test + p_w*R2_test) - (q_w*asy_res$mu1 + p_w*asy_res$mu2)) / 
    (sqrt(q_w^2*asy_res$sig11 + p_w^2*asy_res$sig22 + 2*p_w*q_w*asy_res$sig12))
  max_t = max(abs(max_diff_t), max_w_t)
  p_max = 1-2*(pnorm(max_t)-0.5)*pnorm(max_t)
  test_statistic = max_t
  return(list(test_statistic = test_statistic, p_value = p_max))
}
