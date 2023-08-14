#' get the approximate test statistic and p-value based on asymptotic theory using robust generalized edge-count test 
#' 
#' @importFrom stats pchisq
#' @param asy_res analytic expressions of expectations, variances and covariances
#' @param R1_test weighted within-sample edge-counts of sample 1
#' @param R2_test weighted within-sample edge-counts of sample 2
#' 
#' @return A list containing the following components:
#' \item{test_statistic}{the asymptotic test statistic using robust generalized graph-based test.}
#' \item{p_value}{the asymptotic p-value using robust generalized graph-based test.}
#' 
#' @keywords internal
#' @export
#' 
asy_gen <- function(asy_res, R1_test, R2_test){
  sigma = matrix(c(asy_res$sig11, asy_res$sig12, asy_res$sig12, asy_res$sig22), nrow = 2)
  temp_asy = c(R1_test - asy_res$mu1, R2_test - asy_res$mu2) %*% solve(sigma)
  z_gen = temp_asy %*% c(R1_test - asy_res$mu1, R2_test - asy_res$mu2)
  test_statistic = z_gen[1, 1]
  return(list(test_statistic = test_statistic, p_value = 1-pchisq(test_statistic, df = 2)))
}
