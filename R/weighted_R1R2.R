#' get weighted within-sample edge-counts and between-sample edge-counts
#' 
#' @param E an edge matrix representing a similarity graph. Each row represents an edge and records the indices of two ends of an edge in two columns. The indices of observations in sample 1 are from 1 to n1 and indices of observations in sample 2 are from 1+n1 to n1+n2. 
#' @param n1 number of observations in sample 1.
#' @param wei a vector of weights of each edge.
#' 
#' @return 
#' \item{R1}{the weighted within-sample edge-count for sample 1.}
#' \item{R2}{the weighted within-sample edge-count for sample 2.}
#' \item{R}{the weighted between-sample edge-count.}
#' 
#' @keywords internal
#' @export
#' 
weighted_R1R2 <- function(E, n1, wei){
  G1 = 1:n1
  
  R1 = R2 = R = 0
  e1 = E[, 1] %in% G1
  e2 = E[, 2] %in% G1
  
  R1 = sum(wei[e1 + e2 == 2])
  R2 = sum(wei[e1 + e2 == 0])
  R = sum(wei[e1 + e2 == 1])
  
  return(list(R1=R1, R2=R2, R=R))
}
