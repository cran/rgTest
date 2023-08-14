#' Get distance matrix
#'
#' This function returns the distance matrix using L2 distance.
#'
#' @param y dataset of the pooled data
#' 
#' @return A distance matrix based on the L2 distance.
#' 
#' @export
#' 
#' @examples
#' data(example0)
#' data = as.matrix(example0$data)     # pooled dataset
#' getdis(data)
#' 
getdis <- function(y){ #L2 distance
  n = dim(y)[1]
  G = y%*%t(y)
  g = diag(G)
  dis = sqrt( matrix(rep(g,n),n) + matrix(rep(g,n),n,byrow=T) - 2*G )
  dis
}
