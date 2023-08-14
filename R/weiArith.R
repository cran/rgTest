#' Weighted function
#' 
#' This weight function returns the inverse of the arithmetic average of the node degrees of an edge.
#' 
#' @param a node degree of one end of an edge
#' @param b node degree of another end of an edge
#' 
#' @return The weight uses the arithmetic average of the node degrees of an edge.
#' 
#' @export
#' 
#' @examples
#' # For an edge where one end has a node degree of 5
#' # another end has a node degree of 6
#'  weiArith(6, 5)
#' 
weiArith <- function(a, b){
  1/sqrt(a*b)
}