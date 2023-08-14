#' Weighted function
#' 
#' This weight function returns the inverse of the max node degree of an edge.
#' 
#' @param a node degree of one end of an edge
#' @param b node degree of another end of an edge
#' 
#' @return The weight uses the max node degrees of an edge.
#' 
#' @export
#' 
#' @examples
#' # For an edge where one end has a node degree of 5
#' # another end has a node degree of 6
#' weiMax(6, 5)
#' 
weiMax <- function(a, b){
  1/max(a, b)
}