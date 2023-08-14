#' Weighted function
#' 
#' This weight function returns the inverse of the geometric average of the node degrees of an edge.
#' 
#' @param a node degree of one end of an edge
#' @param b node degree of another end of an edge
#' 
#' @return The weight uses the geometric average of the node degrees of an edge.
#' 
#' @export
#' 
#' @examples
#' # For an edge where one end has a node degree of 5
#' # another end has a node degree of 6
#' weiGeo(6, 5)
#' 
weiGeo <- function(a, b){
  2/(a + b)
}