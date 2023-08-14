#' construct k-mst
#' 
#' @importFrom ade4 mstree
#' @importFrom stats as.dist
#' @param y data
#' @param dis distance matrix
#' @param k parameter in K-MST, with default 1
#' 
#' @return An edge matrix representing a similarity graph. Each row represents an edge and records the indices of two ends of an edge in two columns. 
#' 
#' @keywords internal
#' @export
#' 
kmst <- function(y=NULL, dis=NULL, k=1){
  if (is.null(dis) && is.null(y)){
    cat("Please input data or the distance matrix!\n")
    return(0)
  }
  if (is.null(dis)) dis = getdis(y)
  mymst = mstree(as.dist(dis),k)
  cbind(mymst[,1], mymst[,2])
}
