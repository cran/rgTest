#' get lists of permuted weighted within-sample edge-counts and between-sample edge-counts
#' @param n_per number of permutations.
#' @param E an edge matrix representing a similarity graph. Each row represents an edge and records the indices of two ends of an edge in two columns. The indices of observations in sample 1 are from 1 to n1 and indices of observations in sample 2 are from 1+n1 to n1+n2. 
#' @param n1 number of observations in sample 1.
#' @param n2 number of observations in sample 2.
#' @param wei a vector of weights of each edge.
#' @param progress_bar a logical evaluating to TRUE or FALSE indicating whether a progress bar of the permutation should be printed.
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' 
#' @return 
#' \item{R1}{the permuted weighted within-sample edge-counts for sample 1.}
#' \item{R2}{the permuted weighted within-sample edge-counts for sample 2.}
#' \item{R}{the permuted weighted between-sample edge-counts.}
#' 
#' @keywords internal
#' @export
#' 
permu_edge <- function(n_per, E, n1, n2, wei, progress_bar = FALSE){
  obs = n1+n2
  if(progress_bar){
    pb = txtProgressBar(min = 0, max = n_per, initial = 0) 
    temp = sapply(1:n_per, function(peri){
      setTxtProgressBar(pb,peri)
      per = sample(obs)
      new_E = matrix(per[E], ncol = 2)
      edgeinfo = weighted_R1R2(new_E, n1, wei)
      return(unlist(edgeinfo))
    })
    close(pb)
    distri_R1 = temp['R1', ]
    distri_R2 = temp['R2', ]
    distri_R = temp['R', ]
  }else{
    temp = sapply(1:n_per, function(peri){
      per = sample(obs)
      new_E = matrix(per[E], ncol = 2)
      edgeinfo = weighted_R1R2(new_E, n1, wei)
      unlist(edgeinfo)
    })
    
    distri_R1 = temp['R1', ]
    distri_R2 = temp['R2', ]
    distri_R = temp['R', ]
  }
  return(list(R1=distri_R1, R2=distri_R2, R=distri_R))
}
