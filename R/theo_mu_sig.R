#' get analytic expressions of expectations, variances and covariances 
#' 
#' @param E an edge matrix representing a similarity graph. Each row represents an edge and records the indices of two ends of an edge in two columns. The indices of observations in sample 1 are from 1 to n1 and indices of observations in sample 2 are from 1+n1 to n1+n2. 
#' @param n1 number of observations in sample 1
#' @param n2 number of observations in sample 2
#' @param weights weights assigned to each edges
#' 
#' @return 
#' \item{mu}{the expectation of the between-sample edge-count.}
#' \item{mu1}{the expectation of the within-sample edge-count for sample 1.}
#' \item{mu2}{the expectation of the within-sample edge-count for sample 2.}
#' \item{sig}{the variance of the between-sample edge-count.}
#' \item{sig11}{the variance of the within-sample edge-count for sample 1.}
#' \item{sig22}{the variance of the within-sample edge-count for sample 2.}
#' \item{sig12}{the covariance of the within-sample edge-counts.}
#' 
#' @keywords internal
#' @export
#' 
theo_mu_sig <- function(E, n1, n2, weights){
  N = n1+n2
  Ebynode = vector("list", N)
  Ebynode_deg = vector("list", N)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
    Ebynode_deg[[E[i,1]]] = c(Ebynode_deg[[E[i,1]]],weights[i])
    Ebynode_deg[[E[i,2]]] = c(Ebynode_deg[[E[i,2]]],weights[i])
  }
  
  mu_1 = sum(weights)*n1*(n1-1)/N/(N-1)
  mu_2 = sum(weights)*n2*(n2-1)/N/(N-1)
  mu = 2*sum(weights)*n1*n2/N/(N-1)
  
  part11 = sum(weights^2)*n1*(n1-1)/N/(N-1)
  part21 = sum(weights^2)*n2*(n2-1)/N/(N-1)
  part1 = 2*sum(weights^2)*n1*n2/N/(N-1)
  #edge pair
  ss = 0
  for (i in 1:N) {
    if(length(Ebynode[[i]]) == 1){
      next
    }
    ss = ss + sum(sapply(1:(length(Ebynode[[i]])-1), function(ii) sum(Ebynode_deg[[i]][ii]*Ebynode_deg[[i]][(ii+1):length(Ebynode[[i]])])))
  }
  part12 = 2*ss*n1*(n1-1)*(n1-2)/N/(N-1)/(N-2)
  part22 = 2*ss*n2*(n2-1)*(n2-2)/N/(N-1)/(N-2)
  part2 = 2*ss*(n1*n2*(n2-1)/N/(N-1)/(N-2) + n1*n2*(n1-1)/N/(N-1)/(N-2))
  #different edges
  s_4 = sum(sapply(1:nrow(E), function(ii) sum(weights[ii]*weights)))
  s_4 = s_4 - 2*ss - sum(weights^2)
  
  part13 = s_4*n1*(n1-1)*(n1-2)*(n1-3)/N/(N-1)/(N-2)/(N-3)
  part23 = s_4*n2*(n2-1)*(n2-2)*(n2-3)/N/(N-1)/(N-2)/(N-3)
  part3 = s_4*4*n1*n2*(n1-1)*(n2-1)/N/(N-1)/(N-2)/(N-3)
  
  sig1 = part11+part12+part13 - mu_1^2
  sig2 = part21+part22+part23 - mu_2^2
  sig = part1+part2+part3 - mu^2
  
  part_cov = s_4*n1*(n1-1)*n2*(n2-1)/N/(N-1)/(N-2)/(N-3)
  sig_12 = part_cov - mu_1*mu_2 
  return(list(mu = mu, mu1 = mu_1, mu2 = mu_2, sig = sig, sig11 = sig1, sig22 = sig2, sig12 = sig_12))
}
