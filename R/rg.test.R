#' Robust graph-based two sample test
#' 
#' @description Performs robust graph-based two sample test. 
#' 
#' @usage 
#' rg.test(data.X, data.Y, dis = NULL, E = NULL, n1, n2, k = 5, weigh.fun, perm.num = 0, 
#' test.type = list("ori", "gen", "wei", "max"), progress_bar = FALSE)
#' 
#' @param data.X a numeric matrix for observations in sample 1.
#' @param data.Y a numeric matrix for observations in sample 2.
#' @param dis a distance matrix of the pooled dataset of sample 1 and sample 2. The indices of observations in sample 1 are from 1 to n1 and indices of observations in sample 2 are from 1+n1 to n1+n2 in the pooled dataset.
#' @param E an edge matrix representing a similarity graph. Each row represents an edge and records the indices of two ends of an edge in two columns. The indices of observations in sample 1 are from 1 to n1 and indices of observations in sample 2 are from 1+n1 to n1+n2. 
#' @param n1 number of observations in sample 1.
#' @param n2 number of observations in sample 2.
#' @param k parameter in K-MST, with default 5.
#' @param weigh.fun weighted function which returns weights of each edge and is a function of node degrees.
#' @param perm.num number of permutations used to calculate the p-value (default=1000). Use 0 for getting only the approximate p-value based on asymptotic theory. 
#' @param test.type type of graph-based test. This must be a list containing elements chosen from "ori", "gen", "wei", and "max", with default 'list("ori", "gen", "wei", "max")'. "ori" refers to robust orignial edge-count test, "gen" refers to robust generalized edge-count test, "wei" refers to robust weighted edge-count test and "max" refers to robust max-type edge-count tests. 
#' @param progress_bar a logical evaluating to TRUE or FALSE indicating whether a progress bar of the permutation should be printed.
#' 
#' @details The input should be one of the following: 
#' 1. datasets of the two samples; 
#' 2. the distance matrix of the pooled dataset; 
#' 3. the edge matrix generated from a similarity graph. 
#' 
#' Typical usages are:
#' 
#' \preformatted{rg.test(data.X, data.Y, n1, n2, weigh.fun, ...)}
#' \preformatted{rg.test(dis, n1, n2, weigh.fun, ...)}
#' \preformatted{rg.test(E, n1, n2, weigh.fun, ...)}
#' If the data matrices or the distance matrix are used, the similarity graph is generated using K-MST.
#' 
#' @return A list containing the following components:
#' \item{asy.ori.statistic}{the asymptotic test statistic using robust original graph-based test.}
#' \item{asy.ori.pval}{the asymptotic p-value using robust original graph-based test.}
#' \item{asy.gen.statistic}{the asymptotic test statistic using robust generalized graph-based test.}
#' \item{asy.gen.pval}{the asymptotic p-value using robust generalized graph-based test.}
#' \item{asy.wei.statistic}{the asymptotic test statistic using robust weighted graph-based test.}
#' \item{asy.wei.pval}{the asymptotic p-value using robust weighted graph-based test.}
#' \item{asy.max.statistic}{the asymptotic test statistic using robust max-type graph-based test.}
#' \item{asy.max.pval}{the asymptotic p-value using robust max-type graph-based test.}
#' \item{perm.ori.pval}{the p-value based on permutation using robust original graph-based test.}
#' \item{perm.gen.pval}{the p-value based on permutation using robust generalized graph-based test.}
#' \item{perm.wei.pval}{the p-value based on permutation using robust weighted graph-based test.}
#' \item{perm.max.pval}{the p-value based on permutation using robust max-type graph-based test.}
#' 
#' 
#' @export
#' 
#' @examples
#' ## Simulated from Student's t-distribution. 
#' ## Observations for the two samples are from different distributions.
#' data(example0)
#' data = as.matrix(example0$data)     # pooled dataset
#' label = example0$label              # label of observations
#' s1 = data[label == 'sample 1', ]    # sample 1
#' s2 = data[label == 'sample 2', ]    # sample 2
#' num1 = nrow(s1)                     # number of observations in sample 1
#' num2 = nrow(s2)                     # number of observations in sample 2
#' 
#' ## Graph-based two sample test using data as input
#' rg.test(data.X = s1, data.Y = s2, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 0)
#' 
#' ## Graph-based two sample test using distance matrix as input
#' dist = example0$distance
#' rg.test(dis = dist, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 0)
#' 
#' ## Graph-based two sample test using edge matrix of the similarity graph as input
#' E = example0$edge
#' rg.test(E = E, n1 = num1, n2 = num2, weigh.fun = weiMax, perm.num = 0)
#' 
rg.test <- function(data.X, data.Y, dis = NULL, E = NULL, n1, n2, k = 5, weigh.fun, perm.num = 0, test.type = list("ori", "gen", "wei", "max"), progress_bar = FALSE){
  #check inputs
  if (missing(data.X) && is.null(dis)  && missing(E) || (missing(data.Y) && is.null(dis)  && is.null(E))){
    stop("Please input data, the distance matrix or the edge matrix!\n")
  }
  if(missing(n1) || missing(n2)){
    stop("Please enter sample size!\n")
  }
  
  if(!all(test.type %in% list("ori", "gen", "wei", "max"))){
    stop("Please choose correct test type!\n")
  }

  ### get edge matrix
  n = n1+n2
  if(missing(E)){
    if (missing(dis)) {
      y = rbind(data.X, data.Y)
      E = kmst(y=y, k=k)
    }else{
     E = kmst(dis=dis, k=k)
    }
  }

  #asymptotic results
  Ebynode = vector("list", n)
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }
  
  nodedeg = sapply(Ebynode, length)
  
  weis = rep(0, nrow(E))
  for(i in 1:nrow(E)){
    weis[i] = weigh.fun(nodedeg[E[i, 1]], nodedeg[E[i, 2]])
  }
  
  weighted_edge_count = weighted_R1R2(E, n1, weis)
  # browser()
  asy_res = theo_mu_sig(E, n1, n2, weis)
  
  res_list = list()
  if("ori" %in% test.type){
    #original edge count test
    ori_test = (weighted_edge_count$R - asy_res$mu)/sqrt(asy_res$sig)
    ori_asy_p = pnorm(ori_test)
    res_list = c(res_list, asy.ori.statistic = ori_test, asy.ori.pval = ori_asy_p)
  }
  if("gen" %in% test.type){
    #generalized edge count test
    gen_asy_p = asy_gen(asy_res, weighted_edge_count$R1, weighted_edge_count$R2)
    res_list = c(res_list, asy.gen.statistic = gen_asy_p$test_statistic, asy.gen.pval = gen_asy_p$p_value)
  }
  if("wei" %in% test.type){
    #weighted edge count test
    wei_asy_p = asy_wei(asy_res, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2)
    res_list = c(res_list, asy.wei.statistic = wei_asy_p$test_statistic, asy.wei.pval = wei_asy_p$p_value)
  }
  if("max" %in% test.type){
    #max-type edge count test
    max_asy_p = asy_max(asy_res, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2)
    res_list = c(res_list, asy.max.statistic = max_asy_p$test_statistic, asy.max.pval = max_asy_p$p_value)
  }

  if(perm.num != 0 && is.numeric(perm.num)){
    if(perm.num > 10000){
      message("Doing more than 10,000 permutations could be time consuming.\n")
    }
    n_per = perm.num
    # temp_list = permu_edge(n_per, weis, progress_bar)
    temp_list = permu_edge(n_per, E, n1, n2, weis, progress_bar)
    if("ori" %in% test.type){
      #original edge count test
      ori = sum(temp_list$R < weighted_edge_count$R)
      ori = ori/n_per
      res_list = c(res_list, perm.ori.pval =ori)
    }
    if("gen" %in% test.type){
      #generalized edge count test
      gen_p = permu_gen(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n_per)
      gen_p = gen_p/n_per
      res_list = c(res_list, perm.gen.pval = gen_p)
    }
    if("wei" %in% test.type){
      #weighted edge count test
      wei_p = permu_wei(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2, n_per)
      wei_p = wei_p/n_per
      res_list = c(res_list, perm.wei.pval = wei_p)
    }
    if("max" %in% test.type){
      #max-type edge count test
      max_p = permu_max(temp_list$R1, temp_list$R2, weighted_edge_count$R1, weighted_edge_count$R2, n1, n2, n_per)
      max_p = max_p/n_per
      res_list = c(res_list, perm.max.pval = max_p)
    }
    
  }

  if(perm.num != 0 && !is.numeric(perm.num)){
    stop('incorrect perm_num')
  }
  
  return(res_list)
}

