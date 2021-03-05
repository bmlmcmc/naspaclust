#' Noise distributions 
#' @description Generating various distributions for noise in optimisation algorithm
#' @param distr the distribution between \code{"normal"}, \code{"uniform"}, \code{"levy"}, \code{logchaotic}, \code{kentchaotic}
#' @param n number of observations.
#' @param randomN random seed for initialisation.
#' @param r weight in logistic chaotic between [0,4]. Can be used when \code{ei.distr='logchaotic'}. Can be defined as \code{chaos} in \code{opt_param}. Default is 4.
#' @param m mapping parameter in kent chaotic between [0,1]. Can be used when \code{ei.distr='kentchaotic'}. Can be defined as \code{map} in \code{opt_param}. Default is 0.7.
#' @param ind Levy distribution index for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{ind} in \code{opt_param}. Default is 1.
#' @param skew Levy distribution skewness for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{skew} in \code{opt_param}. Default is 0.
#' @param sca Levy distribution scale for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{sca} in \code{opt_param}. Default is 1.

eiDist <- function(distr='normal',n,randomN=40,r=4,m=0.7,ind=1,skew=0,sca=1) {
  set.seed(randomN)
  if(distr=='uniform') {
    return(runif(n,-1,1))
  }
  else if(distr=='normal') {
    return(rnorm(n,0,1))
  }
  else if (distr=="levy") {
    # require(stabledist)
    return(rstable(n,ind,skew,sca))
  }
  else if (distr=="logchaotic") {
    return(logchaotic(n,r,randomN))
  }
  else if (distr=="kentchaotic") {
    return(kentchaotic(n,m,randomN))
  }
}

#' @rdname eiDist
#' @param seed the random number

logchaotic <- function(n,r=4,seed=1) {
  set.seed(seed)
  x0 <- runif(1)
  x <- c()
  for(i in 1:n) {
    x0 <- r*x0*(1-x0)
    x <- c(x,x0)
  }
  return(x)
}

#' @rdname eiDist
#' @param seed the random number

kentchaotic <- function(n,m=0.7,seed) {
  set.seed(seed)
  x0 <- runif(1)
  x <- c(x0)
  for(i in 2:n) {
    if(x0>0 && x0<=m) {
      x0 <- x0/m
    }
    else {
      x0 <- (1-x0)/(1-m)
    }
    x <- c(x,x0)
  }
  return(x)
}

#' @rdname eiDist
#' @param alpha the current position of alpha
#' @param iter the current iteration
#' @param maxiter the maximum iteration
#' @param type the update type

update_alpha <- function(alpha, iter, maxiter, type) {
  if(type==1) {
    return(1e-5+(alpha-(1e-5))*exp(-iter))
  }
  else if(type==2) {
    return(alpha*runif(1,0.95,0.99)^iter)
  }
  else if(type==3) {
    delta <- 1-(10^(-4)/9^(1/maxiter))
    return(1-delta*alpha)
  }
  else if (type==4) {
    return((1.11*10^(-4))^(5/maxiter)*alpha)
  }
  else if (type==5){
  	return(alpha*(1-(iter/maxiter)))
  }
}