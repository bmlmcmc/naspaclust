#' Fuzzy Geographicaly Weighted Clustering with Artificial Bee Colony Optimization
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Artificial Bee Colony
#' @param data an object of data with d>1. Can be \code{matrix} or \code{data.frame}. If your data is univariate, bind it with \code{1} to get a 2 columns.
#' @param pop an n*1 vector contains population.
#' @param distmat an n*n distance matrix between regions.
#' @param ncluster an integer. The number of clusters.
#' @param m degree of fuzziness or fuzzifier. Default is 2.
#' @param distance the distance metric between data and centroid, the default is euclidean, see \code{\link{cdist}} for details.
#' @param order, minkowski order. default is 2.
#' @param alpha the old membership effect with [0,1], if \code{alpha} equals 1, it will be same as fuzzy C-Means, if 0, it equals to neighborhood effect.
#' @param a spatial magnitude of distance. Default is 1.
#' @param b spatial magnitude of population. Default is 1.
#' @param max.iter maximum iteration. Default is 500.
#' @param error error tolerance. Default is 1e-5.
#' @param randomN random seed for initialisation (if uij or vi is NA). Default is 0.
#' @param vi.dist a string of centroid population distribution between \code{"uniform"} (default) and \code{"normal"}. Can be defined as \code{vi.dist=} in \code{opt_param}.
#' @param nant number of ants population. Can be defined as \code{npar=} in \code{opt_param}.
#' @param n.onlooker number of onlooker bees, Can be defined as \code{n.onlooker} in \code{opt_param}.
#' @param limit number of turns to eliminate ant with no solutions. Can be defined as \code{limit} in \code{opt_param}.
#' @param pso whether to add PSO term in bee's movement. Either \code{TRUE} or \code{FALSE}. Can be defined as \code{pso} in \code{opt_param}.
#' @param aco.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.

#' @return an object of class \code{"fgwc"}.\cr
#' An \code{"fgwc"} object contains as follows:
#' \itemize{
#' \item \code{converg} - the process convergence of objective function
#' \item \code{f_obj} - objective function value
#' \item \code{membership} - membership matrix
#' \item \code{centroid} - centroid matrix
#' \item \code{validation} - validation indices (there are partition coefficient (\code{PC}), classification entropy (\code{CE}), 
#' SC index (\code{SC}), separation index (\code{SI}), Xie and Beni's index (\code{XB}), IFV index (\code{IFV}), and Kwon index (Kwon))
#' \item \code{max.iter} - Maximum iteration
#' \item \code{cluster} - the cluster of the data
#' \item \code{finaldata} - The final data (with the cluster)
#' \item \code{call} - the syntax called previously
#' \item \code{time} - computational time.
#' }
#' @details Fuzzy Geographically Weighted Clustering (FGWC) was developed by \insertCite{fgwc;textual}{naspaclust} by adding 
#' neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. Furthermore,
#' the Artificial Bee Colony (ACO) algorithm was developed by \insertCite{Karaboga2007;textual}{naspaclust} in order to get a more optimal
#' solution of a certain complex function. FGWC using ACO has been implemented previously by \insertCite{fgwcaco1;textual}{naspaclust} and \insertCite{fgwcaco2;textual}{naspaclust}.

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- acofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',4,nant=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the ACO parameter
#' aco_param <- c(vi.dist='normal',npar=5,pso=FALSE,same=15,n.onlooker=5,limit=5) 
#' ##FGWC with ACO optimization algorithm
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'aco',param_fgwc,aco_param) 

#' @export

acofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100, randomN=0, vi.dist="uniform", nant=10, nintervals=10, aco.same=10, algorithm='as', ant_alpha=2, ant_beta=2, q0=0.9, 
          rho=0.2, Q=1, p0=0.5){
  # require(beepr)
  randomnn <- randomN
  ptm<-proc.time()
  n <- nrow(data)
  d <- ncol(data)
  iter=0
  beta <- 1-alpha
  same=0
  data <- as.matrix(data)
  if (alpha ==1) {
    pop <- rep(1,n)
    distmat <- matrix(1,n,n)
  }
  datax <- data
  pop <- matrix(pop,ncol=1)
  mi.mj <- pop%*%t(pop)
  minmaxdata <- rbind(apply(data,2,min),apply(data,2,max))
  ant <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, nant)
  ant.swarm <- ant$centroid
  ant.other <- ant$membership
  ant.fit <- ant$I

  ant.finalpos <- ant$centroid[[which.min(ant.fit)]]
  ant.finalpos.other <- ant$membership[[which.min(ant.fit)]]
  ant.fit.finalbest <- ant$I[[which.min(ant.fit)]]
  conv <- c(ant.fit[which.min(ant.fit)])
  t <- rep(0,nant)
  repeat{
  	minmax <- c(which.min(ant.fit)[1],which.max(ant.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
    ant.swarm <- lapply(1:nant, function(x) global_update(data,ant.swarm[[x]],m,distance,order))
  
    ant.other <- lapply(1:nant, function(x) uij(data,ant.swarm[[x]],m,distance,order))
  	ant.other <- lapply(1:nant, function(x) renew_uij(data,ant.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	ant.swarm <- lapply(1:nant, function(x) vi(data,ant.other[[x]],m))
  	ant.fit <- sapply(1:nant, function(x) jfgwcv(data,ant.swarm[[x]],m,distance,order))
    best <- which(ant.fit==min(ant.fit))[1]
    ant.curbest <- ant.swarm[[best]]
    ant.curbest.other <- ant.other[[best]]
    ant.fit.curbest <- ant.fit[best]
    conv <- c(conv,ant.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (ant.fit.curbest<=ant.fit.finalbest) {
      ant.finalpos <- ant.curbest
      ant.finalpos.other <- ant.curbest.other
      ant.fit.finalbest <- ant.fit.curbest
    }
    randomN <- randomN+nant
    if (iter==max.iter || same==aco.same) break
  }
  finaldata=determine_cluster(datax,ant.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  aco <- list("converg"=conv,"f_obj"=jfgwcv(data,ant.finalpos,m,distance,order),"membership"=ant.finalpos.other,"centroid"=ant.finalpos,
              "validation"=index_fgwc(data,cluster,ant.finalpos.other,ant.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(aco)
}

divide_data <- function(data,nintervals){
  mindata <- apply(data,2,min)
  maxdata <- apply(data,2,max)
  res <- c()
  for(i in 1:ncol(data)){
    res <- cbind(res,seq(mindata[i],maxdata[i],(maxdata[i]-mindata[i])/nintervals))
  }
  return(res)
}

state_trans <- function(ant, fitness, algorithm, rho, best, worst, p0){

}

global_up <- function(){

}

local_up <- function(){

}