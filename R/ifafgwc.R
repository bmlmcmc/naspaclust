########################################################
#############INTELLIGENT FIREFLY ALGORITHM##############
########################################################
ifafgwc <- function (data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1,
						error=1e-5, max.iter=100,randomN=0,vi.dist="uniform", ei.distr="normal",
						fa.same=10, nfly=10, ffly.no=2, ffly.dist='euclidean', ffly.order=2, gamma=1, ffly.beta=1,
            ffly.alpha=1, r.chaotic=4,m.chaotic=0.7,ind.levy=1,skew.levy=0,scale.levy=1,ffly.alpha.type=4) {
  #require(beepr)
  randomnn <- randomN
  ptm<-proc.time()
  n <- nrow(data)
  d <- ncol(data)
  iter=0
  gen=1
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
  ffly.finalbest <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster,
                                m, alpha, a, b, randomN, 1)
  inten.finalbest <- ffly.finalbest$I
  conv <- c(inten.finalbest)
  ffly.finalpos <- ffly.finalbest$centroid
  ffly.finalpos.other <- ffly.finalbest$membership
  ffly.new <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster,
                        m, alpha, a, b, randomN, nfly)
  ffly.swarm <- ffly.new$centroid
  ffly.other <- ffly.new$membership
  inten <- ffly.new$I
  repeat {
    set.seed(randomN)
    ffly.alpha <- update_alpha(ffly.alpha,gen,max.iter,ffly.alpha.type)
    set.seed(randomN)
    ffly.swarm <- ffly.new$centroid <- moving(ffly.new,ffly.no,ffly.beta,gamma,ffly.alpha,ffly.dist,ffly.order,ei.distr,
      r.chaotic,m.chaotic,ind.levy,skew.levy,sca.levy,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b,randomN)
    ffly.other <- ffly.new$membership <- lapply(1:nfly, function(x) uij(data,ffly.swarm[[x]],m,distance,order))
    ffly.other <- ffly.new$membership <- lapply(1:nfly, function(x) renew_uij(data,ffly.new$membership[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
    ffly.swarm <- ffly.new$centroid <- lapply(1:nfly, function(x) vi(data,ffly.new$membership[[x]],m))
    inten <- ffly.new$I <- sapply(1:nfly, function(x) jfgwcv(data,ffly.new$centroid[[x]],m,distance,order))
    best <- which(inten==min(inten))[1]
    ffly.curbest <- ffly.swarm[[best]]
    ffly.curbest.other <- ffly.other[[best]]
    inten.curbest <- inten[best]
    conv <- c(conv,inten.finalbest)
    if (abs(conv[gen+1]-conv[gen])<error) {
      same <- same+1
    }
    else {
      same <- 0
    }
    if (inten.curbest<=inten.finalbest) {
      ffly.finalpos <- ffly.curbest
      ffly.finalpos.other <- ffly.curbest.other
      inten.finalbest <- inten.curbest
    }
    gen <- gen+1
    randomN <- randomN+nfly
    if (gen==max.iter || same==fa.same) break
  }##end repeat
  print(class(ffly.finalpos.other))
  if (class(ffly.finalpos.other)=="list") {
  	new_uij <- ffly.finalpos.other[[1]]
  	vi <- ffly.finalpos[[1]]
  }
  else {
  	new_uij <- ffly.finalpos.other
  	vi <- ffly.finalpos
  }
  finaldata=determine_cluster(datax,new_uij)
  cluster=finaldata[,ncol(finaldata)]
  ifa <- list("converg"=conv,"f_obj"=jfgwcv(data,vi,m,distance,order),"membership"=new_uij,"centroid"=vi,
              "validation"=index_fgwc(data,cluster,new_uij,vi,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=gen,"same"=same,"time"=proc.time()-ptm)
  print(c(order, ncluster,m, randomN))
  #result <- list(ifa=ifa,fgwc=fgwc)
  class(ifa) <- 'fgwc'
  return (ifa)
}

init.swarm <- function(data, pop, distmat, distance, order, vi.dist, ncluster,
                        m, alpha, a, b, randomN, nfly) {
  inten <- rep(0,nfly)
  beta <- 1-alpha
  start.uij <- lapply(1:nfly, function(x) gen_uij(data,ncluster,nrow(data),randomN+x))
	start.uij <- lapply(1:nfly, function (x) renew_uij(data,start.uij[[x]],pop,distmat,alpha,beta,a,b))
	start.vi <- lapply(1:nfly, function (x) vi(data,start.uij[[x]],m))
  for(i in 1:nfly) {
    inten[i] <- jfgwcv(data,start.vi[[i]],m,distance,order)
  }
  result <- list("membership"=start.uij,"centroid"=start.vi,"I"=inten)
  return(result)
}

intel.ffly <- function(ffly.list,no) {
  best <- order(ffly.list$I,decreasing = F)[1:no]
  intel.uij <- lapply(best, function(x) ffly.list$membership[[x]])
  intel.vi <- lapply(best, function(x) ffly.list$centroid[[x]])
  inten <- ffly.list$I[best]
  result <- list("membership"=intel.uij,"centroid"=intel.vi,"I"=inten)
  return(result)
}

swarm_dist <- function (swarm1,swarm2,distance,order) {
  # jarak<-rep(0,nrow(swarm1))
  # for (i in 1:nrow(swarm1)) {
  #     diff <- abs(swarm1[i,]-swarm2[i,])^order
  #     jarak[i] <- sum(diff)^(1/order)
  # }
  return(diag(dist(swarm1,swarm2,distance,order)))
}

moving <- function(ffly.all,no,ff.beta,gamma,ff.alpha,ffly.dist,ffly.order,ei.distr,r.chaotic,m.chaotic,ind.levy,skew.levy,sca.levy,
                  data,m,distance,order,mi.mj,dist,alpha,beta,a,b,randomN){##menggerakkan firefly
  times <- 0
  intel.ffly <- intel.ffly(ffly.all,no)
  dd <- dim(ffly.all$centroid[[1]])
  ffly <- ffly.all$centroid
  fit <- ffly.all$I
  for(i in 1:length(intel.ffly$centroid)){
    for(j in 1:length(ffly)){
      r <- diag(cdist(ffly[[j]],intel.ffly$centroid[[i]],ffly.dist,ffly.order))
      ei <- matrix(eiDist(ei.distr,dd[1]*dd[2],randomN+i+j,r.chaotic,m.chaotic,ind.levy,skew.levy,sca.levy),ncol=dd[2])
      if (fit[j] > intel.ffly$I[i]){
        ffly[[j]]+beta*exp(-gamma*r^2)*(intel.ffly$centroid[[i]]-ffly[[j]])+(ff.alpha*ei)
      }
      else{
        times <- times+1
        if(times==no){
          ffly[[j]]+(ff.alpha*ei)
        }
      }
      fit[j] <- jfgwcv2(data,ffly[[j]],m,distance,order,mi.mj,dist,alpha,beta,a,b)
    }
  }
  return(ffly)
}
