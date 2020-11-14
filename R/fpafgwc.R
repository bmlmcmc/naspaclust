fpafgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1,
					error=1e-5, max.iter=100, randomN=0, vi.dist="uniform", nflow=10, p=0.8, gamma=1, lambda=1.5, delta=0,
          ei.distr='normal', flow.same=10,flow.maxiter=100,r=4,m.chaotic=0.7,skew=0,sca=1){
  require(stabledist)
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

  flow <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster,
                                m, alpha, a, b, randomN, nflow)
  flow.swarm <- flow$centroid
  flow.other <- flow$membership
  flow.fit <- flow$I

  flow.finalpos <- flow$centroid[[which.min(flow.fit)]]
  flow.finalpos.other <- flow$membership[[which.min(flow.fit)]]
  flow.fit.finalbest <- flow$I[[which.min(flow.fit)]]
  conv <- c(flow.fit[which.min(flow.fit)])
  repeat{
  	minmax <- c(which.min(flow.fit)[1],which.max(flow.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
  	pollen <- flow.finalpos
    flow.swarm <- pollination(flow.swarm,p,pollen,gamma,lambda,delta,randomN,ei.distr,r,m.chaotic,ind,skew,sca)
    flow.other <- lapply(1:nflow, function(x) uij(data,flow.swarm[[x]],m,distance,order))
  	flow.other <- lapply(1:nflow, function(x) renew_uij(data,flow.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	flow.swarm <- lapply(1:nflow, function(x) vi(data,flow.other[[x]],m))
  	flow.fit <- sapply(1:nflow, function(x) jfgwcv(data,flow.swarm[[x]],m,distance,order))
    best <- which(flow.fit==min(flow.fit))[1]
    flow.curbest <- flow.swarm[[best]]
    flow.curbest.other <- flow.other[[best]]
    flow.fit.curbest <- flow.fit[best]
    conv <- c(conv,flow.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (flow.fit.curbest<=flow.fit.finalbest) {
      flow.finalpos <- flow.curbest
      flow.finalpos.other <- flow.curbest.other
      flow.fit.finalbest <- flow.fit.curbest
    }
    randomN <- randomN+nflow
    if (iter==max.iter || same==flow.same) break
  }
  finaldata=determine_cluster(datax,flow.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  fpa <- list("converg"=conv,"f_obj"=jfgwcv(data,flow.finalpos,m,distance,order),"membership"=flow.finalpos.other,"centroid"=flow.finalpos,
              "validasi"=index_fgwc(data,cluster,flow.finalpos.other,flow.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"maxgeneration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(fpa)
}

pollination <- function(flow,p,pollen,gamma,lambda,delta,seed,ei.distr,r,m,ind,skew,sca){
  set.seed(seed<-seed+10)
  rand <- runif(length(flow))
  dd <- dim(pollen)
  return(lapply(1:length(flow),function(x){
    if(rand[x]<p){ ##global pollination
      flow[[x]]+gamma*matrix(rstable(dd[1]*dd[2],lambda,skew,sca,delta),ncol=dd[2])*(pollen-flow[[x]])
    }
    else{ ##local pollination
      ei <- matrix(eiDist(ei.distr,dd[1]*dd[2],seed+x,r,m,lambda,skew,sca),ncol=dd[2])
      no <- 1:length(flow)
      sample <- sample(no[-x],2)
      a = sample[1]
      b = sample[2]
      flow[[x]]+ei*(flow[[a]]-flow[[b]])
    }
  }))
}
