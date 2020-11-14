gsafgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",npar=10,par.no=2,par.dist='euclidean', par.order=2,
          gsa.same=10, G=1, g.type='sim.annealing',vmax=0.7, pso=F,
          wmax=0.9,wmin=0.4,chaos=4,x0="F",map=0.7,ind=1,skew=0,sca=1){
  require(beepr)
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

  par <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, npar)
  par.swarm <- par$centroid
  par.other <- par$membership
  par.fit <- par$I

  par.finalpos <- par$centroid[[which.min(par.fit)]]
  par.finalpos.other <- par$membership[[which.min(par.fit)]]
  par.fit.finalbest <- par$I[[which.min(par.fit)]]
  v <- lapply(1:npar, function(x) matrix(0, ncluster, d))
  #{set.seed(randomN+x+100); matrix(runif(ncluster*d, 0,1), ncluster, d)}
  pbest <- par$centroid
  pfit <- par$I
  conv <- c(par.fit[which.min(par.fit)])
  repeat{
  	minmax <- c(which.min(par.fit)[1],which.max(par.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
    G <- G*runif(1,0.95,1)
    v <- force_v(par,par.no,G,v,vmax,par.dist,par.order,randomN)
    par.swarm <- lapply(1:npar, function (x) v[[x]] + par.swarm[[x]])
    if(pso==T){
      par.swarm <- lapply(1:npar,function(x) new.move(par.swarm[[x]],pbest[[x]],par.finalpos,randomN+x))
    }
    par.other <- lapply(1:npar, function(x) uij(data,par.swarm[[x]],m,distance,order))
  	par.other <- par$membership <- lapply(1:npar, function(x) renew_uij(data,par.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	par.swarm <- par$centroid <- lapply(1:npar, function(x) vi(data,par.other[[x]],m))
    par.fit <- par$I <- sapply(1:npar, function(x) jfgwcv(data,par.swarm[[x]],m,distance,order))

    if(pso==T){ 
      pbest.ind <- which(par.fit<pfit)
      if(length(pbest.ind)>0){
        for(i in pbest.ind){
          pbest[[i]] <- par.swarm[[i]]
          pfit[i] <- par.fit[i]
        }
      }
    }

    best <- which(par.fit==min(par.fit))[1]
    par.curbest <- par.swarm[[best]]
    par.curbest.other <- par.other[[best]]
    par.fit.curbest <- par.fit[best]
    conv <- c(conv,par.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (par.fit.curbest<=par.fit.finalbest) {
      par.finalpos <- par.curbest
      par.finalpos.other <- par.curbest.other
      par.fit.finalbest <- par.fit.curbest
    }
    randomN <- randomN+npar
    if (iter==max.iter || same==gsa.same) break
  }
  finaldata=determine_cluster(datax,par.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  gsa <- list("converg"=conv,"f_obj"=jfgwcv(data,par.finalpos,m,distance,order),"membership"=par.finalpos.other,"centroid"=par.finalpos,
              "validasi"=index_fgwc(data,cluster,par.finalpos.other,par.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"maxgeneration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(gsa)
}

force_v <- function(par,no,G,v,vmax,par.dist,par.order,randomN){
  dd <- dim(par$centroid[[1]])
  intel.par <- intel.ffly(par,no)
  mass <- (par$I-max(par$I))/(min(par$I)-max(par$I))
  Mass <- mass/sum(mass)
  Mass.intel <- sort(Mass,decreasing=T)[1:no]
  v1 <- v
  for(i in 1:length(par$centroid)){
    Fij <- lapply(1:no,c)
    for(j in 1:no){
      r <- diag(cdist(par$centroid[[i]],intel.par$centroid[[j]],par.dist,par.order))
      set.seed(randomN <- randomN+1)
      eps <- runif(length(r),0,1e-6)
      set.seed(randomN <- randomN+1)
      rand <- matrix(runif(dd[1]*dd[2]),ncol=dd[2])
      Fij[[j]] <- rand*G*Mass[i]*Mass.intel[j]*(intel.par$centroid[[j]]-par$centroid[[i]])/(r+eps)
    }
    Fi <- Reduce("+",Fij)
    a <- Fi/Mass[i]
    set.seed(randomN <- randomN+1)
    rand <- matrix(runif(dd[1]*dd[2]),ncol=dd[2])
    v1[[i]] <- rand*v1[[i]]+a
    # for(i in 1:nrow(v1[[i]])){
    #   v1[[i]][v1[[i]]< -vmax,] <- -vmax
    #   v1[[i]][v1[[i]]> vmax,] <- vmax
    # }
  }
  return(v)
}

new.move <- function(par,pbest,gbest,randomN){ ##Li dan Dong, 2017 GSA new technique
  dd <- dim(par)
  mu <- (par+pbest+gbest)/3
  sigma <- sqrt(((par-mu)^2+(pbest-mu)^2+(gbest-mu)^2)/3)
  set.seed(randomN+100)
  c1 <- matrix(runif(dd[1]*dd[2], 0,1), ncol=dd[2])
  set.seed(randomN+101)
  c2 <- matrix(runif(dd[1]*dd[2], 0,1), ncol=dd[2])
  z <- sqrt(-2*log(c1))*cos(2*pi*c2)
  return(mu+sigma*z)
}

# weight candidate
# constant, chaotic, simulated annealing, natural exp 1 2, exp decreasing
update_inertia <- function(w.inert, wmax, wmin, z, iter, maxiter) {
  if(w.inert=="constant") {
    return(wmax)
  }
  else if(w.inert=="chaotic") {
    return((wmax-wmin)*(1-iter/maxiter)+(wmax*z))
  }
  else if(w.inert=="sim.annealing") {
    return(wmin+((wmax-wmin)*0.95^(iter-1)))
  }
  else if(w.inert=="nat.exponent1") {
    return(wmin+((wmax-wmin)*exp(iter/(maxiter/10))))
  }
  else if(w.inert=="nat.exponent2") {
    return(wmin+((wmax-wmin)*exp((iter/(maxiter/10)^2))))
  }
}