psofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",npar=10,
          vmax=0.7, pso.same=10, c1=0.49, c2=0.49, w.inert='sim.annealing',
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
  v <- lapply(1:npar, function(x) {set.seed(randomN+x+100); matrix(rnorm(ncluster*d, 0,1), ncluster, d)})
  conv <- c(par.fit[which.min(par.fit)])
  pbest <- par$centroid
  pfit <- par$I
  repeat{
  	minmax <- c(which.min(par.fit)[1],which.max(par.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
    theta <- update_inertia(w.inert,wmax,wmin,map,iter,max.iter)
    v <- lapply(1:npar,function(x) update_v(theta,v[[x]],vmax,c1,c2,par.finalpos,pbest[[x]],par.swarm[[x]],randomN+x))
    par.swarm <- lapply(1:npar, function (x) v[[x]] + par.swarm[[x]])
    par.other <- lapply(1:npar, function(x) uij(data,par.swarm[[x]],m,distance,order))
  	par.other <- lapply(1:npar, function(x) renew_uij(data,par.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	par.swarm <- lapply(1:npar, function(x) vi(data,par.other[[x]],m))
  	par.fit <- sapply(1:npar, function(x) jfgwcv(data,par.swarm[[x]],m,distance,order))
    pbest.ind <- which(par.fit<pfit)
    if(length(pbest.ind)>0){
      for(i in pbest.ind){
        pbest[[i]] <- par.swarm[[i]]
        pfit[i] <- par.fit[i]
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
    if (iter==max.iter || same==pso.same) break
  }
  finaldata=determine_cluster(datax,par.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  pso <- list("converg"=conv,"f_obj"=jfgwcv(data,par.finalpos,m,distance,order),"membership"=par.finalpos.other,"centroid"=par.finalpos,
              "validasi"=index_fgwc(data,cluster,par.finalpos.other,par.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"maxgeneration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(pso)
}

update_v <- function(theta,v0,vmax,c1,c2,gbest,pbest,particle,randomN) {
  n <- nrow(particle)
  d <- ncol(particle)
  set.seed(randomN <- randomN+1)
  e1 <- matrix(runif(n*d),ncol=d,nrow=n)
  set.seed(randomN <- randomN+1)
  e2 <- matrix(runif(n*d),ncol=d,nrow=n)
  v_new <- theta*v0+c1*e1*(gbest-particle)+c2*e2*(pbest-particle)
  for(i in 1:ncol(v_new)) {
    x <- which(v_new[,i]<(-vmax))
    v_new[x,i] <- -vmax
    x <- which(v_new[,i]>(vmax))
    v_new[x,i] <- vmax
  }
  return(v_new)
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