abcfgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100, randomN=0, vi.dist="uniform", nfood=10, n.onlooker=5, limit=4, pso=F, abc.same=10){
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
  minmaxdata <- rbind(apply(data,2,min),apply(data,2,max))
  food <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, nfood)
  food.swarm <- food$centroid
  food.other <- food$membership
  food.fit <- food$I

  food.finalpos <- food$centroid[[which.min(food.fit)]]
  food.finalpos.other <- food$membership[[which.min(food.fit)]]
  food.fit.finalbest <- food$I[[which.min(food.fit)]]
  conv <- c(food.fit[which.min(food.fit)])
  t <- rep(0,nfood)
  repeat{
  	minmax <- c(which.min(food.fit)[1],which.max(food.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
    candfood <- employed.bee(food.swarm,food.fit,pso,food.finalpos,randomN,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)
    candfood.fit <- sapply(1:nfood, function(x) jfgwcv2(data,candfood[[x]],m,distance,order,mi.mj,distmat,alpha,beta,a,b))
    newfood <- compare(candfood,food.swarm,candfood.fit,food.fit,t,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)
    onlook.food <- onlooker.bee(newfood$swarm,newfood$fit,newfood$t,newfood$prob,n.onlooker,pso,food.finalpos,
      randomN+1,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)
    t <- onlook.food$t
    food.swarm <- scout.bee(onlook.food$swarm,t,limit,minmaxdata,randomN+2)
    food.other <- lapply(1:nfood, function(x) uij(data,food.swarm[[x]],m,distance,order))
  	food.other <- lapply(1:nfood, function(x) renew_uij(data,food.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	food.swarm <- lapply(1:nfood, function(x) vi(data,food.other[[x]],m))
  	food.fit <- sapply(1:nfood, function(x) jfgwcv(data,food.swarm[[x]],m,distance,order))
    best <- which(food.fit==min(food.fit))[1]
    food.curbest <- food.swarm[[best]]
    food.curbest.other <- food.other[[best]]
    food.fit.curbest <- food.fit[best]
    conv <- c(conv,food.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (food.fit.curbest<=food.fit.finalbest) {
      food.finalpos <- food.curbest
      food.finalpos.other <- food.curbest.other
      food.fit.finalbest <- food.fit.curbest
    }
    randomN <- randomN+nfood
    if (iter==max.iter || same==abc.same) break
  }
  finaldata=determine_cluster(datax,food.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  abc <- list("converg"=conv,"f_obj"=jfgwcv(data,food.finalpos,m,distance,order),"membership"=food.finalpos.other,"centroid"=food.finalpos,
              "validasi"=index_fgwc(data,cluster,food.finalpos.other,food.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"maxgeneration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(abc)
}


employed.bee <- function(swarm,fitness,pso,gbest,seed,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
  real <- 1:length(swarm)
  set.seed(seed <- seed+2)
  sample1 <- sample(1:length(swarm),length(swarm))
  while(sum(real==sample1)!=0){
    ind <- which(real==sample1)
    set.seed(seed <- seed+1)
    sample1[ind] <- sample(1:length(swarm),length(ind))
  }
  swarm <- lapply(1:length(swarm), function(x){
    set.seed(seed+10+x)
    phi <- matrix(runif(ncol(swarm[[x]])*nrow(swarm[[x]]),-1,1),ncol=ncol(swarm[[x]]))
    a <- sample1[x]
    psi <- 0
    if(pso==T){
      set.seed(seed+10+x)
      psi <- matrix(runif(ncol(swarm[[x]])*nrow(swarm[[x]]),0,1.5),ncol=ncol(swarm[[x]]))
    }
    swarm[[x]]+phi*(swarm[[x]]-swarm[[a]])+psi*(gbest-swarm[[x]])
  })
  return(swarm)
}

onlooker.bee <- function(swarm,fit,t,prob,n.onlooker,pso,gbest,seed,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
  real <- order(prob,decreasing=T)[1:n.onlooker]
  set.seed(seed <- seed+2)
  sample1 <- sample(1:length(swarm),n.onlooker)
  while(sum(real==sample1)!=0){
    ind <- which(real==sample1)
    set.seed(seed <- seed+1)
    sample1[ind] <- sample(1:length(swarm),length(ind))
  }
  oldswarm <- lapply(real, function(x) swarm[[x]])
  oldfit <- fit[real]
  t2 <- t[real]
  newswarm <- lapply(1:length(real), function(x){
    set.seed(seed+10+x)
    phi <- matrix(runif(ncol(swarm[[x]])*nrow(swarm[[x]]),-1,1),ncol=ncol(swarm[[x]]))
    a <- sample1[x]
    psi <- 0
    if(pso==T){
      set.seed(seed+10+x)
      psi <- matrix(runif(ncol(oldswarm[[x]])*nrow(oldswarm[[x]]),0,1.5),ncol=ncol(oldswarm[[x]]))
    }
    oldswarm[[x]]+phi*(oldswarm[[x]]-swarm[[a]])+psi*(gbest-oldswarm[[x]])
  })
  newfit <- sapply(1:length(newswarm),function(x) jfgwcv2(data,newswarm[[x]],m,distance,order,mi.mj,dist,alpha,beta,a,b))
  swarmlast <- compare(newswarm,oldswarm,newfit,oldfit,t2,data,m,distance,order,mi.mj,dist,alpha,beta,a,b)
  j = 0
  for(i in real){
    j <- j+1
    swarm[[i]] <- swarmlast$swarm[[j]]
    fit[i] <- swarmlast$fit[j]
    t[i] <- swarmlast$t[j]
  }
  return(list(swarm=swarm,fit=fit,t=t))
}

scout.bee <- function(swarm,t,limit,minmaxdata,seed){
  ind <- which(t==limit)
  for(i in ind){
    set.seed(seed <- seed+5+i)
    r <- matrix(runif(ncol(swarm[[i]])*nrow(swarm[[i]])),ncol=ncol(swarm[[i]]))
    swarm[[i]] <- minmaxdata[1,]+r*(minmaxdata[2,]-minmaxdata[1,])
  }
  return(swarm)
}

compare <- function(newswarm,oldswarm,newfit,oldfit,t,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
  ind <- which(newfit<oldfit)
  for(i in ind){
    oldswarm[[i]] <- newswarm[[i]]
    oldfit[i] <- newfit[i]
    t[i] <- 0
  }
  t[-ind] <- t[-ind]+1
  obj <- sapply(1:length(oldswarm),function(x) jfgwcv2(data,oldswarm[[x]],m,distance,order,mi.mj,dist,alpha,beta,a,b))
  prob <- (1/obj)/sum(1/obj)
  return(list(swarm=oldswarm,fit=obj,prob=prob,t=t))
}