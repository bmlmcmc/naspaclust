tlbofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",nstud=10,vmax=0.4, tlbo.same=10,
          nselection=10,elitism=F,n.elite=2){
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

  stud <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, nstud)
  stud.swarm <- stud$centroid
  stud.other <- stud$membership
  stud.fit <- stud$I

  stud.finalpos <- stud$centroid[[which.min(stud.fit)]]
  stud.finalpos.other <- stud$membership[[which.min(stud.fit)]]
  stud.fit.finalbest <- stud$I[[which.min(stud.fit)]]
  conv <- c(stud.fit[which.min(stud.fit)])
  repeat{
  	minmax <- c(which.min(stud.fit)[1],which.max(stud.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
  	teacher <- stud$centroid[[best]]
    class.ave <- Reduce('+',stud.swarm)/nstud
    studs.new <- teacher.phase(stud.swarm,stud.fit,teacher,class.ave,randomN,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)
    studs.new <- learner.phase(studs.new$studs,studs.new$fit,randomN+5,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)

    if(elitism==T){
    	stud.swarm <- elitism(studs.new$studs,studs.new$fit,stud.finalpos,stud.fit.finalbest,n.elite,randomN+6)
    }
    else{
    	stud.swarm <- studs.new$studs
    }
    stud.other <- lapply(1:nstud, function(x) uij(data,stud.swarm[[x]],m,distance,order))
  	stud.other <- lapply(1:nstud, function(x) renew_uij(data,stud.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	stud.swarm <- lapply(1:nstud, function(x) vi(data,stud.other[[x]],m))
  	stud.fit <- sapply(1:nstud, function(x) jfgwcv(data,stud.swarm[[x]],m,distance,order))
    best <- which(stud.fit==min(stud.fit))[1]
    stud.curbest <- stud.swarm[[best]]
    stud.curbest.other <- stud.other[[best]]
    stud.fit.curbest <- stud.fit[best]
    conv <- c(conv,stud.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (stud.fit.curbest<=stud.fit.finalbest) {
      stud.finalpos <- stud.curbest
      stud.finalpos.other <- stud.curbest.other
      stud.fit.finalbest <- stud.fit.curbest
    }
    randomN <- randomN+nstud
    if (iter==max.iter || same==tlbo.same) break
  }
  finaldata=determine_cluster(datax,stud.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  tlbo <- list("converg"=conv,"f_obj"=jfgwcv(data,stud.finalpos,m,distance,order),"membership"=stud.finalpos.other,"centroid"=stud.finalpos,
              "validasi"=index_fgwc(data,cluster,stud.finalpos.other,stud.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"maxgeneration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(tlbo)
}

teacher.phase <- function(studs,studs.fit,teacher,average,seed,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
  set.seed(seed <- seed+2)
  tf <- matrix(1+runif(ncol(teacher)*nrow(teacher)),ncol=ncol(teacher))
  set.seed(seed <- seed+1)
  r <- matrix(runif(ncol(teacher)*nrow(teacher)),ncol=ncol(teacher))
  diff <- r*(teacher-tf*average)
  nstud <- length(studs)
  stud2 <- lapply(1:nstud, function (x) studs[[x]]+diff)
  stud.fit2 <- sapply(1:nstud, function (x) jfgwcv2(data,stud2[[x]],m,distance,order,mi.mj,dist,alpha,beta,a,b))
  whichone <- which(stud.fit2<studs.fit)
  for(i in whichone){
    studs[[i]] <- stud2[[i]]
    studs.fit[i] <- stud.fit2[i]
  }
  return(list(studs=studs,fit=studs.fit))
}

learner.phase <- function(studs,studs.fit,seed,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
	set.seed(seed <- seed+2)
	sample1 <- sample(1:length(studs),length(studs))
	set.seed(seed <- seed+2)
	sample2 <- sample(1:length(studs),length(studs))
	for(i in 1:length(studs)){
		set.seed(seed <- seed+1)
		r <- matrix(runif(ncol(studs[[1]])*nrow(studs[[1]])),ncol=ncol(studs[[1]]))	
		if(sample1[i]!=sample2[i]){
	  	a <- sample1[i]
	    b <- sample2[i]
	    if(studs.fit[a]<studs.fit[b]) stud2 <- studs[[a]]+r*(studs[[a]]-studs[[b]])
		  else stud2 <- studs[[a]]+r*(studs[[b]]-studs[[a]])
		  stud.fit2 <- jfgwcv2(data,stud2,m,distance,order,mi.mj,dist,alpha,beta,a,b)
			if(stud.fit2<studs.fit[a]){
				studs[[a]] <- stud2
    		studs.fit[a] <- stud.fit2
			}
	  }
	}
	return(list(studs=studs,fit=studs.fit))
}

elitism <- function(swarm,fit,gbest,gfit,n.elite,seed){
	worst <- order(fit,decreasing = T)[1:n.elite]
	for(i in worst){
		swarm[[i]] <- gbest
		fit[i] <- gfit
	}
	dup <- which(duplicated(fit))
	for(i in dup){
		set.seed(seed <- seed+i)
		r <- matrix(rnorm(ncol(swarm[[1]])*nrow(swarm[[1]])),ncol=ncol(swarm[[1]]))
		swarm[[i]] <- swarm[[i]]+r
	}
	return(swarm)
}