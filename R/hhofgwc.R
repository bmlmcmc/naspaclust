hhofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",nhh=10,hh.alg='heidari',
          A=c(2,1,0.5),p=0.5,hh.maxiter=200,hh.same=10,levy.beta=1.5,update.type=5){
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

  hh <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, nhh)
  hh.swarm <- hh$centroid
  hh.other <- hh$membership
  hh.fit <- hh$I

  hh.finalpos <- hh$centroid[[which.min(hh.fit)]]
  hh.finalpos.other <- hh$membership[[which.min(hh.fit)]]
  hh.fit.finalbest <- hh$I[[which.min(hh.fit)]]

  conv <- c(hh.fit[which.min(hh.fit)])
  repeat{
  	minmax <- c(which.min(hh.fit)[1],which.max(hh.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
  	rabbit <- hh.finalpos
    sample <- sample(1:nhh,nhh,T)
    set.seed(randomN+iter)
    if(tolower(hh.alg) == 'bairathi'){
      A[1] <- update_alpha(A[1],iter,max.iter,update.type)
      hh.swarm <- lapply(1:nhh, function(x) hh.attack.bairathi(hh.swarm[[x]],hh.swarm,hh.finalpos,A,p,sample[x],
        randomN+x,best))
    }
  	else {
      E0 <- A[1]*(2*runif(nhh)-1)
      E <- update_alpha(E0,iter,max.iter,update.type)
      hh.swarm <- lapply(1:nhh, function(x) hh.attack.heidari(hh.swarm[[x]],hh.swarm,rabbit,E[x],A,p,
      	sample[x],levy.beta,randomN+x,best,worst,hh.fit[x],data,m,distance,order,mi.mj,distmat,alpha,beta,a,b))
    }
    set.seed(randomN+iter)
  	hh.other <- lapply(1:nhh, function(x) uij(data,hh.swarm[[x]],m,distance,order))
  	hh.other <- lapply(1:nhh, function(x) renew_uij(data,hh.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	hh.swarm <- lapply(1:nhh, function(x) vi(data,hh.other[[x]],m))
  	hh.fit <- sapply(1:nhh, function(x) jfgwcv(data,hh.swarm[[x]],m,distance,order))
    best <- which(hh.fit==min(hh.fit))[1]
    hh.curbest <- hh.swarm[[best]]
    hh.curbest.other <- hh.other[[best]]
    hh.fit.curbest <- hh.fit[best]
    conv <- c(conv,hh.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (hh.fit.curbest<=hh.fit.finalbest) {
      hh.finalpos <- hh.curbest
      hh.finalpos.other <- hh.curbest.other
      hh.fit.finalbest <- hh.fit.curbest
    }
    randomN <- randomN+nhh
    if (iter==max.iter || same==hh.same) break
  }
  finaldata=determine_cluster(datax,hh.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  hho <- list("converg"=conv,"f_obj"=jfgwcv(data,hh.finalpos,m,distance,order),"membership"=hh.finalpos.other,"centroid"=hh.finalpos,
              "validasi"=index_fgwc(data,cluster,hh.finalpos.other,hh.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"maxgeneration"=iter,"same"=same,"time"=proc.time()-ptm)
  return(hho)
}

hh.attack.bairathi <- function(hawk,hawks,gbest,A,p,rand,seed,best){
  dd <- dim(hawk)
  set.seed(seed<-seed+10)
  rs <- runif(1)
  set.seed(seed<-seed+10)
  c1 <- A[1]*(2*runif(1)-1)
  set.seed(seed)
  c2 <- matrix(2*(1-runif(dd[1]*dd[2])),ncol=dd[2])
  if(rs<p && c1>=A[2]) return(hawks[[rand]]-c1*(c2*hawks[[rand]]-hawk)) ##exploration
  else if(c1>=A[3]) return(hawks[[best]]-c1*(c2*hawks[[best]]-hawk)) ##local exploitation
  else return((gbest-hawk)-c1*(c2*gbest-hawk)) ##global exploitation
}


hh.attack.heidari <- function(hawk,hawks,rabbit,E,A,p,rand,levy.beta,seed,best,worst,fithawk,data,m,
	distance,order,mi.mj,dist,alpha,beta,a,b){
	set.seed(seed<-seed+10)
	rs <- runif(3)
	dd <- dim(hawk)
	hawks.m <- Reduce('+',hawks)/length(hawks)
	if(E>=A[2]){ ##exploration phase
		if(rs[1]>=p){ ##q
			return(hawks[[rand]]-rs[2]*(hawks[[rand]]-2*rs[3]*hawk))
		}
		else{
			return(rabbit-hawks.m-rs[2]*(hawks[[worst]]+rs[3]*(hawks[[best]]-hawks[[worst]])))
		}
	}
	else{ ##exploitation phase
		set.seed(seed<-seed+10)
		J <- matrix(2*(1-runif(dd[1]*dd[2])),ncol=dd[2])
		if(rs[1]>=p){ 
			if(E>=A[3]) return((rabbit-hawk)-E*(J*rabbit-hawk)) ##soft besiege
			else return(rabbit-E*(rabbit-hawk)) ##hard besiege
		}
		else{ 
			S <- matrix(rnorm(dd[1]*dd[2]),ncol=dd[2])
			LF <- matrix(rlevy(dd[1]*dd[2],levy.beta,seed),ncol=dd[2])
			if(E>=A[3]){ ##soft besiege with progressive rapid dives
				y <- rabbit-E*(J*rabbit-hawk)
			}
			else{ ##hard besiege with progressive rapid dives
				y <- rabbit-E*(J*rabbit-hawks.m)
			}
			z <- y+S*LF
			fity <- jfgwcv2(data,y,m,distance,order,mi.mj,dist,alpha,beta,a,b)
			fitz <- jfgwcv2(data,z,m,distance,order,mi.mj,dist,alpha,beta,a,b)
			if(fity<fithawk) return(y)
			else if(fitz<fithawk) return(z)
			else return(hawk)
		}
	}	
}

rlevy <- function(n,beta,seed){
  set.seed(seed+100)
  u <- runif(n)
  set.seed(seed+101)
  v <- runif(n)
  sigma1 <- gamma(1+beta)*sin(pi*beta/2)
  sigma2 <- gamma((1+beta)/2)*beta*2^((beta-1)/2)
  return(u*((sigma1/sigma2)/v)^(1/beta))
}