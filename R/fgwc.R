########################################################
#####################CLASSICAL FGWC#####################
########################################################
library(rdist)
fgwcuv <- function(data, pop, distmat, kind=NA,ncluster=2, m=2, distance='euclidean', order=2,
                      alpha=0.7, a=1, b=1, max.iter=500, error=1e-5,
                     randomN=0, uij=NA, vi=NA) {
  ##populasi berupa matriks n x 1
  ##jarak berupa matriks jarak antar daerah
  ##alpha + beta = 1
  ##m = fuzzifier
  ptm <- proc.time()
  stopifnot(kind=="v" || kind=="u" || any(is.na(kind))==T)
  n <- nrow(data)
  d <- ncol(data)
  beta = 1-alpha
  iter = 0
  conv <- c(0)
  if (is.matrix(data)==F) {
    data <- as.matrix(data)
  }
  ##jika alfa =1, akan menjadi fuzzy c-means,
  ##populasi dan matriks jarak dianggap matriks 1
  if (alpha==1) {
    pop <- rep(1,n)
    distmat <- matrix(1,n,n)
  }
  ##membaca matriks populasi dan menjadikannya matriks dengan 1 kolom
  pop <- matrix(pop,ncol=1)
  mi.mj <- pop%*%t(pop)
  ##membaca pendekatan FGWC jika dikosongkan atau u,
  ##pendekatannya akan menjadi matriks keanggotaan, inisialisasi awal matriks keanggotaan
  if (any(is.na(kind))==T || kind=="u"){ ##fgwc biasa = fgwc u
    if (any(is.na(uij))==T) {
      set.seed(randomN)
      uij <- matrix(runif(n*ncluster,0,1),ncol=ncluster)
      new_uij  <- uij/rowSums(uij)
    }
    else {
      new_uij <- uij
    }
    old_uij <- new_uij+1
    while (max(abs(new_uij-old_uij))>error && iter<max.iter) {
      old_uij <- new_uij
      vi <- vi (data,old_uij,m) ##mengubah matriks keanggotaan menjadi centroid
      uij <- uij(data,vi,m,distance,order) ##centroid yang didapat diubah menjadi matriks keanggotaan
      new_uij <- renew_uij(data,uij$u,mi.mj,distmat,alpha,beta,a,b) ##memodifikasi matriks keanggotaan dengan memanfaatkan matriks jarak dan populasi
      iter = iter+1
      conv <- c(conv,sum(new_uij^m*uij$d))
    }
  }
  else { ## fgwc.v
    ##stopifnot(any(is.na(vi))==F)
    if (is.na(vi)) { ##jika centroid tidak diinisialisasi, akan dilakukan generate dengan menggunakan distribusi uniform
      vi <- gen_vi(data,ncluster,"uniform",randomN)
    }
    v_new <- vi
    uij <- uij(data,v_new,m,'euclidean',2) ##menghitung matriks keanggotaan dengan menggunakan jarak euclidean untuk awal
    new_uij <- uij$u ##matriks keanggotaan yang disimpan dimasukkan ke yang baru
    v_old <- vi+1
    while (abs(sum(v_new-v_old))>error && iter<max.iter) {
      v_old <- v_new
      uij <- uij(data,v_old,m,distance,order)##memperbarui matriks keanggotaan dengan menggunakan jarak minkowski sesuai orde
      new_uij <- renew_uij(data,uij$u,mi.mj,distmat,alpha,beta,a,b) ##memodifikasi matriks keanggotaan dengan memanfaatkan matriks jarak dan populasi
      v_new <- vi(data,new_uij,m) ##menghitung centroid yang baru
      vi <- v_new ##centroid yang lama diperbarui
      conv <- c(conv,jfgwcv(data,vi,m,distance,order)) ##menghitung fungsi objektif dengan menggunakan pendekatan v
      iter = iter+1
    }
  }
  fgwc_obj <- sum(new_uij^m*uij$d) ##fungsi objektif akhir
  finaldata <- determine_cluster(data,new_uij)
  cluster <- finaldata[,ncol(finaldata)]
  result <- list("converg"=conv[-1],"f_obj"=fgwc_obj, "membership"=new_uij,"centroid"=vi, "validasi"=index_fgwc(data,cluster,new_uij,vi,m,exp(1)) ,"max.iter" = iter,
                   "cluster"=cluster, "finaldata"=finaldata, "call"=match.call(), "time" = proc.time()-ptm)
  print(c(order, ncluster,m, randomN))
  return (result)
}

##membentuk distance matrix dari lattitude longitude
formdistmat <- function (datalonglat,p=3) {
  jarak <- matrix(0,nrow(datalonglat),nrow(datalonglat))
  for (i in 1:nrow(datalonglat)) {
    for (j in 1:nrow(datalonglat)) {
      jarak[i,j] <- sum((matrix(datalonglat[i,],ncol=1)-matrix(datalonglat[j,],ncol=1))^p)^(1/p)
    }
  }
  return(jarak)
}

##vi dari nilai keanggotaan
vi <- function(data,uij,m) {
  return (t(uij^m)%*%data/colSums(uij^m))
}

##uij dari centroid
uij <- function(data,vi,m,distance,order=2) {
  u <- matrix(0,nrow(data),nrow(vi))
  d <- cdist(data,vi,distance,order)^2
  x <- (d)^(1/(m - 1))
  u <-  (1/x)/rowSums(1/x)
  res <- list("d"=d, "u"=u)
  return(res)
}

##menentukan cluster dari data
determine_cluster <- function(data,uij) {
	clust = apply(uij,1,which.max)
	return(cbind.data.frame(data,cluster=clust))
}

##memodifikasi matriks keanggotaan dengan memanfaatkan matriks jarak dan populasi
renew_uij <- function(data,old_uij,mi.mj,dist,alpha,beta,a,b) {
  diag(dist) <- Inf
  wij <- mi.mj^b/dist^a
  ##new uij
  wijmuj <- wij %*% old_uij
  A <- rowSums(wijmuj)
  new_uij <-  alpha*old_uij + (beta/A)*wijmuj
  return(new_uij)
}

##generate matrik keanggotaan
gen_uij <- function(data,ncluster,n,randomN) {
  set.seed(randomN)
  uij <- matrix(runif(ncluster*n,0,1),n,ncluster)
  return(uij/rowSums(uij))
}

##generate pusat cluster
gen_vi <- function(data,ncluster,gendist,randomN) {##generate centroid
  p <- ncol(data)
  piclass <- matrix(0,ncluster,p)
  for (i in 1:p) {
    set.seed(randomN)
    if (gendist=="normal"){
      piclass[,i] <- rnorm(ncluster,mean(data[,i]),sd(data[,i]))
    }
    else if (gendist=="uniform"){
      piclass[,i] <- runif(ncluster,min(data[,i]),max(data[,i]))
    }
  }
  return(piclass)
}

##fungsi objektif
jfgwcu <- function(data,uij,m,distance,order) { ##fungsi objektif fgwc-u
  vi <- (t(uij^m)%*%data/colSums(uij^m))
  d <- cdist(data,vi,distance,order)^2
  return(sum((uij^m)*d))
}

##fungsi objektif
jfgwcu2 <- function(data,uij,m,distance,order) { ##fungsi objektif fgwc-u
  u <- renew_uij(data,u$u,mi.mj,dist,alpha,beta,a,b)
  vi <- (t(uij^m)%*%data/colSums(uij^m))
  d <- cdist(data,vi,distance,order)^2
  return(sum((uij^m)*d))
}

jfgwcv  <- function(data,vi,m,distance,order) { ##fungsi objektif fgwc-v
  u <- matrix(0,nrow(data),nrow(vi))
  d <- cdist(data,vi,distance,order)^2
  x <- (d)^(1/(m - 1))
  for (i in 1:nrow(d)) {
    for (j in 1:ncol(d)) {
      temp <- (d[i,j]/d[i,])^(1/(m-1))
      u[i,j] <- 1/sum(temp)
    }
  }
  return(sum((u^m)*d))
}

jfgwcv2  <- function(data,vi,m,distance,order,mi.mj,dist,alpha,beta,a,b) { ##fungsi objektif fgwc-v
  u <- uij(data,vi,m,distance,order)
  u <- renew_uij(data,u$u,mi.mj,dist,alpha,beta,a,b)
  vi <- vi(data,u,m)
  d <- cdist(data,vi,distance,order)^2
  for (i in 1:nrow(d)) {
    for (j in 1:ncol(d)) {
      temp <- (d[i,j]/d[i,])^(1/(m-1))
      u[i,j] <- 1/sum(temp)
    }
  }
  return(sum((u^m)*d))
}

##penghitungan index
index_fgwc <- function(data,cluster,uij,vi,m,a=exp(1)) {
  result<-list()
  result$PC <- PC1(uij)
  result$CE <- CE1(uij,a)
  result$SC <- SC1(data,cluster,uij,vi,m)
  result$SI <- SI1(data,uij,vi)
  result$XB <- XB1(data,uij,vi,m)
  result$IFV <- IFV1(data,uij,vi,m)
  result$Kwon <- Kwon1(data,uij,vi,m)
  return(result)
}

########################################################
#################VALIDATION MEASUREMENT#################
########################################################

##kelompok yang optimum dinyatakan dengan nilai PC yang maksimum.
PC1 <- function(uij) {
  return(sum(uij^2)/nrow(uij))
}

##kelompok yang optimum dinyatakan dengan nilai indeks CE yang minimum.
CE1 <- function(uij,a=exp(1)) {##
  return(sum(uij*log(uij,a))/(-nrow(uij)))
}

##Partisi yang optimum dinyatakan dengan nilai indeks SC yang minimum.
SC1 <- function(data,cluster,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  pt1 <- colSums((uij^m)*d)
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[i,]-vi[k,])^2)
    }
    Ni <- length(which(cluster==i))
    vkvi[i,]<-Ni*vkvi[i,]
  }
  pt2 <- colSums(vkvi)
  return(sum(pt1/pt2))
}

##Jumlah kelompok yang optimum dinyatakan dengan nilai indeks S yang minimum.
SI1 <- function(data,uij,vi) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[k,]-vi[i,])^2)
    }
  }
  diag(vkvi) <- Inf
  pt1 <- sum((uij^2)*d)
  pt2 <- nrow(data)*min(vkvi)
  return(sum(pt1/pt2))
}

##Jumlah kelompok yang optimal dinyatakan dengan nilai XB yang minimum.
XB1 <- function(data,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
      if (d[i,j]==0) {
        d[i,j]==Inf
      }
    }
  }
  pt1 <- sum((uij^m)*d)
  pt2 <- nrow(data)*min(d)
  return(pt1/pt2)
}

##Ketika nilai IFV maksimum maka kualitas cluster semakin baik.
IFV1 <- function(data,uij,vi,m) {
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[k,]-vi[i,])^2)
    }
  }
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  sigmaD <- sum(d)/(nrow(data)*nrow(vi))
  SDmax <- max(vkvi)
  log2u <- colSums(log(uij,2))/nrow(data)
  u2ij <- colSums(uij^2)
  inside <- sum(u2ij*(log(nrow(vi),2)-log2u))
  return(sum(u2ij*((log(nrow(vi),2)-log2u)^2)/nrow(data)*(SDmax/sigmaD)))
}

##Ketika nilai Kwon minimum maka kualitas cluster semakin baik.
Kwon1 <- function(data,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  s <- matrix(0,nrow(vi))
  vivj <- matrix(0,nrow(vi),nrow(vi))
  for (j in 1:nrow(vi)) {
    s[j,] <- (sum(vi[j,]-colMeans(data))^2)
  }
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  for (i in 1:nrow(vi)) {
    for (j in 1:nrow(vi)) {
      vivj[i,j] <- sum((vi[i,]-vi[j,])^2)
    }
  }
  diag(vivj) <- Inf
  pt1 <- colSums((uij^m)*d)
  pt2 <- sum(s)/nrow(vi)
  pt3 <- min(vivj)
  return(sum((pt1+pt2)/pt3))
}
