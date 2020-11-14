
########################################################
#################PLOTTING AND MAPPING###################
########################################################

##pembentukan biplot
bi.plot <- function(data,cluster) {
  library(factoextra)
  library(FactoMineR)
  res.pca <- prcomp(data, scale. = T)
  bipl <- fviz_pca_biplot(X=res.pca, axes = c(1, 2), geom = c("point", "text"), col.ind = "black", col.var = "steelblue", label = "all",
                          invisible = "none", repel = FALSE, habillage = cluster, palette = NULL,
                          addEllipses = FALSE, title = "Biplot Cluster")+theme_minimal()
  return(bipl)
}

##pembentukan radarplot
radarplot <- function(data, cluster){
  library(fmsb)
  hsl <- unique(cluster)
  c <- length(hsl)
  n.data <- nrow(data)
  m <- ncol(data)
  s <- matrix(NA, nrow = c,ncol = m)
  kel <- 1
  cluster <- matrix(c(1:c), nrow = c, ncol = 1)
  colnames(cluster) <- "cluster"
  #mulai

  maximum_value_data<-matrix(NA, ncol = ncol(data))
  minimum_value_data<-matrix(NA, ncol = ncol(data))
  for (d in 1:ncol(data)) {
    minimum_value_data[1,d]<-data[which.min(data[,d]),d]
    maximum_value_data[1,d]<-data[which.max(data[,d]),d]
  }
  minimum_value_data <- c(NA, minimum_value_data)
  maximum_value_data <- c(NA, maximum_value_data)

  for(i in 1:c){
    aa<-matrix(data[which(cluster==hsl[i]),],ncol=m)
    s[i,]<-colMeans(aa)
  }
  print(s)
  colnames(s) <- colnames(data)
  s <- cbind(cluster, s)
  sgab<-rbind(minimum_value_data,maximum_value_data,s)
  datachart <- as.data.frame(sgab)
  COL<-colorRampPalette(c("red", "blue"))(nrow(datachart)-2)
  chart <- radarchart(datachart[,-1], pcol = COL, cglcol = "grey80", seg = 10,
                      title = "Radar Chart Cluster")
  legend(2, 1, legend = levels(as.factor(datachart$cluster)), title = "cluster",
         col = COL, seg.len = 2, border = "transparent", pch = 16, lty = 1)

  #selesai

  return(chart)
}

##pembentukan pemetaan
##pengurutan data dengan memanfaatkan file dbf dari shp
fitmap <- function (data,data_key,dbf) {
  if (is.character(dbf)) {
    require(foreign)
    dbf <- read.dbf(dbf)
  }
  else {
    dbf <- dbf
  }
  if (is.factor(data_key)) {
    data_key <- as.character(data_key)
  }
  sor <- rep(0,nrow(data))
  for (i in 1:ncol(dbf)) {
    testing <- dbf[order(dbf[,i],decreasing = F),i]==data_key
    bool=any(testing)
    ##print(bool)
    if (bool==T) {
      key <- dbf[,i]
    }
  }
  for (j in 1:nrow(data)) {
    sor[j] <- which(data_key==key[j])
  }
  return(sor)
}

##membentuk data peta sekaligus matriks jarak
checkmap <- function(dist) {
  options(warn=-1)
  if (!is.character(dist)) {
    dist <- as.matrix(dist)
  }
  else {
    a <- dist
    if (class(dist)!="SpatialPolygonsDataFrame") {
      map<-readOGR(dist)
    }
    centroid <- gCentroid(map, byid = T)
    distance <- as.matrix(spDists(centroid, longlat = T))
  }
  res <- list("map"=map,"dist"=distance)
  return(res)
}

##membuat peta dengan menggunakan data peta yang didapat sebelumnya
createmap <- function (map,cluster) {
  datamap <- map_data(map)
  clst <- rep(0,nrow(datamap))
  for (i in 1:length(cluster)) {
    b <- which(datamap$region==i-1)
    clst[b] <- as.character(cluster[i])
  }
  cluster <- clst
  peta <- ggplot(data = datamap) +
    geom_polygon(aes(x = long, y = lat, fill = cluster, group = group), color = "black", ) +
    coord_fixed(1.3)+ggtitle("Java, Bali, and Nusa Tenggara Map of Social Vulnerability")+theme(panel.background = element_blank(),
                                                                                                axis.line=element_blank(),
                                                                                                axis.text.x=element_blank(),
                                                                                                axis.text.y=element_blank(),
                                                                                                axis.ticks=element_blank(),
                                                                                                axis.title.x=element_blank(),
                                                                                                axis.title.y=element_blank())

  return(peta)
}

##pemetaan secara keseluruhan
mapping <- function(data,cluster,shapefile,datakey2) {
  require(rgeos)
  require(rgdal)
  require(foreign)
  require(ggplot2)
  require(ggmap)
  require(maps)
  require(mapdata)
  stopifnot(any(!is.na(datakey2) | !is.na(shapefile)))
  datakey <- gsub("\r","",datakey2)
  dir <- gsub(".shp",".dbf",shapefile)
  urutan <- fitmap(data,datakey,dir)
  data <- data[urutan,]
  cluster <- cluster[urutan]
  print(cluster)
  distance <- checkmap(shapefile)
  map <- distance$map
  return(createmap(map,cluster))
}
