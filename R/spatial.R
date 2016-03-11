#' Calculate wave fetch
#'
#' @param x coord of point at wich to calc fetch.
#' @param y coord of point at wich to calc fetch.
#' @param rad Radius of fetch distance.
#' @param dir Number of compass directions to calc fetch.
#' @param coast SpatialPolgons of coastline to use.
#' @return The fetch distance at each compass point.
#' @examples
#' fetch(0, 0,1000,32,coast)
#' @export
fetch <- function(x,y,rad,dir=32,coast){
  
  int <- 360/dir
  t <- seq(0,359,int)
  xr = (rad * cos(deg.t.rad(t)))+x
  yr = (rad * sin(deg.t.rad(t)))+y
  
  lns <- lapply(1:length(t),function(t) Lines(Line(cbind(c(x,xr[t]),c(y,yr[t]))),ID=t))
  Sl = SpatialLines(lns)
  
  ld <- gDifference(Sl,coast,byid = T)
  
  lns.1 <- lapply(1:length(t),function(i) Lines(ld@lines[[i]]@Lines[[1]],ID=i))
  Sl = SpatialLines(lns.1)
  l <- gLength(Sl,byid = T)
  names(l) <- t
  l
}


###Leafet##########################################

geoserverWMS <- function(.,layer,group,z){
  addWMSTiles(.,"http://148.252.96.22:8081/geoserver/geodata/wms?",
              layers = layer,group = group,
              options = WMSTileOptions(format = "image/png", transparent = TRUE,zIndex=z))
}

#' WGS84
#'
#' @param x Spatial* object
#'
#' @return Spatial object transformed to WGS 1984
#' @export
#'
#' @examples 
#' WGS84(coast)
WGS84 <- function(x){spTransform(x,CRS("+init=epsg:4326"))}
#' UTM29
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
UTM29 <- function(x){spTransform(x,CRS("+init=epsg:32629"))}
#' UTM30
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
UTM30 <- function(x){spTransform(x,CRS("+init=epsg:32630"))}
#' UTM31
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
UTM31 <- function(x){spTransform(x,CRS("+init=epsg:32631"))}
#' BNG
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
BNG <- function(x){spTransform(x,CRS("+init=epsg:27700"))}


#' max_extent
#'
#' @param rlist 
#'
#' @return
#' @export
#'
#' @examples
max_extent <- function(rlist){
  # given list of rasters
  # returns union of extent
  xmin=min(sapply(rl,FUN=function(x){extent(x)@xmin}))
  xmax=max(sapply(rl,FUN=function(x){extent(x)@xmax}))
  ymin=min(sapply(rl,FUN=function(x){extent(x)@ymin}))
  ymax=max(sapply(rl,FUN=function(x){extent(x)@ymax}))
  
  extent(c(xmin,xmax,ymin,ymax))
}




#' long2UTM
#'
#' @param long 
#'
#' @return
#' @export
#'
#' @examples
long2UTM <- function(long) {
  (floor((long + 180)/6) %% 60) + 1
}

#' aspect.rad
#'
#' @param by 
#'
#' @return
#' @export
#'
#' @examples
aspect.rad = function(by){
  # creates aspect in radians with edges fixed
  aspect.rad <- raster::terrain(by,opt='aspect',unit='radians')    
  f <- matrix(1, nrow=5, ncol=5)
  v <- raster::focal(aspect.rad, w=f,
                     fun=function(x, ...) mean(x,na.rm=T), pad=T, padValue=NA,na.rm=T)#,
  v[is.na(v)] = median(aspect.rad,na.rm=T)
  vm = raster::merge(aspect.rad,v)  
  raster::mask(vm,by)
  
}

#' slope.rad
#'
#' @param by 
#'
#' @return
#' @export
#'
#' @examples
slope.rad = function(by){
  # creates aspect in radians with edges fixed
  slope.rad <- raster::terrain(by,opt='slope',unit='radians')    
  f <- matrix(1, nrow=5, ncol=5)
  v <- raster::focal(slope.rad, w=f, fun=function(x, ...) mean(x,na.rm=T), pad=T, padValue=NA,na.rm=T)#,
  v[is.na(v)] = median(slope.rad,na.rm=T)
  vm = raster::merge(slope.rad,v)  
  raster::mask(vm,by) 
}

#' Title
#'
#' @param byc 
#' @param w 
#'
#' @return
#' @export
#'
#' @examples
roughness = function(byc,w){
  f <- matrix(1, nrow=w, ncol=w)
  rough <- focal(byc, w=f, fun=function(x, ...) max(x,na.rm=T) - min(x,na.rm=T), pad=T, padValue=NA)#, na.rm=TUE)
  mask(rough,byc)}



############################################
#' VRM
#'
#' @param by 
#' @param cs 
#'
#' @return
#' @export
#'
#' @examples
VRM <- function(by,cs){
  #require(raster)
  #by = byc
  # Description: This tool measures terrain ruggedness by calculating the vector ruggedness measure
  #              described in Sappington, J.M., K.M. Longshore, and D.B. Thomson. 2007. Quantifiying
  #              Landscape Ruggedness for Animal Habitat Anaysis: A case Study Using Bighorn Sheep in
  #              the Mojave Desert. Journal of Wildlife Management. 71(5): 1419 -1426.
  # Requirements: Spatial Analyst 
  # Author: Mark Sappington
  # Date: 2/1/2008
  
  #by <- raster("bathyN2_UTM29.img")
  #projection(by) <- projection(gp)
  #cs <- 5 
  
  # Create Slope and Aspect rasters
  asp.rad <- raster::slopeAspect(by,out='aspect')
  slp.rad <- raster::slopeAspect(by,out='slope')
  
  # Calculate x, y, and z rasters
  
  xy <- sin(slp.rad)
  z <- cos(slp.rad)
  
  x <- sin(asp.rad) * xy
  y <- cos(asp.rad) * xy
  
  # Calculate sums of x, y, and z rasters for selected neighborhood size
  
  
  #writeRaster(xsum,'xsum.img')
  
  f <- matrix(1, nrow=cs, ncol=cs)
  
  xb <- focal(x, w=f, fun=function(xi, ...) mean(xi,na.rm=T), pad=T, padValue=NA,na.rm=T)#,
  xb[is.na(xb)] = median(xb,na.rm=T)
  xm  = merge(x,xb)
  xsum <- focal(xm,w=cs,fun=sum)^2
  xsum = mask(xsum,by)
  
  yb <- focal(y, w=f, fun=function(yi, ...) mean(yi,na.rm=T), pad=T, padValue=NA,na.rm=T)#,
  yb[is.na(yb)] = median(yb,na.rm=T)
  ym  = merge(y,yb)
  ysum <- focal(ym,w=cs,fun=sum)^2
  ysum = mask(ysum,by)
  
  zb <- focal(z, w=f, fun=function(zi, ...) mean(zi,na.rm=T), pad=T, padValue=NA,na.rm=T)#,
  zb[is.na(zb)] = median(zb,na.rm=T)
  zm  = merge(z,zb)
  zsum <- focal(zm,w=cs,fun=sum)^2
  zsum = mask(zsum,by)
  
  #xsum <- focal(x,w=cs,fun=sum)^2
  #ysum <- focal(y,w=cs,fun=sum)^2
  #zsum <- focal(z,w=cs,fun=sum)^2
  
  # Calculate the resultant vector
  
  R <- sqrt(xsum + ysum + zsum)
  
  plot(R)
  
  max <- cs^2
  
  VRM <- 1-(R/ max)
  VRM
  
}


#' RGBimage.seg
#'
#' @param r 
#' @param g 
#' @param b 
#' @param factor 
#' @param output 
#'
#' @return
#' @export
#'
#' @examples
RGBimage.seg <- function(r,g,b,factor=25,output){
  #factor=25
  #output='img12_25'
  
  # segemnt at different scales
  writeRaster(r,filename='imgR.sgrd',overwrite=T)
  writeRaster(g,filename='imgG.sgrd',overwrite=T)
  writeRaster(b,filename='imgB.sgrd',overwrite=T)
  
  # region growing
  SAGA.seed("imgR.sgrd;imgG.sgrd;imgB.sgrd",factor=factor,output)
  
  SAGA.SRG(paste(output,'_seed.sgrd',sep=''),"imgR.sgrd;imgG.sgrd;imgB.sgrd",paste(output,'SRG',sep=''))
  
  # extract object features
  Rf <- getfeatures(raster(paste(output,'SRG.sdat',sep=''))+2,r)
  Gf <- getfeatures(raster(paste(output,'SRG.sdat',sep=''))+2,g)
  Bf <- getfeatures(raster(paste(output,'SRG.sdat',sep=''))+2,b)
  
  
  dimnames(Rf)[[2]] <- paste(dimnames(Rf)[[2]],'R',sep='')
  dimnames(Gf)[[2]] <- paste(dimnames(Gf)[[2]],'G',sep='')
  dimnames(Bf)[[2]] <- paste(dimnames(Bf)[[2]],'B',sep='')
  
  cbind(Rf,Gf,Bf)}


#' getfeatures
#'
#' @param seg 
#' @param r 
#'
#' @return
#' @export
#'
#' @examples
getfeatures <- function(seg,r){
  
  segm <- as.matrix(seg)
  rm <- scale(as.matrix(r))
  
  rm[is.na(rm)] <- median(rm,na.rm=T)
  
  rm.fB <- computeFeatures.basic(segm,rm)
  dimnames(rm.fB)[[2]]
  
  rm.fM <- computeFeatures.moment(segm,rm)
  dimnames(rm.fM)[[2]]
  
  rm.fS <- computeFeatures.shape(segm)
  dimnames(rm.fS)[[2]]
  
  rm.fH <- computeFeatures.haralick(segm,rm,haralick.scales=1,haralick.nbins=32)
  head(rm.fH)
  dimnames(rm.fH)[[2]]
  
  feats <- cbind(rm.fB[,1:3],rm.fM[,3:5],rm.fS[,1:3],rm.fH)
  
  dimnames(feats)[[2]] <- c("mean","sd","mad","majax","eccen","thet",
                            "area","perim","radi","asm","con",
                            "cor","var","idm","sav","sva","sen",
                            "ent","dva","den","f12","f13")
  
  feats}


#' greyscale
#'
#' @param RGB 
#'
#' @return
#' @export
#'
#' @examples
greyscale <- function(RGB){
  gs <- 0.2989 * RGB[[1]] + 0.5870 * RGB[[2]] + 0.1140 * RGB[[3]]
  return(gs)
}


#' filter.by.dist
#'
#' @param data1 
#' @param data2 
#' @param datefield 
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
filter.by.dist <- function(data1,data2,datefield,d){
  require(fields)
  #data1 <- 
  #data2 <- 
  #d=10
  #datefield='datecmb'
  # given 2 point datasets
  # will filter remove points from data1 one that fall within d 
  #of data2 and with matching values in datefield
  # assumes d is km and coords are lat lon
  
  
  # for each point in data1 get vector of ids for points in data2 that fall within d
  
  filt.dist <- function(i){
    #i <- 4
    #get points in data2 that are < d from p
    p <- data1[i,]
    dm <- as.vector(rdist.earth(p@coords,data2@coords))<d
    if(is.null(datefield)){T%in%dm} else {
      p$datecmb%in%data2@data[dm,datefield]}
    
  }
  
  l <- sapply(1:nrow(data1),FUN=filt.dist)
  
  data1[!l,]
  
  
}



#' CHtoGrid
#'
#' @param data 
#' @param res 
#'
#' @return
#' @export
#'
#' @examples
CHtoGrid <- function(data,res){
  
  #data <- cruise.list[[1]]
  #res = 0.1
  
  ext <- extent(data)
  
  xmin <- ext@xmin
  xmax <- ext@xmax
  ymin <- ext@ymin
  ymax <- ext@ymax
  
  pg <- Spatial.Point.Grid(xlim=c(ext@xmin,ext@xmax),
                           ylim=c(ext@ymin,ext@ymax),
                           xres=res,yres=res)
  
  ch <- convexhull(data)    
  coast.poly <- get.coast(data)
  
  ext1 <- gDifference(ch,coast.poly)
  grid <- gIntersection(pg,ext1)
  return(grid)
  
}




#' get.lim
#'
#' @param x 
#' @param exp 
#'
#' @return
#' @export
#'
#' @examples
get.lim <- function(x,exp=0){
  #x <- t$Lon_M
  #exp <- 0.1
  mar <- (max(x)-min(x))*exp
  c(min(x)-mar,max(x)+mar)
}








#' get.extents
#'
#' @param shps 
#'
#' @return
#' @export
#'
#' @examples
get.extents <- function(shps){
  # given a list of spatialdataframes
  # will return their merged extents
  
  
  merge <- merge.spdf(shps)  
  merge.cln <- FilterErrCoords(merge)
  return(extent(merge.cln))
}



#' extent.expand
#'
#' @param e 
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
extent.expand <- function(e,f=0.1){
  
  xi <- (e@xmax-e@xmin)*(f/2)
  yi <- (e@ymax-e@ymin)*(f/2)
  
  xmin <- e@xmin-xi
  xmax <- e@xmax+xi
  ymin <- e@ymin-yi
  ymax <- e@ymax+yi
  
  return(extent(c(xmin,xmax,ymin,ymax)))  
}



#' lst.shpfls
#'
#' @param list 
#'
#' @return
#' @export
#'
#' @examples
lst.shpfls <- function(list=dir()){
  FbE(list,'.shp')
}

#' stack.list
#'
#' @param l 
#' @param stk 
#' @param fldr 
#'
#' @return
#' @export
#'
#' @examples
stack.list <- function(l,stk=NULL,fldr=NULL){
  #stk=NULL
  #fldr='derivatives/'
  # given a list of rasters and a dir
  # returns stack of listed rasters
  if (!is.null(fldr)){wd <- getwd()
  setwd(fldr)}  
  if(is.null(stk)){stk <- raster(l[1])
  l <- l[-1]}
  
  for (r in l) {
    stk <- stack(stk,raster(r))
  }
  setwd(wd)
  return(stk)}


#' aspect.2.EN
#'
#' @param aspect 
#'
#' @return
#' @export
#'
#' @examples
aspect.2.EN <- function(aspect){
  
  eastf <- function(x){sin((x*pi)/180)}
  northf <- function(x){cos((x*pi)/180)}
  east <- calc(aspect,eastf)
  north <- calc(aspect,northf)
  
  return(stack(east,north))
}


#' set.nodata
#'
#' @param lyr 
#' @param nodata 
#'
#' @return
#' @export
#'
#' @examples
set.nodata <- function(lyr,nodata=-9999){
  #lyr <- bs.class
  lyr@file@nodatavalue <- nodata  
  lyr <- trim(lyr)
  setMinMax(lyr)
  return(lyr)
}



#' standardise.origins
#'
#' @param rasters 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
standardise.origins <- function(rasters,method='bilinear'){
  # given a list of rasters this function will make all 
  # rasters conform to the origin of the first in the list
  ol <- sapply(rasters,origin)
  i <- ol[1,] != origin(rasters[[1]])[1]
  
  ii <- 1
  for(r in rasters[i]){
    #r <- rasters[i][1]
    r <- resample(r,rasters[[1]],method=method)
    rasters[i][ii] <- r
    ii <- ii+1
  }
  return(rasters)
}

#' duplicate.mean
#'
#' @param spdf 
#' @param X 
#'
#' @return
#' @export
#'
#' @examples
duplicate.mean <- function(spdf,X){
  #X <- 'X'
  
  x <- spdf@data[,X]
  zd <- zerodist(spdf)  #all prow pairs with identical coords
  zdu <- unique(zd[,1]) #unique parirs from c1
  zdu <- zdu[!zdu%in%zd[,2]]
  
  
  # for each unique point
  for (p in zdu){
    #p = zdu[3]
    #p
    pairi <- zd[zd[,1]==p,2]
    #pairi
    
    spdf@data[p,X] <- mean(c(x[p],x[pairi]))
    
  }
  
  
  i <- (!1:nrow(spdf) %in%(zd)) | (1:nrow(spdf)%in%zdu)
  
  return(spdf[i,])}




#' duplicate.max
#'
#' @param spdf 
#' @param X 
#'
#' @return
#' @export
#'
#' @examples
duplicate.max <- function(spdf,X){
  #X <- 'X'
  #str(spdf)
  x <- spdf@data[,X]
  zd <- zerodist(spdf)  #all prow pairs with identical coords
  zdu <- unique(zd[,1]) #unique parirs from c1
  zdu <- zdu[!zdu%in%zd[,2]]
  
  
  # for each unique point
  for (p in zdu){
    #p = zdu[3]
    #p
    pairi <- zd[zd[,1]==p,2]
    #pairi
    
    spdf@data[p,X] <- max(c(x[p],x[pairi]))
    
  }
  
  
  i <- (!1:nrow(spdf) %in%(zd)) | (1:nrow(spdf)%in%zdu)
  
  return(spdf[i,])}



#' jitter
#'
#' @param spdf 
#' @param jiter 
#'
#' @return
#' @export
#'
#' @examples
jitter <- function(spdf,jiter){
  # given a spatial point dataframe this will 
  # apply a specified amount of jitter to avoid points having 
  # idetical coords
  
  rnge <- seq(0+jiter/10,jiter,jiter/10)
  while (nrow(spdf) != nrow(remove.duplicates(spdf))){
    zd <- zerodist(spdf,unique.ID=T)
    for (i in 1:length(zd)){
      #i <- 1
      if(i != zd[i]){
        xy <- sample(c(1,2),1)
        dist <- sample(rnge,1)
        spdf@coords[zd[i],xy] <- spdf@coords[zd[i],1]+dist
      }
    }
  }
  return(spdf)}


#' lon.cor
#'
#' @param lon 
#' @param lat 
#'
#' @return
#' @export
#'
#' @examples
lon.cor <- function(lon,lat){
  return(lon * cos((lat*pi)/180))
  
  #longitude * cos((latitude*pi)/180)
  #return(lon / cos(lat*pi/180))
}

#' lon.cori
#'
#' @param lon 
#' @param lat 
#'
#' @return
#' @export
#'
#' @examples
lon.cori <- function(lon,lat){
  #return(lon * cos((lat*pi)/180))
  #longitude * cos((latitude*pi)/180)
  return(lon / cos((lat*pi)/180))
}

#' stack.kmeans
#'
#' @param stack 
#' @param cols 
#' @param cents 
#' @param plot 
#' @param scale 
#'
#' @return
#' @export
#'
#' @examples
stack.kmeans <- function(stack,cols=NA,cents=3:5,plot=F,scale=F){
  # given a raster stack object this function will 
  # apply a multivariate kmeans classification
  
  i <- complete.cases(stack[])
  data <- stack[i]
  
  if (!is.na(cols)){data <- data[,cols]}
  
  if (scale==T){data <- scale(data,center=T)}
  
  if (plot==T){kmeans.plot(data)}  
  
  clst.lst <- list()
  for (c in cents){
    clust <- kmeans(data,centers=c,nstart=99)
    clst.lst <- append(clst.lst,list(clust))
    names(clst.lst)[length(clst.lst)] <- paste('c',c,sep='')
  }
  
  stack <- stack(stack)  
  rlst <- list()
  for (c in clst.lst){
    cr <- stack@layers[[1]]
    cr[!is.na(cr)] <- c['cluster'][[1]]
    rlst <- append(rlst,list(cr))
  }
  clst.stack <- stack(rlst)
  layerNames(clst.stack) <- names(clst.lst)  
  clst.lst <- append(clst.lst,list(stack=clst.stack))    
  return(clst.lst)
}

#' stack.pca
#'
#' @param stack 
#' @param scale 
#'
#' @return
#' @export
#'
#' @examples
stack.pca <- function(stack,scale=T){
  ## performs pca on raster stack
  ## returns results of pca and new stack of pca scores
  
  i <- complete.cases(stack[])
  mydata <- stack[i]
  if (scale==T){mydata <- scale(mydata,center=T)}
  pca <- princomp(mydata,cor=T)  
  pca.stack <- stack
  values(pca.stack)[i] <- pca$scores
  dimnames(pca.stack@data@values)[[2]] <- names(pca$sdev)
  pca.stack@layernames <- names(pca$sdev)
  print(summary(pca))
  plot(pca)
  
  return(list(pca=pca,stack=pca.stack))
  
}


#' stack.cor.mat
#'
#' @param stack 
#' @param smple.size 
#'
#' @return
#' @export
#'
#' @examples
stack.cor.mat <- function(stack,smple.size=1000){
  # creates corrleation matrix from random samples of raster stack 
  
  smple <- sampleRandom(x=stack,size=smple.size)
  
  pairs(smple,lower.panel=panel.smooth, upper.panel=panel.cor)
  
}



#' fillGaps
#'
#' @param data 
#' @param grid 
#'
#' @return
#' @export
#'
#' @examples
fillGaps <- function(data,grid){
  require(gstat)
  # interpolates to fill the gaps in a raster
  
  #grid <- extent.fin
  #data <- wave.msk
  
  grid[!is.na(grid[])] <- -9999
  
  data.gaps <- merge(data,grid)
  data.gaps[data.gaps!=-9999] <- NA
  data.gaps.sp <- rasterToPoints(data.gaps,spatial=T)
  data.sp <- rasterToPoints(data,spatial=T)
  gaps.idw <- idw(formula=layer~1,data.sp,newdata=data.gaps.sp,nmax=5,idp=1)
  data.gaps.sp$pred <- gaps.idw$var1.pred
  gaps.r <- rasterize(data.gaps.sp,data.gaps,field='pred')
  
  filled <- merge(gaps.r,data)
  filled.flt <- focalFilter(filled,filter=filt.gau(),fun=sum,na.rm=FALSE,pad=TRUE,padValue=NA)
  final <- merge(data,filled.flt,filled)
  
  return(final)
}





#' extent.polygon
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
extent.polygon <- function(data){
  e<-extent(data)
  crds <- matrix(c(e@xmin,e@xmin,e@xmax,e@xmax,e@xmin,e@ymin,e@ymax,e@ymax,e@ymin,e@ymin),ncol=2)
  polys <- Polygons(list(Polygon(crds)),ID=c(1))
  sp <-SpatialPolygons(list(polys))    
  projection(sp) <- projection(data)
  return(sp)
}


#' clip.to.extent
#'
#' @param data 
#' @param extent 
#'
#' @return
#' @export
#'
#' @examples
clip.to.extent <- function(data,extent){
  #clips vector to extent of another vector
  ext <- extent.polygon(extent)
  projection(ext) <- projection(data)
  data.clip <- gIntersection(data,ext)
  return(data.clip)
}


#' get.polygon.utm
#'
#' @param polygon 
#'
#' @return
#' @export
#'
#' @examples
get.polygon.utm <- function(polygon){
  # gets utm of polygons based upon their centroid
  load('C:/geodata/R/utm.rdata')  
  poly.cent <- gCentroid(polygon,byid=T)
  poly.ov <- overlay(poly.cent,utm)
  poly.utms <- sapply(poly.ov,function(x){utm@data[utm@data$ZONE_ID==x,'ZONE']})
  return(poly.utms)
}


#' Spatial.Point.Grid
#'
#' @param xlim 
#' @param ylim 
#' @param xres 
#' @param yres 
#'
#' @return
#' @export
#'
#' @examples
Spatial.Point.Grid <- function(xlim,ylim,xres,yres){
  
  #,c(48.6,53.8),0.1
  #xlim <- c(-2,1.8)
  #ylim <- c(48.6,53.8)
  #xres <- 0.025
  
  
  if (xlim[1] > xlim[2]){xres <- -xres}
  if (ylim[1] > ylim[2]){yres <- -yres}
  
  xvals <- seq(xlim[1],xlim[2],xres)
  yvals <- seq(ylim[1],ylim[2],yres)
  
  df <- data.frame()
  i <- 1
  ii <- length(xvals)
  
  for (y in yvals){
    #y <- yvals[4]    
    df[i:ii,'i'] <- i:ii  
    df[i:ii,'y'] <- y
    df[i:ii,'x'] <- xvals
    i <- i + length(xvals)
    ii <- ii + length(xvals)
  }
  
  spdf <- SpatialPointsDataFrame(cbind(df$x,df$y),df)
  return(spdf)
  
}



#' dd
#'
#' @param deg 
#' @param min 
#' @param hem 
#'
#' @return
#' @export
#'
#' @examples
dd <- function(deg,min,hem=NULL){
  dd <- deg+(min/60)
  if(!is.null(hem)){dd[tolower(hem)=='w'] <- -dd[tolower(hem)=='w'] }
  
  dd  
}


#' knn
#'
#' @param point 
#' @param data 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
knn <- function(point,data,k){
  #returns k nearest neighbours
  require(spatstat)
  n.ppp <- as.ppp(point)
  data.ppp <- as.ppp(data)
  n.nn <- nncross(data.ppp,n.ppp)
  n.nn$id <- c(1:nrow(n.nn))
  nn.sort <- n.nn[order(n.nn$dist),] 
  nn.i <- nn.sort[c(1:k),'id']
  return(data[nn.i,])
}


#' merge.spdf
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
merge.spdf <- function(data){
  # merge spatial dataframes retaining common fields
  # must be in the same coord sys
  varlist= c()
  for (n in names(data[[1]])){
    exist=F
    for (d in data[-1]){
      if (n %in% names(d)){exist=T}
      else {exist = F
      break}
    }
    if (exist == T){varlist <- union(varlist,c(n))}
  }
  
  merge <- data[[1]][varlist]
  for (d in data[-1]){merge <- rbind(merge, d[varlist])}
  return(merge)
}



#' convexhull
#'
#' @param points 
#'
#' @return
#' @export
#'
#' @examples
convexhull <- function(points){
  require(spatstat)
  ## convex hull
  #library(rgdal)
  cv <- convexhull.xy(points@coords[,1],points@coords[,2])
  cvdf <- SpatialPolygonsDataFrame(owin2SP(cv),data.frame('id'=1))
  cvdf@proj4string <- points@proj4string
  return(cvdf)}



#' points.to.lines
#'
#' @param points 
#' @param multi 
#'
#' @return
#' @export
#'
#' @examples
points.to.lines <- function(points,multi){
  ## points to lines
  library(rgdal)
  if (multi == FALSE){
    Sl1 = Line(points@coords)
    S1 = Lines(list(Sl1))
    Sl = SpatialLines(list(S1))
    Sldf <- SpatialLinesDataFrame(Sl, data.frame('id'=1),match.ID=F)
  } else {
    
    llist <- NULL
    for (i in 1:(dim(points@coords)[1]-1)){
      lstart <- points@coords[i,]
      lend <- points@coords[i+1,]
      Sl1 <- Line(rbind(lstart,lend))
      S1 <- Lines(list(Sl1),ID=as.character(i))
      llist <- union(llist, list(S1))
    }
    Sl = SpatialLines(llist)
    Sldf <- SpatialLinesDataFrame(Sl, data.frame(c(1:length(llist))),match.ID=F)
  }
  Sldf@proj4string <- points@proj4string
  return(Sldf)}


## reproject

#' projectVector
#'
#' @param data 
#' @param cs 
#'
#' @return
#' @export
#'
#' @examples
projectVector <- function(data,cs){
  library(rgdal)
  crs <- crs.list()[[1]]
  codes <- crs.list()[[2]]
  pdata <- spTransform(data,CRS(paste('+init=epsg:',codes[which(crs==cs)],sep='')))
  return(pdata)
}

## define project

#' defineProjection
#'
#' @param data 
#' @param cs 
#'
#' @return
#' @export
#'
#' @examples
defineProjection <- function(data,cs){
  library(sp)
  crs <- crs.list()[[1]]
  codes <- crs.list()[[2]]
  proj4string(data) <- CRS(paste('+init=epsg:',codes[which(crs==cs)],sep=''))
  return(data)
}

#' EPSG
#' EPSDG codes
#' @return
#' @export
#'
#' @examples
EPSG <- function(){
  crs <- c('WGS84','UTM29','UTM30','UTM31','UTM32','UTM33','BNG')
  codes <- c(4326,32629,32630,32631,32632,32633,27700)
  return(list(crs,codes))
}


# 

#' owin2Polygons
#' convert spatstat objects to sp classes
#'
#' @param x 
#' @param id 
#'
#' @return
#' @export
#'
#' @examples
owin2Polygons <- function(x, id="1") {
  library(sp)
  stopifnot(is.owin(x))
  x <- as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

#' tess2SP
#' convert spatstat objects to sp classes
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
tess2SP <- function(x) {
  library(sp)
  stopifnot(is.tess(x))
  y <- tiles(x)
  nam <- names(y)
  z <- list()
  for(i in seq(y))
    z[[i]] <- owin2Polygons(y[[i]], nam[i])
  return(SpatialPolygons(z))
}

#' owin2SP
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
owin2SP <- function(x) {
  library(sp)
  stopifnot(is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

#' classify.psa.vector
#'
#' @param MUD 
#' @param SAND 
#' @param GRAV 
#' @param max 
#'
#' @return
#' @export
#'
#' @examples
classify.psa.vector <- function(MUD,SAND,GRAV,max=100){
  #MUD=bgs$mud[1:100]
  #SAND=bgs$sand[1:100]
  #GRAV=bgs$grav[1:100]
  #max=100
  #require(epicalc)
  
  if (max==1){
    MUD <- MUD*100
    SAND <- SAND*100
    GRAV <- GRAV*100
  }
  
  MUD[GRAV==100] <- 0.01
  SAND[GRAV==100] <- 0.01
  GRAV[GRAV==100] <- 99.98
  
  smr <- SAND/MUD
  
  folk.smr<- cut(smr,c(-Inf,psa.rc$FOLK1$brks,Inf))
  levels(folk.smr) <- c(1,10,100,1000)
  
  eunis1.smr <- cut(smr,  c(-Inf,psa.rc$EUNIS1$brks,Inf))
  levels(eunis1.smr) <- c(1,10,100)
  eunis2.smr <- cut(smr,c(-Inf,psa.rc$EUNIS2$brks,Inf))
  levels(eunis2.smr) <- c(1,10,100,1000)
  eunis3.smr <- cut(smr,c(-Inf,psa.rc$EUNIS3$brks,Inf))
  levels(eunis3.smr) <- c(1,10,100,1000,10000)
  
  #reclass Grav
  #if (max == 1){grav.brks <- c(0,0.01,0.05,0.3,0.8,1)}
  #if (max == 100){grav.brks <- c(-1,1,5,30,80,101)}  
  grav.brks <- c(-1,1,5,30,80,101)
  
  grav.rc <- cut(GRAV,grav.brks)
  levels(grav.rc) <- c(1,2,3,4,5)
  
  f <- as.numeric(as.vector(folk.smr))*as.numeric(as.vector(grav.rc))
  e1 <- as.numeric(as.vector(eunis1.smr))*as.numeric(as.vector(grav.rc))
  e2 <- as.numeric(as.vector(eunis2.smr))*as.numeric(as.vector(grav.rc))
  e3 <- as.numeric(as.vector(eunis3.smr))*as.numeric(as.vector(grav.rc))
  
  df <- data.frame(folk1=f,folk2=f,eunis1=e1,eunis2=e2,eunis3=e3)
  #df
  
  df$folk1 = factor(df$folk1,levels=psa.rc$FOLK1$fromto)
  levels(df$folk1) = psa.rc$FOLK1$class
  levels(df$folk1) = psa.rc$FOLK1$names  
  
  df$folk2 = factor(df$folk2,levels=psa.rc$FOLK2$fromto)
  levels(df$folk2) = psa.rc$FOLK2$class
  levels(df$folk2) = psa.rc$FOLK2$names  
  
  df$eunis1 = factor(df$eunis1,levels=psa.rc$EUNIS1$fromto)
  levels(df$eunis1) = psa.rc$EUNIS1$class
  levels(df$eunis1) = psa.rc$EUNIS1$names  
  
  df$eunis2 = factor(df$eunis2,levels=psa.rc$EUNIS2$fromto)
  levels(df$eunis2) = psa.rc$EUNIS2$class
  levels(df$eunis2) = psa.rc$EUNIS2$names  
  
  df$eunis3 = factor(df$eunis3,levels=psa.rc$EUNIS3$fromto)
  levels(df$eunis3) = psa.rc$EUNIS3$class
  levels(df$eunis3) = psa.rc$EUNIS3$names  
  
  
  return(df)
  
}

#' classify.psa.raster
#'
#' @param MUD 
#' @param SAND 
#' @param GRAV 
#' @param max 
#'
#' @return
#' @export
#'
#' @examples
classify.psa.raster <- function(MUD,SAND,GRAV,max=1){
  require(raster)
  # given MSG layers
  # returns raster stack of all major classifications
  #max=1
  #MUD <- MSGr[[3]]#raster(paste('CH_mud.img',sep=''))
  #SAND <- MSGr[[2]]#raster(paste('CH_sand.img',sep=''))
  #GRAV <- MSGr[[1]]#raster(paste('CH_grav.img',sep=''))
  #calc_smr
  smr <- SAND/MUD
  #smr <- overlay(SAND,MUD,fun=function(x,y){return(x/y)})
  
  reclass_smr <- function(brks,cls){
    #brks <- c(0.1111,1,9)
    #cls <- c(1,10,100,1000)
    from <- union(c(0),brks)
    to <- union(brks,c(Inf))
    return(reclassify(smr,cbind(from,to,cls)))
  }
  
  folk.smr <- reclass_smr(psa.rc$FOLK1$brks,c(1,10,100,1000))
  eunis1.smr <- reclass_smr(psa.rc$EUNIS1$brks,c(1,10,100))
  eunis2.smr <- reclass_smr(psa.rc$EUNIS2$brks,c(1,10,100,1000))
  eunis3.smr <- reclass_smr(psa.rc$EUNIS3$brks,c(1,10,100,1000,10000))
  
  #reclass Grav
  if (max == 1){grav.brks <- c(0.01,0.05,0.3,0.8)}
  if (max == 100){grav.brks <- c(1,5,30,80)}  
  from <- union(c(0),grav.brks)
  to <- union(grav.brks,c(Inf))
  grav.cls <- c(1,2,3,4,5)
  
  grav.rc <- reclassify(GRAV,rcl=cbind(from,to,grav.cls))
  grav.rc
  # generate and rc folk
  
  folk1 <- reclassify(folk.smr*grav.rc,
                      cbind(psa.rc$FOLK1$fromto-1,psa.rc$FOLK1$fromto+1,psa.rc$FOLK1$class))
  folk2 <- reclassify(folk.smr*grav.rc,
                      cbind(psa.rc$FOLK2$fromto-1,psa.rc$FOLK2$fromto+1,psa.rc$FOLK2$class))
  eunis1 <- reclassify(eunis1.smr*grav.rc,
                       cbind(psa.rc$EUNIS1$fromto-1,psa.rc$EUNIS1$fromto+1,psa.rc$EUNIS1$class))
  eunis2 <- reclassify(eunis2.smr*grav.rc,
                       cbind(psa.rc$EUNIS2$fromto-1,psa.rc$EUNIS2$fromto+1,psa.rc$EUNIS2$class))
  eunis3 <- reclassify(eunis3.smr*grav.rc,
                       cbind(psa.rc$EUNIS3$fromto-1,psa.rc$EUNIS3$fromto+1,psa.rc$EUNIS3$class))
  
  stk <- stack(folk1,folk2,eunis1,eunis2,eunis3)
  #sampleRandom(stk,size=100)
  
  names(stk) <- c('folk1','folk2','EUNIS1','EUNIS2','EUNIS3')
  return(stk)
}



#' poly_coords
#'
#' @param shapefile 
#'
#' @return
#' @export
#'
#' @examples
poly_coords<- function(shapefile){
  #http://spatialanalysis.co.uk/2010/09/27/maps-with-ggplot2/
  # accesses shapefile geometry for ggplot2
  
  if (nrow(data.frame(shapefile$ID))< 1) 
    
  {
    
    print ("No ID field in SpatialPolygon")
    
  }else{
    
    Order<-0 
    YX3<- as.numeric("XX", "XX", "XX", "XX")
    num_polys<- nrow(shapefile@data)+1
    YX3<- as.numeric("XX", "XX", "XX")
    
    curr_poly<- shapefile@data[1,]
    curr_poly_start_row <- 1
    poly_old= F
    
    for(curr_row in curr_poly_start_row:num_polys)
    {
      curr_poly_row<-shapefile@data[curr_row,]
      curr_poly_end_row = curr_row - 1	
      Poly_n= shapefile@data[curr_poly_start_row:curr_poly_end_row,]
      curr_poly_start_row = curr_row
      Poly_Name<-as.vector(Poly_n$ID)
      Poly<-shapefile[shapefile$ID==Poly_Name,]
      PolyCoords<-lapply(slot(Poly, "polygons"), function(x) lapply(slot(x,
                                                                         "Polygons"), function(y) slot(y, "coords")))
      PolyCoordsY<-PolyCoords[[1]][[1]][,1]
      PolyCoordsX<-PolyCoords[[1]][[1]][,2]
      Order<- 1:nrow(data.frame(PolyCoordsX)) + max(Order)
      if (poly_old != Poly_n$ID)
      {
        YX1<- data.frame(Poly_Name, Order, PolyCoordsY, PolyCoordsX)
        YX2<-rbind(YX3,YX1)
        YX3<-YX2
      }
      poly_old<-Poly_n$ID
    }
    
    join<-merge(YX3, shapefile@data, by.x="Poly_Name", by.y= "ID", all=T)
    join[order(join$Order),][1:nrow(join)-1,]
  }
}

#' ncdfraster
#'
#' @param f 
#' @param var 
#' @param crs 
#' @param extent 
#'
#' @return
#' @export
#'
#' @examples
ncdfraster <- function(f,var,crs,extent){
  
  chlsd <- get.var.ncdf(f,var)
  chlsd <- t(chlsd)
  chlsd <- chlsd[nrow(chlsd):1,]
  chlsd <- raster(chlsd)
  extent(chlsd) <- extent#c(min(lon),max(lon),min(lat),max(lat))
  #chlsd <- crop(chlsd,r)
  proj4string(chlsd) <- crs
  chlsd}


#' sdf
#'
#' @param spdf 
#'
#' @return
#' @export
#'
#' @examples
sdf <- function(spdf){
  cbind(spdf@coords,spdf@data)  
}




