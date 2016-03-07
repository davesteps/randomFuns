


fetch <- function(x,y,rad,dir=32,coast){
  #x,y: coords of point at wich to calc fecth
  #rad: radius for max fetch
  # dir: number of compass directions to calc
  # coast: vecotr coast line
  # x=xy$x
  # y=xy$y
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






# projectSetup <- function(name='',fld=getwd()){
# 
#   #name='sdf'
#   #fld='C:/Users/ds10/Desktop/test'
#   
#   folderStructure(fld)
#   if(name!=''){name=paste(name,'_',sep='')}
#   
# 
#   zz <- file(file.path(getwd(),paste(name,"initial_script.R",sep='')), "w")  # open an output file connection
#   cat("#########################",file = zz, sep = "\n")
#   cat("# Author: David Stephens", file = zz, sep = "\n")
#   cat("# Title:",paste(name,"initial_script.R",sep=''), "\n", file = zz)
#   cat('# Created automatically on:',as.character(Sys.Date()),"\n", file = zz)
#   cat("#########################", file = zz)
#   close(zz)  
#   
# }



folderStructure<-function(fld=getwd()){
  #fld=getwd()
  dir.create(file.path(fld, 'R'), showWarnings = FALSE)
  dir.create(file.path(fld, 'data'), showWarnings = FALSE)
  dir.create(file.path(fld, 'output'), showWarnings = FALSE)
  dir.create(file.path(fld, 'docs'), showWarnings = FALSE)
  dir.create(file.path(fld, 'figs'), showWarnings = FALSE)
  
  zz <- file("data/_README.txt", "w")  # open an output file connection
  cat("This folder is a data repository and should be read-only", file = zz, sep = "\n")
  cat('Created automatically on',as.character(Sys.Date()), file = zz)
  close(zz)  
  zz <- file("R/_README.txt", "w")  # open an output file connection
  cat("This folder is a repository for R code",file = zz, sep = "\n")
  cat('Created automatically on',as.character(Sys.Date()), file = zz)
  close(zz)  
  zz <- file("docs/_README.txt", "w")  # open an output file connection
  cat("This folder is for project documents",file = zz, sep = "\n")
  cat('Created automatically on',as.character(Sys.Date()), file = zz)
  close(zz)  
  zz <- file("figs/_README.txt", "w")  # open an output file connection
  cat("This folder is for figures generated automatically",file = zz, sep = "\n")
  cat('Created automatically on',as.character(Sys.Date()), file = zz)
  close(zz)  
  zz <- file("output/_README.txt", "w")  # open an output file connection
  cat("This folder is for automatically generated outputs",
      "eg simuation output, processed datasets, logs, or other processed things.",file = zz, sep = "\n")
  cat('Created automatically on',as.character(Sys.Date()), file = zz)
  close(zz)  
  
  print('Project directories setup')
  
}



WGS84 <- function(x){spTransform(x,CRS("+init=epsg:4326"))}
UTM29 <- function(x){spTransform(x,CRS("+init=epsg:32629"))}
UTM30 <- function(x){spTransform(x,CRS("+init=epsg:32630"))}
UTM31 <- function(x){spTransform(x,CRS("+init=epsg:32631"))}
BNG <- function(x){spTransform(x,CRS("+init=epsg:27700"))}



get.lim <- function(x,exp=0){
  #x <- t$Lon_M
  #exp <- 0.1
  mar <- (max(x)-min(x))*exp
  c(min(x)-mar,max(x)+mar)
}


max_extent <- function(rlist){
  # given list of rasters
  # returns union of extent
  xmin=min(sapply(rl,FUN=function(x){extent(x)@xmin}))
  xmax=max(sapply(rl,FUN=function(x){extent(x)@xmax}))
  ymin=min(sapply(rl,FUN=function(x){extent(x)@ymin}))
  ymax=max(sapply(rl,FUN=function(x){extent(x)@ymax}))
  
  extent(c(xmin,xmax,ymin,ymax))
}

long2UTM <- function(long) {
  (floor((long + 180)/6) %% 60) + 1
}

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

slope.rad = function(by){
  # creates aspect in radians with edges fixed
  slope.rad <- raster::terrain(by,opt='slope',unit='radians')    
  f <- matrix(1, nrow=5, ncol=5)
  v <- raster::focal(slope.rad, w=f, fun=function(x, ...) mean(x,na.rm=T), pad=T, padValue=NA,na.rm=T)#,
  v[is.na(v)] = median(slope.rad,na.rm=T)
  vm = raster::merge(slope.rad,v)  
  raster::mask(vm,by) 
}

roughness = function(byc,w){
  f <- matrix(1, nrow=w, ncol=w)
  rough <- focal(byc, w=f, fun=function(x, ...) max(x,na.rm=T) - min(x,na.rm=T), pad=T, padValue=NA)#, na.rm=TUE)
  mask(rough,byc)}

############################################
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


greyscale <- function(RGB){
  gs <- 0.2989 * RGB[[1]] + 0.5870 * RGB[[2]] + 0.1140 * RGB[[3]]
  return(gs)
}


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



rasterDerivatives <- function(inR, out, tri=F,tpi=F,rgh=F,mor=F,sob=F,gau=F,med=F){
  #inR <- bsgau
  #out <- 'Derivatives/NSGC_bsgau'
  
  if(tri==T){terrain(inR,'TRI',filename=paste(out,'_TRI.img',sep=''))}
  if(tpi==T){terrain(inR,'TPI',filename=paste(out,'_TPI.img',sep=''))}
  
  if(rgh==T){
    f <- matrix(1, nrow=3, ncol=3)
    focal(inR, w=f, fun=function(x, ...) max(x) - min(x), pad=TRUE,
          padValue=NA, na.rm=TRUE,filename=paste(out,'_RGH.img',sep=''))}
  
  if(sob==T){
    focal(inR, w=f.sobelx,filename=paste(out,'_SOBx.img',sep=''))
    focal(inR, w=f.sobely,filename=paste(out,'_SOBy.img',sep=''))
  }
  if(mor==T){
    MOR = MoranLocal(inR,w=5)
    writeRaster(MOR,paste(out,'_MOR.img',sep=''))
  }
  
  if(gau==T){focal(inR,w=f.gau(1.5,5),na.rm=T,filename=paste(out,'_gau.img',sep=''))}#,,pad=T,padValue=NA)}
  if(med==T){focal(inR,w=3,fun=median,na.rm=T,filename=paste(out,'_med.img',sep=''))}#,na.rm=T,pad=T,padValue=NA)}
  
}


# writeSHP <- function(x,d = getwd(),n){
#   
#   writeOGR(obj=x,dsn=d,layer=n,driver='ESRI Shapefile')
#   
# }


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
# 
# quick.krige <- function(data,var,res=0.01,grid=NULL,nmax=10){
#   #given a spdf and var name
#   # automaticaly kriges var to CH grid
#   
#   #data <- t1.2010
#   #var <- 'littert'
#   #res <- 0.1
#   #grid = NULL
#   
#   require(automap)
#   df <- data.frame(cbind(lon=data@coords[,1],lat=data@coords[,2],var=data@data[,var]))
#   df <- df[complete.cases(df),]
#   
#   df$loncor = lon.cor(df$lon,df$lat)
#   
#   coordinates(df) <- ~loncor+lat
#   
#   d <- autofitVariogram(var~1,input_data=df)
#   
#   print(plot(d))
#   
#   if (is.null(grid)){
#     grid <- CHtoGrid(data,res)
#   }
#   
#   
#   
#   grid.loncor <- data.frame(cbind(lon=grid@coords[,1],lat=grid@coords[,2]))
#   grid.loncor$loncor <- lon.cor(grid.loncor$lon,grid.loncor$lat)
#   coordinates(grid.loncor) <- ~loncor+lat 
#   
#   kr <- autoKrige(var~1,input_data=df,new_data=grid.loncor,block=res,nmax=nmax)
#   
#   grid.df <- SpatialPointsDataFrame(grid,data=data.frame(pred=kr$krige_output$var1.pred,var=kr$krige_output$var1.var))
#   
#   grid.r <- raster(extent(grid.df),nrow=length(unique(grid@coords[,2])),ncol=length(unique(grid@coords[,1])))
#   grid.r <- rasterize(grid.df,y=grid.r,field='pred')
#   grid.r1 <- rasterize(grid.df,y=grid.r,field='var')
#   
#   return(list(d,stack(grid.r,grid.r1))) 
# }
# 
# 

read.shps <- function(){
  #given list of shapefiles in folder
  # opens them into a list
  s <-   sapply(shps,FUN=read.shp)
  return(s)
}


get.extents <- function(shps){
  # given a list of spatialdataframes
  # will return their merged extents
  
  
  merge <- merge.spdf(shps)  
  merge.cln <- FilterErrCoords(merge)
  return(extent(merge.cln))
}

# get.coast <- function(data){
#   ext <- extent(data)
#   
#   xmin <- ext@xmin
#   xmax <- ext@xmax
#   ymin <- ext@ymin
#   ymax <- ext@ymax
#   
#   
#   coast <- map("worldHires", fill=T, plot = FALSE, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
#   IDs <- sapply(strsplit(coast$names, ":"), function(x) x[1])
#   coast.poly <- map2SpatialPolygons(coast, IDs=IDs)#, proj4string=CRS("+proj=longlat +datum=wgs84"))
#   coast.poly <- defineProjection(coast.poly,'WGS84')
#   slot(coast.poly, "polygons") <- lapply(slot(coast.poly, "polygons"), checkPolygonsHoles)
#   coast.poly <- unionSpatialPolygons(coast.poly, as.character(c(1:length(coast.poly))))
#   return(coast.poly)
# }


#geodata <- function(){
#  return(source('http://dl.dropbox.com/u/4683316/R/zspatial_data.r'))
#  }

extent.expand <- function(e,f=0.1){
  
  xi <- (e@xmax-e@xmin)*(f/2)
  yi <- (e@ymax-e@ymin)*(f/2)
  
  xmin <- e@xmin-xi
  xmax <- e@xmax+xi
  ymin <- e@ymin-yi
  ymax <- e@ymax+yi
  
  return(extent(c(xmin,xmax,ymin,ymax)))  
}


FilterErrCoords <- function(x){
  x[!(x@coords[,2] < 0.1 & x@coords[,1] < 0.1),]
}

read.shp <- function(shp){
  #require(rgdal)
  #shp <- shps[1]
  readOGR(shp,left(shp,nchar(shp)-4))
  
}


list.shapefiles <- function(list=dir()){
  FbE(list,'.shp')
}

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


aspect.2.EN <- function(aspect){
  eastf <- function(x){sin((x*pi)/180)}
  northf <- function(x){cos((x*pi)/180)}
  east <- calc(aspect,eastf)
  north <- calc(aspect,northf)
  
  return(stack(east,north))
}


set.nodata <- function(lyr,nodata=-9999){
  #lyr <- bs.class
  lyr@file@nodatavalue <- nodata  
  lyr <- trim(lyr)
  setMinMax(lyr)
  return(lyr)
}



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


lon.cor <- function(lon,lat){
  return(lon * cos((lat*pi)/180))
  
  #longitude * cos((latitude*pi)/180)
  #return(lon / cos(lat*pi/180))
}

lon.cori <- function(lon,lat){
  #return(lon * cos((lat*pi)/180))
  #longitude * cos((latitude*pi)/180)
  return(lon / cos((lat*pi)/180))
}

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


stack.cor.mat <- function(stack,smple.size=1000){
  # creates corrleation matrix from random samples of raster stack 
  
  smple <- sampleRandom(x=stack,size=smple.size)
  
  pairs(smple,lower.panel=panel.smooth, upper.panel=panel.cor)
  
}



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


tab.to.mat <- function(data,x,y,z){
  
  lats <- sort(x=unique(data[,y]),decreasing=T)
  lons <- sort(unique(data[,x]))  
  nrows <- length(lats)
  ncols <- length(lons)
  lons
  
  mat <- matrix(nrow=nrows,ncol=ncols,dimnames= list(lats,lons))
  
  for (i in 1:nrow(data)){
    #print(paste(i ,'of', nrow(data),sep=' '))
    lat <- data[i,y]
    lon <- data[i,x]  
    v <- data[i,z]
    
    mat[as.character(lat),as.character(lon)] <- v
  }
  
  return(mat)
  
}


mat.to.raster <- function(mat){
  #converts mat to raster
  ymn <- min(as.numeric(dimnames(mat)[[1]]))
  ymx <- max(as.numeric(dimnames(mat)[[1]]))
  xmn <- min(as.numeric(dimnames(mat)[[2]]))
  xmx<- max(as.numeric(dimnames(mat)[[2]]))
  
  
  return(raster(mat,xmn,xmx,ymn,ymx,proj4string(NA)))
}


extent.polygon <- function(data){
  e<-extent(data)
  crds <- matrix(c(e@xmin,e@xmin,e@xmax,e@xmax,e@xmin,e@ymin,e@ymax,e@ymax,e@ymin,e@ymin),ncol=2)
  polys <- Polygons(list(Polygon(crds)),ID=c(1))
  sp <-SpatialPolygons(list(polys))    
  projection(sp) <- projection(data)
  return(sp)
}


clip.to.extent <- function(data,extent){
  #clips vector to extent of another vector
  ext <- extent.polygon(extent)
  projection(ext) <- projection(data)
  data.clip <- gIntersection(data,ext)
  return(data.clip)
}


get.polygon.utm <- function(polygon){
  # gets utm of polygons based upon their centroid
  load('C:/geodata/R/utm.rdata')  
  poly.cent <- gCentroid(polygon,byid=T)
  poly.ov <- overlay(poly.cent,utm)
  poly.utms <- sapply(poly.ov,function(x){utm@data[utm@data$ZONE_ID==x,'ZONE']})
  return(poly.utms)
}


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



dd <- function(deg,min,hem=NULL){
  dd <- deg+(min/60)
  if(!is.null(hem)){dd[tolower(hem)=='w'] <- -dd[tolower(hem)=='w'] }
  
  dd  
}


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



convexhull <- function(points){
  require(spatstat)
  ## convex hull
  #library(rgdal)
  cv <- convexhull.xy(points@coords[,1],points@coords[,2])
  cvdf <- SpatialPolygonsDataFrame(owin2SP(cv),data.frame('id'=1))
  cvdf@proj4string <- points@proj4string
  return(cvdf)}



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

projectVector <- function(data,cs){
  library(rgdal)
  crs <- crs.list()[[1]]
  codes <- crs.list()[[2]]
  pdata <- spTransform(data,CRS(paste('+init=epsg:',codes[which(crs==cs)],sep='')))
  return(pdata)
}

## define project

defineProjection <- function(data,cs){
  library(sp)
  crs <- crs.list()[[1]]
  codes <- crs.list()[[2]]
  proj4string(data) <- CRS(paste('+init=epsg:',codes[which(crs==cs)],sep=''))
  return(data)
}

crs.list <- function(){
  crs <- c('WGS84','UTM29','UTM30','UTM31','UTM32','UTM33','BNG')
  codes <- c(4326,32629,32630,32631,32632,32633,27700)
  return(list(crs,codes))
}


# convert spatstat objects to sp classes

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

owin2SP <- function(x) {
  library(sp)
  stopifnot(is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

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

psa.rc <- list(
  
  FOLK1 = list( fromto= c(1,10,100,1000,2,20,200,2000,3,30,300,3000,4,40,400,4000,5,50,500,5000),
                class=  c(1,2,3,4,5,6,7,8,9,9,10,11,12,12,13,14,15,15,15,15),
                names=  c('M','sM','mS','S','(g)M','(g)sM','(g)mS','(g)S','gM','gmS','gS','mG','msG','sG','G'),
                colors= c("#90C8A1","#E6E91C","#E3DC46","#FDFCA0","#92CBDF","#ACC9B5","#E5E0B6","#FED022","#ACCCE5","#DDB198","#FBD493","#CDC6E5","#DCB2C0","#F7D1DE","#F499BA"),
                brks =  c(0.1111,1,9)),
  
  FOLK2 = list( fromto= c(1,10,100,1000,2,20,200,2000,3,30,300,3000,4,40,400,4000,5,50,500,5000),
                class=  c(1,2,3,4,1,2,3,4,5,5,6,7,8,8,9,10,11,11,11,11),
                names=  c('M','sM','mS','S','gM','gmS','gS','mG','msG','sG','G'),
                colors= c("#90C8A1","#E6E91C","#E3DC46","#FDFCA0","#ACCCE5","#DDB198","#FBD493","#CDC6E5","#DCB2C0","#F7D1DE","#F499BA"),
              
                brks =  c(0.1111,1,9)),  
  
  EUNIS1 = list( fromto= c(1,10,100,2,20,200,3,30,300,4,40,400,5,50,500),
                 class=  c(1,2,2,1,2,2,3,3,4,3,3,4,4,4,4),
                 names=  c('Mud and sandy mud','Sand and muddy sand','Mixed sediment','Coarse sediment'),
                 brks =  c(4,9)),
  
  EUNIS2 = list( fromto= c(1,10,100,1000,2,20,200,2000,3,30,300,3000,4,40,400,4000,5,50,500,5000),
                 class=  c(1,2 ,3  ,4   ,1,2 ,3  ,4   ,5,5 ,3  ,4   ,5,5 ,5  ,5   ,6,6 ,6  ,6),
                 names=  c('Fine mud','Sandy mud','Muddy sand','Fine sand','Mixed sediment','Coarse sediment'),
                 brks =  c(1,6,9)),
  
  EUNIS3 = list( fromto= c(1,10,100,1000,10000,2,20,200,2000,20000,3,30,300,3000,30000,4,40,400,
                           4000,40000,5,50,500,5000,50000),
                 class=  c(1,2,3,3,4,   1,2,3,3,4,   5,5,5,6,6,   5,5,5,6,6,  6,6,6,6,6),
                 names=  c('Fine mud','Sandy mud','Muddy sand','Fine sand','Mixed sediment','Coarse sediment'),
                 brks =  c(0.1111,4,9,19))
)


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


### ##////// Obselete functions that call gp.geoprocessor



grid.to.kml <- function(data,layername,outputdir){
  require(graphics)
  library(RSAGA)
  library(maptools)
  library(rgdal)
  
  tf <- tempfile()
  tf.asc <- paste(tf,'.asc',sep='')
  tf.sgrd <- paste(tf,'.sgrd',sep='')
  tf.sgrd.ll <- paste(tf,'ll','.sgrd',sep='')
  tf.asc.ll <- paste(tf,'ll','.asc',sep='')
  f.png <- paste(outputdir, layername,'.png',sep='')
  f.kml <- paste(outputdir, layername,'.kml',sep='')
  
  writeAsciiGrid(data,tf.asc)
  rsaga.esri.to.sgrd(tf.asc,tf.sgrd)
  
  rsaga.geoprocessor(lib="pj_proj4", 2, param=list(SOURCE_PROJ=paste('"',
                                                                     proj4string(data), '"',sep=''), TARGET_PROJ="\"+proj=longlat
                                                   +datum=WGS84\"",SOURCE=tf.sgrd, TARGET=tf.sgrd.ll,
                                                   TARGET_TYPE=0,INTERPOLATION=1))
  
  rsaga.sgrd.to.esri(tf.sgrd.ll,tf.asc.ll)
  
  grd.ll<-readGDAL(tf.asc.ll)
  proj4string(grd.ll) <- CRS("+proj=longlat + datum=WGS84")
  grd.kml <- GE_SpatialGrid(grd.ll)
  
  png(file=f.png,width=grd.kml$width,height=grd.kml$height,bg="transparent")
  par(mar=c(0,0,0,0),xaxs='i',yaxs='i')
  image(as.image.SpatialGridDataFrame(grd.ll[1]),col=rainbow(100,alpha=.6),
        xlim=grd.kml$xlim,ylim=grd.kml$ylim)
  
  ?as.image.SpatialGridDataFrame(grd.ll[1])
  kmlOverlay(grd.kml,kmlfile=f.kml, imagefile=f.png,name = layername)
  dev.off()
  unlink(paste(tf,'*',sep=''))
}


ascii.to.raster <- function(infile,outfile,type){
  library(RPyGeo)
  env <- rpygeo.build.env(python.path='C:/Python25/',overwriteoutput=1)
  rpygeo.geoprocessor('ASCIIToRaster_conversion',args=c(infile,outfile,type),working.directory='C:/workspace/',env=env)
}



merge.raster <- function(inlist,outfile){
  library(RPyGeo)
  list <- paste(inlist,sep='',collapse=',')
  inexp <- paste('MERGE(',list,')',sep='')
  env <- rpygeo.build.env(python.path='C:/Python25/',overwriteoutput=1,extent = "MAXOF")
  rpygeo.geoprocessor('SingleOutputMapAlgebra_sa',args=c(inexp,outfile),working.directory='C:/workspace/',env=env)
}


Repair.geometry <- function(infile){
  library(RPyGeo)
  env <- rpygeo.build.env(python.path='C:/Python25/',overwriteoutput=1)
  rpygeo.geoprocessor('RepairGeometry',args=c(infile),working.directory='C:/workspace/',env=env)
}


######MISC###################

didYouMean=function(input){
  require(RCurl)
  
  input=gsub(" ", "+", input)
  doc=getURL(paste("https://www.google.com/search?q=",input,"/", sep=""))
  
  
  dym=gregexpr(pattern ='Did you mean',doc)
  srf=gregexpr(pattern ='Showing results for',doc)
  
  
  if(length(dym[[1]])>1){
    doc2=substring(doc,dym[[1]][1],dym[[1]][1]+1000)
    s1=gregexpr("?q=",doc2)
    s2=gregexpr("/&amp;",doc2)
    new.text=substring(doc2,s1[[1]][1]+2,s2[[1]][1]-1)
    return(gsub("[+]"," ",new.text))
    break
  }
  
  else if(srf[[1]][1]!=-1){
    doc2=substring(doc,srf[[1]][1],srf[[1]][1]+1000)
    s1=gregexpr("?q=",doc2)
    s2=gregexpr("/&amp;",doc2)
    new.text=substring(doc2,s1[[1]][1]+2,s2[[1]][1]-1)
    return(gsub("[+]"," ",new.text))
    break
  }
  else(return(gsub("[+]"," ",input)))
}  


permute <- function(v1,v2,n=1e3){
  obs.diff <- abs(mean(v1)-mean(v2))
  
  vc <- c(v1,v2)
  l <- length(vc)
  c=0
  for(p in 1:n){
    si1 <- sample(1:l,size=length(v1),replace=F)
    v1p <- vc[si1]
    v2p <- vc[!((1:l)%in%si1)]
    if((mean(v1p)-mean(v2p))>obs.diff){c=c+1}}
  c/n}

permuteJB = function(v1, v2, nreps=10000) {
  #Jon Barry's permute function
  #*************************************************************
  obs.diff = abs(mean(v1) - mean(v2))
  v12 = c(v1,v2)
  l12 = length(v12)
  l1 = length(v1) 
  sim.diff = rep(0,nreps)
  for (j in 1:nreps) {
    perm = sample(v12)
    sim.diff[j] = mean(perm[1:l1]) - mean(perm[(l1+1):l12])
  }
  bigger = sim.diff[abs(sim.diff) >= obs.diff]
  pvalue = length(bigger)/(nreps)
  pvalue
}

## 

ncdfraster <- function(f,var,crs,extent){
  
  chlsd <- get.var.ncdf(f,var)
  chlsd <- t(chlsd)
  chlsd <- chlsd[nrow(chlsd):1,]
  chlsd <- raster(chlsd)
  extent(chlsd) <- extent#c(min(lon),max(lon),min(lat),max(lat))
  #chlsd <- crop(chlsd,r)
  proj4string(chlsd) <- crs
  chlsd}

I0 <- function(IR){
  ### to prepare irradiance:
  # daily value expressed as Wm-2
  # x4.15 to convert to UE m-2s-1 
  # x60 to convert to UE m-2m-1
  # x60 to convert to UE m-2h-1
  # x24 to convert to UE m-2d-1
  # /1e6 to convert to E m-2d-1
  # x0.45 to convert to par E m-2d-1
  # x0.94 to correct for surface reflection
  #4.15 * 60 * 60 * 24 * 0.94  * 0.45 / 1e6 = 0.1516709
  IR*0.1516709
}

lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

Zeu <- function(kd){
  ### to prepare Kd:
  # calculation of depth of photic zone as depth of 1% of I0
  # Zeu = 4.61/Kd
  4.61/kd}


PP <- function(IR,kd,chl){
  ## given vectors of kf,chl and IR
  # returns PP 
  #IR <- abs(PPdf$IR[i])
  #kd <- exp(PPdf$Kd[i])
  #chl <- exp(PPdf$chl[i])-0.01
  ##### calc PP #################
  ## calculate PP;
  # according to Cloern 1987
  # units: I0 E m-2 d-1; chl mg m-3; Zeu m 
  # sum by year to get annual PP
  1.46 + 0.73*chl*Zeu(kd)*I0(IR)
}

EBIMage_install = function(){
  
  #http://joelgranados.wordpress.com/2012/04/05/very-painful-ebimage-install-on-windows/
  source("http://bioconductor.org/biocLite.R")
  biocLite("EBImage")
}

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


margin <- function(mat){
  # given a matrix of probabilities
  # returns the confidence margon
  c <- ncol(mat)
  apply(mat,1,function(i) sort((i))[c]-sort((i))[c-1])}

swap.quadrant <- function(ft){
  fts <- ft
  Y <- nrow(ft) 
  y <- Y/2
  X <- ncol(ft)
  x <- X/2
  fts[1:y,1:x] <- ft[(y+1):Y,(x+1):X] 
  fts[1:y,(x+1):X] <- ft[(y+1):Y,1:x]
  fts[(y+1):Y,1:x] <- ft[1:y,(x+1):X]
  fts[(y+1):Y,(x+1):X] <- ft[1:y,1:x]
  return(fts)}

plotFF <- function(r,returnFF=F){
  #r <- bs
  if(class(r)=='RasterLayer'){r <- as.matrix(r)}
  rm <- r
  
  rm[is.na(r)] <- median(rm,na.rm=T)
  
  ft <- fft(rm)
  
  mag <- Mod(ft)
  phase <- Arg(ft)
  
  fts <- swap.quadrant(mag)
  
  ftf <- log(fts)
  ftr <- raster(ftf)
  
  #plot(r)
  plot(ftr)
  
  if(returnFF==T){return(ftr)}
  
}

imageFF <- function(r,d,l,t,lwl_buf=F,bd=0.1){

  #l <- 0.3#length of slice
  #d <- 45#major direction
  #t <- 6#taper window
  #lwl_buf=T
  #bd=0.01
  #r <- raster(nrows=241, ncols=241)
  #r[] <- 1:ncell(r)
  
  if(class(r)== "RasterLayer"){rm <- as.matrix(r) 
  } else{rm <- r}
  
  #if(nrow(r)%%2==1){
  #  rm <- rbind(rm,rep(NA,ncol(rm)))}
  #if(ncol(r)%%2==1){
  #  rm <- cbind(rm,rep(NA,nrow(rm)))}
  nai <- is.na(rm)
  rm.m <- mean(rm,na.rm=T)
  rm.sd <- sd(as.vector(rm),na.rm=T)
  
  rm[nai] <- rnorm(length(rm[nai]),rm.m,sd=rm.sd)
  
  #rm[nai] <- median(rm,na.rm=T)
  #rm[nai] <- mean(rm,na.rm=T)
  
  
  ft <- fft(rm)
  
  mag <- Mod(ft)
  phase <- Arg(ft)
  
  fts <- swap.quadrant(mag)
  
  ftf <- log(fts)
  ftr <- raster(ftf)
  
  extent(ftr)<- extent(c(-0.5,0.5,-0.5,0.5))
  #pie_slice <- function(l,d,t)
  
  d1 <- deg.t.rad(d-(t/2))
  d2 <- deg.t.rad(d+(t/2))
  
  x1 = l*cos(d1)
  y1 = l*sin(d1)
  x2 = l*cos(d2)
  y2 = l*sin(d2)
  
  plot(ftr)
  #points(x1,y1)
  #points(-x1,-y1)
  #points(x2,y2)
  #points(-x2,-y2)
  
  crds <- cbind(c(0,x1,x2,0,-x1,-x2,0),c(0,y1,y2,0,-y1,-y2,0))
  
  p <- Polygons(list(Polygon(crds)),0)
  sp <- SpatialPolygons(list(p))
  plot(sp,add=T)
  
  if(lwl_buf==T){
    require(rgeos)
    p <- SpatialPoints(cbind(0,0))
    pb <- gBuffer(p,width=bd)
    plot(pb,add=T)
    sp <- gDifference(spgeom1=sp,spgeom2=pb)
    plot(sp,add=T)  
  }
  
  PS <- rasterize(sp,ftr)
  APS <- merge(PS,ftr)
  plot(APS)
  APSm <- as.matrix(APS)
  APSms <- swap.quadrant(APSm)
  mp <- matrix(complex(modulus=exp(APSms),argument=phase),nrow=nrow(APSm),ncol=ncol(APSm))
  fftinv <- fft(mp,inverse=T)/length(mp)
  #APSr <- raster(Mod(fftinv))
  modfft <- Mod(fftinv)
  modfft[nai] <- NA
  #modfft
  
  #if(nrow(r)%%2==1){modfft <- modfft[1:nrow(r),]}
  #if(ncol(r)%%2==1){modfft <- modfft[,1:ncol(r)]}
  
  if(class(r)== "RasterLayer"){
    modfft <- raster(modfft) 
    extent(modfft) <- extent(r)
    projection(modfft) <- projection(r)}
  
  modfft
  
}


average.by.station <- function(df,stnfld,var){
  #df <- ti
  #stnfld <- 'latlontime'
  #var <- 'chl'
  # given dataframe and station field
  # returns data frame with var averaged
  require(plyr)
  
  stns <-   unique(df[,stnfld])
  
  #identify stations with duplicate records  
  dupl <- sapply(stns,FUN=function(i)length(df[df[,stnfld]==i,var])>1)
  
  summary(dupl)
  
  stn.mean <- function(s){
    tdf <- df[df[,stnfld]==s,][1,]
    tdf[,var] <- mean(df[df[,stnfld]==s,][,var])
    tdf}
  
  
  dfmean <- ldply(stns[dupl],.fun=stn.mean)
  
  nondupi <- df[,stnfld] %in% stns[!dupl]
  
  newdf <- rbind(dfmean,df[nondupi,])
  
  newdf}







kappa <- function(x,y){
  require(vcd)
  1-Kappa(table(x,y))$Unweighted[1]
}

class.error <- function(xi,yi){
  #xi <- knn.cv(train=x[,BCi],cl=yC,k=8)
  #yi <- yC
  
  t <- table(xi,yi)
  
  1-(sum(diag(t))/sum(t))
  
  
}


class.error.mean <- function(xi,yi){
  #xi <- knn.cv(train=x[,BCi],cl=yC,k=100)
  #yi <- yC
  
  t <- table(xi,yi)
  ii <- sapply(1:nrow(t),FUN=function(i) diag(t)[i]/sum(t[i,]))
  ii[is.na(ii)] <- 0
  mean(1-ii,na.rm=T)
  
  
}

class.error.cmb <- function(xi,yi){
  #xi <- knn.cv(train=x[,BCi],cl=yC,k=100)
  #yi <- yC
  
  t <- table(xi,yi)
  t
  ii <- sapply(1:nrow(t),FUN=function(i) diag(t)[i]/sum(t[i,]))
  ii[is.na(ii)] <- 0
  (max(1-ii)-min(1-ii))*class.error(xi,yi)
  
  
}


knn.cv2 <- function(train,cl,cost,k,n){
  #train <- x[,BCi]
  #cl=yC
  #cost='kappa'
  #k=2
  #n = 100
  cost.call <- call(cost)
  l <- c(1:n)
  ll <- sapply(l,FUN= function(xi) eval(call(cost,x=knn.cv(cl=cl,train=train,k=k),y=cl)))
  as.vector(ll)
}


alr1 <- function(x,y){
  y <- log(x/y)
  y[y==-Inf] <- min(y[y!=-Inf])/2
  y[y==Inf] <- max(y[y!=Inf])*2
  y
}

hclust.centre <- function (data, hfit,nclust) { 
  # returns cluster centers from hclust
  
  clust <- cutree(hfit, k=nclust)
  
  nvars=length(data[1,]) 
  ntypes=max(clust) 
  centroids<-matrix(0,ncol=nvars,nrow=ntypes) 
  for(i in 1:ntypes) { 
    c<-rep(0,nvars) 
    n<-0 
    for(j in names(clust[clust==i])) { 
      n<-n+1 
      c<-c+data[j,] 
    } 
    centroids[i,]<-c/n 
  } 
  rownames(centroids)<-c(1:ntypes) 
  colnames(centroids)<-colnames(data) 
  centroids 
}

kmeans_CH <- function(centers){
  # returns C-H stat from clusters
  km <- kmeans(pcs,centers=centers)
  return(calinhara(pcs,km$cluster))
}

cor.mat <- function(x){
  pairs(x,
        lower.panel=panel.smooth, upper.panel=panel.cor)  
}




alrInv1 <- function(x1,x2){
  
  x1e <- exp(x1)
  x2e <- exp(x2)
  
  rat <- x1e + x2e + 1
  
  c1 <- x1e/rat
  c2 <- x2e/rat
  c3 <- 1 - (c1+c2)
  
  if (class(c1)[1] =="RasterLayer"){return(stack(c1,c2,c3))
  } else {return(cbind(c1,c2,c3))}
  
}

alrInv2 <- function(y1,y2,y3){
  
  #M <- 0.04
  #G <- 0.16
  #S <- 0.8
  #  
  #y1 <- log(M/G) + sample(seq(-0.1,0.1,0.01),1)
  #y2 <- log(S/G) + sample(seq(-0.1,0.1,0.01),1)
  #y3 <- log(M/S) + sample(seq(-0.1,0.1,0.01),1)
  
  r1 <- alrInv1(y1,y2)
  r2 <- alrInv1(-y2,y3)
  r3 <- alrInv1(-y1,-y3)
  
  f1 <- (r1[,1]+r2[,2]+r3[,3])/3
  f2 <- (r1[,2]+r2[,3]+r3[,2])/3
  f3 <- (r1[,3]+r2[,1]+r3[,1])/3
  
  cbind(f1,f2,f3)
}

mcapply <- function(x,fun,cores=NULL){
  require(parallel)
  # clusterApply() for Windows
  if (is.null(cores)){cl <- makeCluster(detectCores()+1)
  } else {cl <- makeCluster(cores)}
  
  runtime <- system.time({
    out <- clusterApply(cl=cl, x=x, fun=fun)
  })[3]
  print(runtime)
  stopCluster(cl) # Don't forget to do this--I frequently do
  return(out)
}


partialPlot2 <- function(rf,y,x,v2,class,n.pt=5){
  #given RF, class and varname
  #returns bivartie partial depeendancy plot
  #(v1 has to be the first var in x )
  
  df <- cbind(y,x)
  v1 <- dimnames(df)[[2]][2]
  
  v2range <- seq(min(df[,v2]),max(df[,v2]),(max(df[,v2])-min(df[,v2]))/(n.pt-1))
  
  v1c<- NULL
  v2c<- NULL
  pc <- c()
  
  for (v in v2range){
    #v <- v2range[1]
    
    df[,v2] <- v
    
    p <- partialPlot(rf,df,x.var=2,class,n.pt=n.pt,plot=F,rug=F)
    
    v1c <- c(v1c,p$x)
    v2c <- c(v2c,rep(v,n.pt))
    pc <- c(pc,p$y)
    
  }
  
  df <- data.frame(cbind(v1c,v2c,pc))
  
  qplot(v1c, v2c, data=df, geom="tile", fill=pc,xlab=v1,ylab=v2,main=class) +
    scale_fill_gradient2(high='#33CC33',mid='#FFFF66',low='#FF3300')
  
}




# Gaussian filter for square cells
f.gau <- function(sigma, n=5) {
  m <- matrix(nc=n, nr=n)
  col <- rep(1:n, n)
  row <- rep(1:n, each=n)
  x <- col - ceiling(n/2)
  y <- row - ceiling(n/2)
  # according to http://en.wikipedia.org/wiki/Gaussian_filter
  m[cbind(row, col)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
  # sum of weights should add up to 1  
  m / sum(m)
}

f.laplace <- matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3)
f.sobelx <- matrix(c(1,2,1,0,0,0,-1,-2,-1), nrow=3)
f.sobely <- matrix(c(-1,0,1,-2,0,2,-1,0,1), nrow=3)
f.mean <- matrix(1/9, nrow=3, ncol=3)




logit <- function(x){log((x/1) / (1-(x/1)))}
rlogit <- function(x){exp(x)/(1+exp(x))}



tplot <- function(text,main){
  #plots text outputs
  
  textplot(capture.output(text))
  title(main)
}


ggplotfunctions <- function(){
  return(source('http://dl.dropbox.com/u/4683316/R/ggpot_functions.r'))}

spatialfunctions <- function(){
  source('http://dl.dropbox.com/u/4683316/R/zspatial_functions.r')
  
}

postgisfunctions <- function(){
  source('http://dl.dropbox.com/u/4683316/R/zpostgis.r')
  
}


FbE <- function(list,ext){
  #list <- dir()
  #ext <- '.shp'
  #given a list of filenames filter
  #returns only matching .ext
  return(grep(paste('*',ext,sep=''),list,value=T))
  
}


right <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

left <- function(x, n){
  substr(x, 0, n)
}


split.dataframe.quantile <- function(data,field,quantiles=10){
  #data <- mb
  #field <-'depth'
  #quantiles=5
  # given data and number of quantiles
  # divides data into specified quatiles
  f <- quantile.factor(data[,field],quantiles)
  s <- split.data.frame(data,f=f)
  return(s)
}



#
quantile.factor <- function(data,quantiles){
  #given a numeric vector
  #returns factor from quantiles
  q <- c(1:quantiles)/quantiles
  
  brks <- c(min(data),quantile(data,q))
  
  b <- cut(data,brks,include.lowest=T,labels=names(brks)[2:length(brks)])
  
  return(b)
  
}


eastness <- function(x){sin((x*pi)/180)}
northness <- function(x){cos((x*pi)/180)}

cfplot <- function(data,xlab='x',ylab='Cumulative distribution'){
  require(ggplot2)
  #data <- data$hrs_km
  
  data <- data.frame(data)
  
  data$ecdf <- ecdf(data[[1]])(data[[1]])
  
  names(data) <- c('x','y')
  
  ggplot(data, aes(x,y)) + geom_step()+
    scale_x_continuous(xlab)+
    scale_y_continuous(ylab)
}



splot <- function(x,y,xlab='x',ylab='y'){
  require(ggplot2)
  df <- data.frame(x=x,y=y)
  p <- ggplot(df,aes(x,y))+geom_point()+
    scale_x_continuous(xlab)+scale_y_continuous(ylab)
  return(p)
}




CleanColNames <- function(df){
  #df <- new.sp
  # given sp/df
  # removes spaces from colnames reduces length to 11
  nmes <- names(df)
  clean <- function(n){
    #n <- names(df)[1]
    n <- substr(gsub(' ','_',n),1,10)
  }  
  names(df) <- sapply(nmes,FUN=clean)
  return(df)
}


scrub.names <- function(data){
  # given a dataframe or spatial df scrubs names 
  # and replaces them with v1,v2....
  # for easy export to shapefile
  names(data) <- sapply(1:length(names(data)),function(x)paste('V',x,sep=''))  
  return(data)
}


return.first <- function(data,field){
  #field <- 'Station Number'
  
  # given a data frame and field name
  # will return the first observation for each unique item in field
  
  uv <- unique(data[,field])
  df <- data.frame()
  for (v in uv[!is.na(uv)]){
    #v <- uv[2]
    print(v)
    df <- rbind(df, data[data[,field]==v,][1,])
  }
  return(df)  
  
}

OpenExcel <- function(file,sheet='Sheet1'){
  require(RODBC)
  channel <- odbcConnectExcel(file)
  sh1 <- sqlFetch(channel, sheet)
  close(channel)
  return(sh1)
}
OpenExcel2007 <- function(file,sheet='Sheet1'){
  require(RODBC)
  channel <- odbcConnectExcel2007(file)
  sh1 <- sqlFetch(channel, sheet)
  close(channel)
  return(sh1)
}



panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  
  text(0.5, 0.5, txt, cex = cex * r)
  text(.8, .8, Signif, cex=cex, col=2)
}







#######################################


rand.part <- function(data,nfold){
  #randomly partitions the dataset into and nfold partitions
  #data <- data1
  #nfold <- 100
  
  rdm <- sample(c(1:nrow(data)))
  psize <- nrow(data)%/%nfold
  
  i=1
  ii = psize
  for (p in 1:nfold){
    data@data[rdm[i:ii],'nfold'] <- p
    i = i + psize
    ii = ii + psize
    if (p == nfold-1){ii = nrow(data)}
  } 
  return(data)
}


rad.t.deg <- function(x){(180/pi)*x}
deg.t.rad <- function(x){(pi/180)*x}

polar.t.cart <- function(r,deg){
  x <- r*cos(deg.t.rad(deg))
  y <- r*sin(deg.t.rad(deg))
  return(cbind(x,y))
}

GD.t.AD <- function(deg){
  con.deg <- function(a){ 
    d = 0
    if (a >= 0 & a <= 90){
      d = 90 - a
    } else if ( a > 90 & a <= 360){ 
      d = 450 - a}
    return(d)
  }
  return(sapply(deg,con.deg))}

positive.degrees<- function(deg){
  pos.deg <- function(d){ 
    if (d < 0) {return(360 +d)} 
    else {return(d)}
  }
  return(sapply(deg,pos.deg))}

negative.degrees<- function(deg){
  pos.deg <- function(d){ 
    if (d > 180) {return(-(360-d))} 
    else {return(d)}
  }
  return(sapply(deg,pos.deg))}



fit.ellipse <- function (x, y = NULL)
{
  # Least squares fitting of an ellipse to point data
  # using the algorithm described in:
  #   Radim Halir & Jan Flusser. 1998.
  #   Numerically stable direct least squares fitting of ellipses.
  #   Proceedings of the 6th International Conference in Central Europe
  #   on Computer Graphics and Visualization. WSCG '98, p. 125-132
  #
  #https://stat.ethz.ch/pipermail/r-help/2010-September/254560.html
  # Adapted from the original Matlab code by Michael Bedward
  # <michael.bedward at gmail.com>
  #
  # Arguments:
  # x, y - the x and y coordinates of the data points
  #
  # Returns a list with the following elements:
  #
  # coef - coefficients of the ellipse as described by the general
  #        quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0
  #
  # center - center x and y
  #
  # major - major semi-axis length
  #
  # minor - minor semi-axis length
  #
  
  EPS <- 1.0e-8
  dat <- xy.coords(x, y)
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y)
  D2<- cbind(dat$x, dat$y, 1)
  
  S1 <- t(D1) %*% D1
  S2 <- t(D1) %*% D2
  S3 <- t(D2) %*% D2
  T <- -solve(S3) %*% t(S2)
  M <- S1 + S2 %*% T
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2)
  evec <- eigen(M)$vec
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2
  a1 <- evec[, which(cond > 0)]
  f <- c(a1, T %*% a1)
  names(f) <- letters[1:6]
  
  # calculate the center and lengths of the semi-axes
  b2 <- f[2]^2 / 4
  center <- c((f[3] * f[4] / 2 - b2 * f[5]), (f[1] * f[5] / 2 - f[2] *
                                                f[4] / 4)) / (b2 - f[1] * f[3])
  names(center) <- c("x", "y")
  
  num <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 -
                f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6])
  den1 <- (b2 - f[1]*f[3])
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2)
  den3 <- f[1] + f[3]
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 *
                                                               (-den2 - den3)) ))
  
  # calculate the angle of rotation
  term <- (f[1] - f[3]) / f[2]
  angle <- atan(1 / term) / 2
  
  list(coef=f, center = center, major = max(semi.axes), minor =
         min(semi.axes), angle = unname(angle))
}



get.ellipse <- function ( fit, n=360 )
{
  # Calculate points on an ellipse described by
  # the fit argument as returned by fit.ellipse
  #
  # n is the number of points to render
  
  tt <- seq(0, 2*pi, length=n)
  sa <- sin(fit$angle)
  ca <- cos(fit$angle)
  ct <- cos(tt)
  st <- sin(tt)
  
  x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa
  y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca
  cbind(x=x, y=y)
}


col.ramp <- function(name='Spectral'){
  # get a nice color ramp
  library(RColorBrewer)
  pal <- brewer.pal(10,name)
  ramp <- colorRampPalette(pal)
  return(ramp(1000))
}

kmeans.plot <- function(data,n=10,main=''){
  # scree plot for wss
  wss <- lapply(1:n, function(i) sum(kmeans(data,centers=i,nstart=5)$withinss)) %>% 
    unlist()
  plot(1:n, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares",main=main)
}

cut.quantile <- function(x,n) {
  cut(x, breaks=quantile(x, probs = seq(0, 1, by = (1/n))),
      labels=c(1:n), include.lowest=TRUE)
}

sdf <- function(spdf){
  cbind(spdf@coords,spdf@data)  
}


#logit <- function(data){
# logit transform, replaces 0 and 100 
#	i = data == 0							#
#	data <- replace(data,i,0.01)			#
#	i = data == 100
#	data <- replace(data,i,99.99)

#	return(log((data/100) / (1-(data/100))))	
#}

#bt <- function(data) {
#back transform logits
# 	return((exp(data)/(1+exp(data)))*100)   
#}

# stats <- function(x){
#   #summary stats including skewness and kurt
#   library(moments)
#   sum <- summary(x)
#   sum[7] <- round(sd(x),3)
#   sum[8] <- round(skewness(x),3)
#   sum[9] <- round(kurtosis(x),3)
#   names(sum)[7:9] <- c('Std.','Skew.','Kurt.')
#   return(sum)
# }
# 
# factor.stats <- function(factors,vars,data){	
#   # returns summary for vars based on factors
#   slist <- list()
#   for (v in vars){
#     print(v)
#     t<- tapply(data[,v],data[,factors],stats)
#     tble <- NULL
#     for(e in t){tble <- rbind(tble,e)}
#     row.names(tble)<- names(t)
#     slist <- append(slist,list(tble))
#   }
#   names(slist) <- vars	
#   return(slist)
# }


day.of.year <- function(data,format)
{return(strptime(data, format)$yday+1)}

moving.average<-function(x){
  y <- numeric(length(x)-2)
  for (i in 2:(length(x)-1)){
    y[i] <-(x[i-1]+x[i]+x[i+1])/3
  }
  return(y)
}


loess.aic <- function (x) {
  if (!(inherits(x,"loess"))) stop("Error: argument must 
                                   be a loess object")
  # extract values from loess object
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum( x$residuals^2 ) / (n-1)
  delta1 <- x$one.delta
  delta2 <- x$two.delta
  enp <- x$enp
  
  aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2)
  #aicc1<- n*log(sigma2) + n* ( (delta1/(delta2*(n+enp)))/(delta1^2/delta2)-2 )
  aicc1<- n*log(sigma2) + n* ( (delta1/delta2)*(n+enp) / (delta1^2 / delta2) -2 )
  gcv  <- n*sigma2 / (n-traceL)^2
  
  result <- list(span=span, aicc=aicc, aicc1=aicc1, gcv=gcv)
  return(result)
}

bestLoess <- function(model, criterion=c("aicc", "aicc1", "gcv"),
                      spans=c(.05, .95)){
  criterion <- match.arg(criterion)
  f <- function(span) {
    mod <- update(model, span=span)
    loess.aic(mod)[[criterion]]
  }
  result <- optimize(f, spans)
  list(span=result$minimum, criterion=result$objective)
}

# RMSSE <- function(obs,pred){
#   rmse <- RMSE(obs-pred)
#   return(round(rmse/sd(obs),3))}
# 
# RMSE <- function(data){
#   return(round(sqrt(mean((data)^2)),3))}
# 
# ME <- function(data){
#   return(round(mean(data),3))}
# 
# MAE <- function(data){
#   return(round(mean(abs(data)),3))}
# 
# MSE <- function(data){
#   return(round(mean(data^2),3))}
# 
# Rsqr <- function(obs,pred){
#   return(round(cor(obs,pred)^2,3))}
# 
# R <- function(obs,pred){
#   return(round(cor(obs,pred),3))}
# 

# split.dataset <- function(dataset,vp){
#   vp_ <- nrow(dataset)%/%100 * vp
#   ii <- seq_len(nrow(dataset))
#   ind1 <- sample(ii, vp_)
#   ind2 <- ii[!ii %in% ind1]
#   valid <- dataset[ind1, ]
#   train <- dataset[ind2, ]
#   rl <- list(valid,train)
#   names(rl) <- c('valid','train')
#   return(rl)
# }

# compare <- function(observed,models,data){
#   # returns data frame with errors calc for each model
#   
#   obs <- data[,observed]
#   #models <- colnames(data[models,])[-1]
#   table <- NULL
#   
#   for (m in models){
#     pred <- data[,m]
#     table <- cbind(table,calc.errs(obs,pred))
#   }
#   
#   rownames(table) <- c('ME','MAE','MSE','r','r2','R.res','RMSE','RMSSE')	
#   colnames(table) <- models
#   
#   return(table)
# }
# 
# calc.errs <- function(obs,pred){
#   res <- obs-pred
#   tble <- rbind(ME(res),MAE(res),MSE(res),R(pred,obs),Rsqr(pred,obs),R(pred,res),RMSE(res),RMSSE(obs,pred))
#   rownames(tble) <- c('ME','MAE','MSE','r','r2','R.res','RMSE','RMSSE')
#   return(tble)
# }




install.packages2 <- function(){
  libs <- c(
    "rpart",
    "nnet",
    "neuralnet",
    "gam",
    "relaimpo",
    "RSAGA",
    "RPyGeo", 
    "maptools",
    "rgdal",
    "sp",
    "ncdf",
    "raster",
    "RCurl",
    "rgeos",
    "spatstat",
    "stringr",
    "rgl",
    "moments",
    "RColorBrewer",
    "classInt",
    "RODBC",
    "googleVis",
    "ggplot2",
    "gpclib",
    "maps",
    "mapproj",
    "vcd",
    "mda",
    "mgcv",
    "Ellipse",
    "papply",
    "igraph",
    "Hmisc",
    "TSP",
    "EBImage",
    "abind",
    "sm",
    "randomForest",
    "lattice",
    "twitteR",
    "ROAuth",
    "Biostrings",
    "bootstrap",
    "car",
    "RNiftyReg",
    "oro.nifti",
    "fpc",
    "plotGoogleMaps",
    "adimpro",
    "compositions",
    "epicalc",
    "geoR",
    "mapdata",
    "spmaps",
    "geosphere",
    "maxent",
    "quantmod",
    "tree",
    "prp",
    "spdep",
    "gridExtra",
    "reshape2",
    "lubridate" ,
    "shiny",
    "test",
    "ROCR", 
    "randomFields",
    "dplyr",
    "ggmap",
    "ade4",
    "ggtern",
    "leaflet",
    "shinyRGL",
    "shinyIncubator",
    "devtools",
    "shinyAce",
    "leaflet",
    "animation",
    'shinyBS',
    'shinyjs',
    'shinythemes',
    'shinybootstrap2',
    'ggvis')

  install.packages(libs)


}




##R --arch x86_64 CMD INSTALL --configure-args="-with-gdal-config=/Library/Frameworks/GDAL.framework/unix/bin/gdal-config -with-proj-lib=/Library/Frameworks/PROJ.framework/unix/lib/ -with-proj-include=/Library/Frameworks/PROJ.framework/unix/include/" rgdal_0.6-33.tar.gz
##R --arch x86_64 CMD INSTALL --configure-args="-with-geos-config=/Library/Frameworks/GEOS.framework/unix/bin/geos-config" rgeos_0.1-4.tar.gz


####Gggplot#############



ggROC = function(Yp,Yv,n=100){
  
  ROCp = performance(prediction(Yp,(Yv)), "tpr","fpr")
  plot(ROCp,colorize=T)
  
  ROCdf =data.frame('FP'=ROCp@x.values[[1]],'TP'=ROCp@y.values[[1]]) 
  names(ROCdf) = c('FP','TP')
  qdf = data.frame(Th=c(0.25,0.5,0.75),
                   FP=approx(x=ROCp@alpha.values[[1]],y=ROCp@x.values[[1]],xout=c(0.25,0.5,0.75))$y,
                   TP=approx(x=ROCp@alpha.values[[1]],y=ROCp@y.values[[1]],xout=c(0.25,0.5,0.75))$y)
  
  simdf = NULL
  for(i in 1:n){
    ROCp = performance(prediction(Yp,sample(Yv)), "tpr","fpr")
    simdf =rbind(simdf,data.frame(sim=i,'FP'=ROCp@x.values[[1]],'TP'=ROCp@y.values[[1]])) 
  }
  
  print(ggplot(simdf,aes(x=FP,y=TP,group=sim))+
          geom_line(alpha=0.1)+
          geom_line(data=ROCdf,aes(x=FP,y=TP,group=NULL),col=2)+
          geom_point(data=qdf,aes(x=FP,y=TP,group=NULL),col=2))}


fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- format(l, scientific = TRUE) 
  # quote the part before the exponent to keep all the digits 
  l <- gsub("^(.*)e", "'\\1'e", l) 
  # turn the 'e+' into plotmath format 
  l <- gsub("e", "%*%10^", l) 
  # return this as an expression 
  parse(text=l) 
} 

Boruta.plot <- function(BC){
  require(ggplot2)
  require(plyr)
  
#   require('Boruta')
#   iris.extended<-data.frame(iris,apply(iris[,-5],2,sample));
#   names(iris.extended)[6:9]<-paste("Nonsense",1:4,sep="");
#   #Run Boruta on this data
#   Boruta(Species~.,data=iris.extended,doTrace=2)->BC
#   
  BCdf <- BC$ImpHistory
  
  Zscores <- as.vector(BCdf)
  
  feature <- as.vector(sapply(dimnames(BCdf)[[2]],FUN=rep,times=nrow(BCdf)))
  
  df <- data.frame(Zscores,feature)
  
  ggplot(df,aes(x=feature,y=Zscores))+geom_boxplot()
  
  as.vector(BCdf)
  
}
# 
# ggbase <- function(data,zoom=0.2,coast=T,ext=NULL){
#   #given data or extent
#   # returns a base plot with coast
#   
#   if(is.null(ext)) {ext <- extent(data)}
#   ext.exp <- extent.expand(ext,zoom)
#   xmin <- ext.exp@xmin
#   xmax <- ext.exp@xmax
#   ymin <- ext.exp@ymin
#   ymax <- ext.exp@ymax
#   
#   x <- data@coords[,1]
#   y <- data@coords[,2]
#   
#   
#   xbreaks <- xybreaks(x,xmin,xmax,3)
#   ybreaks <- xybreaks(y,ymin,ymax,3)
#   
#   #ybreaks <- classIntervals(y,n=5,style='pretty')$brks
#   
#   #xbreaks <- seq(ceiling(xmin),floor(xmax),0.1)
#   #ybreaks <- seq(ceiling(ymin),floor(ymax),0.1)
#   
#   c <- NULL
#   if (coast==T){c <- coast.geom(xmin,xmax,ymin,ymax)}              
#   
#   qplot(data@data)+
#     coord_map(ylim=c(ymin,ymax),xlim=c(xmin,xmax),fast=T)+ 
#     scale_x_continuous(breaks=xbreaks, labels=xbreaks,name='x')+
#     scale_y_continuous(breaks=ybreaks, labels=ybreaks,name='y')+c+
#     theme_bw()
#   
#   
# }
# 
# ggClass <- function(data,groupvar=NULL,ext=NULL,label=NULL,
#                     coast=T,main=NULL,pal='Spectral',zoom=0.2){
#   
#   #  ext=NULL
#   #  label=NULL
#   #  main =NULL
#   #  coast = F
#   #  zoom=0.2
#   #  data <- SURVEY NAM
#   #  groupvar <- '5#Biotope (EUNIScode)'
#   #  pal='Spectral'
#   
#   if (is.null(main)){main <- groupvar}
#   
#   label.geom <- NULL
#   if(!is.null(label)){l <- as.character(data@data[,label])
#                       label.geom = geom_text(aes(x+0.08,y+0.08,label=l))}
#   
#   p <- ggbase(data=data,zoom=zoom,coast=coast)  
#   
#   p+opts(title=main)+
#     geom_point(aes(x=x,y=y,colour = factor(data@data[,groupvar])))+
#     scale_color_brewer(name=groupvar,palette=pal)+
#     label.geom
#   
# }

# ggBubble <- function(data,sizevar,sizebrks=NULL,colvar=NULL,colbrks=NULL,ext=NULL,label=NULL,
#                      coast=T,main=NULL,lowcol='#66CC66',highcol='#CC0000',zoom=0.2){
#   
#   #label=NULL
#   #main =NULL
#   #sizebrks=NULL
#   #colvar=NULL
#   #colbrks = NULL
#   #coast = F
#   #zoom=0.2
#   #data <- stns
#   #sizevar <- 'test'
#   
#   
#   if(is.null(sizebrks)){sizebrks <- classIntervals(data@data[,sizevar],n=5,style='pretty')$brks}
#   
#   if(is.null(colvar)){colvar <- sizevar
#                       colbrks <- sizebrks}
#   
#   if(is.null(colbrks)){colbrks <- classIntervals(data@data[,colvar],n=5,style='pretty')$brks}
#   
#   if (is.null(main)){main <- sizevar}
#   
#   label.geom <- NULL
#   if(!is.null(label)){l <- as.character(data@data[,label])
#                       label.geom = geom_text(aes(x+0.08,y+0.08,label=l))}
#   
#   p <- ggbase(data=data,zoom=zoom,coast=coast)  
#   
#   collim <- c(min(colbrks),max(colbrks))
#   szelim <- c(min(sizebrks),max(sizebrks))
#   
#   p+opts(title=main)+
#     geom_point(aes(x=x,y=y,size =data@data[,sizevar],colour = data@data[,colvar]))+
#     scale_area(name=sizevar,to=c(2,10),breaks=sizebrks,limits=szelim)+  
#     scale_color_continuous(name=colvar,low=lowcol,high=highcol,breaks=colbrks,limits=collim)+
#     label.geom
#   
# }

# ggkrige <- function(data,var,brks=NULL,ext=NULL,pred=NULL,res=0.01,kv.alpha =F,
#                     coast=T,main=NULL,lowcol='#66CC66',highcol='#CC0000',zoom=0.2){
#   # Given a SPDF and value field
#   # returns ggplot with OK estimations
#   
#   #brks <- y1v1.c$brks
#   #data <- cruise.list[[1]]
#   #var <- var.names[1]
#   
#   
#   if(is.null(pred)){pred <- quick.krige(data,var,res=res)}
#   
#   grid.df <- as.data.frame(rasterToPoints(pred[[2]]))
#   grid.geom <- geom_tile(data=grid.df,aes(x=x,y=y,fill=layer.1))
#   
#   kv = NULL
#   if (kv.alpha == T){
#     grid.df$kv <- 1-(grid.df$layer.2 / var(data@data[,var],na.rm=T))
#   }
#   
#   
#   if(is.null(ext)) {ext <- extent(data)}
#   ext.exp <- extent.expand(ext,zoom)
#   xmin <- ext.exp@xmin
#   xmax <- ext.exp@xmax
#   ymin <- ext.exp@ymin
#   ymax <- ext.exp@ymax
#   
#   if (is.null(main)){main = paste('OK interpolation of ',var,sep='')}
#   
#   x <- data@coords[,1]
#   y <- data@coords[,2]
#   
#   xbreaks <- seq(ceiling(xmin),floor(xmax),1)
#   ybreaks <- seq(ceiling(ymin),floor(ymax),1)
#   
#   c=NULL
#   if (coast==T){c <- coast.geom(xmin,xmax,ymin,ymax)}
#   if (is.null(brks)){limits=NULL
#   } else {limits=c(min(brks),max(brks))}
#   
#   qplot(x,y, data=grid.df,geom='tile',fill=layer.1,alpha=kv)+
#     coord_map(ylim=c(ymin,ymax),xlim=c(xmin,xmax))+
#     #scale_fill_brewer(name=var,breaks=y1v1.c$brks)+
#     scale_fill_gradient(name=var,low=lowcol,high=highcol,breaks=brks,limits=limits)+
#     #coord_map(ylim=c(ymin,ymax),xlim=c(xmin,xmax),fast=T)+ 
#     scale_x_continuous(breaks=xbreaks, labels=xbreaks)+
#     scale_y_continuous(breaks=ybreaks, labels=ybreaks)+
#     c+
#     opts(title=main) + theme_bw()
# }

coast.geom <- function(xlim,ylim){
  require(mapdata)
  coast <- map_data("worldHires", xlim = xlim,ylim=ylim)
  #coast.fort <- fortify(coast)
  coast.geom <- geom_polygon(data=coast.fort, aes(x=long, y=lat, group=group), colour= "#999999", fill="#999999", lwd=0.2)
  
  
  return(coast.geom)
  
  
}

# multiplot <- function(..., plotlist=NULL, cols) {
#   require(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # Make the panel
#   plotCols = cols                          # Number of columns of plots
#   plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
#   
#   # Set up the page
#   grid.newpage()
#   pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
#   vplayout <- function(x, y)
#     viewport(layout.pos.row = x, layout.pos.col = y)
#   
#   # Make each plot, in the correct location
#   for (i in 1:numPlots) {
#     curRow = ceiling(i/plotCols)
#     curCol = (i-1) %% plotCols + 1
#     print(plots[[i]], vp = vplayout(curRow, curCol ))
#   }
#   
# }



xybreaks <- function(xy,min,max,n=3){
  xy <- seq(min,max,1)
  breaks <- classIntervals(xy,n=n+2,style='pretty')$brks
  
  if (tail(breaks,1) > max){breaks <- breaks[-length(breaks)]}
  if (head(breaks,1) < min){breaks <- breaks[-1]}
  
  return(breaks)
}
