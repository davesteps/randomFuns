

#' folderStructure
#'
#' Create basic folder structure
#' @param fld Directory
#' @examples
#' folderStructure()
#' @export
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



#' didYouMean
#'
#' returns google suggestion for a misspelled word - doesn't seem to work though
#' @param input string
#'
#' @return sting
#' @export
#'
#' @examples
#' didYouMean('arse')
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

# permute <- function(v1,v2,n=1e3){
#   obs.diff <- abs(mean(v1)-mean(v2))
#   
#   vc <- c(v1,v2)
#   l <- length(vc)
#   c=0
#   for(p in 1:n){
#     si1 <- sample(1:l,size=length(v1),replace=F)
#     v1p <- vc[si1]
#     v2p <- vc[!((1:l)%in%si1)]
#     if((mean(v1p)-mean(v2p))>obs.diff){c=c+1}}
#   c/n}


#' permuteJB
#'
#' Jon Barrys original function for non-parametric sig test between means
#'
#' @param v1 numeric vector
#' @param v2 numeric vector
#' @param nreps number of imputations
#'
#' @return p-value
#' @export
#'
#' @examples
#' permuteJB(rnorm(100,0,sd = 1),rnorm(100,1,sd = 1))
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





#' I0
#'
#' convert Wm-2 to mE(?)
#' @param IR Surface Irradiance
#'
#' @return converted value
#' @export
#'
#' @examples
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

#' lm_eqn
#'
#' @param df data frame with x and y vars
#'
#' @return
#' @export
#'
#' @examples
lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

#' Zeu
#'
#' Convert Kd to Zeu(1% light depth)
#' @param kd 
#'
#' @return
#' @export
#'
#' @examples
Zeu <- function(kd){
  ### to prepare Kd:
  # calculation of depth of photic zone as depth of 1% of I0
  # Zeu = 4.61/Kd
  4.61/kd}


#' PP
#' given vectors of kf,chl and IR, calculates primary Pproduction
#' @param IR 
#' @param kd 
#' @param chl 
#'
#' @return
#' @export
#' @seealso Cloern 1987
#' @examples
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




#' lmp
#' No Idea what this does
#'
#' @param modelobject 
#'
#' @return
#' @export
#'
#' @examples
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


#' margin
#' Returns the decision margin from a majority vote calssifier
#'
#' @param mat matrix of votes from classifier
#'
#' @return
#' @export
#'
#' @examples
margin <- function(mat){
  # given a matrix of probabilities
  # returns the confidence margon
  c <- ncol(mat)
  apply(mat,1,function(i) sort((i))[c]-sort((i))[c-1])}




#' swap.quadrant
#' inverts the four quadrants of a matrix
#'
#' @param ft matrix
#'
#' @return
#' @export
#'
#' @examples
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



#' plotFF
#' plots the fourier transformation of an image
#'
#' @param r raster or matrix
#' @param returnFF 
#'
#' @return
#' @export
#'
#' @examples
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

#' imageFF
#' Apply a 'pie slice' filter to fourier transform of image
#'
#' @param r raster
#' @param d direction
#' @param l radius
#' @param t taper
#' @param lwl_buf internal buffer radius 
#' @param bd buffer width
#'
#' @return
#' @export
#'
#' @examples
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







#' kappa
#' Unweighted kappa coefficient
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
kappa <- function(x,y){
  
  1-vcd::Kappa(table(x,y))$Unweighted[1]
}

#' CE
#'
#' @param xi 
#' @param yi 
#'
#' @return
#' @export
#'
#' @examples
CE <- function(xi,yi){
  #xi <- knn.cv(train=x[,BCi],cl=yC,k=8)
  #yi <- yC
  t <- table(xi,yi)
  1-(sum(diag(t))/sum(t))
}


#' MCE 
#'
#' @param xi 
#' @param yi 
#'
#' @return
#' @export
#'
#' @examples
MCE <- function(xi,yi){
  #xi <- knn.cv(train=x[,BCi],cl=yC,k=100)
  #yi <- yC
  
  t <- table(xi,yi)
  ii <- sapply(1:nrow(t),FUN=function(i) diag(t)[i]/sum(t[i,]))
  ii[is.na(ii)] <- 0
  mean(1-ii,na.rm=T)
  
  
}

#' MCE2
#'
#' @param xi 
#' @param yi 
#'
#' @return
#' @export
#'
#' @examples
MCE2 <- function(xi,yi){
  #xi <- knn.cv(train=x[,BCi],cl=yC,k=100)
  #yi <- yC
  
  t <- table(xi,yi)
  t
  ii <- sapply(1:nrow(t),FUN=function(i) diag(t)[i]/sum(t[i,]))
  ii[is.na(ii)] <- 0
  (max(1-ii)-min(1-ii))*CE(xi,yi)
  
  
}


#' knn.cv2
#'
#' @param train 
#' @param cl 
#' @param cost 
#' @param k 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
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


#' alr
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
alr <- function(x,y){
  y <- log(x/y)
  y[y==-Inf] <- min(y[y!=-Inf])/2
  y[y==Inf] <- max(y[y!=Inf])*2
  y
}

#' hclust.centre
#'
#' @param data 
#' @param hfit 
#' @param nclust 
#'
#' @return
#' @export
#'
#' @examples
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

#' kmeans_CH
#'
#' @param centers 
#'
#' @return
#' @export
#'
#' @examples
kmeans_CH <- function(centers){
  # returns C-H stat from clusters
  km <- kmeans(pcs,centers=centers)
  return(calinhara(pcs,km$cluster))
}

#' cor.mat
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
cor.mat <- function(x){
  pairs(x,
        lower.panel=panel.smooth, upper.panel=panel.cor)  
}




#' alrInv
#' Inverse alr
#'
#' @param x1 
#' @param x2 
#'
#' @return
#' @export
#'
#' @examples
alrInv <- function(x1,x2){
  
  x1e <- exp(x1)
  x2e <- exp(x2)
  
  rat <- x1e + x2e + 1
  
  c1 <- x1e/rat
  c2 <- x2e/rat
  c3 <- 1 - (c1+c2)
  
  if (class(c1)[1] =="RasterLayer"){return(stack(c1,c2,c3))
  } else {return(cbind(c1,c2,c3))}
  
}

#' mcapply
#'
#' @param x 
#' @param fun 
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
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


#' partialPlotBi
#'
#' @param rf 
#' @param y 
#' @param x 
#' @param v2 
#' @param class 
#' @param n.pt 
#'
#' @return
#' @export
#'
#' @examples
partialPlotBi <- function(rf,y,x,v2,class,n.pt=5){
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





#' f.gau
#'
#' Gaussian filter for square cells
#' @param sigma 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
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




#' logit
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
logit <- function(x){log((x/1) / (1-(x/1)))}
rlogit <- function(x){exp(x)/(1+exp(x))}



#' fltByExt
#'
#' @param list 
#' @param ext 
#'
#' @return
#' @export
#'
#' @examples
fltByExt <- function(list,ext){
  #list <- dir()
  #ext <- '.shp'
  #given a list of filenames filter
  #returns only matching .ext
  return(grep(paste('*',ext,sep=''),list,value=T))
  
}


#' rght
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
rght <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#' lft
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
lft <- function(x, n){
  substr(x, 0, n)
}



#' eastness
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
eastness <- function(x){
  sin((x*pi)/180)}

#' northness
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
northness <- function(x){
  cos((x*pi)/180)}




#######################################



#' radToDeg
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
radToDeg <- function(x){(180/pi)*x}
#' degToRad
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
degToRad <- function(x){(pi/180)*x}

#' polar.t.cart
#'
#' @param r 
#' @param deg 
#'
#' @return
#' @export
#'
#' @examples
polar.t.cart <- function(r,deg){
  x <- r*cos(deg.t.rad(deg))
  y <- r*sin(deg.t.rad(deg))
  return(cbind(x,y))
}




#' fit.ellipse
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
fit.ellipse <- function (x, y = NULL){
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



#' get.ellipse
#'
#' @param fit 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
get.ellipse <- function( fit, n=360 ){
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


#' col.ramp
#'
#' @param name 
#'
#' @return
#' @export
#'
#' @examples
col.ramp <- function(name='Spectral'){
  # get a nice color ramp
  library(RColorBrewer)
  pal <- brewer.pal(10,name)
  ramp <- colorRampPalette(pal)
  return(ramp(1000))
}

#' kmeans.plot
#'
#' @param data 
#' @param n 
#' @param main 
#'
#' @return
#' @export
#'
#' @examples
kmeans.plot <- function(data,n=10,main=''){
  # scree plot for wss
  wss <- lapply(1:n, function(i) sum(kmeans(data,centers=i,nstart=5)$withinss)) %>% 
    unlist()
  plot(1:n, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares",main=main)
}

#' cut.quantile
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
cut.quantile <- function(x,n) {
  cut(x, breaks=quantile(x, probs = seq(0, 1, by = (1/n))),
      labels=c(1:n), include.lowest=TRUE)
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


#' day.of.year
#'
#' @param data 
#' @param format 
#'
#' @return
#' @export
#'
#' @examples
day.of.year <- function(data,format)
{return(strptime(data, format)$yday+1)}

#' moving.average
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
moving.average<-function(x){
  y <- numeric(length(x)-2)
  for (i in 2:(length(x)-1)){
    y[i] <-(x[i-1]+x[i]+x[i+1])/3
  }
  return(y)
}


#' loess.aic
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
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

#' bestLoess
#'
#' @param model 
#' @param criterion 
#' @param spans 
#'
#' @return
#' @export
#'
#' @examples
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









#' ggROC
#'
#' @param Yp 
#' @param Yv 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
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


#' fancy_scientific
#'
#' @param l 
#'
#' @return
#' @export
#'
#' @examples
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

#' ggBoruta
#'
#' @param BC 
#'
#' @return
#' @export
#'
#' @examples
ggBoruta <- function(BC){

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



#' coast.geom 
#'
#' @param xlim 
#' @param ylim 
#'
#' @return
#' @export
#'
#' @examples
coast.geom <- function(xlim,ylim){
  require(mapdata)
  coast <- map_data("worldHires", xlim = xlim,ylim=ylim)
  #coast.fort <- fortify(coast)
  coast.geom <- geom_polygon(data=coast.fort, aes(x=long, y=lat, group=group), colour= "#999999", fill="#999999", lwd=0.2)
  
  
  return(coast.geom)
  
  
}





#' xybreaks
#'
#' @param xy 
#' @param min 
#' @param max 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
xybreaks <- function(xy,min,max,n=3){
  xy <- seq(min,max,1)
  breaks <- classIntervals(xy,n=n+2,style='pretty')$brks
  
  if (tail(breaks,1) > max){breaks <- breaks[-length(breaks)]}
  if (head(breaks,1) < min){breaks <- breaks[-1]}
  
  return(breaks)
}
