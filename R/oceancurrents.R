FrontStats <- function(sst,plot=FALSE,plot.int=100,maxhar=6,neg.gradient=TRUE) {
    ny <- length(sst$lat); nx <- length(sst$lon); nyy <- ny-1; nt <- length(sst$tim)
    lats <- 0.5*(sst$lat[1:nyy]+sst$lat[2:ny]); x <- seq(-1,1,length=length(sst$lon)); y <- rep(NA,nx)
    coefs <- matrix(rep(NA,12*nt),12,nt)
    rownames(coefs) <- c("Year","Month","Day",paste(rep("beta",7),1:7,sep="_"),"R-squared","q_0.95 (dt/dy)")
    attr(coefs,'longitudes') <- range(sst$lon)
    attr(coefs,'latitudes') <- range(sst$lat)
    if (plot) {oldpar <- par(); x11(); par(mfrow=c(2,1)) }

    for (ii in 1:length(sst$tim)) {
      sstmap <- t(sst$dat[ii,,])
      dsstdy1 <- dY(sst$lon,sst$lat,sstmap,maxhar=maxhar)
      dsstdy2 <- dY(sst$lon,sst$lat,dsstdy1$dZ,maxhar=maxhar)
      dy11 <- dsstdy1$dZ[,1:nyy]; dy12 <- dsstdy1$dZ[,2:ny]
      dy21 <- dsstdy2$dZ[,1:nyy]; dy22 <- dsstdy2$dZ[,2:ny]
      i.low <-  (dy21*dy22 < 0) 
      i.low[,c(1:3,(nyy-2):nyy)] <- FALSE
      dsstdy1$dZ[,c(1:4,(ny-3):ny)] <- NA
      for (iii in 1:nx) {
        srt <- order(0.5*(dy11[iii,i.low[iii,]]+dy12[iii,i.low[iii,]]),decreasing=!neg.gradient)
        y[iii] <-  lats[i.low[iii,]][srt][1]
      }
      fit <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6))
      if ( (plot) & (mod(ii,plot.int)==0) ) {
        image(sst$lon,sst$lat,sstmap)
        contour(sst$lon,sst$lat,dsstdy1$Z.fit,add=TRUE)
        points(sst$lon,y,col="grey95",cex=0.5)
        lines(sst$lon,predict(fit),lwd=3,col="grey70")

        plot(sst$lat,sstmap[round(nx/2),],type="l",lwd=3,col="grey",main=sst$lon[ii],sub=ii)
        lines(sst$lat,dsstdy1$Z.fit[round(nx/2),],col="red")
        lines(rep(y[round(nx/2)],2),range(sstmap[round(nx/2),],na.rm=TRUE),lty=2,col="darkred")
      }
      stats <- summary(fit)
      coefs[,ii] <- c(sst$yy[ii],sst$mm[ii],sst$dd[ii],
                      as.numeric(fit$coefficients),stats$r.squared,quantile(dsstdy1$dZ,0.95,na.rm=TRUE))
      print(round(as.numeric(coefs[,ii],2)))
    }
    if (plot) par(oldpar)
    invisible(coefs)
}


Gulfstream <- function(fname="sst.wkmean",plot=FALSE,maxhar=12,plot.int=100) {
  list <- list.files(pattern=fname)
  list <- list[grep(".nc",list)]
  for (i in 1:length(list)) {
    sst <- retrieve.nc(list[i],x.rng=c(-80,-30),y.rng=c(25,55),v.nam="sst")
    coefs <- FrontStats(sst,plot=plot,maxhar=maxhar,plot.int=plot.int)
    if (i==1) {
      Coefs <- coefs
    } else {
      Coefs <- cbind(Coefs,coefs)
    }
  }
  invisible(Coefs)
}


Kuroshio<- function(fname="sst.wkmean",plot=FALSE,maxhar=12,plot.int=100) {
  list <- list.files(pattern=fname)
  list <- list[grep(".nc",list)]
  for (i in 1:length(list)) {
    sst <- retrieve.nc(list[i],x.rng=c(110,180),y.rng=c(25,55),v.nam="sst")
    coefs <- FrontStats(sst,plot=plot,maxhar=maxhar)
    if (i==1) {
      Coefs <- coefs
    } else {
      Coefs <- cbind(Coefs,coefs)
    }
  }
  invisible(Coefs)
}

CurrentClimatology <- function(current="gulfstream",parameter=1) {
  Currentname <- eval(parse(text=names(formals(CurrentHistory))[1]))
  if (is.character(current))  {
    eval(parse(text=paste("data(",current,",envir=environment())")))
    current <- eval(parse(text=current))
  }
  main <- switch(as.character(parameter),
                 "1"="latitude (reg.coef. 0)","2"="dy/dx (reg.coef. 1)", "3"="d2y/dx2 (reg.coef. 2)",
                 "4"="d3y/dx3 (reg.coef. 3)","5"="d4y/dx4 (reg.coef. 4)","6"="d5y/dx5 (reg.coef. 5)",
                 "7"="d6y/dx6 (reg.coef. 6)","8"="R-squared","9"="quantile(dT/dY,0.95)")
  main=paste(Currentname,main,sep=": ")
  plot((current[2,]-1)/12+current[3,]/365.25,current[parameter+3,],pch=19,cex=0.5,
       main=main,xlab="month",ylab="",sub="cyclones: statistics from FrontStats")
  grid()

}
  
CurrentHistory <- function(current="gulfstream",parameter=1) {
#  print(formals(CurrentHistory))
#  print(eval(parse(text=names(formals(CurrentHistory))[1])))
  Currentname <- eval(parse(text=names(formals(CurrentHistory))[1]))
  if (is.character(current))  {
    eval(parse(text=paste("data(",current,",envir=environment())")))
    current <- eval(parse(text=current))
  }
  main <- switch(as.character(parameter),
                 "1"="latitude (reg.coef. 0)","2"="dy/dx (reg.coef. 1)", "3"="d2y/dx2 (reg.coef. 2)",
                 "4"="d3y/dx3 (reg.coef. 3)","5"="d4y/dx4 (reg.coef. 4)","6"="d5y/dx5 (reg.coef. 5)",
                 "7"="d6y/dx6 (reg.coef. 6)","8"="R-squared","9"="quantile(dT/dY,0.95)")
  main=paste(Currentname,main,sep=": ")
  plot(current[1,]+(current[2,]-1)/12+current[3,]/365.25,current[parameter+3,],type="l",
       main=main,xlab="years",ylab="",sub="cyclones: statistics from FrontStats",col="grey40")
  grid()
  points(current[1,]+(current[2,]-1)/12+current[3,]/365.25,current[parameter+3,],pch=19,cex=0.3)
}
  
  
CurrentSpectrum <- function(current="gulfstream",parameter=1,frequency=FALSE) {
  Currentname <- eval(parse(text=names(formals(CurrentHistory))[1]))
  if (is.character(current))  {
    eval(parse(text=paste("data(",current,",envir=environment())")))
    current <- eval(parse(text=current))
  }
  main <- switch(as.character(parameter),
                 "1"="latitude (reg.coef. 0)","2"="dy/dx (reg.coef. 1)", "3"="d2y/dx2 (reg.coef. 2)",
                 "4"="d3y/dx3 (reg.coef. 3)","5"="d4y/dx4 (reg.coef. 4)","6"="d5y/dx5 (reg.coef. 5)",
                 "7"="d6y/dx6 (reg.coef. 6)","8"="R-squared","9"="quantile(dT/dY,0.95)")
  main=paste(Currentname,main,sep=": ")
  t <- current[1,]+(current[2,]-1)/12+current[3,]/365.25
  S <- spectrum(current[parameter+3,])
  if (frequency) {S$freq <- S$freq/diff(t)[1]; xylog <- "y"; xlab="1/years"} else
                 {S$freq <- diff(t)[1]/S$freq; xylog <- "xy"; xlab="years"}
  plot(S$freq,S$spec,type="l",
       main=main,xlab=xlab,ylab="",sub="cyclones: statistics from FrontStats",col="grey40",log=xylog)
  grid()
}


Current2monthly <- function(current="gulfstream",parameter=1) {
  Currentname <- eval(parse(text=names(formals(CurrentHistory))[1]))
  if (is.character(current))  {
    eval(parse(text=paste("data(",current,",envir=environment())")))
    current <- eval(parse(text=current))
  }
  main <- switch(as.character(parameter),
                 "1"="latitude (reg.coef. 0)","2"="dy/dx (reg.coef. 1)", "3"="d2y/dx2 (reg.coef. 2)",
                 "4"="d3y/dx3 (reg.coef. 3)","5"="d4y/dx4 (reg.coef. 4)","6"="d5y/dx5 (reg.coef. 5)",
                 "7"="d6y/dx6 (reg.coef. 6)","8"="R-squared","9"="quantile(dT/dY,0.95)")

  yy <- as.numeric(rownames(table(current[1,]))); ny <- length(yy)
  
  dat <- matrix(rep(NA,12*ny),ny,12)
  
  for (iy in 1:ny) {
    for (im in 1:12) {
      iym <- is.element(current[1,],yy[iy]) & is.element(current[2,],im)
      if (sum(iym)>0) dat[iy,im] <- mean(current[parameter+3,iym],na.rm=TRUE)
    }
  }
  lon <- switch(Currentname,"gulfstream"=-55,"kuroshio"=145,otherwise=NULL)
  lat <- switch(Currentname,"gulfstream"=40,"kuroshio"=40,otherwise=NULL)
  S <- station.obj(x=dat,yy=yy,obs.name="Current statistics",unit="dimensionless",
                   ele=-1,location=Currentname,lon=lon,lat=lat,alt=0,ref="cyclones: statistics from FrontStats")
  
  plotStation(S,type="b",col="grey40",pch=19,
       main=main,xlab=xlab,ylab="",sub="cyclones: statistics from FrontStats")
  grid()
  invisible(S)
}

CurrentLatitude <- function(current="gulfstream",lon=NULL,nx=50,plot=TRUE) {
  Currentname <- eval(parse(text=names(formals(CurrentHistory))[1]))
  if (is.character(current))  {
    eval(parse(text=paste("data(",current,",envir=environment())")))
    current <- eval(parse(text=current))
  }
  Lon <- switch(Currentname,"gulfstream"=c(-80,-30),"kuroshio"=c(110,180),otherwise=NULL)
  x <- seq(-1,1,length=nx)
  if (!is.null(Lon)) lons <- seq(min(Lon),max(Lon),length=nx) else lons <- x
  if (!is.null(lon)) {
    if (length(lon)==1) {
     ii <- min( (1:50)[lons >= lon] )
     x <- x[ii]; lons <- lons[ii]
  } else if (length(lon)==2) {
     ii <-  max( (1:50)[lons <= lon[1]] ):min( (1:50)[lons >= lon[2]] )
     x <- x[ii]; lons <- lons[ii]
   }
  }
  nt <- length(current[1,])
  y <- matrix(rep(NA,nt*length(x)),length(x),nt)
  for (i in 1:nt) {
    y[,i] <- current[4,i]+ current[5,i]*x   + current[6,i]*x^2 + current[7,i]*x^3 +
                           current[8,i]*x^4 + current[9,i]*x^5 + current[10,i]*x^6
  }
  y[c(1:5,45:50),] <-  NA
  breaks <- seq(min(y,na.rm=TRUE),max(y,na.rm=TRUE),length=31)
  if (plot) {
    image(lons,current[1,]+(current[2,]-1)/12+current[3,]/365.25,y,breaks=breaks,col=rainbow(30),
                  main=paste(Currentname,"latitude excursions"),xlab="longitudes",ylab="Time")
    grid()
    x11()
    plot(lons,rowMeans(y),xlab="longitude",ylab="mean latitude")
    grid()
    x11()
    plot(current[1,]+(current[2,]-1)/12+current[3,]/365.25,colMeans(y),xlab="year",ylab="zonal mean latitude")
    grid()
  }
  results <- list(current=Currentname,lons=lons,year=current[1,]+(current[2,]-1)/12+current[3,]/365.25,y=y)
  invisible(results)
}

CurrentValidation <-  function(current="gulfstream",fname="sst.wkmean",plot=TRUE,las=1) {
  if (Sys.info()[1]!="Linux") print("This works completely only on Linux platforms")
  ok.convert <- system("which convert"); ok.gifmerge <- system("which gifmerge")
  if (ok.convert>0) print("ImageMagic should be installed for animating the graphics") 
  if (ok.gifmerge>0) print("gifmerge should be installed for animating the graphics") 
  Currentname <- eval(parse(text=names(formals(CurrentHistory))[1]))
  y <- CurrentLatitude(Currentname,plot=FALSE)
  x11()
  
  lats <- switch(Currentname,"gulfstream"=c(25,55),"kuroshio"=c(0,70),otherwise=NULL)
  list <- list.files(pattern=fname)
  list <- list[grep(".nc",list)]
  iii <- 1
  x11(); oldpar <- par()
  for (i in 1:length(list)) {
    sst <- retrieve.nc(list[i],x.rng=range(y$lons),y.rng=lats,v.nam="sst")

    for (ii in 1:length(sst$tim)) {
      map <- t(sst$dat[ii,,])
      
      filled.contour(sst$lon,sst$lat,map,main=paste(ii,sst$yy[ii],sst$mm[ii],sst$dd[ii],sep="-"),sub=Currentname)
      mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
      on.exit(par(par.orig))
      w <- (3 + mar.orig[2]) * par('csi') * 2.54
      layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w))) 
      par(las = las)
      mar <- mar.orig
      mar[4] <- 1
      par(mar=mar)

      addland()
      lines(y$lons,y$y[,iii],lwd=2,lty=2,col="grey80")
      #print(rbind(y$lons,y$y[,iii]))
      if (iii < 10) ciii <- paste("00",iii,sep="") else
      if (iii < 100) ciii <- paste("0",iii,sep="") else ciii <- iii
      filname <- paste(Currentname,"_validation-",ciii,".eps",sep="")
      if (Sys.info()[1]=="Linux") {
        if (ok.convert==0) {
          dev.copy2eps(file=filname)
          system(paste("convert ",filname," ",substr(filname,1,nchar(filname)-3),"gif",sep=""))
          system(paste("rm -f",filname))
        }
      }
      iii <- iii+1
    }
  }
  if (Sys.info()[1]=="Linux") {
    if (ok.gifmerge==0) {
      system(paste("gifmerge -l0 -30  ",Currentname,"_validation-*.gif > ",Currentname,".gif",sep=""))
      system(paste("rm -f ",Currentname,"_validation-*.gif",sep=""))
    }
  }
  par(oldpar)
}
