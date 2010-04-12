
windstat <- function(fname=NULL,ws0=20,topo="etopo60.Rdata",cyclone=TRUE,
                     x.rng=c(5,35),y.rng=c(55,72),cmp=FALSE,mon=c(12,1,2),ERA40=TRUE) {
require(clim.pact)
print(paste("Latest verion   ",topo))
env <- environment()
  cmon <- c("Jan","Feb","Mar","Apr","May","Jun",
            "Jul","Aug","Sep","Oct","Nov","Dec")

if (is.null(fname)) cmp <- FALSE
if (!is.null(fname)) {
  print("Using the given object")
  if (class(fname)=="list") {results <- fname; class(results) <- c("CCI.object"); print("list--->CCI.object")}
  else if (class(fname)[1]=="CCI.object") {results <- fname}
  else if (is.character(fname)) {
    load(fname)
    print(paste("loading",fname))
    if (exists("results")){
      if (class(results)=="list") class(results) <- "CCI.object" 
    }
 } else stop('Cannot determine the object')
  print(range(results$yy,na.rm=TRUE))
} else {
  if (ERA40) {
    data(Storms.ERA40,envir = env)
    results <- Storms.ERA40
  } else {
  data(Storms.NMC,envir = env)
  results <- Storms.NMC
  }
}

results <- discard.bad(results)
results <- gradient.wind(results)

subtit<-paste("region: ",x.rng[1],"E...",x.rng[2],"E / ",y.rng[1],"N...",y.rng[2],"N.",
              " Threshold= ",ws0,sep="")
print(subtit)
dsh0 <- instring("-",subtit); Es <- instring("E",subtit); Ns <- instring("N",subtit); dsh <- dsh0
if (dsh0[1] !=0) {
  for (i in 1:length(dsh0)) {
    ln <- nchar(subtit)
    pE <- min(Es[Es >= dsh[1]]); pN <- min(Ns[Ns >= dsh[1]])
    #print(c(ln,pE,pN))
    if (dsh[1] <= max(Es)) subtit <- paste(substr(subtit,1,dsh[1]-1),
                                           substr(subtit,dsh[1]+1,pE-1),"W",
                                           substr(subtit,pE+1,ln),sep="") else
    if (dsh[1] > max(Es))  subtit <- paste(substr(subtit,1,dsh[1]-1),
                                           substr(subtit,dsh[1]+1,pN-1),"S",
                                           substr(subtit,pE+1,ln),sep="") 
    print(subtit)
    dsh <- instring("-",subtit); Es <- instring("E",subtit); Ns <- instring("N",subtit)
  }
}

if (cmp) {
  #print("Compare two results")
  result2 <- results
  data(Storms.NMC,envir = env); results <- Storms.NMC
  nobs2 <- length(result2$yy)
  t.scal2 <- 365.25/nobs2
  print(range(results$yy))
  print(range(result2$yy))
# Weed out any double entries:
  date <- result2$yy*10000+result2$mm*100+result2$dd
  chk.ovrlp <- table(date)
  n.obs <- as.numeric(chk.ovrlp)
  t.ovrlp <- as.numeric(rownames(chk.ovrlp))
  date.ovlp <- date[is.element(date,t.ovrlp[n.obs > median(n.obs)])]
  t.scl2 <- 1/median(n.obs)
  print(paste("Double entries: where more than",median(n.obs),"observations"))
  #print(date.ovlp)
  for (i in 1:length(date.ovlp)) {
    ii <- is.element(result2$yy*10000+result2$mm*100+result2$dd,date.ovlp)
    result2$psl[ii][(median(n.obs)+1):sum(ii)] <- NA
  }
}

# Weed out any double entries:
date <- results$yy*10000+results$mm*100+results$dd
chk.ovrlp <- table(date)
n.obs <- as.numeric(chk.ovrlp)
t.ovrlp <- as.numeric(rownames(chk.ovrlp))
date.ovlp <- date[is.element(date,t.ovrlp[n.obs > median(n.obs)])]
t.scl <- 1/median(n.obs)
print(paste("Double entries: where more than",median(n.obs),"observations"))
#print(date.ovlp)
for (i in 1:length(date.ovlp)) {
  ii <- is.element(results$yy*10000+results$mm*100+results$dd,date.ovlp)
  results$psl[ii][(median(n.obs)+1):sum(ii)] <- NA
}
print(paste("Time resolution is",t.scl,"days"))
            
nobs <- length(results$yy)
t.scal <- 365.25/nobs
nlon <- min(results$lon,na.rm=TRUE); xlon <- max(results$lon,na.rm=TRUE);
nlat <- min(results$lat,na.rm=TRUE); xlat <- max(results$lat,na.rm=TRUE);
#results$psl[,3:10] <- NA; results$lon[,3:10] <- NA; results$lat[,3:10]

print(topo)
if (file.exists(topo)) load(topo) else data(etopo60,envir=environment())
lons <- etopo60[[1]]
lons[lons>180] <- lons[lons>180]-360
srtx <- order(lons)
lons <- lons[srtx]
mask <- etopo60[[3]][srtx,]
lats <- etopo60[[2]]

nx <- length(lons)
ny <- length(lats)
mask[mask<=0] <- 0
mask[mask>0] <- 1
mask[lons < x.rng[1]] <- 0
mask[lons > x.rng[2]] <- 0
mask[,lats< y.rng[1]] <- 0
mask[,lats> y.rng[2]] <- 0

ix <- (lons >= nlon) & (lons <= xlon)
iy <- (lats >= nlat) & (lats <= xlat)

lons <- lons[ix]
lats <- lats[iy]
mask <- mask[ix,iy]

nx <- 15; ny=15
clons <- seq(nlon,xlon,length=nx)
clats <- seq(nlat,xlat,length=ny)
dx <- 0.5*(clons[2]-clons[1])
dy <- 0.5*(clats[2]-clats[1])

print("Mapping cyclone counts")
max.wind <- matrix(rep(0,nx*ny),nx,ny)
if (cmp) max.wind2 <- matrix(rep(0,nx*ny),nx,ny)
for (j in 1:ny) {
  for (i in 1:nx) {
    inside <-  (results$lon >= clons[i]-dx) & (results$lon < clons[i]+dx) &
               (results$lat >= clats[j]-dy) & (results$lat < clats[j]+dy) &
               (results$max.wind > ws0)
    max.wind[i,j] <- max(results$v.grad[inside],na.rm=TRUE)
    if (cmp) {
      inside <- (result2$lon >= clons[i]-dx) & (result2$lon < clons[i]+dx) &
                (result2$lat >= clats[j]-dy) & (result2$lat < clats[j]+dy) &
                (result2$max.wind > ws0)
      max.wind2[i,j] <- max(result2$v.grad[inside],na.rm=TRUE)
    }
  }
}
my.col <- rgb(c(1,0.7),c(1,0.7),c(1,0.7))
print(clons); print(clats)
print(paste("Grid box size= ",2*dx,"x",2*dy))

if (cmp) {
  postscript(file = "windstat_1a.eps",onefile=FALSE,horizontal=FALSE)

  par(las=1,mfcol=c(1,2))
  levels=seq(round(floor(min(c(max.wind)))),round(ceiling(max(c(max.wind)))),length=11)
  image(lons,lats,mask,main="Cyclone count per year",lwd=2,
      xlab="Longitude (deg E)",ylab="Latitude (deg N)",
      sub=paste("Period:",min(results$yy,na.rm=TRUE),"-",
      max(results$yy,na.rm=TRUE)," ws0=",ws0),col=my.col,levels=c(0,1))
  polygon(c(x.rng[1],rep(x.rng[2],2),rep(x.rng[1],2)),
        c(rep(y.rng[1],2),rep(y.rng[2],2),y.rng[1]),border="grey70",lty=2)
  contour(clons,clats,max.wind,lwd=2,add=TRUE,levels=levels)
  grid()
  addland()
  levels=seq(round(floor(min(c(max.wind-max.wind2)))),round(ceiling(max(c(max.wind-max.wind2)))),length=11)
  image(lons,lats,mask,main="Cyclone count difference per year",lwd=2,
      xlab="Longitude (deg E)",ylab="Latitude (deg N)",
      sub=paste("Period:",min(results$yy,na.rm=TRUE),"-",
      max(results$yy,na.rm=TRUE)," ws0=",ws0),col=my.col,levels=c(0,1))
  polygon(c(x.rng[1],rep(x.rng[2],2),rep(x.rng[1],2)),
        c(rep(y.rng[1],2),rep(y.rng[2],2),y.rng[1]),border="grey70",lty=2)
  contour(clons,clats,max.wind-max.wind2,lwd=2,add=TRUE,levels=levels)
  grid()
  addland()
  while (dev.cur() > 1) dev.off()
}


postscript(file = "windstat.scratch.eps",onefile=FALSE,horizontal=FALSE)
years <- as.numeric(rownames(table(results$yy)))
if (cmp) years <- as.numeric(rownames(table(c(results$yy,result2$yy))))
nyrs <- length(years)
ncyymm <- rep(NA,12*nyrs)
ncmm <- rep(NA,12*nyrs); dim(ncmm) <- c(12,nyrs); scmm <- ncmm 
yymm <- sort(rep(years,12)) + (rep(1:12,nyrs)-0.5)/12
ncyymm2 <- ncyymm; ncmm2 <- ncmm; scmm2 <- scmm
print(t.scal)

for (iy in 1:nyrs) {
  for (im in 1:12) {
    immyy <- (results$mm==im) & (results$yy==years[iy]) 
    lon <- results$lon[immyy]
    lat <- results$lat[immyy]
    wind.speed <- results$max.wind[immyy]
    if (cyclone) inside <- (lon >= x.rng[1]) & (lon <= x.rng[2]) &
                             (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                             (wind.speed > ws0) else
                 inside <- (lon >= x.rng[1]) & (lon <= x.rng[2]) &
                             (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                             (wind.speed < ws0)
    imiy <- (iy-1)*12+im
    if (sum(immyy,na.rm=TRUE)>0) ncyymm[imiy] <- max(results$v.grad[inside],na.rm=TRUE)
    ncmm[im,iy]<- ncyymm[imiy]
    if (cmp) {
      immyy <- (result2$mm==im) & (result2$yy==years[iy]) 
      lon <- result2$lon[immyy]
      lat <- result2$lat[immyy]
      wind.speed <- result2$max.wind[immyy]
      if (cyclone) inside <- (lon >= x.rng[1]) & (lon <= x.rng[2]) &
                             (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                             (wind.speed > ws0) else
                   inside <- (lon >= x.rng[1]) & (lon <= x.rng[2]) &
                             (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                             (wind.speed < ws0)
      imiy <- (iy-1)*12+im
      if (sum(immyy,na.rm=TRUE)>0) ncyymm2[imiy] <- max(result2$v.grad[inside],na.rm=TRUE)
      ncmm2[im]<-ncyymm2[imiy]
    }
    if (sum(is.finite(ncyymm))>0) {
      print("2>"); print(length(yymm)); print(summary(yymm)); print(length(ncyymm)); print(summary(ncyymm))
      plot(yymm+0.03,as.vector(ncyymm),type="s",xlab="time",
           ylab="max gradient windspeed",lwd=2,
           main=paste("Cyclone windspeed",
           results$yy[1],"-",max(results$yy,na.rm=TRUE)))
      lines(yymm,as.vector(ncyymm),type="s",lwd=2)
      lines(yymm,as.vector(ncyymm2),type="s",lwd=1,col="grey40")
     }
  }
}


nt <- length(ncyymm); ind <- 1:nt
ift <- ind[!is.na(ncyymm)][1]-1; itf <- reverse(ind[!is.na(ncyymm)])[1]
ifix <- is.na(ncyymm) & c(rep(FALSE,ift),rep(TRUE,nt-ift)) & 
                        c(rep(TRUE,itf),rep(FALSE,nt-itf))
ncyymm[ifix] <- 0

if (cmp) {
ift <- ind[!is.na(ncyymm2)][1]-1; itf <- reverse(ind[!is.na(ncyymm2)])[1]
ifix <- is.na(ncyymm2) & c(rep(FALSE,ift),rep(TRUE,nt-ift)) & 
                         c(rep(TRUE,itf),rep(FALSE,nt-itf))
ncyymm2[ifix] <- 0
}
dev.off()

# Final figure

# Fix - REB 15.11.2005:
print(summary(ncyymm))
ncyymm[!is.finite(ncyymm)] <- 0
print("HERE")
print(summary(ncyymm))
print(range(ncyymm,na.rm=TRUE))
x11(); plot(ncyymm); lines(ma.filt(ncyymm,12))
print(ncyymm,12)
print(ma.filt(ncyymm,12))

postscript(file = "windstat2.eps",onefile=FALSE,horizontal=FALSE)
par(las=1)
plot(range(yymm,na.rm=TRUE),range(ma.filt(ncyymm,12),na.rm=TRUE)*c(0,1.75),type="n",
         xlab="time",ylab="Gradient wind speed (m/s)",sub=subtit,
         main=paste("Cyclone wind",min(years),"-",max(years)))
lines(yymm,ma.filt(ncyymm,12),lwd=4)
if (cmp) lines(yymm,ma.filt(ncyymm2,12),type="l",lwd=3,col="grey")
lines(yymm,ma.filt(ncyymm,12))
grid()

ncmm <- ncmm/nyrs
trend1 <- lm(ncyymm ~ yymm)
abline(trend1,lty=2,col="grey20",lwd=3)
p.trend <- summary(trend1)
if (cmp) {
  ncmm2 <- ncmm2/nyrs
  trend2 <- lm(ncyymm2 ~ yymm)
  abline(trend2,lty=2,col="grey60",lwd=3)
}
print(summary(lm(ncyymm ~ yymm)))
if (cmp) print(summary(lm(ncyymm2 ~ yymm)))
print("test 1")
print(dev.cur())
dev.off()
print("test 2")


for (iy in 1:nyrs) {
  for (im in 1:12) {
    immyy <- (results$mm==im) & (results$yy==years[iy]) 
    lon <- results$lon[immyy]
    lat <- results$lat[immyy]
    wind.speed <- results$max.wind[immyy]
    if (cyclone) icyclone <- t.scal*sum((lon >= x.rng[1]) & (lon <= x.rng[2]) &
                                        (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                                        (wind.speed > ws0),na.rm=TRUE ) else
                  icyclone <- t.scal*sum((lon >= x.rng[1]) & (lon <= x.rng[2]) &
                                         (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                                         (wind.speed < ws0),na.rm=TRUE )
    scmm[im]<-scmm[im]+(icyclone - ncmm[im])^2
  }
}
scmm <- sqrt(scmm/(nyrs-1) )

if (cmp) {
for (iy in 1:nyrs) {
  for (im in 1:12) {
    immyy <- (result2$mm==im) & (result2$yy==years[iy]) 
    lon <- result2$lon[immyy]
    lat <- result2$lat[immyy]
    wind.speed <- result2$max.wind[immyy]
    if (cyclone) icyclone <- t.scal2*sum((lon >= x.rng[1]) & (lon <= x.rng[2]) &
                                         (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                                         (wind.speed > ws0),na.rm=TRUE ) else
                  icyclone <- t.scal2*sum((lon >= x.rng[1]) & (lon <= x.rng[2]) &
                         (lat >= y.rng[1]) & (lat <= y.rng[2]) &
                         (wind.speed < ws0),na.rm=TRUE )
    scmm2[im]<-scmm2[im]+(icyclone - ncmm2[im])^2
  }
}
scmm2 <- sqrt(scmm2/(nyrs-1) )
}

print("test 3")
print(summary(ncmm+scmm))

postscript(file = "windstat3.eps",onefile=FALSE,horizontal=FALSE)
par(col.axis="white",las=1)
plot(c(1,24),range(ncmm+scmm,na.rm=TRUE)*c(0,1.75),main="seasonal cyclone variability",
     type="n",xlab="Month",ylab="Mean storm count/month",sub=subtit)
sdv1 <- c(rep(ncmm+scmm,2),reverse(rep(ncmm-scmm,2))); sdv1[sdv1 < 0] <- 0
polygon(c(1:24,seq(24,1,by=-1)),sdv1,col="grey70")
if (cmp) {
  sdv2 <- c(rep(ncmm2+scmm2,2),reverse(rep(ncmm2-scmm2,2))); sdv2[sdv2 < 0] <- 0
  polygon(c(1:24,seq(24,1,by=-1)),sdv2,col="grey90",density=16,lwd=2)
}
lines(1:24,rep(ncmm,2),lwd=6)
if (cmp) lines(1:24,rep(ncmm2,2),lwd=4,col="grey50")

par(col.axis="black",las=2)
axis(1,at=1:24,labels=rep(c("Jan","Feb","Mar","Apr","May","Jun",
                            "Jul","Aug","Sep","Oct","Nov","Dec"),2))
axis(2)
grid()
dev.off()

mm <- rep(1:12,nyrs); yy <- sort(rep(years,12))
if (!is.null(mon)) imon <- is.element(mm,mon) else imon <- rep(TRUE,length(mm))

postscript(file = "windstat4.eps",onefile=FALSE,horizontal=FALSE)
par(las=1)
x <- seq(0,max(ncyymm[imon],na.rm=TRUE)+10,by=1)
mu <- mean(ncyymm[imon],na.rm=TRUE)
#poisson <- exp(-mu+x*log(mu)-cumsum(x))
poisson <- dpois(x,lambda=mu)
if (cmp) brks <- seq(0,max(c(ncyymm[imon],ncyymm2[imon]),na.rm=TRUE)+10,by=5) else
         brks <- seq(0,max(ncyymm[imon],na.rm=TRUE)+10,by=5)
subtit <- paste(subtit," (",cmon[mon[1]],"-",cmon[mon[length(mon)]],")",sep="")
if (cmp) h2 <- hist(ncyymm2[imon],breaks=brks,freq=FALSE)
h1 <- hist(ncyymm[imon],breaks=brks,freq=FALSE)
plot(h1$mids,h1$density,lwd=3,main="cyclone semiday count distribution",
     ylim=c(0,max(poisson)),xlab="Mean storm count/month",type="l",sub=subtit)
if (cmp) lines(h2$mids,h2$density,lwd=3,col="grey60")
grid()
lines(x,poisson,lwd=1,lty=2,col="grey30")
dev.off()


stat <- summary(trend1)
p.value1 <- round(100*(1-pf(stat$fstatistic[1],
                            stat$fstatistic[2],
                            stat$fstatistic[3])),1)

cyclones <- station.obj(x=ncyymm,yy=yy,mm=mm,
                        obs.name=paste("Cyclone (ws0= ",ws0,") count",sep=""),
                        location=paste(x.rng[1],"E - ",x.rng[2],"E / ",y.rng[1],"N - ",y.rng[2],sep=""),
                        unit="cyclone-days",station=NA,lat=mean(y.rng),lon=mean(x.rng),alt=NA,
                        wmo.no=NA,
                        ref=paste("Linear trend:",as.numeric(stat$coefficients[2])),
                                   "; p-value=",p.value1,"%")
if (stat$coefficients[2] > 0) cyclones$p.value.1 <- p.value1 else cyclones$p.value.1 <- -p.value1
if (cmp) {
  cyclones$ncyymm2 <- ncyymm2
  stat <- summary(trend2)
  p.value2 <- round(100*(1-pf(stat$fstatistic[1],
                              stat$fstatistic[2],
                              stat$fstatistic[3])),1)

  if (stat$coefficients[2] > 0) cyclones$p.value.2 <- p.value2 else cyclones$p.value.2 <- -p.value2
}


save(file="windstat.Rdata",cyclones)
invisible(cyclones)
}

