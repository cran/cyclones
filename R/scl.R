# R.E. Benestad, Heathrow 25.02.2005
# 10.7 cm flux: http://www.drao.nrc.ca/icarus/www/maver.txt

f10.7cm <- function(url="ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/Penticton_Adjusted/daily/DAILYPLT.ADJ",plot=TRUE,na.rm=TRUE) {
  # old: "ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/DAILYPLT.OBS"
  
  asciidata <- readLines(con=url)
  N <- length(asciidata); keep <- rep(TRUE,N)
  for (i in 1:N) if (nchar(asciidata[i]) <= 9) keep[i] <- FALSE
  asciidata <- asciidata[keep]
  writeLines(asciidata,"f10.7cm.txt")
  aa <- read.fwf(file="f10.7cm.txt",,widths=c(4,2,2,8),
                 col.names=c("Year","Month","Day","f10.7cm"))
  years <- as.numeric(row.names(table(aa$Year)))
  N <- length(years)
  aa.annual <- rep(NA,N)
  for (i in 1:N) {
    ii <- is.element(aa$Year,years[i])
    aa.annual[i] <- mean(aa$f10.7cm[ii],na.rm=na.rm)
  }
  if (plot) {
    plot(years,aa.annual,type="b",lwd=3,pch=19,main="Geomagnetic aa-index",col="grey",cex=0.6,
         sub=url,ylab="annual mean f10.7cm flux",xlab="Year")
  }
  attr(aa,"annual mean") <- rbind(years,aa.annual)
  invisible(aa)
}

sunspots <- function(url="http://sidc.oma.be/DATA/monthssn.dat",plot=TRUE) {
  require(cyclones)
  asciidata <- readLines(url)
  for (i in 1:6) asciidata[i] <- paste(asciidata[i],"NA")
  for (i in (length(asciidata)-20):length(asciidata)) {
    asterisk <- instring("*",asciidata[i])
    print(paste("i=",i,"astersk=",asterisk,asciidata[i]))
    if (length(asterisk)>0) {
      if (asterisk[1]>0) {
        if (length(asterisk)==1) asciidata[i] <- substr(asciidata[i],1,asterisk[1]-1) else
        if (length(asterisk)==2) asciidata[i] <- paste(substr(asciidata[i],1,asterisk[1]-1),
                                                 substr(asciidata[i],asterisk[1]+1,asterisk[2]-1))
      }
    } 
    if (i >= length(asciidata)-5) asciidata[i] <- paste(asciidata[i],"NA")
    print(asciidata[i])
  }
  writeLines(asciidata,"SIDC-sunspotnumber.txt")
  Rs <- read.table("SIDC-sunspotnumber.txt",header=FALSE,
                   col.names=c("yyyymm","year","sunspotnumber","smoothed"),
                   fileEncoding="latin1")
  attr(Rs,"description") <- list(src="monthly sunspot number from ROYAL OBSERVATORY OF BELGIUM",url=url)
  if (plot) {
    plot(Rs$year,Rs$sunspotnumber,type="b",lwd=3,pch=19,main="sunspots",col="grey",cex=0.6,
         sub=paste("ROYAL OBSERVATORY OF BELGIUM",url),ylab="monthly sunspot number",xlab="Year")
    polygon(c(2000,1857,1933,2010),c(0,209,209,0),col="grey90",border="grey90")
    polygon(c(2010,1933,1933,2010),c(0,209,258,0),col="grey80",border="grey80")
    lines(Rs$year,Rs$sunspotnumber,type="b",lwd=3,pch=19,col="grey40")
    lines(Rs$year,Rs$smoothed,lwd=2,col="red")
    grid()

    nh <- 90; N <- length(Rs$sunspotnumber)
    dt1.Rs <- dT(Rs$sunspotnumber,maxhar=nh)
    dt2.Rs <- dT(dt1.Rs$dy,maxhar=nh)
    lines(Rs$year,dt1.Rs$y.fit,lty=2,col="pink")

    imin <- (dt1.Rs$dy[2:N]*dt1.Rs$dy[1:(N-1)] < 0) &
            (dt2.Rs$dy[2:N]+dt2.Rs$dy[1:(N-1)] > 0) &
            (dt1.Rs$y.fit[2:N] < 30)
    minima <- round( (1:(N-1))[imin] * N/(N-1) )

    # check:
    check.scl <- diff(minima)
    print(check.scl)
    short <- (1:length(check.scl))[check.scl < 36]
    print(short); print(Rs$year[minima[short]])
    for (i in 1:length(short)) {
      ii <- minima[short[i]]; iii <- minima[short[i]+1]
      if ( (is.finite(ii)) & is.finite(iii) ) {
        print(c(ii,iii,Rs$sunspotnumber[ii],Rs$sunspotnumber[iii]))
        if (Rs$sunspotnumber[ii] > Rs$sunspotnumber[iii]) minima[short[i]] <- NA else
                                                          minima[short[i]+1] <- NA
      }
    }
    minima <- minima[is.finite(minima)]
    print("after qcheck:"); print(short); print(Rs$year[minima])
    
    points(Rs$year[minima],rep(0,length(minima)),pch=21,col="blue",cex=0.4)
    
    par(fig=c(0.45,0.70,0.72,0.89),new=TRUE,mar=rep(0.5,4),xaxt="n",yaxt="n",cex.axis=0.5)
    plot(Rs$year,Rs$sunspotnumber,xlim=c(2000,2010),pch=19,type="b",cex=0.4)
    polygon(c(1990,2020,2020,1990,1990),c(-10,-10,300,300,-10),col="grey95")
    points(Rs$year,Rs$sunspotnumber,pch=19,type="b",lty=1,cex=0.4)
    grid
    
    x11()
    plot(c(1,max(diff(Rs$year[minima]))),c(0,max(Rs$sunspotnumber,na.rm=TRUE)),type="n",
         main="solar cycle character",xlab="year from cycle start",ylab="monthly sunspot number",
         sub=paste("ROYAL OBSERVATORY OF BELGIUM",url))
    grid()
    for (i in 1:(length(minima)-1)) {
      col <- paste("grey",90-(3*i),sep="")
      if (i==(length(minima)-2)) col="red"
      lines(Rs$year[minima[i]:(minima[i+1]-1)] - Rs$year[minima[i]],
            Rs$sunspotnumber[minima[i]:(minima[i+1]-1)],pch=19,cex=0.6,col=col)
    }
    lines(Rs$year[minima[length(minima)]:length(Rs$sunspotnumber)] - Rs$year[minima[length(minima)]],
          Rs$sunspotnumber[minima[length(minima)]:length(Rs$sunspotnumber)],pch=19,cex=0.6,col="blue")
    legend(10,max(Rs$sunspotnumber,na.rm=TRUE),
           c("previous cycle   ","this cycle   "),lty=1,col=c("red","blue"),bg="grey95",cex=0.6)
  }
  invisible(Rs)
}



Benestad2005 <- function() {
# Figures in Benestad, R.E. (2005) A review of the solar cycle length estimates GRL 32 L15714, doi:10.1029/2005GL023621, August 13

stand <- function(x) {
 x <- (x - mean(c(x),na.rm=TRUE))/sd(c(x),na.rm=TRUE)
}
  
R <- sunspots(plot=FALSE)
R$year <- trunc(R$year)
R$month <- R$yyyymm - R$year*100
R$number <- R$sunspotnumber
#R <- read.table("~/data/indices/sunspot_num.dat",
#                col.names=c("year","month","number"))
yymm <- R$year+(R$month-0.5)/12

N <- length(R$number)
plot(c(1,1),c(1,1),type="n",ylab="Sunspot number",
     xlab="Time",
     main="sunspot number",
     ylim=c(0,300),xlim=range(yymm))
grid()
polygon(c(yymm,yymm[N],yymm[1]),c(R$number,0,0),
        col="grey80",border="grey70",density=40)

nharm <-c(30,50,80)
cols <- c("red","blue","darkgreen")
ltys <- c(1,2,3)
ih <- 0
scl.max <- rep(NA,50*length(nharm)); dim(scl.max) <- c(50,length(nharm))
scl.min <- scl.max; yymm.max <- scl.max; yymm.min <- scl.min

for (nh in nharm) {
  ih <- ih +1
  hfit <- dT(R$number,maxhar=nh)
  h2fit <- dT(hfit$dy,maxhar=nh)
  R.fit <- 0.5*(hfit$y.fit[2:N]+hfit$y.fit[1:(N-1)])
  YYMM <- 0.5*(yymm[2:N] + yymm[1:(N-1)])
  lines(YYMM,R.fit,lwd=1,col=cols[ih],lty=ltys[ih])

  imax <- (hfit$dy[2:N]*hfit$dy[1:(N-1)] < 0) & 
          (h2fit$dy[2:N]+h2fit$dy[1:(N-1)] < 0)
  imx <- (1:(N-1))
  skip.ends <- (imx > 12) & (imx < length(imax)-11)
  imax[!skip.ends] <- FALSE
  imx <- (1:(N-1))[imax]

  max.est <- 0.5*(hfit$y.fit[2:N][imax]+ hfit$y.fit[1:(N-1)][imax])
# Need to remove the end points causing spurious minima/maxima and small wiggles.
  keep <-  (R.fit[imax] > 35) & (R.fit[imx-11] < max.est) & (R.fit[imx+12] < max.est) 
  #print(c(sum(keep),sum(R.fit[imx-11]<max.est),sum(R.fit[imx+12]<max.est)))
  if (sum(keep)> 0) imax[!keep] <- FALSE
  nmax <- sum(imax)

  imin <- (hfit$dy[2:N]*hfit$dy[1:(N-1)] < 0) & 
          (h2fit$dy[2:N]+h2fit$dy[1:(N-1)] > 0) &
          (hfit$y.fit[2:N] < 25)
  imn <- (1:(N-1))
  skip.ends <- ((imn > 12) & (imn < length(imin)-11))
  imin[!skip.ends] <- FALSE
  imn <- (1:(N-1))[imin]

#print(c(length(imn),sum(imin)))
  min.est <- 0.5*(hfit$y.fit[2:N][imin]+ hfit$y.fit[1:(N-1)][imin])
# Need to remove the end points causing spurious minima/maxima and small wiggles.
  keep <- (R.fit[imin] < 20) & (R.fit[imn-5] > min.est) & (R.fit[imn+6] > min.est)     
  if (sum(keep)> 0) {
    #print("Remove")
    #print(rbind(YYMM[!keep],min.est[!keep]))
    if (ih >1) imin[!keep] <- FALSE
    print(c(sum(keep),sum(R.fit[imin] < 15),
            sum(R.fit[imn-12] > min.est),sum(R.fit[imn+12] > min.est)))
  }
  nmin <- sum(imin)

  #print(c(nmax,nmin))
  yymm.max[1:(nmax-1),ih] <- 0.5*(YYMM[imax][1:(nmax-1)]+YYMM[imax][2:nmax])
  scl.max[1:(nmax-1),ih] <- diff(YYMM[imax])
  yymm.min[1:(nmin-1),ih] <- 0.5*(YYMM[imin][1:(nmin-1)]+YYMM[imin][2:nmin])
  scl.min[1:(nmin-1),ih] <- diff(YYMM[imin])
  points(YYMM[imax],R.fit[imax],col=cols[ih],pch=2)
  points(YYMM[imin],R.fit[imin],col=cols[ih],pch=6)
  #print(YYMM[imin])
}
legend(1750,300,c(paste("N=",nharm[1]),paste("N=",nharm[2]),paste("N=",nharm[3])),
       lty=ltys,lwd=1,col=cols,bg="grey95",cex=0.8)
#dev.copy2eps(file="paper30_fig1.eps")


newFig()
plot(c(1,1),c(1,1),type="n",ylab="SCL (years)",
     xlab="Time",
     main="Solar Cycle Length",
     ylim=c(0,20),xlim=range(yymm))
grid()
points(yymm.max,scl.max,pch=19,col="black",cex=1.5)
points(yymm.min,scl.min,pch=19,col="grey35",cex=1.5)
points(yymm.max,scl.max,pch="+",col="grey95",cex=0.9)
points(yymm.min,scl.min,pch="-",col="black",cex=0.9)

# More conventional way...
# Low-pass filter the sunspot curve in order to find max and min
print("Finding SCL the conventional way")
Rz <- R$number; yy.R <- R$year; mm.R <- R$month
lpRz<-filter(Rz,rep(1,27)/27)
dRz1<-filter(diff(lpRz,differences=1),rep(1,27)/27)
dRz2<-filter(diff(lpRz,differences=2),rep(1,27)/27)
ntz<-length(dRz1)

mx<-dRz1[2:ntz]*dRz1[1:(ntz-1)] <= 0 & dRz2 < 0
mn<-dRz1[2:ntz]*dRz1[1:(ntz-1)] <= 0 & dRz2 > 0
ix<-seq(1,length(mx),by=1)
mx<-ix[mx] 
mn<-ix[mn]

mx<-c(mx[!is.na(mx)],length(Rz))
mn<-mn[!is.na(mn)]
ixx<-seq(1,length(Rz),by=1)
      
for (i in seq(1,length(mx),by=1)) {
  ix1<-max(c(1,mx[i]-12))
  ix2<-min(c(mx[i]+12,length(Rz)))
  ix.mask <- ixx >= ix1 & ixx <= ix2
  Rx<-max( Rz[ix1:ix2],na.rm=TRUE )
  mx[i]<- ixx[is.element(Rz, Rx) & ix.mask]
}

for (i in seq(1,length(mn),by=1)) {
  ix1<-max(c(1,mn[i]-12))
  ix2<-min(c(mn[i]+12,length(Rz)))  
  ix.mask <- ixx >= ix1 & ixx <= ix2
  Rn<-min( Rz[ix1:ix2], na.rm=TRUE )
  mn[i]<- ixx[is.element(Rz, Rn) & ix.mask]
}
print("Maxima:")
print(cbind(yy.R[mx],mm.R[mx],round(Rz[mx])))
print("Minima:")
print(cbind(yy.R[mn],mm.R[mn],round(Rz[mn])))


# Compute the Solar Cycle Length (SCL) from R:
ind<-seq(1,ntz,by=1)
scl.xx<-diff(ind[mx])
tim.xx<-yymm[mx]
yy.xx<-0.5*( R$year[mx[1:(length(mx)-1)]] + R$year[mx[2:length(mx)]] )
mm.xx<-0.5*( R$month[mx[1:(length(mx)-1)]] + R$month[mx[2:length(mx)]] )
scl.nn<-diff(ind[mn])
tim.nn<-yymm[mn]
yy.nn<-0.5*( R$year[mn[1:(length(mn)-1)]] + R$year[mn[2:length(mn)]] )
mm.nn<-0.5*( R$month[mn[1:(length(mn)-1)]] + R$month[mn[2:length(mn)]] )
scl<-c(scl.xx,scl.nn)/12
tim.scl<-c(tim.xx,tim.nn)
yy.scl<-c(yy.xx,yy.nn)
mm.scl<-c(mm.xx,mm.nn)
srt<-order(tim.scl)
scl<-scl[srt]
yy.scl<-yy.scl[srt]
mm.scl<-mm.scl[srt]
yymm.scl <- yy.scl + (mm.scl - 0.5)/12
tim.scl<-tim.scl[srt]

#scl.wl<-read.table("/home/rasmus/data/scl_wavelet.dat",as.is=T)
#
#scl.Rz<-data.frame(scl=scl,yy=yy.scl,mm=mm.scl,
#                   yymm=yy.scl + (mm.scl-0.5)/12,
#                   tim=tim.scl,
#                   Rz.max=Rz[mx],yy.max=yy.R[mx],mm.max=mm.R[mx],
#                   Rz.min=Rz[mn],yy.min=yy.R[mn],mm.min=mm.R[mn],
#                   i.max=mx,i.min=mn)

scl[yy.scl > 2000]<-NA
points(yy.scl + (mm.scl-0.5)/12,scl,col="black",pch=5)
legend(1900,3,c("RFC(max-max)    ","RFC(min-min)     ","Rz      "),
       pch=c(20,20,5),col=c("black","grey","black"),cex=0.9)
points(1907.5,2.20,pch="+",col="white",cex=0.8)
points(1907.5,1.40,pch="-",col="black",cex=0.8)
#dev.copy2eps(file="paper30_fig2.eps")


dev.new()
plot(c(1,1),c(1,1),type="n",ylab="SCL (year)",
     xlab="Date",
     main="Solar cycle length",
     ylim=c(0,20),xlim=range(yymm))
grid()
points(yymm.max,scl.max,pch=19,col="darkblue",cex=1.5)
points(yymm.min,scl.min,pch=19,col="steelblue",cex=1.5)
points(yymm.max,scl.max,pch="^",col="lightblue",cex=0.9)
points(yymm.min,scl.min,pch="v",col="black",cex=0.9)
points(yy.scl + (mm.scl-0.5)/12,scl,col="black",pch=5)
legend(1900,3,c("RFC(maks.-maks.)    ","RFC(min.-min.)     ","N      "),
       pch=c(20,20,5),col=c("darkblue","steelblue","black"),cex=0.9)
points(1907.5,2.20,pch="+",col="white",cex=0.8)
points(1907.5,1.40,pch="-",col="black",cex=0.8)
#dev.copy2eps(file="Cicerone_SCL-fig1.eps")

# SCL: visually reproduced from graphs in Thejll and Lassen DMI 99-9 using transparency and millimeter paper.

fk91 <- c(11.63, 11.53, 11.63, 11.53, 11.73, 11.43, 11.43, 11.13, 11.03, 10.83, 10.50, 10.40,
          10.13, 10.20, 10.10, 10.50, 10.40, 10.70, 10.63, 10.63, 10.63, 10.50, 10.30, 10.10)
tl99 <- c(11.52, 11.40, 11.22, 11.64, 11.16, 11.46, 11.52, 11.46, 11.34, 11.10, 10.50, 10.32,
          10.20,  9.84, 10.26,  9.96, 10.26, 10.44, 10.62, 10.74, 10.92, 10.62, 10.50, 10.32,
          10.20, 10.50)
yy.fk <- c(1865.6, 1872.7, 1876.9, 1884.0, 1889.6, 1895.2, 1900.9, 1907.9, 1912.2, 1917.8,
           1923.5, 1929.1, 1931.9, 1937.5, 1941.8, 1948.8, 1951.7, 1958.7, 1962.9, 1971.4,
           1974.2, 1981.3, 1982.7, 1985.5)
yy.tl <- c(1860.9, 1865.6, 1873.4, 1878.0, 1884.2, 1888.9, 1896.7, 1899.8, 1907.6, 1912.3,
           1918.5, 1923.2, 1929.4, 1934.1, 1940.3, 1943.4, 1949.6, 1954.3, 1960.5, 1965.2,
           1973.0, 1976.1, 1982.3, 1985.4, 1993.2, 1996.3)

# Filtered values: unfilter
#"FL91 12221","TL99 121"
# y <- F %*% x
# x <- solve(F) %*% y
n.fk91 <- length(fk91)
n.tl99 <- length(tl99)

F1 <- rep(0,n.fk91^2); dim(F1) <- c(n.fk91,n.fk91)
for (i in 1:(n.fk91-4)) F1[i,i:(i+4)] <- c(1,2,2,2,1)/8
i <- n.fk91 - 3;  F1[i,i:(i+3)] <- c(1,2,2,1)/6
i <- n.fk91 - 2;  F1[i,i:(i+2)] <- c(1,2,1)/4
i <- n.fk91 - 1;  F1[i,i:(i+1)] <- c(1,1)/2
i <- n.fk91;  F1[i,i] <- 1
scl.fk91 <- solve(F1) %*% fk91
scl.fk91[1:3] <- NA

F2 <- rep(0,n.tl99^2); dim(F2) <- c(n.tl99,n.tl99)
for (i in 1:(n.tl99-2)) F2[i,i:(i+2)] <- c(1,2,1)/4
i <- n.tl99 - 1;  F2[i,i:(i+1)] <- c(1,1)/2
i <- n.tl99;  F2[i,i] <- 1
scl.tl99 <- solve(F2) %*% tl99
scl.tl99[1] <- NA
#points(yy.fk,scl.fk91,col="red",pch=7)
#points(yy.tl,scl.tl99,col="blue",pch=8)
#lines(yy.fk,fk91,col="red",lty=2,lwd=2)
#lines(yy.tl,tl99,col="blue",lty=2,lwd=2)


# 10.7 cm radio flux.
# aa-index
#Lassen & Friis-Chirtensen


dev.new()
plot(c(1,1),c(1,1),type="n",ylab="Standardised",
     xlab="Time",
     main="SCL, GCR, 10.7cm flux & sunspot number",
     ylim=c(-2,5),xlim=c(1945,max(yymm)))
grid()
polygon(c(yymm,yymm[N],yymm[1]),stand(c(R$number,0,0)),
        col="grey70",border="grey70",density=30)

data(gcr,env=environment())
#gcr <- read.table("neutron2.dat",header=T, 
#                  col.names=c("Date","SecsOf1904","Climax","Huancayo.NoGC",
#                             "Huancayo.GC","Haleakala.IGY","Haleakala.S/M"))
climax <- stand(filter(as.numeric(gcr$Climax),rep(1,30)/30)) * 3 + 10
nt.gcr <- length(climax)
climax <- -climax[seq(15,nt.gcr,by=30)]
gcr.yy <- as.numeric(substr(as.character(gcr$Date),7,10))
gcr.mm <- as.numeric(substr(as.character(gcr$Date),1,2))
gcr.dd <- as.numeric(substr(as.character(gcr$Date),4,5))
jday.gcr  <- gcr.yy + (gcr.mm-1)/12 + (gcr.dd-0.5)/365.25
jday.gcr <- jday.gcr[seq(15,nt.gcr,by=30)]

F10.7cm <- f10.7cm(plot=FALSE)
#f10.7cm <- read.table("solar10-7cm.txt",header=TRUE,
#           col.names=c("Year","Month","Obsflux","Adjflux","Absflux"))
yymm.10.7 <- F10.7cm$Year + (F10.7cm$Month-0.5)/12

scl.min[length(scl.min)] <- NA
points(yymm.max,stand(scl.max),pch=19,col="black",cex=1.5)
points(yymm.min,stand(scl.min),pch=19,col="grey35",cex=1.5)
points(yymm.max,stand(scl.max),pch="+",col="grey95",cex=0.9)
points(yymm.min,stand(scl.min),pch="-",col="black",cex=0.9)

lines(jday.gcr,stand(climax),lty=1,lwd=2,col="steelblue")
lines(yymm.10.7, stand(gauss.filt(F10.7cm$f10.7cm,25)),lty=2,col="red")
legend(1990,5,c("Rz   ","GCR    ","F10.7cm   "),col=c("grey70","steelblue","red"),lty=c(1,1,2),cex=0.8)
#dev.copy2eps(file="paper30_fig3.eps")

#newFig()
#plot(yymm.10.7,f10.7cm$Obsplux)
#lines(yymm.10.7, f10.7cm$Adjflux)
#lines(yymm.10.7, f10.7cm$Absflux)

results <- list(yymm.max=yymm.max,scl.max=scl.max,
                yymm.min=yymm.min,scl.min=scl.min,
                yy.tl=yy.tl,scl.tl99=scl.tl99,
                yy.fk=yy.fk,scl.fk91=scl.fk91)
invisible(results)
}

