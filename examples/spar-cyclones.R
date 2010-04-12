library(cyclones)
source("cyclones/R/CCI.R")
source("clim.pact/R/distAB.R")
today <- datestr2num(format(Sys.time(), "%d-%b-%Y"))
this.year <- today[1]

path <- "/disk1"
list <-  list.files(path=path,pattern="b0mq-daily-slp-",full.names=TRUE)
list <- list[grep(".nc",list)]
n <-  length(list)

for (i in 1:n) {
  current.file <- list[i]

  result.name <- paste("data",substr(current.file,nchar(path)+1,nchar(current.file)-2),"Rdata",sep="")
  print(result.name)
                       
  CCI(fielddata=current.file,vname="msl",maxhar=16,fname=result.name,tslice=366,EPS=FALSE,
    plot.interval=10,graph.dir="CCI.SPAR",force365.25=TRUE,x.rng=c(-180,180),y.rng=c(0,90))
}
