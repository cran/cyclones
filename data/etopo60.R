load("etopo60.Rdata")
etopo60 <- list(x=top$ETOPO60X,y=top$ETOPO60Y,z=t(top$ROSE))
rm(top); gc(reset=TRUE)


