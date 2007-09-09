gradient.wind <- function(storms=NULL,icyclone=1) {
# Gradient wind: Fleagle and Businger (1980) p. 163. (eq 4.27):
  # print(class(storms))
  if (is.null(storms)) {data(Storms.ERA40,envir =environment()); storms <- Storms.ERA40; rm(Storms.ERA40)} else
                        if (class(storms)[1]=="CCI.object") {print("Use given object")} else
                        if (class(storms)[1]=="character") {
                         if (file.exists(storms)) load(storms) else
                         if (storms=="NMC") {data(storms,envir =environment()); storms <- Storms.NMC}
                       } else stop("Cannnot find the CCI data")
  f <- 2*0.0000729212*sin(pi*storms$lat[,icyclone]/180)
  vg <- storms$max.speed[,icyclone]    
  v.grad <- -0.5*f*pi*storms$radius[,icyclone]*1000*(1 - sqrt(1 + 4*vg/(f*storms$radius[,icyclone]*1000)))
#  print(summary(vg))
#  print(summary(v.grad))
  attr(v.grad,"unit") <- "m/s"
  attr(v.grad,"description") <- "gradient wind"
  storms$v.grad <- v.grad
  invisible(storms)
}
  
