gradient.wind <- function(storms=NULL,icyclone=1) {
# Gradient wind: Fleagle and Businger (1980) p. 163. (eq 4.27):
  
  if (is.null(storms)) data(storms,envir =environment())
  f <- 2*0.0000729212*sin(pi*storms$lat[,icyclone]/180)
  vg <- storms$max.speed[,icyclone]    
  v.grad <- -0.5*f*pi*storms$radius[,icyclone]*1000*(1 - sqrt(1 + 4*vg/(f*storms$radius[,icyclone]*1000)))
  attribute(v.grad,"unit") <- "m/s"
  attribute(v.grad,"description") <- "gradient wind"
  storms$v.grad <- v.grad
  invisible(storms)
}
  
