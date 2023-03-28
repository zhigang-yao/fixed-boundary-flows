## Simulation Demo on Cone
rm(list=ls())
library(rgl)


source("./add_functions.R")
source("./RFBF_smoothing.R")
source("./RFBF_interpolation.R")
source("./RFBF_fitting_function.R")


type = "cone"

case = "s_long" #"s_short" #"c"   

sd = 0.015

showcase = paste(type,"_",case,sep="")

if(showcase=="cone_c"){
  pure_curve = data.matrix(read.csv("./data sets/cone_c.csv", header=FALSE))
  # set parameters 
  resolution <- 30
  h <- 0.2
}else if(showcase=="cone_s_short"){
  pure_curve = data.matrix(read.csv("./data sets/cone_s_short.csv", header=FALSE))
  # set parameters 
  resolution <- 30
  h <- 0.12
}else if(showcase=="cone_s_long"){
  pure_curve = data.matrix(read.csv("./data sets/cone_s_long.csv", header=FALSE))
  # set parameters 
  resolution <- 20
  h <- 0.1
}

cone = data.matrix(read.csv("./data sets/cone.csv", header=FALSE))



n = ncol(pure_curve)
## generate noisy data 
R=1
H=1


manifoldata <- matrix(0,nrow=3,ncol=n)

for (i in 1:n){
  theta_ns <- pi+atan(pure_curve[2,i]/pure_curve[1,i])+rnorm(1,0,sd)
  r_ns <- pure_curve[3,i]*R/H+rnorm(1,0,sd)
  manifoldata[,i] <- c(r_ns*cos(theta_ns),r_ns*sin(theta_ns),H*r_ns/R)
}


## set boundary points
y0 = pure_curve[,20]
y1 = pure_curve[,ncol(pure_curve)-20]

## initial flow
gamma_ini <- gamma_given(resolution,y0,y1,1)

## RFBF algorithm
sol <- RFBF_fitting(gamma_ini, manifoldata, y0,y1,h,0.95)


## smoothing
sol_smoothing = sol


## interpolation 
dist_min = min(sapply(2:ncol(sol_smoothing),function(ii){norm2(sol_smoothing[,ii]-sol_smoothing[,ii-1])}))
dist_ini = min(norm2(gamma_ini[,2]-gamma_ini[,3]),dist_min)

curve_fbf = FBF_interpolation(sol_smoothing,dist_ini,sphere=FALSE)


### plot the result
open3d()
points3d(cone[1,],cone[2,],cone[3,],color="yellow",alpha=0.2)
points3d(manifoldata[1,],manifoldata[2,],manifoldata[3,],alpha=0.9)
rgl.spheres(curve_fbf[1,],curve_fbf[2,],curve_fbf[3,], r = 0.005, color = "red")


