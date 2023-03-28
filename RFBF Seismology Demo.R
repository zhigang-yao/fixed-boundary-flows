rm(list=ls())
### loading data
### lying on the unit sphere
library(rgl)


source("./add_functions.R")
source("./RFBF_smoothing.R")
source("./RFBF_interpolation.R")
source("./RFBF_fitting_function.R")

### load data set 
manifoldata_origin <- data.matrix(read.csv("./data sets/seismology cartesian.csv", header=TRUE,sep=","))
manifoldata_origin <- t(manifoldata_origin[,-1])
dim(manifoldata_origin)

### endpoints
set_idx = 3 # 1 2 3
h_idx = 2 # 1 2 3

if (set_idx==1){
  y0 <- c(-0.4825831,  0.7905980,  0.3769195)
  y1 <- c(-0.7494472,  0.5485727,  0.3706709)
}else if(set_idx==2){
  y0 <- c(-0.5,  0.8,  0.4)
  y0 <- y0/norm2(y0)
  y1 <- c(-0.66,  0.5,  0.4)
  y1 <- y1/norm2(y1)
}else if(set_idx==3){
  y0 <- c(-0.7494472,  0.5485727,  0.3706709)
  y1 <- c(-0.4914957,  0.7690331,  0.4086809)
}

if (h_idx==1){
  h <- 0.075 #300*1.6/6371
}else if(h_idx==2){
  h <- 0.15 #600*1.6/6371
}else if(h_idx==3){
  h <-  0.226
}

# select a subset of the given data to fit RFBF
manifoldata = sub_data(0.23, y0,y1,150,manifoldata_origin)$data


points3d(manifoldata[1,],manifoldata[2,],manifoldata[3,],color="green")
points3d(manifoldata_origin[1,],manifoldata_origin[2,],manifoldata_origin[3,])

rgl.spheres(y0[1],y0[2],y0[3],r = 0.008, color="red")
rgl.spheres(y1[1],y1[2],y1[3],r = 0.008, color="red")


# RFBF algorithm

## initial curve
gamma_ini <- gamma_given(resolution=30,y0,y1,2)

sol <- RFBF_fitting(gamma_ini, manifoldata, y0,y1,h,1)

sol <- apply(sol,2,function(x)x/norm2(x))



## smoothing
sol_smoothing = sol



## interpolation
dist_ini = norm2(gamma_ini[,2]-gamma_ini[,3])

curve_fbf = FBF_interpolation(sol_smoothing,dist_ini)


## plot the final results
rgl.open()

spheres3d(c(0,0,0),radius = 1,color="yellow",alpha=1)

points3d(manifoldata_origin[1,],manifoldata_origin[2,],manifoldata_origin[3,],color="black",alpha=0.9)

rgl.spheres(curve_fbf[1,],curve_fbf[2,],curve_fbf[3,], r = 0.005, color = "red")

