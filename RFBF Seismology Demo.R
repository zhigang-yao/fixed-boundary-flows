#' ANALYSIS SEISMOLOGICAL DATA ON THE UNIT SPHERE
#' 
#' @description The earthquake data contains epicenters of significant 
#' earthquakes (magnitude 5:5 and above) between 1904 and 2015.
#' 
#' parameters used to obtain FBFs 
#' @param set_idx The index corresponding to 3 sets of boundary points
#'       @example set_idx=1 Figure 9 a-c
#'       @example set_idx=2 Figure 9 d-f
#'       @example set_idx=3 Figure 9 g-i
#' @param h_idx he index corresponding to 3 values of scale parameters
#'       @example h_idx = 1 The scale parameter h = 0.075 (300miles)
#'       @example h_idx = 2 The scale parameter h = 0.15 (600miles)
#'       @example h_idx = 3 The scale parameter h = 0.225 (900miles)
#' @param resolution The number of points on the initial flow connecting 
#' y0 and y1
#' @param rho The shrinkage constant used during iterations 
#' @param eps The stopping criterion constant 
#' 
#' 
#' @export curve_fbf A RFBF determined by the given parameters and 
#' visualization in 3D 



rm(list=ls())
library(rgl)

# import functions
source("./RFBF functions/add_functions.R")
source("./RFBF functions/RFBF_smoothing.R")
source("./RFBF functions/RFBF_interpolation.R")
source("./RFBF functions/RFBF_fitting_function.R")

# load data set 
manifoldata_origin <- data.matrix(read.csv("./data sets/seismology cartesian.csv", 
                                           header=TRUE,sep=","))
manifoldata_origin <- t(manifoldata_origin[,-1])
dim(manifoldata_origin)

# determine boundary points 
set_idx = 3 # choose one from 1 2 3
# determine scale parameter 
h_idx = 2 # choose one from 1 2 3

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

# visualize the data set with the boundary points
points3d(manifoldata[1,],manifoldata[2,],manifoldata[3,],color="green")
points3d(manifoldata_origin[1,],manifoldata_origin[2,],manifoldata_origin[3,])

rgl.spheres(y0[1],y0[2],y0[3],r = 0.008, color="red")
rgl.spheres(y1[1],y1[2],y1[3],r = 0.008, color="red")


# RFBF algorithm begins here
# set parameters 
resolution <- 30
rho <- 1
eps <- 1e-2

# initial curve
gamma_ini <- gamma_given(resolution,y0,y1,2)

# Fitting FBF
sol <- RFBF_fitting(gamma_ini, manifoldata, y0,y1,h,rho)

sol <- apply(sol,2,function(x)x/norm2(x))


# smoothing the flow
sol_smoothing = sol


# interpolation for the flow
dist_ini = norm2(gamma_ini[,2]-gamma_ini[,3])

curve_fbf = FBF_interpolation(sol_smoothing,dist_ini)


# plot the final results in 3D
rgl.open()

spheres3d(c(0,0,0),radius = 1,color="yellow",alpha=1)

points3d(manifoldata_origin[1,],manifoldata_origin[2,],manifoldata_origin[3,],color="black",alpha=0.9)

rgl.spheres(curve_fbf[1,],curve_fbf[2,],curve_fbf[3,], r = 0.005, color = "red")

# Remark: visualization on the flat map is done in Python
