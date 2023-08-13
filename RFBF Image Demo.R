#' ANALYSIS IMAGE DATA
#' 
#' @description The image data set contains 264 face images of 66 people, with four 
#' images of each person. Each image has been resized to 50*37 pixels and 
#' becomes a vector in R^{1850}. 
#' 
#' parameters used to obtain FBFs 
#' @param y0 The starting (boundary) point
#' @param y1 The ending (boundary) point
#' @param resolution The number of points on the initial flow connecting 
#' y0 and y1
#' @param h The scale parameter to capture the local variation
#' @param rho The shrinkage constant used during iterations 
#' @param eps The stopping criterion constant 
#' 
#' 
#' @export curve_fbf A RFBF determined by the given parameters and 
#' visualization is done in Python


rm(list=ls())

# import functions
source("./RFBF functions/add_functions.R")
source("./RFBF functions/RFBF_smoothing.R")
source("./RFBF functions/RFBF_interpolation.R")
source("./RFBF functions/RFBF_fitting_function.R")

## LOAD DATA SET
image_data <- data.matrix(read.csv("./data sets/image_faces_33by2by4.csv", 
                                   header=TRUE,sep=","))
dim(image_data)
# Samples: 264
# Features: 1850 each image dimension (50,37)

# scale pixels to [0,1]
manifoldata <- t(image_data)/255

# scale data to the unit sphere 
manifoldata <- apply(manifoldata,2,function(x) x/norm2(x))


# FBF algorithm begins here 

# set boundary points 
y0 <- manifoldata[,10]
y1 <- manifoldata[,89]

# set scale parameter (number of nearest data points to capture local variation)
h <-  10
# set parameters 
resolution <- 20
rho <-1
eps <- 1e-2

# set initial flow
gamma_ini <- gamma_given(resolution,y0,y1,2,dimension=nrow(manifoldata))

# Fitting FBF
sol <- RFBF_fitting(gamma_ini, manifoldata, y0,y1,h,rho,fixed_num_points=TRUE)

sol <- apply(sol,2,function(x)x/norm2(x))

# smoothing the flow 
sol_smoothing = sol


# interpolation for the flow 
gamma_ini_dis <- sapply(1:(ncol(gamma_ini)-1),function(i){norm2(gamma_ini[,i+1]-gamma_ini[,i])})
dist_ini = min(gamma_ini_dis)

curve_fbf = FBF_interpolation(sol_smoothing,dist_ini)

# get the dimension of the FBF obtained 
dim(curve_fbf)

# Remark: visualisation is implemented in Python