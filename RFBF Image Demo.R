## ANALYSIS IMAGE ON SPHERE
rm(list=ls())

source("./add_functions.R")
source("./RFBF_smoothing.R")
source("./RFBF_interpolation.R")
source("./RFBF_fitting_function.R")

## LOAD DATA SET

image_data <- data.matrix(read.csv("./data sets/image_faces_33by2by4.csv", header=TRUE,sep=","))
dim(image_data)
# Samples: 264
# Features: 1850 each image dimension (50,37)


manifoldata <- t(image_data)/255

# scale data to the unit sphere 
manifoldata <- apply(manifoldata,2,function(x) x/norm2(x))


# FBF algorithm

y0 <- manifoldata[,10]
y1 <- manifoldata[,89]

gamma_ini <- gamma_given(resolution=20,y0,y1,2,dimension=nrow(manifoldata))

h <-  10

sol <- RFBF_fitting(gamma_ini, manifoldata, y0,y1,h,1,fixed_num_points=TRUE)

sol <- apply(sol,2,function(x)x/norm2(x))



## smoothing
sol_smoothing = sol


## interpolation
gamma_ini_dis <- sapply(1:(ncol(gamma_ini)-1),function(i){norm2(gamma_ini[,i+1]-gamma_ini[,i])})
dist_ini = min(gamma_ini_dis)

curve_fbf = FBF_interpolation(sol_smoothing,dist_ini)

dim(curve_fbf)

