
### Find the norm of vector x
norm2 <- function(x) {
  return(sqrt(sum (x^2)))
}


### The xyangle function used to find the angle 
### between two vectors x and y

xyangle <- function(x,y) {
  normx <- norm2(x)
  normy <- norm2(y)
  normxy <- sum(x*y)
  
  xy <- normxy/(normx * normy)
  if (abs(xy) > 1) {xy <- round(xy)}
  return(acos(xy)) 
}


# find local points with fixed number 

local_fixed_pts <-function(p,range,manifoldata){
  
  dimension <- nrow(manifoldata)
  
  num_point <- ncol(manifoldata)
  
  local_dist <- rep(0,num_point)
  
  for(i in 1:num_point)
  {
    local_dist[i] <- dist_pt(p,manifoldata[,i])
  }
  local_dist_sorted <- sort(local_dist)
  
  dist_selected <- local_dist_sorted[1:range]
  
  local_pts_idx <- unlist(sapply(1:range, function(i){which(local_dist==dist_selected[i])}))
  
  local_pts_idx <- sort(local_pts_idx)
  
  local_manifoldata <- manifoldata[,local_pts_idx]
  
  return(local_manifoldata)
}

# find local points with fixed distance 

local_pts <-function(p,range,manifoldata){
  
  dimension <- nrow(manifoldata)
  num_point <-ncol(manifoldata)
  
  local_manifoldata <- matrix(NA,dimension,1)
  
  for(i in 1:num_point)
  {
    if(dist_pt(p,manifoldata[,i])<range)
      local_manifoldata <- cbind(local_manifoldata,manifoldata[,i])
  }
  
  local_manifoldata <- local_manifoldata[,2:ncol(local_manifoldata)]
  
  return(local_manifoldata)
}



### Given n points on the sphere and a point X. 
### Rescale n points to the sphere of radius |X|.
### Find the covariance matrix at point X

cov_mat <- function(conedata, X) {
  
  n <- ncol(conedata)
  
  tanvec <- matrix(0,nrow(conedata),1)
  
  avg_mat <- matrix(0,nrow(conedata), nrow(conedata))
  
  for (i in 1: n) {
    tanvec[,1] = X-conedata[,i]
    normtanvec <- norm2(tanvec)
    if(normtanvec > 0) {
      avg_mat = avg_mat + tanvec %*% t(tanvec)}
  }
  
  avg_mat <-  avg_mat/n
  return(avg_mat)
}

## vector field not using the eigenvalues
kth_vector_field <- function(manifoldata, X, vec_direct,k) {
  
  matX <- cov_mat(manifoldata,X)
  
  ktheigen <- kth_eigen(matX,k)
  
  if((ktheigen$kvec)%*%vec_direct>0) {
    
    return(ktheigen$kvec)
    # * ktheigen$kval
    
  } else {
    
    return(-ktheigen$kvec)
    #
  }
  
}


## Function that solves the BVP-DAE given a particular value of delta.

gamma_given<-function(resolution,y0,y1,method,dimension=3,...){
  if(method==1){
    # straight line not in the manifold
    xguess = seq(from=0, to=1, by=(1/resolution))
    gamma_given <- matrix(NA,dimension,resolution+1)
    for(i in 1:(1+resolution)){
      gamma_given[,i]=y0+(xguess[i]*(y1-y0))
    }
  }  
  else if(method==2){
    # curve lies in the manifold
    angle_diff = xyangle(y0,y1)
    str_dist = norm2(y1-y0)
    
    phi_guess = seq(from=0, to=angle_diff, by=(angle_diff/resolution))
    fake_xguess = c(1)
    for(i in phi_guess){
      d = sin(i)/sin((pi/2)+(angle_diff/2)-i)
      fake_xguess = c(fake_xguess, (d/str_dist))
    }
    fake_xguess = fake_xguess[-1]
    gamma_given <- matrix(NA,dimension,resolution+1)
    for(i in 1:(1+resolution)){
      pos_y = (y0+(fake_xguess[i]*(y1-y0)))
      gamma_given[,i] = normalise(pos_y)
    }
  }
  
  return(gamma_given)
}


# select a subset of the given data set
sub_data <- function(h_sub, y0,y1,resolution,manifoldata){
  
  dimension <- nrow(manifoldata)
  
  num_point <-ncol(manifoldata)
  
  y0 <- matrix(y0,ncol=1)
  
  y1 <- matrix(y1,ncol=1)
  
  gamma_ini <- gamma_given(resolution,y0,y1,1)
  
  # sub_manifold data 
  sub_manifold_tmp = NULL
  
  for (i in 1:(resolution+1)){
    data_tmp <- local_pts(gamma_ini[,i],h_sub,manifoldata)
    sub_manifold_tmp <- cbind(sub_manifold_tmp,data_tmp)
  }
  
  data1 <- unique(sub_manifold_tmp[1,])
  
  sub_manifold= NULL
  for (i in 1:length(data1)){
    idx = which(sub_manifold_tmp[1,]==data1[i])
    sub_manifold = cbind(sub_manifold,sub_manifold_tmp[,idx[1]])
  }
  
  return(list(data=sub_manifold,gamma=gamma_ini))
  
}


# map a given point back to the data set
proj_updated <- function(point,rr,data){
  num_data <- ncol(data)
  #step 1: find the nearest point
  dist <- sapply(1:num_data, function(i){norm2(data[,i]-point)})
  nst_pt = data[,which(dist==min(dist))[1]]
  #step 2: find local points around the nearest point 
  pts = local_pts(nst_pt,rr,manifoldata)
  pts_mean = apply(pts,1,mean)
  #step 3: find the projected point
  dist <- sapply(1:num_data, function(i){norm2(data[,i]-pts_mean)})
  proj_pt = data[,which(dist==min(dist))[1]]
  
  return(proj_pt)
}
