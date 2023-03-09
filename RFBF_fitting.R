## RFBF algorithm 

dimension <- nrow(manifoldata)

num_point <-ncol(manifoldata)

y0 <- matrix(y0,ncol=1)

y1 <- matrix(y1,ncol=1)

vec_direct = (y1-y0)/norm2(y1-y0)

## initial curve
gamma_ini <- gamma_given(resolution,y0,y1,1)


#points3d(gamma_ini[1,],gamma_ini[2,],gamma_ini[3,],color="purple")
#points3d(manifoldata[1,],manifoldata[2,],manifoldata[3,])

#################################################################
## FBF algorithm begins 


finish <- FALSE

count <- 0

obj_val <- 0

stop_cri <- 1e-2

while(finish==FALSE){
  
  if (count==0){
    gamma <- gamma_ini
    
  }else{
    h <- roh*h
    gamma <- gamma_proj
    
  }
  
  
  ###  algorithm begins here
  
  gamma_proj_tmp <- NULL
  
  vec_field_tmp <- NULL 
  
  for (i in 1:resolution){
    
    if ((i%%2)==0){
      # find the point in the dataset 
      center_point = proj_updated(gamma[,i],h,manifoldata)
    
      ## local points 
      data_local <- local_pts(center_point,h,manifoldata)
      ## local mean
      local_mean <- matrix(apply(data_local,1,mean),ncol=1)
      
      cov_data_local <- cov_mat(data_local,local_mean)
      
      p_vec <- kth_vector_field(data_local, local_mean, vec_direct,1)
      
      p_vec <- matrix(p_vec,ncol=1)
      ## plot 
      #rgl.open()
      #points3d(center_point[1],center_point[2],center_point[3],color="red")
      #points3d(gamma_ini[1,],gamma_ini[2,],gamma_ini[3,],color="purple")
      
      #points3d(local_mean[1],local_mean[2],local_mean[3],color="blue")
      #points3d(data_local[1,],data_local[2,],data_local[3,],color="green")
      #points3d(manifoldata[1,],manifoldata[2,],manifoldata[3,])
      #arrow3d(p0 = c(local_mean), p1 =c(center_point+0.02*p_vec),width = 0.1, thickness = 0.005,col="blue")
      
      ## projections
      l_proj = (p_vec)%*%t(p_vec)%*%(gamma[,i-1]-local_mean)+local_mean
      r_proj = (p_vec)%*%t(p_vec)%*%(gamma[,i+1]-local_mean)+local_mean
      gamma_proj_tmp = cbind(gamma_proj_tmp,l_proj,r_proj)
      vec_field_tmp <- cbind(vec_field_tmp,p_vec)
    }
    
    
  } 
  
  gamma_proj = y0
  
  for(i in 1:(ncol(gamma_proj_tmp)-1)){
    avg_point = apply(gamma_proj_tmp[,i:(i+1)],1,mean)
    gamma_proj = cbind(gamma_proj,avg_point)
  }
  
  gamma_proj = cbind(gamma_proj,y1)
  
  
  #rgl.open() 
  #points3d(gamma[1,],gamma[2,],gamma[3,],color="purple")
  #points3d(gamma_proj[1,],gamma_proj[2,],gamma_proj[3,],color="red")
  
  
  
  
  vec_field <- vec_field_tmp
  
  
  obj_val_upd <- 0
  
  for (i in 1:resolution){
    
    a <- vec_field[,ceiling(i/2)]
    
    ab_angle <- xyangle(gamma_proj[,i],gamma_proj[,i+1])
    
    obj_val_upd <- obj_val_upd+ab_angle
  }
  
  iter_change <- abs((obj_val_upd-obj_val)/obj_val)
  if(abs(iter_change)<=stop_cri) finish <- TRUE
  print(iter_change)
  ## for next itertion 
  obj_val<- obj_val_upd
  count <- count +1
  
}




