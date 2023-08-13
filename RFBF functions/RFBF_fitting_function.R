### FBF FITTING ALGORITHM 
#' @description An iterative algorithm to obtain FBF using the following 
#' parameters. The algorithm starts from the initial flow and update the flow 
#' to capture the maximum local variation during iterations. Finally, the 
#' algorithm returns a discrete flow until the stopping criterion is met. 
#' 
#' @param gamma_ini The initial flow started FBF fitting
#' @param manifoldata The given data set 
#' @param y0 The starting (boundary) point
#' @param y1 The ending (boundary) point
#' @param h The scale parameter to capture the local variation
#' @param rho The shrinkage constant used during iterations 
#' @param fixed_num_points Whether use fixed number of local points to determine
#' local variation (TRUE/FALSE)
#'      @example fixed_num_points=TRUE The parameter h is to set integer values 
#'      of the number of local points
#'      @example fixed_num_points=FALSE The parameter h is to set the radius of 
#'      local neighborhood
#' @param interpolation Whether implement interpolation during iterations 
#' (TRUE/FALSE)
#' 
#' 
#' @export gamma_proj The updated flow from y0 to y1


RFBF_fitting <- function(gamma_ini, manifoldata, y0,y1,h,rho,fixed_num_points = FALSE, interpolation=FALSE){
  ## settings 
  dimension <- nrow(manifoldata)
  
  num_point <- ncol(manifoldata)
  
  y0 <- matrix(y0,ncol=1)
  
  y1 <- matrix(y1,ncol=1)
  
  vec_direct = (y1-y0)/norm2(y1-y0)
  
  ## FBF algorithm iteration begins here
  
  finish <- FALSE
  
  count <- 0
  
  obj_val <- 0
  
  stop_cri <- eps
  
  while(finish==FALSE){
    
    if (count==0){
      gamma <- gamma_ini
      
    }else{
      h <- rho*h
      gamma <- gamma_proj
      
    }
    
    
    ###  vector field begins
    
    gamma_proj_tmp <- NULL
    
    vec_field_tmp <- NULL 
    
    for (i in 1:(ncol(gamma)-1)){
      
      if ((i%%2)==0){
        # find the point in the data set 
        center_point = proj_updated(gamma[,i],h,manifoldata)
        
        ## local points 
        if (fixed_num_points == TRUE){
          data_local <- local_fixed_pts(center_point,h,manifoldata)
        }else{
          data_local <- local_pts(center_point,h,manifoldata)
        }
        
        ## local mean
        local_mean <- matrix(apply(data_local,1,mean),ncol=1)
        
        cov_data_local <- cov_mat(data_local,local_mean)
        
        p_vec <- kth_vector_field(data_local, local_mean, vec_direct,1)
        
        p_vec <- matrix(p_vec,ncol=1)
        
        
        ## vector fields
        l_proj = (p_vec)%*%t(p_vec)%*%(gamma[,i-1]-local_mean)+local_mean
        r_proj = (p_vec)%*%t(p_vec)%*%(gamma[,i+1]-local_mean)+local_mean
        gamma_proj_tmp = cbind(gamma_proj_tmp,l_proj,r_proj)
        vec_field_tmp <- cbind(vec_field_tmp,p_vec)
      }
      
      
      
    } 
    
    ## update
    gamma_proj = y0
    
    
    for(i in 1:(ncol(gamma_proj_tmp)-1)){
      avg_point = apply(gamma_proj_tmp[,i:(i+1)],1,mean)
      gamma_proj = cbind(gamma_proj,avg_point)
    }
    
    gamma_proj = cbind(gamma_proj,y1)
    
    vec_field <- vec_field_tmp
    
    ## change of value
    obj_val_upd <- 0
    
    for (i in 1:(ncol(gamma)-1)){
      
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
  
  if(interpolation==TRUE){
    gamma_proj = FBF_interpolation(gamma_proj,dist_ini)
    
  }
  
  return(gamma_proj) 
}



