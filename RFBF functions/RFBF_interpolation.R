######################################################################
## linear interpolation 

FBF_interpolation <- function(curve_input,dist_ini,sphere=TRUE){
  
  pt1 = curve_input[,1]
  
  ptn = curve_input[,ncol(curve_input)]
  
  curve_output = matrix(curve_input[,1],nrow=nrow(curve_input))
  
  i = 2
  
  while (i<=ncol(curve_input)){
    
    test_dist = norm2(curve_input[,i]-curve_output[,ncol(curve_output)])
    
    
    if(test_dist>dist_ini){
      ## calculate the number of points 
      npts = ceiling(test_dist/dist_ini)
      
      ## insert points 
      inner_pts = sapply(c(1:npts),function(x){curve_output[,ncol(curve_output)]+x/npts*(curve_input[,i]-curve_output[,ncol(curve_output)])})
      
      if (sphere==TRUE){
        inner_pts = apply(inner_pts,2,function(x){x/norm2(x)})
      }
      
      
      test = norm2(inner_pts[,ncol(inner_pts)]-curve_input[,i])
      
      if (test<1e-3){
        curve_output = cbind(curve_output,inner_pts)
      }else{
        curve_output = cbind(curve_output,inner_pts,curve_input[,i])
      }
      
      i = i+1
      
    }else{
      i = i+1
    }
    
  }
  
  if(norm2(curve_output[,ncol(curve_output)]-ptn)>=1e-3){
    curve_output = cbind(curve_output,ptn)
  }
  
  if(norm2(curve_output[,1]-pt1)>=1e-3){
    curve_output = cbind(curve_output,pt1)
  }
  
  return(curve_output)
  
}



