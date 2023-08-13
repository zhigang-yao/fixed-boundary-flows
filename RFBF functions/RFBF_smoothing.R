### LOCAL SMOOTHING ALGORITHM 
#' @description An algorithm to perform local smoothing to the flow by using the
#' given data set.  
#' 
#' @param curve_input The input flow 
#' @param h_smoothing The scale parameter to perform local smoothing
#'        @example h_smoothing to set integer values if fixed_num = TRUE
#'        @example h_smoothing to set radius of local neighborhood if 
#'        fixed_num = FALSE
#' @param wmean Whether weighted average is implemented (TRUE/FALSE)
#' @param fixed_num Whether use fixed number of local points to determine
#' local variation (TRUE/FALSE)
#' @param sphere Whether the testing manifold is unit spheres (TRUE/FALSE)
#' 
#' 
#' @export curve_output The flow after local smoothing

FBF_smoothing <- function(curve_input,h_smoothing,wmean=FALSE,fixed_num=FALSE,sphere=TRUE){
  
  finish = FALSE
  
  test = curve_output = curve_input
  
  while(finish==FALSE){
    
    for (i in 2:(ncol(curve_input)-1)){
      
      ## local points 
      if(fixed_num==FALSE){
        if (h_smoothing%%1==0){
          print("Enter a valid h!")
        }
        data_local <- local_pts(curve_output[,i],h_smoothing,manifoldata)
      }else{
        if (h_smoothing%%1!=0){
          print("Enter a valid h!")
        }
        data_local <- local_fixed_pts(curve_output[,i],h_smoothing,manifoldata)
      }
      
      ## local mean
      if(wmean==TRUE){
        data_local_dist = apply(data_local-matrix(rep(curve_output[,i],h_smoothing),
                                                  ncol=h_smoothing),2,norm2)
        w = data_local_dist/sum(data_local_dist)
        
        wlocal_mean <- apply(data_local*matrix(rep(w,each=3),ncol=h_smoothing),1,sum)
        
      }else{
        wlocal_mean <- apply(data_local,1,mean)
      }
      
      ## map back to the sphere
      if(sphere==TRUE){
        test[,i] = wlocal_mean/norm2(wlocal_mean)
      }else{
        test[,i] = wlocal_mean
      }
      
      
    }
    
    ## stopping criterion
    d = sum(apply(curve_output-test,2,norm2))
    
    if(d<=1e-4){
      finish=TRUE
    }
    curve_output = test
  }
  
  return(curve_output)
}




