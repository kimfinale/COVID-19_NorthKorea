# function to calculate the doubling time
# T indicates the time interval over which doubling time is calculated
# ie, daily doubling time if T = 1 and weekly rolling doubling time if T = 7
# 
doubling_time = function(T, count){
  n <- length(count) - T
  double_time <- rep(NA, length(count))
  for (i in 1: n){
    start <- i
    end <- i + T
    double_time[i + T] <-  
      (end - start)*log(2)/log(count[end]/count[start])
  }
  return(double_time)
}