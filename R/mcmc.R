mcmc <- function(init, kernel, iters=1000, thin=10, verbose=TRUE) {
  np <- length(init)
  ll <- -Inf
  mat <- matrix(0, nrow = iters, ncol = np)
  colnames(mat) <- names(init)
  x <- init
  ll <- logpost(p=x)
  if (verbose){ 
    message(paste(iters, "iterations"))
  }
  for (i in 1:iters) {
    # if (verbose) {
    #   message(paste(i, ""), appendLF = FALSE)
    # }
    for (j in 1:thin) {
      walk <- kernel(logpost, current=x, loglik=ll)
      x <- walk$current
      ll <- walk$loglik
    }
    if (verbose) {
      print(paste0("i =", i, ", ", unlist(lapply(x, function(x) paste0(x, "")))))
    }
    
    mat[i, ] = x
  }
  if (verbose) message("Done.")
  
  return(mat)
}
