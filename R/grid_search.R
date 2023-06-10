grid_search <- function(parm_grids,
                        dat,
                        error_dist="pois") {

  library(parallel)
  library(doParallel)
  ncores <- detectCores()

  cl <- makeCluster(getOption("cl.cores", ncores-2))
  doParallel::registerDoParallel(cl)

  output <-
    foreach(i=1:nrow(parm_grids),
      .packages = c("COVID19NorthKorea","RcppDE"),
      .inorder=F) %dopar% {
        pars = as.double(parm_grids[i,])
        nll = negloglik(pars=pars, dat=dat, error_dist=error_dist)
        return(list(parm=pars, loglik=-nll))}

  parallel::stopCluster(cl)

  return(output)
}
