find_min_DE = function(iter=30,
                       obs,
                       fn,
                       lower = c(2,  1, -10, 1e-3),
                       upper = c(100,20 ,20, 1),
                       control = DEoptim.control(NP=200,
                                               itermax=200,
                                               trace=FALSE),
                       dist = "pois", ...) {


  library(parallel)
  library(doParallel)
  ncores <- detectCores()

  cl <- makeCluster(getOption("cl.cores", ncores-2))
  doParallel::registerDoParallel(cl)

  fits <-
    foreach(i=1:iter,
            .packages = c("COVID19NorthKorea","RcppDE"),
            .inorder=F) %dopar% {
              set.seed(i)
              out <- DEoptim(
                fn = fn,
                lower = lower,
                upper = upper,
                control = control,
                obs = obs,
                dist = dist, ...)
              return(out)}

  parallel::stopCluster(cl)

  vals = sapply(fits, function(x) x$optim$bestval)
  idmin = which.min(vals)

  return(list(fits=fits, min=fits[[idmin]]))
}
