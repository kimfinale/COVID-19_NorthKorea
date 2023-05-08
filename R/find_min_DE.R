find_min_DE = function(dat,
                       iter=30,
                       fn,
                       lower = c(2,  1, -10, 1e-3),
                       upper = c(100,20 ,20, 1),
                       control = DEoptim.control(NP=200,
                                               itermax=200,
                                               trace=FALSE),
                       error_dist = error_dist) {


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
                dat = dat,
                error_dist = error_dist)
              return(out)}

  parallel::stopCluster(cl)

  vals = sapply(fits, function(x) x$optim$bestval)
  idmin = which.min(vals)

  return(list(fits=fits, min=fits[[idmin]]))
}
