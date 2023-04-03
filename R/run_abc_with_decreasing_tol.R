run_abc_with_decreasing_tol <- function(method,
                                        timelimit,
                                        init_tol,
                                        decrement,
                                        niter,
                                        seed,
                                        model,
                                        prior,
                                        nsim,
                                        target,
                                        alpha,
                                        nrep,
                                        save_res=TRUE){
  library(EasyABC)
  tol <- init_tol
  res_to_save <- NULL # store the results and save when it is the best (i.e., when the next res is NULL)
  tstamp <- tstamp(hour=T, minute=T, second=T)
  for(z in 1:niter){
    cat("seed =", seed, ", iter =", z, ", tolerance =", tol, "\n")
    set.seed(seed)
    res <-
      run_with_timelimit(timelimit,
                         ABC_sequential(method=method,
                                        model=model,
                                        prior=prior,
                                        nb_simul=nsim,
                                        summary_stat_target=target,
                                        alpha=alpha,
                                        tolerance_target=tol,
                                        verbose=F,
                                        M=nrep))

    tstamp <- tstamp()
    if (!is.null(res)) {
      if(z == niter & save_res) {
      saveRDS(res_to_save, paste0("outputs/abc_region_", PARAMETERS$region, "_erlang_",
                                  PARAMETERS$erlang, "_alpha_", alpha,
                                  "_tol_", tol, "_seed_", seed, "_", tstamp, ".rds"))
      }
    }
    else {
      if (!is.null(res_to_save) & save_res) {
        saveRDS(res_to_save, paste0("outputs/abc_region_", PARAMETERS$region, "_erlang_",
                                    PARAMETERS$erlang, "_alpha_", alpha,
                                    "_tol_", tol, "_seed_", seed, "_", tstamp, ".rds"))
      }
      break;
    }

    res_to_save <- res
    tol <- tol - decrement
  }
  return(list(fit=res_to_save, seed=seed, final_tolerance=tol,
              alpha=alpha, final_tstamp=tstamp))
}
