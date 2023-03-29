summarize_model_output <- function(model_output,
                                   probs=c(0.025, 0.25, 0.5, 0.75, 0.975),
                                   end_date=NULL){
  sim1 <- as.data.frame(
    t(apply(model_output[[2]], 1, quantile, probs=probs, na.rm=T)))
  names(sim1) <- paste0("symp_", names(sim1))
  sim2 <- as.data.frame(
    t(apply(model_output[[1]], 1, quantile, probs=probs, na.rm=T)))
  names(sim2) <- paste0("inf_", names(sim2))
  sim <- cbind(sim1, sim2)

  if (is.null(end_date)) {
    dat <- readRDS("data/covid_overall_20230122.rds")
    end_date <- dat$date[1] + PARAMETERS$obslength - 1
  }
  start_date <- end_date - nrow(sim) + 1
  sim$date <- rev(seq(end_date, start_date, by=-1))

  return(sim)
}
