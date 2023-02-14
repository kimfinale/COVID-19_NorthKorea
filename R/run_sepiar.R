# # params2 indicate the paramters that were estimated
#   params <- PARAMETERS # global variables
#   params$measure_var <- outvar
#   run_model <- function(pars, model, outvar) {
#     max_intro <- max(round(pars[,1])) # value for the earliest introduction
#     dflist <- vector('list', length(outvar))
#     # to store symptomatic
#     for (i in 1:nrow(pars)) {
#       df <- data.frame(matrix(NA, nrow=(max_intro + params$obslength), ncol=length(outvar))
#       res <- daily_incidence(pars = pars)
#       for(j in 1;length(ourvar)) {
#         df[(eps+1):nrow(incmat), outvar[i]] <- res[, outvar[i]]
#       }
#
#   }
#
#   inc_summary <-
#     as.data.frame(t(apply(incmat, 1, quantile,
#                           probs=c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)))
#   totinc_summary <-
#     as.data.frame(t(apply(totincmat, 1, quantile,
#                           probs=c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)))
#
#   inc_summary$day <- 0:(nrow(inc_summary)-1)
#   inc_summary$date <-
#     seq(as.Date("2022-05-13") - maxd,
#         as.Date("2022-05-13") + params$obslength - 1, by="day")
#
#   totinc_summary$day <- inc_summary$day
#   totinc_summary$date <- inc_summary$date
#
#   return(list(infection=totinc_summary, symptomatic=inc_summary))
# }
#
