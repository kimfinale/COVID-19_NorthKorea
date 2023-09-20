#' Incidence
#' Runs the model to extract difference of the state variables across
#' the unit_time in days
#'
#' @param pars # parameter values to be applied to the model
#' @param model # parameter values to be applied to the model
#' @param unit_time
#' @param state
#'
#' @return
#' @export
#'
#' @examples
#'
incidence <- function(pars, model=NULL, unit_time=1,
                      state="cumul_symptomatic"){

  params = update_params(pars=pars)

  if(!is.null(model)) {
    out <- model(params)
  }
  else {
    out <- params$model(params)
  }
  # extract by time_filter to match incidence across the
  time_filter <- seq(1, by=round(unit_time/params[["tau"]]),
                     length.out=(params[["ndays"]]))

  out <- out[time_filter, state, drop=FALSE]
  # if measured variables is more than one
  df <- data.frame(matrix(NA, nrow=nrow(out)-1, ncol=length(state)))
  names(df) <- state
  for(i in 1:length(state)){
    df[,i] <- diff(out[, state[i]])
  }

  return(df)
}
