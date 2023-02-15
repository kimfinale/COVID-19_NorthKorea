run_with_timelimit <- function(timelimit, expression) {
  setTimeLimit(timelimit)
  res <- tryCatch({expression}, error=function(e) NULL)
  return(res)
}


