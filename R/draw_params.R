draw_params <- function(prior){
  prior <- list(c("unif",5,150), c("unif",1,20), c("unif",-5,10),
                c("unif",0.1,2))
  params <- rep(NA, length(prior))

  for(i in 1:length(params)){
    params[i] <-
      eval(parse(text=paste0("r", prior[[i]][1], "(1,",
                             prior[[i]][2], ",", prior[[i]][3],")")))
  }
  return(params)
}
