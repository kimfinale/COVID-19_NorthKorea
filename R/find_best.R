find_best = function(id, dat, iter=30){
  # cat("b =", b, "\n")
  # OBS = dat[,id]
  outDE = find_min_DE(iter=iter, dat=dat[,id])
  return(c(outDE$min$optim$bestmem,outDE$min$optim$bestval))
}
