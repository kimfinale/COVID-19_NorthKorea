find_min_DE = function(dat,
                       iter=30,
                       lower=c(1,  1, -10, 1e-3),
                       upper=c(100,20 ,20, 1),
                       control=DEoptim.control(NP=200,
                                               itermax=200,
                                               trace=FALSE)) {
  outlist = vector('list', iter)
  for (i in seq_len(iter)) {
    # cat("i =", i, "\n")
    set.seed(i)
    outlist[[i]] <- DEoptim(fn = negloglik,
                            lower = lower,
                            upper = upper,
                            control = control,
                            dat=dat)
  }
  vals = sapply(outlist, function(x) x$optim$bestval)
  idmin = which.min(vals)

  return(list(all=outlist, min=outlist[[idmin]]))
}
