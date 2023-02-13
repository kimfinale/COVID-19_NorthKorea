create_get_sum_stat_sepiar_erlang_province <- function(){
  pop_province <- data.table::fread("data/census_2008.csv")
  clean_province <- janitor::make_clean_names(pop_province$Province)

  for(i in 1:nrow(pop_province)){
    func_str <- sprintf('get_sum_stat_sepiar_erlang_%s <- function(pars) {
      params <- params_init(erlang=TRUE)
      params$obslength <- 70 #
      params$tau <- 0.2 #
      params$init$S <- %s
      p1 <- round(pars[1])
      params$ndays <- p1 + params$obslength + 1
      params$R0 <- pars[2]
      params$day_intervention <- p1 + pars[3] #
      params$R0_2 <- pars[4] #

      out <- sepiar_erlang_stoch(params)
      day_filter <- seq(1, by = round(1/params$tau), length.out = params$ndays)
      out <- out[day_filter, "CI"]
      daily_inc <- diff(out)

      res <- c(sum(daily_inc[1:p1]), daily_inc[(p1+1):(p1+params$obslength)])

      return(res)
    }', clean_province[i], pop_province$Total[i])

    eval(parse(text = func_str))
  }
}
