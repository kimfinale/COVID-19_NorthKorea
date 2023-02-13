plot_model_data <- function(model, data, var="symptomatic"){
  library(ggplot2)
  p <- NULL
  if (var == "symptomatic"){
   p <- ggplot(model, aes(x=date))+
      geom_ribbon(aes(ymax=`97.5%`,ymin=`2.5%`), fill="steelblue", alpha=0.3) +
      geom_ribbon(aes(ymax=`75%`,ymin=`25%`), fill="steelblue", alpha=0.6) +
      geom_line(aes(y=`50%`), color="steelblue") +
      geom_col(data=data, aes(x=date, y=`당일 발생자수`), fill="brown", alpha=0.5,
             inherit.aes = F) +
      theme_bw()
  }

  return(p)
}
