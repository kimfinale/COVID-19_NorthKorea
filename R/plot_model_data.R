plot_model_data <- function(model, data, var="symptomatic"){
  library(ggplot2)
  sb <- alpha(c("steelblue"), alpha = c(0.4, 0.7, 0.9)) # symptomatic
  br <- alpha(c("brown"), alpha = c(0.5)) # data
  gr <- alpha(c("darkgreen"), alpha = c(0.4, 0.7, 0.9)) # infection

  p <- NULL
  if (var == "symptomatic"){
    p <- ggplot(model, aes(x=date))+
      geom_ribbon(aes(ymax=`97.5%`,ymin=`2.5%`, fill="Model 95% CrI")) +
      geom_ribbon(aes(ymax=`75%`,ymin=`25%`, fill="Model 50% CrI")) +
      geom_line(aes(y=`50%`, color="Model median"), linewidth=1) +
      geom_col(data=data, aes(x=newdate, y=symptomatic, fill="Data"),
               inherit.aes = F) +
      scale_fill_manual("", values=c("Model 95% CrI"= sb[1],
                                     "Model 50% CrI"= sb[2], "Data"=br))+
      scale_color_manual("", values=c("Model median"=sb[3]))+
      labs(x="", y="Case") +
      scale_x_date(date_breaks="2 weeks", date_labels="%Y-%m-%d",
                   limits=c(min(model$date), max(data$newdate)))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      theme(legend.position = c(0.2, 0.7))
  }

  return(p)
}
