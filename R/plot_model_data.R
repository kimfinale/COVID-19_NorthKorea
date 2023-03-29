plot_model_data <- function(model, data, var="symp"){
  library(ggplot2)
  sb <- scales::alpha(c("steelblue"), alpha = c(0.2, 0.55, 0.9)) # symptomatic
  br <- scales::alpha(c("brown"), alpha = c(0.5)) # data
  gr <- scales::alpha(c("darkgreen"), alpha = c(0.2, 0.55, 0.9)) # infection

  p <- NULL
  if (var == "symp"){
    p <- ggplot(model, aes(x=date))+
      geom_ribbon(aes(ymax=`symp_97.5%`,ymin=`symp_2.5%`, fill="Model 95% CrI")) +
      geom_ribbon(aes(ymax=`symp_75%`,ymin=`symp_25%`, fill="Model 50% CrI")) +
      geom_line(aes(y=`symp_50%`, color="Model median"), linewidth=1) +
      geom_col(data=data, aes(x=date, y=symptomatic, fill="Data"),
               inherit.aes = F) +
      scale_fill_manual("", values=c("Model 95% CrI"= sb[1],
                                     "Model 50% CrI"= sb[2], "Data"=br))+
      scale_color_manual("", values=c("Model median"=sb[3]))+
      labs(x="", y="Case") +
      scale_x_date(date_breaks="2 weeks", date_labels="%Y-%m-%d",
                   limits=c(min(model$date), max(data$date)))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      theme(legend.position = c(0.2, 0.7))
  }
  else if (var == "inf"){
    p <- ggplot(model, aes(x=date))+
      geom_ribbon(aes(ymax=`symp_97.5%`,ymin=`symp_2.5%`, fill="Model symp 95% CrI")) +
      geom_ribbon(aes(ymax=`symp_75%`,ymin=`symp_25%`, fill="Model symp 50% CrI")) +
      geom_line(aes(y=`symp_50%`, color="Model symp median"), linewidth=1) +
      geom_ribbon(aes(ymax=`inf_97.5%`,ymin=`inf_2.5%`, fill="Model inf 95% CrI")) +
      geom_ribbon(aes(ymax=`inf_75%`,ymin=`inf_25%`, fill="Model inf 50% CrI")) +
      geom_line(aes(y=`inf_50%`, color="Model inf median"), linewidth=1) +
      geom_col(data=data, aes(x=date, y=symptomatic, fill="Data"),
               inherit.aes = F) +
      scale_fill_manual("", values=c("Model symp 95% CrI"=sb[1],
                                     "Model symp 50% CrI"=sb[2], "Data"=br,
                                     "Model inf 95% CrI"=gr[1],
                                     "Model inf 50% CrI"=gr[2]))+
      scale_color_manual("", values=c("Model symp median"=sb[3],
                                      "Model inf median"=gr[3]))+
      labs(x="", y="Case") +
      scale_x_date(date_breaks="2 weeks", date_labels="%Y-%m-%d",
                   limits=c(min(model$date), max(data$date)))+
      theme_bw()+
      theme(axis.text.x=element_text(angle=60, hjust=1))+
      theme(legend.position = c(0.2, 0.7))
  }

  return(p)
}
