plot_annotate <- function(model, data){
  library(ggplot2)
  p <- ggplot(model, aes(x=date)) +
    geom_ribbon(aes(ymax=`97.5%`,ymin=`2.5%`), fill="steelblue", alpha=0.3) +
    geom_ribbon(aes(ymax=`75%`,ymin=`25%`), fill="steelblue", alpha=0.6) +
    geom_line(aes(y=`50%`), color="steelblue") +
    geom_col(data=data, aes(x=newdate, y=symptomatic), fill="brown", alpha=0.5,
             inherit.aes = F) +
    theme_bw() +
    labs(x="", y="Case") +
    scale_x_date(date_breaks="2 weeks", date_labels="%Y-%m-%d",
                 limits=c(min(model$date), as.Date("2022-06-30"))) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=0.55e6, ymax=0.75e6, alpha=0.3, fill="steelblue") +
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=0.6e6, ymax=0.7e6, alpha=0.6, fill="steelblue") +
    annotate("segment", x=as.Date("2022-02-26"), xend=as.Date("2022-02-28"),
             y=0.65e6, yend=0.65e6, color="steelblue") +
    annotate("text", x = as.Date("2022-03-01"), y = 0.65e6, size = 4,
             label = "Symptomatic (model)", hjust=0, vjust = 0.5)+
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
               ymin=0.3e6, ymax=0.5e6, alpha=0.5, fill="brown") +
    annotate("text", x = as.Date("2022-03-01"), y = 0.4e6, size = 4,
             label = "Symptomatic (data)", hjust=0, vjust = 0.5)
  return(p)
}
