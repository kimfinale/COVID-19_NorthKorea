
# this function is created just not to throw any some annotate functions I used.
plot_annotate2 <- function(model1, model2, data){
  library(ggplot2)
  p <- ggplot(model1, aes(x=date)) +
    geom_ribbon(aes(ymax=`97.5%`,ymin=`2.5%`), fill="steelblue", alpha=0.3) +
    geom_ribbon(aes(ymax=`75%`,ymin=`25%`), fill="steelblue", alpha=0.6) +
    geom_line(aes(y=`50%`), color="steelblue") +
    geom_ribbon(data=model2, aes(x=date, ymax=`97.5%`,ymin=`2.5%`), fill="darkgreen", alpha=0.3, inherit.aes=F) +
    geom_ribbon(data=model2, aes(x=date, ymax=`75%`,ymin=`25%`), fill="darkgreen", alpha=0.6, inherit.aes=F) +
    geom_line(data=model2, aes(x=date, y=`50%`), color="darkgreen",
              inherit.aes = F) +
    geom_col(data=data, aes(x=newdate, y=symptomatic), fill="brown", alpha=0.5,
             inherit.aes = F) +
    theme_bw() +
    labs(x="", y="Case") +
    scale_x_date(date_breaks="1 week", date_labels="%Y-%m-%d",
                 limits=c(min(inc_summary$date), as.Date("2022-06-30"))) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=1.8e6, ymax=2e6, alpha=0.3, fill="darkgreen") +
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=1.85e6, ymax=1.95e6, alpha=0.6, fill="darkgreen") +
    annotate("segment", x=as.Date("2022-02-26"), xend=as.Date("2022-02-28"),
             y=1.9e6, yend=1.9e6, color="darkgreen") +
    annotate("text", x = as.Date("2022-03-01"), y = 1.9e6, size = 4,
             label = "Infection (model)", hjust=0, vjust = 0.5) +
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=1.55e6, ymax=1.75e6, alpha=0.3, fill="steelblue") +
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=1.6e6, ymax=1.7e6, alpha=0.6, fill="steelblue") +
    annotate("segment", x=as.Date("2022-02-26"), xend=as.Date("2022-02-28"),
             y=1.65e6, yend=1.65e6, color="steelblue") +
    annotate("text", x = as.Date("2022-03-01"), y = 1.65e6, size = 4,
             label = "Symptomatic (model)", hjust=0, vjust = 0.5)+
    annotate("rect", xmin=as.Date("2022-02-26"), xmax=as.Date("2022-02-28"),
             ymin=1.3e6, ymax=1.5e6, alpha=0.5, fill="brown") +
    annotate("text", x = as.Date("2022-03-01"), y = 1.4e6, size = 4,
             label = "Symptomatic (data)", hjust=0, vjust = 0.5)

  return(p)
}


