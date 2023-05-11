# final epidemic size
# the simplest SIR model
# R = 1 - exp(-R0*R) where R is the final epidemic size (or R(\infty) for the SIR model)
final_epidemic_size <- function(R0 = 2) {
  y = function(x) x - 1 + exp(-R0*x)
  final_size <- uniroot(y, interval=c(1e-6,1-1e-6))$root

  return(final_size)

}

# # print parameter values
# print_params <- function(params){
#   n <- names(params)
#   for(i in seq_along(n)){
#     cat(paste0(n[i], "=", params[n[i]]), ", ")
#   }
# }

# print parameter values
print_params <- function(params){
  n <- names(params)
  str = paste0(n, "=", params[n])
  print(str)
}

my_discrete_colors <-
  c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
    "#FF7F00","black", "gold1", "skyblue2", "palegreen2", "#FDBF6F",
    "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
# tstamp
# create the time stamp yearmonthday by default and hour, minute, and second can be added
tstamp <- function(year=TRUE, month=TRUE, day=TRUE,
                   hour=FALSE, minute=FALSE, second=FALSE) {
  stamp1 <- c()
  stamp2 <- c()
  if (year & !month & !day) {
    stamp <- format(Sys.time(), "%Y")
  } else if (year & month & !day) {
    stamp1 <- format(Sys.time(), "%Y%m")
  } else if (year & month & day) {
    stamp1 <- format(Sys.time(), "%Y%m%d")
  } else if (!year & month & day) {
    stamp1 <- format(Sys.time(), "%m%d")
  } else if (year & !month & day) {
    stamp1 <- format(Sys.time(), "%Y%d")
  } else if (!year & month & !day) {
    stamp1 <- format(Sys.time(), "%m")
  } else if (!year & !month & day) {
    stamp1 <- format(Sys.time(), "%d")
  } else{ stamp1 <- "You'd better select parameters well."}

  if (hour & !minute & !second) {
    stamp2 <- format(Sys.time(), "%H")
  } else if (hour & minute & !second) {
    stamp2 <- format(Sys.time(), "%H%M")
  } else if (hour & minute & second) {
    stamp2 <- format(Sys.time(), "%H%M%S")
  } else if (!hour & minute & !second) {
    stamp2 <- format(Sys.time(), "%M")
  } else if (!hour & !minute & second) {
    stamp2 <- format(Sys.time(), "%S")
  } else if (!hour & minute & second) {
    stamp2 <- format(Sys.time(), "%M%S")
  } else{}

  if (!is.null(stamp2)) {
    stamp1 <- paste0(stamp1, "T", stamp2)
  }
  return (stamp1)
}

