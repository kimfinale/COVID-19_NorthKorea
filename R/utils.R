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


# these colors were derived from the stackoverflow (https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes)
mycolors25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

mycolors16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")

SB_alpha <- scales::alpha(c("steelblue"), alpha = c(0.4, 0.7, 0.9)) # symptomatic
Br_alpha <- scales::alpha(c("brown"), alpha = c(0.5)) # data
Gr_alpha <- scales::alpha(c("darkgreen"), alpha = c(0.4, 0.7, 0.9)) # infection




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

theme_map <- function(legend_position = c(0.2, 0.3)) {
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(), # bg of the panel
    legend.background = element_blank(), # get rid of legend bg
    legend.box.background = element_blank(),
    panel.spacing = unit( c(0,0,0,0), "null" ),
    plot.margin = unit( c(0,0,0,0), "null" ),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(legend_position[1],legend_position[2]))
}


##############################################################################
# tidy raster images
# adjust the resolution, crop and mask based on the reference shapefile or
# raster
tidy_raster <- function(rst, target_res_km = NULL, func = "mean",
                        extent = NULL, ref = NULL){

  deg_20km <- 0.1666667 # approx. 20 km in degrees
  target_res_deg <- round(target_res_km * deg_20km / 20, digits=7)
  rst_res_deg <- round(res(rst)[1], digits=7)

  if (!is.null(target_res_km)) { # aggregate or disaggregate
    if (target_res_deg > rst_res_deg) {
      rst <-
        raster::aggregate(rst, fact = round(target_res_deg / rst_res_deg),
                          fun = eval(parse(text=func)))
    }
    else if(target_res_deg < rst_res_deg){
      rst <-
        raster::disaggregate(rst, fact = round(rst_res_deg / target_res_deg))
    }
  }

  if (!is.null(ref)) {
    rst <- raster::projectRaster(rst, ref)
  }
  if(!is.null(extent)){
    rst <- raster::crop(rst, extent(extent), snap = "out")
    rst <- raster::mask(rst, mask = extent)
  }
  return(rst)
}

###############################################################################
## map_ratio
## calculates the ratio of the raster such that the ratio can be used to
## save the image
map_ratio <- function(rst){
  ext <- extent(rst)
  horiz <- ext[2] - ext[1]
  vert <- ext[4] - ext[3]
  ratio <- vert / horiz
  return(ratio)
}

#############################################################################
## IR_plot
## Holds default settings for plotting incidence rates
# IR_plot <- function(raster=NULL, shape=NULL, afss=NULL){
IR_plot <- function(raster=NULL, shape=NULL, shape2=NULL,
                    color_ramp="YlOrBr", rev=FALSE, level=9){
  # source("util/ggplot2_theme.R")
  library(raster)
  library(ggplot2)
  library(RColorBrewer)
  library(scico)
  if (is.null(shape)){
    shape <- readRDS("data/africa_sub_Sahara_adm0_shp.rds")
  }
  # if (is.null(shape2)){
  #   afss <- readRDS("data/africa_sub_Sahara_adm1_shp.rds")
  # }
  pal <- brewer.pal(level, color_ramp)
  if (rev) pal <- rev(pal)
  myPalette <- colorRampPalette(pal)
  # sc <- scale_colour_gradientn(colours = myPalette(10000), limits=c(1, 1e4))
  mypal <- myPalette(1e4)
  # mypal <- scico(30, palette = 'lajolla')

  rpts <- as.data.frame(raster::rasterToPoints(raster))
  colnames(rpts) <- c("lon","lat","ir")
  rpts$ir <- as.double(rpts$ir + 1) # 1 is added to use log
  maxir <- 1e4 # arbitrary maximum incidence rate
  rpts$ir[rpts$ir > maxir] <- maxir
  # # identify a way to plot the incidence rate on a log scale
  p <- ggplot(rpts) +
    geom_raster(aes(lon,lat,fill=ir)) +
    # scale_fill_gradientn(trans = "log",
    #                      limits = c(1, 10000),
    #                      breaks = c(10, 100, 1000, 10000),
    #                      colors = c("grey93", "steelblue", "gold", "darkred"),
    #                      "Incidence rate") +
    scale_fill_gradientn(trans = "log10",
                         limits = c(1, 1e4),
                         breaks = c(1, 1e1, 1e2, 1e3, 1e4),
                         colors = mypal,
                         "Incidence rate") +
    # geom_polygon(data = shape2, aes(long, lat, group = group),
    #              fill = NA, inherit.aes = FALSE) +
    # geom_path(data = shape2, aes(long, lat, group = group),
    #           color = "black", linewidth = 0.2, inherit.aes = FALSE) +

    geom_polygon(data = shape, aes(long, lat, group = group),
                 fill = NA, inherit.aes = FALSE) +
    geom_path(data = shape, aes(long, lat, group = group),
              color = "black", linewidth = 0.9, inherit.aes = FALSE) +
    coord_equal() +
    theme_map() +
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=12))

  return (p)
}


#############################################################################
## probocc_plot
## Holds default settings for plotting incidence rates
# IR_plot <- function(raster=NULL, shape=NULL, afss=NULL){
probocc_plot <- function(raster=NULL, shape=NULL, shape2=NULL,
                         color_ramp="YlOrBr", rev=FALSE, level=9, title=NULL,
                         theme=theme_map()){
  # source("util/ggplot2_theme.R")
  library(raster)
  library(ggplot2)
  library(RColorBrewer)
  library(scico)
  if (is.null(shape)){
    shape <- readRDS("data/africa_sub_Sahara_adm0_shp.rds")
  }
  if (is.null(shape2)){
    shape2 <- readRDS("data/africa_sub_Sahara_adm1_shp.rds")
  }
  mypal <- brewer.pal(level, color_ramp)
  if (rev) mypal <- rev(mypal)

  rpts <- as.data.frame(raster::rasterToPoints(raster))
  colnames(rpts) <- c("lon", "lat", "probocc")
  # identify a way to plot the incidence rate on a log scale
  p <- ggplot(rpts) +
    geom_raster(aes(lon, lat, fill=probocc)) +
    scale_fill_gradientn(colors=mypal, "Probability of\noccurrence") +
    geom_polygon(data=shape, aes(long, lat, group = group),
                 fill=NA, inherit.aes=FALSE) +
    geom_path(data = shape, aes(long, lat, group = group),
              color = "black", linewidth = 0.8, inherit.aes = FALSE) +
    geom_polygon(data = shape2, aes(long, lat, group = group),
                 fill = NA, inherit.aes = FALSE) +
    geom_path(data = shape2, aes(long, lat, group = group),
              color = "black", linewidth = 0.4, inherit.aes = FALSE) +
    coord_equal() +
    theme +
    labs(title=title) +
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          plot.title = element_text(size=22))

  return (p)
}


# Plot a raster image using ggplot2
polygon_plot <- function(polygon = NULL,
                         title = "", xlab = "", ylab = "",
                         title_size = 16, legend_title = "",
                         legend_position = c(0.2, 0.3)) {
  library(viridis)
  library(raster)
  library(plyr)
  if(!"id" %in% names(polygon)) {
    polygon$id <- 1:nrow(polygon)
  }
  # grid$fitted_pois <- mod3$summary.fitted.values$mean * grid$cellarea
  poly_pts <- broom::tidy(polygon, region = "id")
  poly_df <- plyr::join(poly_pts, polygon@data, by = "id")
  library(ggplot2)
  p <- ggplot(poly_df) +
    geom_polygon(aes(long, lat, group = group), fill = "gray95", color = "black") +
    scale_color_discrete(legend_title) +
    labs(title = title, x = xlab, y = ylab) +
    theme(panel.background = element_blank(), # bg of the panel
          plot.background = element_blank(), # bg of the plot
          legend.background = element_blank(), # get rid of legend bg
          legend.box.background = element_blank(),
          panel.spacing = unit(c(0,0,0,0), "null"),
          plot.margin = unit(c(0,0,0,0), "null"),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          plot.title = element_text(hjust = 0.5, size = title_size))

  return (p)
}
