###########################################################################
## take the occurrence data and prepare them for a brt model to use
prepare_brt_data <- function( data, raster, abs_occ_ratio = 2 ){
    occ_lonlat <- data.frame( lon=data$lon, lat=data$lat )
    occ_lonlat <- unique( occ_lonlat[ c( "lon", "lat" ) ] )
    occ_lonlat <- occ_lonlat[ complete.cases(occ_lonlat), ]
    covars_occ <- raster::extract( raster, occ_lonlat )
    covars_occ <- covars_occ[ complete.cases(covars_occ), ]
    raster_sample <- raster
    if( length( names( raster ) ) > 1  ) { # is raster stack or brick?
        raster_sample <- raster::subset( raster, 1, drop=TRUE ) # then take the first one as a sample
    }
    abs_lonlat <- dismo::randomPoints( raster_sample, nrow(occ_lonlat)*abs_occ_ratio, p=occ_lonlat, prob=FALSE )
    covars_abs <- raster::extract( raster, abs_lonlat )
    covars_abs <- covars_abs[ complete.cases(covars_abs), ]
    occ_abs_dat <- c( rep(1, nrow(covars_occ)), rep(0, nrow(covars_abs)) )  # 1 and 0 representing presence and absence, respectively
    occ_abs_dat <- data.frame( cbind( occ_abs_dat, rbind( covars_occ, covars_abs ) ) )
    occ_abs_dat <- occ_abs_dat[ complete.cases(occ_abs_dat), ]

    return( occ_abs_dat )
}

###########################################################################
##
poly_occ_random <- function( data, raster ){
    occ_poly <- data %>%
        dplyr::filter( point_loc == 0 ) %>%  # not a point location
        dplyr::filter( !is.na( adm_1 ) ) # remove the records where only adm_0 is available

    samples_poly <- sample_from_polygon_data( occ_poly, raster )# prec is used, but any other raster is fine as long as it has the relevent geo-graphical region (in this case, Africa)
    occ_poly_new <- occ_poly
    occ_poly_new$lon <- samples_poly$lon
    occ_poly_new$lat <- samples_poly$lat

    occ_poly_new <- dplyr::filter( occ_poly_new, !is.na(lon), !is.na(lat) )
    return ( occ_poly_new )
}


###########################################################################
##
point_occ_random <- function( data, raster, xlim = 0.1, ylim = 0.1 ){
    if( !require(dplyr) ){
        cat( "Loading dplyr package...\n" )
        library(dplyr)
    }
    if( !require(raster) ){
        cat( "Loading raster package...\n" )
        library(raster)
    }
    if( !is.data.frame(data) ){
        stop( " data must be a data frame" )
    }
    occ_points <- filter( data, point_loc == 1, !is.na(lon), !is.na(lat) )
    occ_points_new <- occ_points
    drift_lon <- runif( nrow(occ_points), -xlim, xlim )
    drift_lat <- runif( nrow(occ_points), -ylim, ylim )
    for( i in 1 : nrow(occ_points) ){
        occ_points_new$lon[ i ] <- occ_points$lon[ i ] + drift_lon[ i ]
        occ_points_new$lat[ i ] <- occ_points$lat[ i ] + drift_lat[ i ]
    }
    # check if the points are out of bounds by evaluating values for the one of the covariates (e.g., precipitation )
    if( length( names( raster ) ) > 1 ){
        raster <- raster::subset( raster, 1, drop=TRUE )
    }
    vals <- raster::extract( raster,
                             data.frame( lon = occ_points_new$lon, lat = occ_points_new$lat ) )
    index <- which( is.na( vals ) )

    occ_points_new <- occ_points_new[ - index, ]
    occ_ob <- occ_points[ index, ]
    occ_ob_new <- occ_ob
    max_iter <- 100
    iter <- 0
    na_index <- c()
    for( i in 1 : nrow(occ_ob) ){
        occ_ob_new$lon[ i ] <- occ_ob$lon[ i ] + runif( 1, -xlim, xlim )
        occ_ob_new$lat[ i ] <- occ_ob$lat[ i ] + runif( 1, -ylim, ylim )
        val <- raster::extract( raster, data.frame( lon=occ_ob_new$lon[ i ], lat=occ_ob_new$lat[ i ] ) )
        while( is.na(val) & iter < max_iter ){
            occ_ob_new$lon[ i ] <- occ_ob$lon[ i ] + runif( 1, -xlim, xlim )
            occ_ob_new$lat[ i ] <- occ_ob$lat[ i ] + runif( 1, -ylim, ylim )
            val <- raster::extract( raster,
                                    data.frame( lon=occ_ob_new$lon[ i ], lat=occ_ob_new$lat[ i ] ) )
            iter <- iter + 1
        }
        if( iter == max_iter ){
            na_index <- c( na_index, i )
            # cat( "i =", i, ", max iteration,", max_iter, ", reached!\n")
        }
        iter <- 0
    }
    occ_ob_new <- occ_ob_new[ - na_index, ]
    occ_points_new <- bind_rows( occ_points_new, occ_ob )

    return ( occ_points_new )
}


###########################################################################
## sample of longitude (x) and latitude (y) from a list (or stack) of masked rasters
sample_from_raster <- function( raster, sample_size=1 ){
    xy <- data.frame()
    if( !is.null(raster) ){
        # cat( "i =", i, "\n")
        sample <- raster::sampleRandom( raster, sample_size, xy=TRUE ) # sampling without replacement, thus sample size is limited by the number of cells available in the raster
        sample_xy <- c( sample[,1], sample[,2] ) # each column becomes an individual sample
        xy <- rbind( xy, sample_xy )
    } else{
        x <- rbind( x, NA ) # each column becomes an individual sample
        y <- rbind( y, NA )
    }
    colnames( xy ) <- c( "lon", "lat" )
    return( xy )
}


###########################################################################
##
## sample of longitude (x) and latitude (y) from a list (or stack) of masked rasters
lonlat_sample_masked_raster <- function( raster, sample_size=1 ){
    rlist <- raster
    if( class( raster ) == "list" ){
        rlist <- ( raster )
    }
    else if( class( raster ) == "RasterStack" ){
        rlist <- unstack( raster )
    }
    x <- data.frame()
    y <- data.frame()
    for( i in 1:length(rlist) ){
        if( !is.null(rlist[[i]]) ){
            # cat( "i =", i, "\n")
            sample <- sampleRandom( rlist[[i]], sample_size, xy=TRUE ) # sampling without replacement, thus sample size is limited by the number of cells available in the raster
            x <- rbind( x, sample[,1] ) # each column becomes an individual sample
            y <- rbind( y, sample[,2] )
        } else{
            x <- rbind( x, NA ) # each column becomes an individual sample
            y <- rbind( y, NA )
        }

    }
    colnames(x) <- 1:ncol(x)
    colnames(y) <- 1:ncol(y)
    xy <- list( X=x, Y=y )
    return( xy )
}
##########################################################################################3
sample_from_polygon_data <- function( data, sample_raster ){
    if( !require(countrycode) ){
        library(countrycode)
    }
    lonlat_sample <- data.frame( lon = rep( NA, length(data$lon) ),
                                 lat = rep( NA, length(data$lat) ),
                                 yr = data$year_begin)
    adm_level <- 3 # default being the smallest unit
    fls_rds <- list.files( path="ADM_Shapefile", pattern= "*.rds$", full.names=TRUE )
    for( i in 1 : nrow( data ) ){
        # cat( "i =", i, "\n" )
        adm_level <- 3 # default being the smallest unit
        if( is.na( data$adm_3[ i ] ) ){ # ADM_3 is NA?
            adm_level <- 2
            if( is.na( data$adm_2[ i ] ) ){ # ADM_2 is NA?
                adm_level <- 1
                if( is.na( data$adm_1[ i ] ) ){ # ADM_1 is NA?
                    adm_level <- 0
                }
            }
        }
        if( adm_level > 0 ){
            code3 <- countrycode( data$adm_0[ i ], "country.name", "iso3c" ) # get the three-letter country code (ISO 3)
            txt <- paste0( code3, "_adm", adm_level )
            index <- sapply( txt, grep, as.character( fls_rds ) ) # ex of fls_rds file = "GADM_2.8_KEN_adm3.rds"
            if( is.integer( index[ 1 ] ) ){
                adm_poly <- readRDS( fls_rds[ index[ 1 ] ] )
            } else {
                adm_poly <- raster::getData( 'GADM', country = code3, level = adm_level )
            }
            if( nrow( data.frame( adm_poly ) ) > 0 ){
                name <- paste0( "adm_poly@data$NAME_", adm_level ) # this is variable name (when evaluated) to select by adm name
                ## devise a way to increase efficiency of finding the pattern
                dat_col <-
                    dplyr::select( data, contains( "adm" ) ) %>%
                    dplyr::select( contains( eval( as.character( adm_level ) ) ) )

                pattern <- as.character( dat_col[i,1] )
                # cat( "i =", i, ", pattern =", pattern, "\n" )
                if( !is.na( pattern ) ){
                    adm_poly_sub <- raster::subset( adm_poly, grepl( pattern, eval( parse( text=name ) ) ) )
                    if( nrow( data.frame( adm_poly_sub ) ) > 0 ){
                        raster_masked <- raster::mask( sample_raster, adm_poly_sub )
                        xy <- sample_from_raster( raster_masked )
                        lonlat_sample$lon[ i ] <- xy$lon
                        lonlat_sample$lat[ i ] <- xy$lat
                        rm( raster_masked )
                    } else {
                        # cat( "i =", i, txt, "The ADM name does not exist in the polygon file\n" )
                    }
                } else {
                    # cat( "i =", i, txt, "polygon files were not found\n" )
                }
            } else {
                # cat( "i =", i, "only adm_0 is available\n" )
            }
        }
    }
    return( lonlat_sample )
}


################################################################
## transform a shapefile into a dataframe fot ggplot2
shape_to_dataframe <- function( shp=NULL ){
  stopifnot( class(shp) == "SpatialPolygonsDataFrame" )
  shp@data$id <- rownames( shp@data )
  shp_points <- broom::tidy( shp, region="id" )
  shp_df <- dplyr::left_join( shp_points, shp@data, by="id" )
}

################################################################
## random_drift()
## at a given point, select a point satisfying the distance between the two points are exponentially distributed with the mean_dist
## mean_dist in terms of a degree, of which one degree is ~111 km on the equator
random_drift <- function( xy=NULL, mean_dist=0.1, rst=NULL ){
  stopifnot( is.data.frame(xy) && !is.null(rst) && mean_dist > 0 )
  n <- nrow( xy )
  r <- rexp( n, 1/mean_dist ) # sample a radius
  t <- 2*pi*runif( n ) # distribute a random sample from a circle with radius r
  newx <- xy$lon + r*cos(t)
  newy <- xy$lat + r*sin(t)
  # identify the points where raster value is not available
  vals <- raster::extract( rst, cbind( newx, newy ) )
  nas <- which( is.na( vals ) ) #
  newxy <- data.frame( lon = newx[ -nas ], lat = newy[ -nas ], yr = xy$year_begin[ - nas ] )
  iter = 0
  # try to select other random points where
  while( length( nas ) > 0 && iter < 100 ){
    n <- length( nas )
    r <- rexp( n, 1/mean_dist ) # sample a radius
    t <- 2*pi*runif( n ) # distribute a random sample from a circle with radius r
    x <- r*cos( t )
    y <- r*sin( t )
    newx <- xy$lon[ nas ] + x # nas has the original ids for na values
    newy <- xy$lat[ nas ] + y
    vals <- raster::extract( rst, cbind( newx, newy ) )
    nas2 = which( is.na( vals ) )
    fixed = cbind( lon = newx[ -nas2 ], lat = newy[ -nas2 ], yr = xy$year_begin[ nas[ -nas2 ] ] )
    newxy <- rbind( newxy, fixed )
    nas = nas[ nas2 ]
    iter <- iter + 1
    # cat( "iter =", iter, ", nas =", length(nas),"\n"  )
  }
  return( newxy )
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

###############################################################################
## poly_from_coord
## create polygons from four coordinates
polygon_from_coord <- function(coord){
  # library(stringr)
  # library(sp)
  val <- as.numeric(stringr::str_extract_all(coord, "[0-9.-]+")[[1]])
  y <- val[c(1,3,5,7)]
  x <- val[-c(1,3,5,7)]
  xy <- cbind(x, y)
  p <- Polygon(xy)
  ps <- Polygons(list(p), 1)
  sps <- SpatialPolygons(list(ps))
  proj4string(sps) <-
    CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  return (sps)
}


raster_plot <- function(raster = NULL, title = "", xlab = "", ylab = "",
                        legend_title = "",
                        title_size = 16,
                        legend_position = c(0.25, 0.3),
                        breaks = c(10, 100, 1000, 10000),
                        better_color = FALSE,
                        afss = NULL,
                        afssadm1 = NULL,
                        overlay_polygon = FALSE) {

  library(raster)
  library(ggplot2)
  source("util/ggplot2_theme.R")
  if (is.null(afssadm1)){
    afssadm1 <- readRDS("data/africa_sub_Sahara_adm1_shp.rds")
  }
  if (is.null(afss)){
    afss <- readRDS("data/africa_sub_Sahara_adm0_shp.rds")
  }

  rpts <- as.data.frame(raster::rasterToPoints(raster))
  colnames(rpts) <- c("lon", "lat", "val")
  plt <- ggplot(data = rpts, aes(x=lon, y=lat, fill=val)) +
    geom_raster() +
    {if (better_color)
      scale_fill_gradientn(trans = "log10", breaks = breaks,
                           colors = c("grey93", "steelblue", "gold", "darkred"),
                         legend_title)} +
    labs(title=title, x=xlab, y=ylab) +
    geom_path(data = afssadm1, aes(long, lat, group = group),
              color = "black", size = 0.2, inherit.aes = FALSE) +
    geom_polygon(data = afss, aes(long, lat, group = group),
                 fill = NA, inherit.aes = FALSE) +
    geom_path(data = afss, aes(long, lat, group = group),
              color = "black", size = 0.9, inherit.aes = FALSE) +

    coord_equal() +
    theme_map() +
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=12))

  return (plt)
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

######################################################################
## tidy_year
## adjust years accounting for that covariates are years 2000-2017
tidy_year <- function(ya, yb) {
  if (is.na(ya) & is.na(yb)){
    ya <- 2010
    yb <- 2010
  }
  else if(!is.na(ya) & is.na(yb)){
    yb <- ya
  }
  else if(is.na(ya) & !is.na(yb)){
    ya <- yb
  }

  if(ya < 2000){
    ya <- 2000
  }
   if(ya > 2017){
    ya <- 2017
  }
  if(yb < 2000){
    yb <- 2000
  }

  if(yb > 2017){
    yb <- 2017
  }

  return(c(ya, yb))
}

#############################################################################
## get_raster_value
## get raster values based on the xy location and catchment area size(ie,
## number of grid cells)
## extracts values of the neighboring cells when catchment are size is larger
## than a single cell (ie, grid = 9) or the particular cell returns NA

get_raster_value <- function(raster, xy, grid, buffer = 30*1000){
  names(xy) <- c("X", "Y")
  val <- NA
  if (sum(is.na(xy)) == 0){
    val <- raster::extract(raster, xy, method = "simple", buffer = NULL)
    if (grid == 9 | is.na(val)) {
      val <- raster::extract(raster, xy, method = "simple", buffer = buffer)
      val <- unlist(val)
    }
    val <- mean(val, na.rm=T)
  }
  return (val)
}


##########################################################################
## get_covariates
## get covariate based on location
get_covariates <- function(data,
                           covariates,
                           lon_name = "lon",
                           lat_name = "lat",
                           year_begin_name = "STUDY_DATE_BEGIN",
                           year_end_name = "STUDY_DATE_END"){


  distw <- covariates[["distance_to_water_20km_af_2010_20220720"]]
  elev <- covariates[["elevation_20km_africa"]]
  travel <- covariates[["travel_time_cities_20km_af_2017_20220802"]]

  ## create a dataframe that will hold covariates and locations
  occ <- data.frame(X = as.numeric(data[, lon_name, drop=T]),
                    Y = as.numeric(data[, lat_name, drop=T]),
                    grid = as.integer(data[, c("grid"), drop=T]))
  nr <- nrow(occ)

  occ$annual_mean_temp <- rep(NA, nr)
  occ$annual_rainfall <- rep(NA, nr)
  occ$distance_to_water <- rep(NA, nr)
  occ$elevation <- rep(NA, nr)
  occ$improved_sanitation <- rep(NA, nr)
  occ$improved_water <- rep(NA, nr)
  occ$prev_HIV <- rep(NA, nr)
  occ$prev_stunting <- rep(NA, nr)
  occ$travel_time_to_cities <- rep(NA, nr)

  occ$open_defecation <- rep(NA, nr)
  occ$piped_sanitation <- rep(NA, nr)
  occ$piped_water <- rep(NA, nr)
  occ$ppp <- rep(NA, nr)
  occ$surface_water <- rep(NA, nr)
  occ$underweight <- rep(NA, nr)
  occ$wasting <- rep(NA, nr)

  for(i in 1:nr) {
    pt <- occ[i, c("X", "Y")]
    grid <- occ[i, c("grid")]
    # for some reason, raster::extract does not recognize these as numeric

    occ$travel_time_to_cities[i] <- get_raster_value(travel, pt, grid)
    occ$elevation[i] <- get_raster_value(elev, pt, grid)
    occ$distance_to_water[i] <- get_raster_value(distw, pt, grid)

    hiv_list <- list()
    stunt_list <- list()
    temp_list <- list()
    rain_list <- list()
    sanitation_list <- list()
    water_list <- list()

    open_defecation_list <- list()
    piped_sanitation_list <- list()
    piped_water_list <- list()
    ppp_list <- list()
    surface_water_list <- list()
    underweight_list <- list()
    wasting_list <- list()


    cnt <- 1

    varnames <- names(data)

    if (year_begin_name %in% varnames && year_end_name %in% varnames) {
      ya <- data[i, year_begin_name, drop=T]
      yb <- data[i, year_end_name, drop=T]
      yayb <- tidy_year(ya, yb)
      ya <- yayb[1]
      yb <- yayb[2]

      yrs <- seq(ya, yb)
      for (yr in yrs){
        hiv <- covariates[[paste0("prev_HIV_adults_20km_af_", yr, "_20220802")]]
        stunt <- covariates[[paste0("prev_stunting_20km_af_", yr, "_20220802")]]
        temp <- covariates[[paste0("annual_mean_temperature_20km_af_", yr, "_20220802")]]
        rain <- covariates[[paste0("annual_rainfall_20km_af_", yr, "_20220802")]]
        sanitation <- covariates[[paste0("improved_sanitation_20km_af_", yr, "_20220802")]]
        water <- covariates[[paste0("improved_water_20km_af_", yr, "_20220802")]]

        open_defecation <- covariates[[paste0("open_defecation_20km_af_", yr, "_20220802")]]
        piped_sanitation <- covariates[[paste0("piped_sanitation_20km_af_", yr, "_20220802")]]
        piped_water <- covariates[[paste0("piped_water_20km_af_", yr, "_20220802")]]
        ppp <- covariates[[paste0("ppp_20km_af_", yr, "_20220929")]]
        surface_water <- covariates[[paste0("surface_water_20km_af_", yr, "_20220802")]]
        underweight <- covariates[[paste0("underweight_20km_af_", yr, "_20220802")]]
        wasting <- covariates[[paste0("wasting_20km_af_", yr, "_20220802")]]

        hiv_list[[ cnt ]] <- get_raster_value(hiv, pt, grid)
        stunt_list[[ cnt ]] <- get_raster_value(stunt, pt, grid)
        temp_list[[ cnt ]] <- get_raster_value(temp, pt, grid)
        rain_list[[ cnt ]] <- get_raster_value(rain, pt, grid)
        sanitation_list[[ cnt ]] <- get_raster_value(sanitation, pt, grid)
        water_list[[ cnt ]] <- get_raster_value(water, pt, grid)

        open_defecation_list[[ cnt ]] <- get_raster_value(open_defecation, pt, grid)
        piped_sanitation_list[[ cnt ]] <- get_raster_value(piped_sanitation, pt, grid)
        piped_water_list[[ cnt ]] <- get_raster_value(piped_water, pt, grid)
        ppp_list[[ cnt ]] <- get_raster_value(ppp, pt, grid)
        surface_water_list[[ cnt ]] <- get_raster_value(surface_water, pt, grid)
        underweight_list[[ cnt ]] <- get_raster_value(underweight, pt, grid)
        wasting_list[[ cnt ]] <- get_raster_value(wasting, pt, grid)


        cnt <- cnt + 1
      }
      occ$prev_HIV[i] <- mean(unlist(hiv_list), na.rm = TRUE)
      occ$prev_stunting[i] <- mean(unlist(stunt_list), na.rm = TRUE)
      occ$annual_mean_temp[i] <- mean(unlist(temp_list), na.rm = TRUE)
      occ$annual_rainfall[i] <- mean(unlist(rain_list), na.rm = TRUE)
      occ$improved_sanitation[i] <- mean(unlist(sanitation_list), na.rm = TRUE)
      occ$improved_water[i] <- mean(unlist(water_list), na.rm = TRUE)

      occ$open_defecation[i] <- mean(unlist(open_defecation_list), na.rm = TRUE)
      occ$piped_sanitation[i] <- mean(unlist(piped_sanitation_list), na.rm = TRUE)
      occ$piped_water[i] <- mean(unlist(piped_water_list), na.rm = TRUE)
      occ$ppp[i] <- mean(unlist(ppp_list), na.rm = TRUE)
      occ$surface_water[i] <- mean(unlist(surface_water_list), na.rm = TRUE)
      occ$underweight[i] <- mean(unlist(underweight_list), na.rm = TRUE)
      occ$wasting[i] <- mean(unlist(wasting_list), na.rm = TRUE)

    }
    else {
      # when there is no information on the year, we take the mid point
      yr <- 2010
      hiv <- covariates[[paste0("prev_HIV_adults_20km_af_", yr, "_20220802")]]
      stunt <- covariates[[paste0("prev_stunting_20km_af_", yr, "_20220802")]]
      temp <- covariates[[paste0("annual_mean_temperature_20km_af_", yr, "_20220802")]]
      rain <- covariates[[paste0("annual_rainfall_20km_af_", yr, "_20220802")]]
      sanitation <- covariates[[paste0("improved_sanitation_20km_af_", yr, "_20220802")]]
      water <- covariates[[paste0("improved_water_20km_af_", yr, "_20220802")]]

      open_defecation <- covariates[[paste0("open_defecation_20km_af_", yr, "_20220802")]]
      piped_sanitation <- covariates[[paste0("piped_sanitation_20km_af_", yr, "_20220802")]]
      piped_water <- covariates[[paste0("piped_water_20km_af_", yr, "_20220802")]]
      ppp <- covariates[[paste0("ppp_20km_af_", yr, "_20220929")]]
      surface_water <- covariates[[paste0("surface_water_20km_af_", yr, "_20220802")]]
      underweight <- covariates[[paste0("underweight_20km_af_", yr, "_20220802")]]
      wasting <- covariates[[paste0("wasting_20km_af_", yr, "_20220802")]]


      occ$prev_HIV[i] <- raster::extract(hiv, pt, method = "simple")
      occ$prev_stunting[i] <- raster::extract(stunt, pt, method = "simple")
      occ$annual_mean_temp[i] <- raster::extract(temp, pt, method = "simple")
      occ$annual_rainfall[i] <- raster::extract(rain, pt, method = "simple")
      occ$improved_sanitation[i] <- raster::extract(sanitation, pt, method = "simple")
      occ$improved_water[i] <- raster::extract(water, pt, method = "simple")

      occ$open_defecation[i] <- raster::extract(open_defecation, pt, method = "simple")
      occ$piped_sanitation[i] <- raster::extract(piped_sanitation, pt, method = "simple")
      occ$piped_water[i] <- raster::extract(piped_water, pt, method = "simple")
      occ$ppp[i] <- raster::extract(ppp, pt, method = "simple")
      occ$surface_water[i] <- raster::extract(surface_water, pt, method = "simple")
      occ$underweight[i] <- raster::extract(underweight, pt, method = "simple")
      occ$wasting[i] <- raster::extract(wasting, pt, method = "simple")
    }
  }

  # leave only covariates
  occ <- subset(occ, select = -c(X, Y, grid))
  return(occ)
}

###############################################################################
#
average_by_polygon <- function(rst = NULL, adj = NULL,
                   poly = NULL,  admlevel = 0){

  pg_names <- c()
  df <- c()
  if (admlevel == 0) {
    poly$NAME_0 <- clean_country_names(poly$NAME_0)
    pg_names <- unique(poly$NAME_0)
    df <- data.frame(matrix(NA, nrow = length(pg_names), ncol = 2))
    names(df) <- c("adm0", "val")
    df$adm0 <- pg_names
  } else if (admlevel == 1) {
    poly$NAME_0 <- clean_country_names(poly$NAME_0)
    poly$NAME_1 <- clean_country_names(poly$NAME_1)
    pg_names <- dplyr::distinct(data.frame(adm0=poly$NAME_0, adm1=poly$NAME_1))
    df <- data.frame(matrix(NA, nrow = nrow(pg_names), ncol = 3))
    names(df) <- c("adm0", "adm1", "val")
    df$adm0 <- pg_names$adm0
    df$adm1 <- pg_names$adm1
  }

  for (i in 1:nrow(df)) {
    pg <- c()
    if (admlevel == 0) {
      pg <- poly[poly$NAME_0 == df$adm0[i], ]
    } else if (admlevel == 1) {
      pg <- poly[poly$NAME_0 == df$adm0[i] & poly$NAME_1 == df$adm1[i], ]
    }

    prob <- raster::extract(rst, pg)[[1]]
    if (!is.null(adj)) {
      pop <- raster::extract(adj, pg)[[1]] # adjust by population size
      casetot <- sum(pop*prob, na.rm=TRUE)
      df$val[i] <- casetot/ sum(pop, na.rm=TRUE)
    } else {
      df$val[i] <- mean(prob, na.rm=TRUE)
    }
  }
  return (df)
}


drawmap <- function(rst, shp, legend_title = "") {
  rmask <- mask(rst, shp)
  rpts <- rasterToPoints(rmask)
  rptsdf <- as.data.frame(rpts)
  colnames(rptsdf) <- c("lon", "lat", "ir")
  rptsdf$ir <- as.double(rptsdf$ir + 1) # 1 is added to use log
  # summary(rptsdf$ir)
  maxir <- 1e4 # arbitrary maximum incidence rate
  rptsdf$ir[rptsdf$ir > maxir] <- maxir
  # identify a way to plot the incidence rate on a log scale
  library(ggplot2)
  p <- ggplot(rptsdf) +
    geom_raster(aes(lon, lat, fill = ir)) +
    scale_fill_gradientn(trans = "log", breaks = c(10, 100, 1000, 10000),
                         colors = c("grey93", "steelblue", "gold", "darkred"),
                         legend_title) +
    geom_polygon(data = shp, aes(long, lat, group = group),
                 fill = NA, inherit.aes = FALSE) +
    geom_path(data = shp, aes(long, lat, group = group),
              color = "black", size = 0.9, inherit.aes = FALSE) +
    coord_equal() +
    theme_map() +
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=12))

  return(p)
}

# extract covariates for incidence data with single row (d) using covariate raster(rst)
ext_cov_single <- function(rst, d){ # extract cov from a single row data.frame
  covariate_val <- NA
  loc <- matrix(c(d$LONGITUDE, d$LATITUDE), nrow=1, ncol=2, byrow=T)
  cells <- NULL
  if(grepl(".shp$", d$SHAPE_FILE)){
    loc <- raster::shapefile(d$SHAPE_FILE)
    cells <- raster::extract(rst, loc, cellnumbers = TRUE)
  }
  else {
    cells <- raster::extract(rst, loc,
                             buffer = sqrt(d$CATCHMENT_AREA_KM2)*1000, cellnumbers = TRUE)
  }
  if(!is.null(cells) & sum(!is.na(cells)) > 0) {
    cells <- as.data.frame(do.call(rbind, cells))
    cells <- cells[!duplicated(cells$cell),]
    covariate_val <- cells[, "value"]
    # val <- val[!is.na(val)]
    # covariate_val <- mean(val)
  }
  return(covariate_val)
}

# extract covariates for incidence data (d) using covariate raster(rst)
extract_covariates <- function(rst, d, func = "mean") {
  covariate_vals <- rep(NA, nrow(d))
  if(nrow(d) > 1){
    for(i in 1:nrow(d)){
      di <- d[i,]
      if(i == 1 | (i > 1 & !(identical(d[i, c(1:6)], d[i-1, c(1:6)])))) {
        if (func == "mean") {
          covariate_vals[i] <- mean(ext_cov_single(rst, di), na.rm=TRUE)
        }
        if (func == "median") {
          covariate_vals[i] <- median(ext_cov_single(rst, di), na.rm=TRUE)
        }
      }
      else{
        covariate_vals[i] <- covariate_vals[i-1]
      }
    }
  }
  return(covariate_vals)
}


# computes the incidence rate by region (subnational region, country, or Africa
# subregion)
# compute the number of cases per region based on the population per grid and
# incidence rate per grid
case_by_region <- function(case = NULL, shape=NULL,
                           region=c("country","subnational", "subregion")) {
  library(raster)
  library(data.table)
  if(is.null(shape)) {
    shape <- readRDS("data/africa_sub_Sahara_adm1_shp.rds")
  }
  afregions <- fread("data/africa_subregion_country_names.csv")

  names(case) <- "value" # change the name for writing convenience
  # results will be returned in raster format
  casenew <- case
  # Create a vector with the same length.
  # This is for presumed efficiency but not tested
  casenew_vec <- rep(NA, length(case[])) # initialize the vector with NA

  areas <- unique(shape$NAME_0)
  if(region == "subnational"){
    areas <- unique(shape$NAME_1)
  }
  else if (region == "subregion") {
    areas <- unique(afregions$Subregion)
  }
  df <- data.frame(matrix(NA, nrow=length(areas), ncol=2))
  names(df) <- c("Area", "Value")
  df$Area <- areas
  for (ar in areas) {
    if (region == "country"){
      poly <- shape[shape$NAME_0 == ar, ] # SpatialPolygon
    }
    else if (region == "subnational"){
      poly <- shape[shape$NAME_1 == ar, ] # SpatialPolygon
    }
    else { # Africa subregion
      cntries <- afregions[Subregion == ar, Country]
      poly <- shape[shape$NAME_0 %in% cntries, ] # SpatialPolygon
    }

    case_grid <- raster::extract(case, poly, df = TRUE, cellnumbers = TRUE)
    case_area <- sum(case_grid$value, na.rm=T)

    casenew_vec[unlist(case_grid$cell)] <- case_area
    df[df$Area == ar, ]$Value <- case_area

  }
  casenew[] <- casenew_vec
  # change the column names for easier interpretation
  names(df) <- c("Area", "Case")
  return(list(raster = casenew, data = df))
}

# computes the incidence rate by region (subnational region, country, or Africa
# subregion)
# compute the number of cases per region based on the population per grid and
# incidence rate per grid
IR_by_region <- function(ppp = NULL, ir = NULL, shape = NULL,
                         region = c("country","subnational","subregion"),
                         PYO = 1e5) {
  library(raster)
  library(data.table)
  if (is.null(ppp)) {
    ppp <- readRDS("data/covariates/prediction/ppp_20km_af_2017_20221208.rds")
  }
  if (is.null(shape)) {
    shape <- readRDS("data/africa_sub_Sahara_adm1_shp.rds")
  }
  afregions <- fread("data/africa_subregion_country_names.csv")

  names(ppp) <- "value" # change the name for writing convenience
  names(ir) <- "value" # change the name for writing convenience

  # results will be returned in raster format
  rst <- ir
  # Create a vector with the same length.
  # This is for presumed efficiency but not tested
  rstvec <- rep(NA, length(rst[])) # initialize the vector with NA

  areas <- unique(shape$NAME_0)
  if (region == "subnational") {
    areas <- unique(shape$NAME_1)
  }
  else if (region == "subregion") {
    areas <- unique(afregions$Subregion)
  }

  df <- data.frame(matrix(NA, nrow=length(areas), ncol=4))
  names(df) <- c("Area", "Pop", "IR", "NumCase")
  df$Area <- areas

  for (ar in areas) {
    if (region == "country"){
      poly <- shape[shape$NAME_0 == ar, ] # SpatialPolygon
    }
    else if (region == "subnational"){
      poly <- shape[shape$NAME_1 == ar, ] # SpatialPolygon
    }
    else { # Africa subregion
      cntries <- afregions[Subregion == ar, Country]
      poly <- shape[shape$NAME_0 %in% cntries, ] # SpatialPolygon
    }

    ppp_grid <- raster::extract(ppp, poly, df = TRUE, cellnumbers = TRUE)
    ir_grid <- raster::extract(ir, poly, df = TRUE, cellnumbers = TRUE)

    # dummy variable is created to indicate where IR is valid and
    # population size (denominator) is only used when IR is valid
    ir_dummy <- as.integer(!is.na(ir_grid$value))
    pop_area <- sum(ppp_grid$value * ir_dummy, na.rm=T)
    # case_area would mean that the expected number of cases
    case_area <- sum(ir_grid$value * ppp_grid$value, na.rm=T)
    ir_area <- case_area / pop_area
    rstvec[unlist(ppp_grid$cell)] <- ir_area

    df[df$Area == ar,]$Pop <- pop_area
    df[df$Area == ar,]$IR <- ir_area
    df[df$Area == ar,]$NumCase <- case_area / PYO
  }
  rst[] <- rstvec
  # change the column names for easier interpretation
  return(list(raster = rst, data = df))
}


#-----------------------------------------------------------------------------
# population by age and country based on the population grid and
# the country shape file and returns a dataframe
#
pop_by_age_country <- function(ppp=NULL, age_prop=NULL, shape=NULL) {
  library(raster)
  library(data.table)
  if(is.null(ppp)) {
    ppp <- readRDS("data/covariates/prediction/ppp_20km_af_2017_20221208.rds")
  }
  if(is.null(age_prop)) {
    age_prop <- readRDS("data/pop_prop_by_age.rds")
  }
  if(is.null(shape)) {
    shape <- readRDS("data/africa_sub_Sahara_adm0_shp.rds")
  }
  names(ppp) <- "value" # change the name for writing convenience
  # afregions <- fread("data/africa_subregion_country_names.csv")

  cntries <- unique(shape$NAME_0)

  age <- c("0_1y", "2_4y", "5_14y", "over14y")
  df <- expand.grid(Country = cntries,
                    Age = c(age, "all"),
                    Pop=NA)

  for (cn in cntries) {
    props <-
      age_prop[country == cn, c(prop_0_1y,prop_2_4y,prop_5_14y,prop_over14y)]
    print(paste0("country = ", cn))
    poly <- shape[shape$NAME_0 == cn, ] # SpatialPolygon
    cells <- raster::extract(ppp, poly, df = TRUE, cellnumbers = TRUE)

    df[(df$Country == cn & df$Age == "all"),]$Pop <-
      sum(cells$value, na.rm=T)
    # df[(df$Country == cn & df$Age == "all"),]$Pop <- format(round(sum(cells$value, na.rm=T), digits=digits), big.mark=",", trim=TRUE)
    for (j in 1:length(age)) {
      df[(df$Country == cn & df$Age == age[j]),]$Pop <-
        sum(cells$value, na.rm=T) * props[j]
      # df[(df$Country == cn & df$Age == age[j]),]$Pop <-
      #  format(round(sum(cells$value, na.rm=T) * props[j], digits=digits), big.mark=",", trim=TRUE)
    }
  }

  return(df)

}

