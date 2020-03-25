#' To determine which points are in the sampling area
#' 
#' @param extent
#' @param sp A SpatialPolygons representation of the section of the lake that can be sampled
#' @param dist Numerical argument to control the within transect set spacing
#' @return Returns a data.frame with the coordinates over which data should be simulated
#' @export
crop.sample.area <- function(extent, sp, dist, units){
  x_min <- extent[1]
  x_max <- extent[2]
  y_min <- extent[3]
  y_max <- extent[4]
  dx <- get.geo.dist(x_min, y_min, x_max, y_min, units)
  dy <- get.dist(x_min, y_min, x_min, y_max, units)
  n_x <- floor(dx/dist)
  n_y <- floor(dy/dist)
  x_grid <- seq(x_min, x_max, length = n_x)
  y_grid <- seq(y_min, y_max, length = n_y)
  xy_grid <- expand.grid(x_grid, y_grid)
  names(xy_grid) <- c("x", "y")
  xy_samp <- xy_grid[point.in.SpatialPolygons(xy_grid$x, xy_grid$y, sp),]
  return(xy_samp)
}


