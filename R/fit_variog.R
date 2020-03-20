#' To estimate the degree of spatial autocoorelation in cpue count data using an exponential variogram
#' 
#' @param count_data A SpatialPointsDataFrame storing the cpue data
#' @return Returns a vector containing estimates of the nugget, partial sill, and range parameter respectively
fit.variog <- function(count_data){
  v <- variogram(count~1, count_data)
  v_fit <- fit.variogram(v, model = vgm("Exp"))
  nug <- v_fit$psill[1]
  p_sill <- v_fit$psill[2]
  range <- v_fit$range[2]
  return(c(nug=nug, p_sill=p_sill, range=range))
}
