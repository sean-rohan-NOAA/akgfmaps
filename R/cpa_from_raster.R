#' Calculate cold pool area from a raster
#' 
#' Calculate area of raster cells less than or equal to a temperature threshold.
#' 
#' @param x A RasterLayer object.
#' @param raster_units Character vector indicating units for x and y dimensions in the raster. Default = "m" for meters.
#' @param temperature_threshold Numeric. Temperature threshold for value.
#' @export

cpa_from_raster <- function(x, 
                                     raster_units = "m", 
                                     temperature_threshold) {
  if(raster_units == "m") {
    cell_area <- raster::res(x)[1]*raster::res(x)[2]/1e6 
  }
  
  n_cells <- sum(x@data@values <= temperature_threshold, 
                 na.rm = TRUE)
  
  total_area <- n_cells * cell_area
  
  return(total_area)
}