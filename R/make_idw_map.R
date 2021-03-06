#' Make IDW maps of CPUE for the EBS/NBS
#' 
#' This function can be used to make inverse-distance-weighted plots for the eastern Bering Sea and northern Bering Sea
#' 
#' @param x Data frame which contains at minimum: CPUE, LATITUDE, and LONGITUDE. Can be passed as vectors instead (see below). Default value: \code{NA}
#' @param region Character vector indicating which plotting region to use. Options: bs.south, bs.north, bs.all
#' @param grid.cell Numeric vector of length two specifying dimensions of grid cells for extrpolation grid, in degrees c(lon, lat). Default c(0.05, 0.05)
#' @param COMMON_NAME Common name
#' @param LATITUDE Latitude (degrees north)
#' @param LONGITUDE Longitude (degrees east; Western hemisphere is negative)
#' @param CPUE_KGHA Catch per unit effort in kilograms per hectare
#' @param extrap.box Optional. Vector specifying the dimensions of the extrapolation grid. Elements of the vector should be named to specify the minimum and maximum x and y values c(xmn, xmx, ymn, ymx). If not provided, region will be used to set the extrapolation area.
#' @param set.breaks Suggested. Vector of break points to use for plotting. Alternatively, a character vector indicating which break method to use. Default = "jenks"
#' @param grid.cell Optional. Numeric vector of length two, specifying the resolution for the extrapolation grid in degrees. Default c(0.05,0.05)
#' @param in.crs Character vector containing the coordinate reference system for projecting the extrapolation grid.
#' @param out.crs Character vector containing the coordinate reference system for projecting the extrapolation grid.
#' @param key.title Character vector which will appear in the legend above CPUE (kg/ha). Default = "auto" tries to pull COMMON_NAME from input.
#' @param log.transform Character vector indicating whether CPUE values should be log-transformed for IDW. Default = FALSE.
#' @param idw.nmax Maximum number of adjacent stations to use for interpolation. Default = 8
#' @param idp Inverse distance weighting power. Default = 2 
#' @param use.survey.bathymetry Logical indicating if historical survey bathymetry should be used instead of continuous regional bathymetry. Default = TRUE
#' @param return.continuous.grid If TRUE, also returns an extrapolation grid on a continuous scale.
#' @return Returns a list containing: 
#' (1) plot: a ggplot IDW map;
#' (2) extrapolation.grid: the extrapolation grid with estimated values on a discrete scale;
#' (3) continuous.grid: extrapolation grid with estimates on a continuous scale;
#' (4) region: the region;
#' (5) n.breaks: the number of level breaks;
#' (6) key.title: title for the legend;
#' (7) crs: coordinate reference system as a PROJ6 (WKT2:2019) string; 
#' @author Sean Rohan \email{sean.rohan@@noaa.gov}
#' @export

make_idw_map <- function(x = NA, 
                         COMMON_NAME = NA, 
                         LATITUDE = NA, 
                         LONGITUDE = NA, 
                         CPUE_KGHA = NA, 
                         region = "bs.south", 
                         extrap.box = NA, 
                         set.breaks = "jenks", 
                         grid.cell = c(0.05, 0.05), 
                         in.crs = "+proj=longlat", 
                         out.crs = "auto", 
                         key.title = "auto", 
                         log.transform = FALSE, 
                         idw.nmax = 4,
                         use.survey.bathymetry = TRUE, 
                         return.continuous.grid = TRUE) {
  
  # Convert vectors to data fram if x is NA---------------------------------------------------------
  if(is.na(x)) {
    x <- data.frame(COMMON_NAME = COMMON_NAME,
                    LATITUDE = LATITUDE,
                    LONGITUDE = LONGITUDE,
                    CPUE_KGHA = CPUE_KGHA)
  }
  
  # Set legend title--------------------------------------------------------------------------------
  if(key.title == "auto") {
    key.title <- x$COMMON_NAME[1]
  }
  
  # Set up mapping region---------------------------------------------------------------------------
  if(is.na(extrap.box)) {
    if(region %in% c("bs.south", "sebs")) {extrap.box = c(xmn = -179.5, xmx = -157, ymn = 54, ymx = 63)}
    if(region %in% c("bs.all", "ebs")) {extrap.box = c(xmn = -179.5, xmx = -157, ymn = 54, ymx = 68)}
  }
  
  # Load map layers---------------------------------------------------------------------------------
  map_layers <- akgfmaps::get_base_layers(select.region = region, set.crs = out.crs)
  
  # Assign CRS to handle automatic selection--------------------------------------------------------
  out.crs <- map_layers$crs
  
  # Use survey bathymetry---------------------------------------------------------------------------
  if(use.survey.bathymetry) {
    map_layers$bathymetry <- get_survey_bathymetry(select.region = region, set.crs = out.crs)
  }
  
  # Assign CRS to input data------------------------------------------------------------------------
  x <- sf::st_as_sf(x, coords = c(x = "LONGITUDE", y = "LATITUDE"), crs = sf::st_crs(in.crs)) %>% 
    sf::st_transform(crs = map_layers$crs)
  
  # Inverse distance weighting----------------------------------------------------------------------
  idw_fit <- gstat::gstat(formula = CPUE_KGHA~1, locations = x, nmax = idw.nmax)
  
  # Predict station points--------------------------------------------------------------------------
  stn.predict <- predict(idw_fit, x)

  # Generate extrapolation grid---------------------------------------------------------------------
  sp_extrap.raster <- raster::raster(xmn = extrap.box['xmn'],
                                     xmx=extrap.box['xmx'],
                                     ymn=extrap.box['ymn'],
                                     ymx=extrap.box['ymx'],
                                     ncol=(extrap.box['xmx']-extrap.box['xmn'])/grid.cell,
                                     nrow=(extrap.box['ymx']-extrap.box['ymn'])/grid.cell,
                                     crs = crs(in.crs)) %>% 
    projectRaster(crs = crs(x))
  
  # Predict, rasterize, mask------------------------------------------------------------------------
  extrap.grid <- predict(idw_fit, as(sp_extrap.raster, "SpatialPoints")) %>% 
    sf::st_as_sf() %>% 
    sf::st_transform(crs = crs(x)) %>% 
    stars::st_rasterize() %>% 
    sf::st_join(map_layers$survey.area, join = st_intersects) 
  
  # Return continuous grid if return.continuous.grid is TRUE ---------------------------------------
  if(return.continuous.grid) {
    continuous.grid <- extrap.grid
  } else {
    continuous.grid <- NA
  }

  # Format breaks for plotting----------------------------------------------------------------------
  # Automatic break selection based on character vector.
  alt.round <- 0 # Set alternative rounding factor to zero based on user-specified breaks
  
  if(is.character(set.breaks[1])) {
    set.breaks <- tolower(set.breaks)
    
    # Set breaks ----
    break.vals <- classInt::classIntervals(x$CPUE_KGHA, n = 5, style = set.breaks)$brks
    
    # Setup rounding for small CPUE ----
    alt.round <- floor(-1*(min((log10(break.vals)-2)[abs(break.vals) > 0])))
    
    set.breaks <- c(-1, round(break.vals, alt.round))
  }
  
  # Ensure breaks go to zero------------------------------------------------------------------------
  if(min(set.breaks) > 0) {
    set.breaks <- c(0, set.breaks)
  }
  
  if(min(set.breaks) == 0) {
    set.breaks <- c(-1, set.breaks)
  }
  
  # Ensure breaks span the full range---------------------------------------------------------------
  if(max(set.breaks) < max(stn.predict$var1.pred)){
    set.breaks[length(set.breaks)] <- max(stn.predict$var1.pred) + 1
  }
  
  
  # Trim breaks to significant digits to account for differences in range among species-------------
  dig.lab <- 7
  set.levels <- cut(stn.predict$var1.pred, set.breaks, right = TRUE, dig.lab = dig.lab)
  
  if(alt.round > 0) {
    while(dig.lab > alt.round) { # Rounding for small CPUE
        dig.lab <- dig.lab - 1
        set.levels <- cut(stn.predict$var1.pred, set.breaks, right = TRUE, dig.lab = dig.lab)
      }
    } else { # Rounding for large CPUE
      while(length(grep("\\.", set.levels)) > 0) {
        dig.lab <- dig.lab - 1
        set.levels <- cut(stn.predict$var1.pred, set.breaks, right = TRUE, dig.lab = dig.lab)
      }
    }
      
  # Cut extrapolation grid to support discrete scale------------------------------------------------
  extrap.grid$var1.pred <- cut(extrap.grid$var1.pred, set.breaks, right = TRUE, dig.lab = dig.lab)
  # extrap.grid$discrete_layer <- cut(extrap.grid$layer, set.breaks, right = TRUE, dig.lab = dig.lab)

  # Which breaks need commas?-----------------------------------------------------------------------
  sig.dig <- round(set.breaks[which(nchar(round(set.breaks)) >= 4)])
  
  # Drop brackets, add commas, create 'No catch' level to legend labels-----------------------------
  make_level_labels <- function(vec) {
    vec <- as.character(vec)
    vec[grep("-1", vec)] <- "No catch"
    vec <- sub("\\(", "\\>", vec)
    vec <- sub("\\,", "–", vec)
    vec <- sub("\\]", "", vec)
    if(length(sig.dig) > 3) {
      for(j in 1:length(sig.dig)) {
        vec <- sub(sig.dig[j], format(sig.dig[j], nsmall=0, big.mark=","), vec)
      }
    }
    return(vec)
  }
  
  # Assign level names to breaks for plotting-------------------------------------------------------
  extrap.grid$var1.pred <- factor(make_level_labels(extrap.grid$var1.pred), levels = make_level_labels(levels(set.levels)))
  
  # Number of breaks for color adjustments----------------------------------------------------------
  n.breaks <- length(levels(set.levels))
  
  # Make plot---------------------------------------------------------------------------------------
    p1 <- ggplot() +
    geom_sf(data = map_layers$survey.area, fill = NA) +
    geom_stars(data = extrap.grid) +
      geom_sf(data = map_layers$survey.area, fill = NA) +
      geom_sf(data = map_layers$akland, fill = "grey80") +
      geom_sf(data = map_layers$bathymetry) +
      geom_sf(data = map_layers$graticule, color = alpha("grey70", 0.3)) +
      scale_fill_manual(name = paste0(key.title, "\nCPUE (kg/ha)"), 
                        values = c("white", RColorBrewer::brewer.pal(9, name = "Blues")[c(2,4,6,8,9)]), 
                        na.translate = FALSE, # Don't use NA
                        drop = FALSE) + # Keep all levels in the plot
      scale_x_continuous(breaks = map_layers$lon.breaks) + 
      scale_y_continuous(breaks = map_layers$lat.breaks) + 
      coord_sf(xlim = map_layers$plot.boundary$x,
               ylim = map_layers$plot.boundary$y) +
      theme(panel.border = element_rect(color = "black", fill = NA),
            panel.background = element_rect(fill = NA, color = "black"),
            legend.key = element_rect(fill = NA, color = "grey70"),
            legend.position = c(0.12, 0.18),
            axis.title = element_blank(),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            plot.background = element_rect(fill = NA, color = NA))
  
  return(list(plot = p1,
              map_layers = map_layers,
              extrapolation.grid = extrap.grid,
              continuous.grid = continuous.grid,
              region = region,
              n.breaks = n.breaks,
              key.title = key.title,
              crs = out.crs))
}
