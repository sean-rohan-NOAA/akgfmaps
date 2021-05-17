#' Calculate cold pool area in the southeastern Bering Sea using varying interpolation methods
#' 
#' Calculate cold pool area using nearest-neighbor interpolation, inverse distance weighting, and ordinary kriging methods.
#' 
#' @param dat Input data frame. Must include columns for latitude, longitde, and variable to be interpolated.
#' @param dat.year Year
#' @param var.col Character vector denoting the name of the variable column.
#' @param lat.col Character vector indicating the name of the latitude column.
#' @param lon.col Character vector indicating the name of the longitude column.
#' @param in.crs Character vector (PROJ4 or WKT2/PROJ6 string) indicating the coordinate reference system for the data.
#' @param interpolation.crs Character vector (PROJ4 or WKT2/PROJ6 string) indicating the coordinate reference system to use for the interpolation. Historically EPSG:3338
#' @param cell.resolution Dimensions of interpolation grid cell in meters.
#' @param nm Maximum number of cells to use for interpolation. 
#' @param pre Prefix for file names in output (in development.)
#' @param write.to.file If true, writes output to GeoTIFF format file. Default = FALSE.
#' @export

calculate_cold_pool_area <- function(dat,
                                dat.year,
                                var.col,
                                lat.col,
                                lon.col,
                                in.crs = "+proj=longlat",
                                interpolation.crs,
                                cell.resolution,
                                nm = Inf,
                                pre = NA,
                                write.to.file = FALSE)
{
  names(dat)[which(names(dat) == var.col)] <- "var.col"
  names(dat)[which(names(dat) == lat.col)] <- "lat.col"
  names(dat)[which(names(dat) == lon.col)] <- "lon.col"
  
  # Load EBS survey exent for masking ----
  sebs_layers <- akgfmaps::get_base_layers(select.region = "sebs", 
                                           set.crs = interpolation.crs)
  
  # Set dimensions for raster cells ----
  n_dim <- floor(abs(-1625000 - -35000))/cell.resolution
  
  # Make raster for interpolation ----
  sp_interp.raster <- raster::raster(xmn = -1625000, 
                                     xmx = -35000, 
                                     ymn = 379500, 
                                     ymx = 1969500, 
                                     nrow = n_dim, 
                                     ncol = n_dim)
  raster::projection(sp_interp.raster) <- interpolation.crs
  
  # Transform data for interpolation ----
  sp_interp.df <- unique(dat)
  sp::coordinates(sp_interp.df) <- c(x = "lon.col", y = "lat.col")
  sp::proj4string(sp_interp.df) <- sp::CRS(in.crs)
  sp_interp.df <- sp::spTransform(sp_interp.df, sp::CRS(interpolation.crs))
  
  # Nearest-neighbor
  nn_fit <- gstat::gstat(formula = var.col ~ 1, 
                         locations = sp_interp.df, 
                         set = list(idp = 0), 
                         nmax = 4)
  nn.predict <- predict(nn_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_nn <- rasterize_and_mask(nn.predict, 
                               sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_nn <- rasterize_and_mask(nn.predict, 
                                    sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Inverse distance weighting with nmax = 4 (Lauth method) ----
  idw_nmax4_fit <- gstat::gstat(formula = var.col ~ 1, 
                          locations = sp_interp.df,
                          set = list(idp = 2),
                          nmax = 4)
  idw_nmax4.predict <- predict(idw_nmax4_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_idw_nmax4 <- rasterize_and_mask(idw_nmax4.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_idw_nmax4 <- rasterize_and_mask(idw_nmax4.predict, 
                                    sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Inverse distance weighting  ----
  idw_fit <- gstat::gstat(formula = var.col ~ 1, 
                          locations = sp_interp.df,
                          set = list(idp = 2),
                          nmax = nm)
  idw.predict <- predict(idw_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_idw <- rasterize_and_mask(idw.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_idw <- rasterize_and_mask(idw.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Set up a new IDW for ordinary kriging ----
  idw_vgm_fit <- gstat::gstat(formula = var.col ~ 1, 
                              locations = sp_interp.df, 
                              nmax = nm)
  
  # Ordinary Kriging: Exponential VGM----
  exp.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Exp")))
  exp_fit <- gstat(formula = var.col ~ 1, 
                   locations = sp_interp.df,
                   model = exp.vgfit, 
                   nmax = nm)
  exp.predict <- predict(exp_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_exp <- rasterize_and_mask(exp.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_exp <- rasterize_and_mask(exp.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Ordinary Kriging: Spherical VGM----
  sph.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Sph")))
  sph_fit <- gstat::gstat(formula = var.col ~ 1,
                          locations = sp_interp.df, 
                          model = sph.vgfit, 
                          nmax = nm)
  sph.predict <- predict(sph_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_sph <- rasterize_and_mask(sph.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_sph <- rasterize_and_mask(sph.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Ordinary Kriging: Bessel VGM----
  bes.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Bes")))
  bes_fit <- gstat::gstat(formula = var.col ~ 1,
                          locations = sp_interp.df, 
                          model = bes.vgfit, 
                          nmax = nm)
  bes.predict <- predict(bes_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_bes <- rasterize_and_mask(bes.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_bes <- rasterize_and_mask(bes.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Ordinary Kriging: Bessel VGM----
  gau.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Gau")))
  gau_fit <- gstat::gstat(formula = var.col ~ 1, 
                          locations = sp_interp.df, 
                          model = gau.vgfit, 
                          nmax = nm)
  gau.predict <- predict(gau_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_gau <- rasterize_and_mask(gau.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_gau <- rasterize_and_mask(gau.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Ordinary Kriging: Circular VGM----
  cir.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Cir")))
  cir_fit <- gstat::gstat(formula = var.col ~ 1, 
                          locations = sp_interp.df, 
                          model = cir.vgfit, 
                          nmax = nm)
  cir.predict <- predict(cir_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_cir <- rasterize_and_mask(cir.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_cir <- rasterize_and_mask(cir.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Ordinary Kriging: Matern VGM----
  mat.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Mat")))
  mat_fit <- gstat::gstat(formula = var.col ~ 1, 
                          locations = sp_interp.df, 
                          model = mat.vgfit, 
                          nmax = nm)
  mat.predict <- predict(mat_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_mat <- rasterize_and_mask(mat.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_mat <- rasterize_and_mask(mat.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Ordinary Kriging: Stein's Matern VGM----
  ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Ste")))
  ste_fit <- gstat::gstat(formula = var.col ~ 1, 
                          locations = sp_interp.df, 
                          model = ste.vgfit, 
                          nmax = nm)
  ste.predict <- predict(ste_fit, as(sp_interp.raster, "SpatialGrid"))
  
  cpe_lte2_ste <- rasterize_and_mask(ste.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_ste <- rasterize_and_mask(ste.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  # Thin-plate spline ----
  tps_fit <- fields::Tps(coordinates(sp_interp.df), 
                         sp_interp.df$var.col)
  tps.predict <- raster::interpolate(sp_interp.raster, 
                                     tps_fit)
  
  cpe_lte2_tps <- rasterize_and_mask(tps.predict, 
                                sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 2)
  
  cpe_lte0_tps <- rasterize_and_mask(tps.predict, 
                                     sebs_layers$survey.area) %>% 
    cpa_from_raster(temperature_threshold = 0)
  
  cpe_df <- data.frame(year = dat.year,
                       nn_lte2 = cpe_lte2_nn, 
                       idw_nmax4_lte2 = cpe_lte2_idw_nmax4, 
                       idw_lte2 = cpe_lte2_idw, 
                       bes_lte2 = cpe_lte2_bes, 
                       exp_lte2 = cpe_lte2_exp, 
                       cir_lte2 = cpe_lte2_cir, 
                       gau_lte2 = cpe_lte2_gau,
                       sph_lte2 = cpe_lte2_sph, 
                       mat_lte2 = cpe_lte2_mat,
                       ste_lte2 = cpe_lte2_ste,
                       tps_lte2 = cpe_lte2_tps,
                       nn_lte0 = cpe_lte0_nn, 
                       idw_nmax4_lte0 = cpe_lte0_idw_nmax4,
                       idw_lte0 = cpe_lte0_idw,
                       bes_lte0 = cpe_lte0_bes, 
                       exp_lte0 = cpe_lte0_exp, 
                       cir_lte0 = cpe_lte0_cir, 
                       gau_lte0 = cpe_lte0_gau,
                       sph_lte0 = cpe_lte0_sph, 
                       mat_lte0 = cpe_lte0_mat,
                       ste_lte0 = cpe_lte0_ste,
                       tps_lte0 = cpe_lte0_tps)
  
  if(write.to.file) {
    print("Writing rasters")
    if(!dir.exists("./output/")) {
      dir.create("./output/")
    }
    if(!dir.exists("./output/raster/")) {
      dir.create("./output/raster/")
    }

    # write unmasked surfaces to raster
    raster::writeRaster(raster::raster(nn.predict), 
                        filename = paste0("./output/raster/nn_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(idw.predict), 
                        filename = paste0("./output/raster/idw_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(idw_nmax4.predict), 
                        filename = paste0("./output/raster/idw_nmax4_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(exp.predict), 
                        filename = paste0("./output/raster/exp_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(sph.predict), 
                        filename = paste0("./output/raster/sph_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(bes.predict), 
                        filename = paste0("./output/raster/bes_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(gau.predict), 
                        filename = paste0("./output/raster/gau_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(cir.predict), 
                        filename = paste0("./output/raster/cir_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(mat.predict), 
                        filename = paste0("./output/raster/mat_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(raster::raster(ste.predict), 
                        filename = paste0("./output/raster/ste_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
    
    raster::writeRaster(tps.predict, 
                        filename = paste0("./output/raster/tps_", dat.year, "_", var.col, ".tif" ),
                        format = "GTiff", 
                        overwrite = TRUE)
  }
  
  return(cpe_df)
  
}