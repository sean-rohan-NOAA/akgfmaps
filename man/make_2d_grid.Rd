% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_2d_grid.R
\name{make_2d_grid}
\alias{make_2d_grid}
\title{Function to interpolate a 2D grid of environmental data set using ordinary kriging}
\usage{
make_2d_grid(
  x,
  extrapolation_2d_grid = NA,
  par_name,
  sel_year,
  vgm_init,
  log_transform = FALSE,
  transform = "none",
  nst_options = list(tails = "none", make_plot = FALSE),
  nmax = 500,
  make_variogram = FALSE,
  save_csv = TRUE,
  suffix = ""
)
}
\arguments{
\item{x}{Data frame with columns for an environmental variable, vessel, cruise, haul, stationid, latitude, longitude, and cdepth.}

\item{extrapolation_2d_grid}{A data frame containing latitude, longitude and cdepth columns. Default NA automatically loads extrapolation grid.}

\item{par_name}{Character vector indicating the name of the variable column.}

\item{sel_year}{Numeric vector indicating the survey year.}

\item{vgm_init}{List of initial parameters to pass to gstat::vgm function}

\item{log_transform}{Should variables be log_transformed prior to subsequent transformation?}

\item{transform}{What transformation method should be applied? Options "none" = No transformation; nst" = normal score transform based on cdf using akgfmaps::normal_transform() with backtransformation using akgfmaps::backtransform_normal(), based on GSLIB nscore and backtr; "boxcox" - Box-Cox transformation and back transformation; "z" = z-score transform to mean zero, variance one without using cdf.}

\item{nst_options}{List with options (tails, make_plot) for backtransformation when transform = "nst". See documentation for akgfmaps::backtransform_normal().}

\item{nmax}{Numeric. Maximum number of data points to use for kriging (default = 500.}

\item{save_csv}{Logical. Should the output be saved to a csv file?}

\item{suffix}{Character vector added to the end of the filename if save_csv is TRUE.}
}
\description{
Function to interpolate a 2D grid of environmental data set using ordinary kriging
}
