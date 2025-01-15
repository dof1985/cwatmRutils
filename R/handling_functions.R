# Data handling functions - helpers

#' na_to_val
#'
#' The function takes a data object, e.g., matrix, RasterLayer, or a 
#' vector and convert NA values to a selected value (default is  0)
#'
#' @param r input data object (e.g., RasterLayer)
#' @param val alternative value to be assigned instead of NA (Default: 0)
#' @return the data object with all NAs replaced by the assigned value

na_to_val <- function(r, val = 0) {
  r[is.na(r)] <- val
  return(r)
}

#' val_to_na
#'
#' The function takes a data object, e.g., matrix, RasterLayer, or a 
#' vector and convert a selected value (default is  0) to NA
#'
#' @param r input data object (e.g., RasterLayer)
#' @param val selected value to be replaced by NA (Default: 0)
#' @return the data object with all selected values replaced by NAs

val_to_na <- function(r, val = 0) {
  r[r == val] <- NA
  return(r)
}

#' getCellsGrid
#'
#' Assigns the cell id into each RasterLayer cell
#'
#' @param r input RasterLayer
#' @return a RasterLayer with cell ids assigned in each cell

getCellsGrid <- function(r) {
  tmp <- r
  tmp[] <- seq_len(raster::ncell(r))
  return(tmp)
}