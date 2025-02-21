#' @title to_tiff_and_stage
#' @description Creates tiff from arco, dumps it out in a folder
#' @param ehydro_root path to location you want to store ehydro data
#' @param process_index index of sheet you want to stage
#' @param domain PARAM_DESCRIPTION
#' @param is_quiet PARAM_DESCRIPTION, Default: FALSE
#' @returns OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  to_tiff_and_stage(ehydro_root,process_index,domain,is_quiet = FALSE)
#'  }
#' }
#' @seealso
#'  \code{\link[sf]{st_read}}, \code{\link[sf]{st_crs}}, \code{\link[sf]{st_zm}}
#'  \code{\link[glue]{glue}}
#'  \code{\link[data.table]{fread}}
#'  \code{\link[sfheaders]{sf_point}}
#'  \code{\link[raster]{rasterize}}, \code{\link[raster]{raster}}
#'  \code{\link[terra]{writeRaster}}
#' @rdname to_tiff_and_stage
#' @export
#' @importFrom sf st_read st_set_crs st_crs st_zm
#' @importFrom glue glue
#' @importFrom data.table fread
#' @importFrom sfheaders sf_point
#' @importFrom raster rasterize raster
#' @importFrom terra writeRaster

to_tiff_and_stage = function(ehydro_root,process_index,domain,is_quiet = FALSE) {
  # sinew::moga(file.path(getwd(),"R/to_tiff_and_stage.R"),overwrite = TRUE)
  #
  # devtools::load_all()
  # ehydro_root = "~/data/ePydRo"
  # is_quiet = FALSE
  # process_index <- 79
  # domain <- file.path(ehydro_root,"CONUS",fsep = .Platform$file.sep)

  ## -- Start --
  fn_time_start <- Sys.time()

  epydro_tracker <- sf::st_read(file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep))
  epydro_tracker_entry <- epydro_tracker[epydro_tracker$sheet_index==process_index,]
  epydro_tracker_sourcedata_name <- gsub('.{4}$', '', basename(epydro_tracker_entry$sourcedata))

  message(glue::glue("Rasterizing {epydro_tracker_sourcedata_name} for domain: {domain}"))

  if(!file.exists(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,"epydro_xyz.csv",fsep = .Platform$file.sep))) {
    print_warning_block()
    message("This survey has not yet been processed, attempting to do that now...")
    scrape_and_transform(ehydro_root = ehydro_root,sheet_index = sheet_index,overwrite=FALSE,vdatum_domain_override=NULL,is_quiet = is_quiet)
  }

  # Open and rasterize (to_tiff)
  dat <- data.table::fread(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,"epydro_xyz.csv",fsep = .Platform$file.sep))
  sfdat <- sfheaders::sf_point(dat,x = 'X',y = 'Y',z = 'Z') %>% sf::st_set_crs(sf::st_crs("EPSG:6349"))
  sfdat$Z_depth <- dat$Z

  # Stage in domain folder (and_stage)
  dir.create(domain)
  surface <- raster::rasterize(sf::st_zm(sfdat, drop = TRUE, what = "ZM"),raster::raster(sfdat,resolution = 0.003),field="Z_depth",fun="min")
  terra::writeRaster(surface,file.path(domain,paste0(epydro_tracker_sourcedata_name,".tif"),fsep = .Platform$file.sep),wopt=list(gdal=c("COMPRESS=LZW","of=COG")))

  # wrap up
  return(TRUE)
}
