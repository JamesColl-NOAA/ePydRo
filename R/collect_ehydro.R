#' @title collect_ehydro
#' @description collectes and formats ehydro SurveyJob.shp for use in the ePydRo workflow
#' @param ehydro_root path to location you want to store ehydro data
#' @param overwrite if TRUE will overwrite any processed ehydro tracker if FALSE will attempt to merge SurveyJob.shp with existing tracker, Default: FALSE
#' @param is_quiet is_quiet flag to determine whether print statements are suppressed, TRUE to suppress messages and FALSE to show them, Default: FALSE
#' @returns OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  ePydRo::collect_ehydro(ehydro_root = "~/data/ePydRo",overwrite=FALSE,is_quiet = FALSE)
#'  }
#' }
#' @seealso
#'  \code{\link[glue]{glue}}
#'  \code{\link[sf]{st_read}}, \code{\link[sf]{st_write}}
#' @rdname collect_ehydro
#' @export
#' @importFrom glue glue
#' @importFrom sf st_read st_write

collect_ehydro = function(ehydro_root,overwrite=FALSE,is_quiet = FALSE) {
  # sinew::moga(file.path(getwd(),"R/collect_ehydro.R"),overwrite = TRUE)
  #
  # devtools::load_all()
  # ehydro_root = "~/data/ePydRo"
  # is_quiet = FALSE

  # Manually scrape and dump at root: https://geospatial-usace.opendata.arcgis.com/datasets/80a394bae6b547f1b5788074261e11f1

  #> names(raw_ehydro_catalog) for posterity
  #[1] "surveyjobi" "sdsid"      "sdsfeature" "sdsmetadat" "surveytype" "channelare" "dateupload" "usacedistr" "surveydate" "surveyda_1" "sourcedata" "sourceproj" "mediaidfk"  "projecteda" "sdsfeatu_1"
  #[16] "dateloaded" "datenotifi" "sourceda_1" "plotsheetl" "sourceagen" "globalid"   "geometry"

  ## -- Start --
  fn_time_start <- Sys.time()
  if (!is_quiet) {
    message(glue::glue("Parsing into {ehydro_root}"))
  }

  ehydro_zip <- list.files(path = ehydro_root, pattern = "eHydro_Survey_Data_.*\\.zip$",full.names = TRUE)
  unzip(ehydro_zip, exdir = ehydro_root)

  raw_ehydro_catalog <- sf::st_read(file.path(ehydro_root,"SurveyJob.shp",fsep = .Platform$file.sep))

  cleaned_ehydro_catalog <- raw_ehydro_catalog[!sf::st_is_empty(raw_ehydro_catalog),,drop=FALSE] |> sf::st_make_valid()

  epydro_catalog <- cleaned_ehydro_catalog
  epydro_catalog$sheet_index <- row.names(epydro_catalog)
  epydro_catalog$horizontal_crs <- epydro_catalog$sourceproj #in
  epydro_catalog$vertical_datum <- "" #in
  epydro_catalog$vertical_datum_offset <- ""
  epydro_catalog$pos_z_convention <- ""
  epydro_catalog$vertical_datum_conversion <- "" # Formatted as 'MSL = {epydro_catalog$vertical_datum} - {epydro_catalog$vertical_datum_offset}'
  epydro_catalog$epydro_entwine <- ""
  epydro_catalog$epydro_url <- ""
  sf::st_write(epydro_catalog,file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep))

  dir.create(file.path(ehydro_root,"survey",fsep = .Platform$file.sep), showWarnings = FALSE)

  return(TRUE)
}
