#' @title debug_vdatum
#' @description debug_vdatum
#' @param ehydro_root path to location you want to store ehydro data
#' @param process_index index of sheet you want to stage
#' @param vdatum_domain_override passthrough for vdatum region override, can be ak, as, contiguous, gcnmi, prvi, westcoast, chesapeak_delaware, if left NULL will use 'contiguous', Default: NULL
#' @returns TRUE
#' @details Puts string together to test the response of https://vdatum.noaa.gov/vdatumweb/
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  ePydRo::debug_vdatum(ehydro_root, process_index = pr_surveys[3,]$sheet_index, vdatum_domain_override = 'prvi')
#'  }
#' }
#' @seealso
#'  \code{\link[sf]{st_read}}, \code{\link[sf]{st_layers}}, \code{\link[sf]{geos_unary}}, \code{\link[sf]{geos_combine}}, \code{\link[sf]{st_as_sf}}, \code{\link[sf]{st_transform}}, \code{\link[sf]{st_crs}}, \code{\link[sf]{st_coordinates}}, \code{\link[sf]{st_write}}
#'  \code{\link[httr]{GET}}, \code{\link[httr]{write_disk}}, \code{\link[httr]{http_error}}, \code{\link[httr]{content}}
#'  \code{\link[data.table]{fread}}, \code{\link[data.table]{fwrite}}
#'  \code{\link[glue]{glue}}
#'  \code{\link[lubridate]{decimal_date}}, \code{\link[lubridate]{ymd}}
#'  \code{\link[lidR]{LAS-class}}, \code{\link[lidR]{las_check}}, \code{\link[lidR]{writeLAS}}
#' @rdname scrape_and_transform
#' @export
#' @importFrom sf st_read st_layers st_centroid st_union st_as_sf st_transform st_crs st_coordinates st_write
#' @importFrom httr GET write_disk http_error content
#' @importFrom data.table fread fwrite
#' @importFrom glue glue
#' @importFrom lubridate decimal_date ymd
#' @importFrom lidR LAS las_check writeLAS

debug_vdatum = function(ehydro_root,process_index,vdatum_domain_override) {
  # sinew::moga(file.path(getwd(),"R/debug_vdatum.R"),overwrite = TRUE)
  #
  # devtools::load_all()
  # ehydro_root = "~/data/ePydRo"
  # process_index <- 23795
  # vdatum_domain_override='prvi'

  ## -- Start --
  fn_time_start <- Sys.time()

  # Light input parsing
  if(is.null(vdatum_domain_override)) {
    vdatum_domain_override <- "contiguous"
  }
  if(vdatum_domain_override == "prvi") {
    t_v_frame_override <- "prvd92"
  }

  # Load sheet
  epydro_tracker <- sf::st_read(file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep))
  epydro_tracker_entry <- epydro_tracker[epydro_tracker$sheet_index==process_index,]
  epydro_tracker_sourcedata_name <- gsub('.{4}$', '', basename(epydro_tracker_entry$sourcedata))

  gdb_path <-file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,paste0(epydro_tracker_sourcedata_name,".gdb"),fsep = .Platform$file.sep)
  gdb_present <- file.exists(gdb_path)

  if(!gdb_present) {
    message("Why!")
    message(glue::glue("it thinks you want: {epydro_tracker_sourcedata_name}"))
    return(FALSE)
  }

  # GDB
  gdb_pts <- NULL
  if(gdb_present) {
    layernames <- sf::st_layers(gdb_path)$name
    if('SurveyPoint_HD' %in% layernames) {
      message('Using HD geodatabase field')
      gdb_pts <- sf::st_read(gdb_path,layer = 'SurveyPoint_HD')
    } else if('SurveyPoint' %in% layernames) {
      message('Using thinned geodatabase field')
      gdb_pts <- sf::st_read(gdb_path,layer = 'SurveyPoint')
    } else {
      message('Missing field')
      gdb_pts <- NULL
    }
  }

  # Get offset: loose get_vdatum_corrrection start...
  ## Geodatabase (in)sanity checks
  if(!is.null(gdb_pts)) {
    if(unique(gdb_pts$elevationDatum) %>% length() > 1 || unique(gdb_pts$elevationUOM) %>% length() > 1 || unique(gdb_pts$SurveyDateStamp) %>% length() > 1) {
      print_waring_block()
      message("Mixed point datums?")
    }
    in_date <- unique(gdb_pts$SurveyDateStamp)[1]
    in_datum_unit <- unique(gdb_pts$elevationUOM)[1]
    in_datum <- unique(gdb_pts$elevationDatum)[1]

    # Parse from a projection
    out <- get_datum_from_crs(gdb_pts)
    proj_datum <- out[1]
    proj_datum_unit <- out[2]

    if(!(proj_datum==in_datum) || !(proj_datum_unit==in_datum_unit)) {
      message('Projection and survey datum mismaatch, defaulting to sheet values')
      message(glue::glue('From sheet:{in_datum} with {in_datum_unit}'))
      message(glue::glue('From projection:{proj_datum} with {proj_datum_unit}'))
    }
    if(!(epydro_tracker_entry$surveydate==in_date)) {
      message('Date mismatch? default to surveyjob')
    }
  }

  # Date manipulations
  this_date_YYYY <- format(epydro_tracker_entry$surveydate, format = "%Y")
  this_date_mm <- format(epydro_tracker_entry$surveydate, format = "%m")
  this_date_dd <- format(epydro_tracker_entry$surveydate, format = "%d")
  this_date_start <- lubridate::decimal_date(lubridate::ymd(paste0(this_date_YYYY, "-", this_date_mm, "-", this_date_dd)))

  # Geography
  in_mean_pt <- sf::st_centroid(sf::st_union(gdb_pts)) %>% sf::st_as_sf()
  in_mean_pt_proj <- sf::st_transform(in_mean_pt,sf::st_crs("EPSG:6349"))

  # Cartography scratch
  # mapview::mapview(gdb_pts) +
  #   mapview::mapview(in_mean_pt,color = "red")

  # Datum unit normalization (to meter)
  if (in_datum_unit == "usSurveyFoot") {
    elev_unit_norm = (1200 / 3937)
  } else if (in_datum_unit == "Foot") {
    elev_unit_norm = 0.3048
  } else {
    elev_unit_norm = 1
  }

  message(glue::glue("Trying {epydro_tracker_sourcedata_name}"))
  message("Sheet row:")
  message(epydro_tracker_entry)
  message('')
  message(glue::glue('get_vdatum_correction(midpoint = {as.character(sf::st_coordinates(in_mean_pt_proj)[1, ][2])} / {as.character(sf::st_coordinates(in_mean_pt_proj)[1, ][1])}'))
  message(glue::glue('in_region_override={vdatum_domain_override},s_v_unit_override="m", t_v_unit_override="m",s_h_frame_override="NAD83_2011",s_v_frame_override={in_datum}'))
  message(glue::glue('t_h_frame_override="NAD83_2011",t_v_frame_override={t_v_frame_override},in_epoch_override={this_date_start},out_epoch_override={lubridate::decimal_date(Sys.time())},is_quiet = is_quiet)'))

  datum_url <- get_vdatum_correction(midpoint = in_mean_pt_proj,
                                     in_region_override=vdatum_domain_override,
                                     s_v_unit_override="m",
                                     t_v_unit_override="m",
                                     s_h_frame_override="NAD83_2011",
                                     s_v_frame_override=in_datum,
                                     t_h_frame_override="NAD83_2011",
                                     t_v_frame_override=t_v_frame_override,
                                     in_epoch_override=this_date_start,
                                     out_epoch_override=lubridate::decimal_date(Sys.time()),
                                     is_quiet = is_quiet)
  message(glue::glue('Attempted: {datum_url}'))
  return(TRUE)
}
