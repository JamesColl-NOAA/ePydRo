#' @title get_vdatum_correction
#' @description Constructs and gets the offset between datums from the [vdatum](https://vdatum.noaa.gov/vdatumweb/) api
#' @param midpoint the point in space you need to transform
#' @param in_region_override vdatum region override, can be ak, as, contiguous, gcnmi, prvi, westcoast, chesapeak_delaware, Default: 'contiguous'
#' @param s_v_unit_override source vertical unit, Default: 'm'
#' @param t_v_unit_override target vertical unit, Default: 'm'
#' @param s_h_frame_override source horozontal projection, Default: 'NAD83_2011'
#' @param s_v_frame_override source vertical datum, Default: 'MLLW'
#' @param t_h_frame_override target horozontal projeection, Default: 'NAD83_2011'
#' @param t_v_frame_override target vertical datum, Default: 'NAVD88'
#' @param in_epoch_override in epoch (survey time) as decimal date, Default: lubridate::decimal_date(Sys.time())
#' @param out_epoch_override out epoc (desired timestamp) as decimapl date, Default: lubridate::decimal_date(Sys.time())
#' @param is_quiet flag to determine whether print statements are suppressed, TRUE to suppress messages and FALSE to show them, Default: FALSE
#' @returns numeric
#' @details [vdatum](https://vdatum.noaa.gov/vdatumweb/) is a complex api.  This function does not help that.  Because we know the unit transformation needed to fix the datum to our normalized standard (merters as spesified in sf::st_crs("EPSG:6349")), I can make a lot of assumptions about how I get the datum normalization value from this utility.  This function expects a reprojected point and auxilary data and calculates the difference between the input datum at 0 and NAVD88.  From there, the correction can be recorded and applied.  Make note of the pretty direct ramp I deploy to get to this point.  Many of these properties could be infered other ways.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  # See API docs at https://vdatum.noaa.gov/docs/services.html#step170
#'  }
#' }
#' @seealso
#'  \code{\link[lubridate]{decimal_date}}
#'  \code{\link[glue]{glue}}
#' @rdname get_vdatum_correction
#' @export
#' @importFrom lubridate decimal_date
#' @importFrom glue glue

get_vdatum_correction = function(midpoint,
                                 in_region_override="contiguous",
                                 s_v_unit_override="m",
                                 t_v_unit_override="m",
                                 s_h_frame_override="NAD83_2011",
                                 s_v_frame_override="MLLW",
                                 t_h_frame_override="NAD83_2011",
                                 t_v_frame_override="NAVD88",
                                 in_epoch_override=lubridate::decimal_date(Sys.time()),
                                 out_epoch_override=lubridate::decimal_date(Sys.time()),
                                 is_quiet = FALSE) {
  # sinew::moga(file.path(getwd(),"R/get_vdatum_correction.R"),overwrite = TRUE)
  #
  #
  # pt = midpoint
  # units = "SI Units"
  # proj_string = "EPSG:26915"
  # in_epoch_override = epydro_tracker_entry$surveydate
  # out_epoch_override = as.integer(as.POSIXct(Sys.time()))
  # is_quiet = FALSE

  ## -- Start --
  fn_time_start <- Sys.time()

  # Awful name normalization: https://vdatum.noaa.gov/docs/services.html#step170
  if(s_v_frame_override=="NAVD_88"){
    s_v_frame_override <- "NAVD88"
  } else if(s_v_frame_override=="North American Datum 1983") {
    s_v_frame_override <- "NAD83_2011"
  }

  if(in_region_override=='chesapeak_delaware') {
    s_h_frame_override <- 'IGS14'
  }

  # rm(resp)
  datum_url <- glue::glue("https://vdatum.noaa.gov/vdatumweb/api/convert?region={in_region_override}&",
                          "s_x={as.character(sf::st_coordinates(midpoint)[1, ][1])}&s_y={as.character(sf::st_coordinates(midpoint)[1, ][2])}&",
                          "s_v_unit={s_v_unit_override}&t_v_unit={t_v_unit_override}&",
                          "s_h_frame={s_h_frame_override}&s_v_frame={s_v_frame_override}&",
                          "t_h_frame={t_h_frame_override}&t_v_frame={t_v_frame_override}&",
                          "&epoch_in={in_epoch_override}&epoch_out={out_epoch_override}")
  return(datum_url)

  # resp <- httr::GET(datum_url)
  # if (httr::http_error(resp)) {
  #   print_warning_block()
  #   message(paste('poorly formed url - Request URL:', datum_url))
  # }
  # jsonRespParsed <- httr::content(resp, as = "parsed")
  # if(!is.numeric(as.numeric(jsonRespParsed$t_z))) {
  #   print_warning_block()
  #   message("Something broke")
  #   jsonRespParsed$t_z <- 0
  # }
  #
  # return(jsonRespParsed$t_z)
}
