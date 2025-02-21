#' @title scrape_and_transform
#' @description download, normalize, and entwine ehydro
#' @param ehydro_root path to location you want to store ehydro data
#' @param process_index index of sheet you want to stage
#' @param overwrite overwrite files on disk if found, Default: FALSE
#' @param vdatum_domain_override passthrough for vdatum region override, can be ak, as, contiguous, gcnmi, prvi, westcoast, chesapeak_delaware, if left NULL will use 'contiguous', Default: NULL
#' @param is_quiet flag to determine whether print statements are suppressed, TRUE to suppress messages and FALSE to show them, Default: FALSE
#' @returns TRUE
#' @details Stages data in ARCO form
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  ePydRo::scrape_and_transform(ehydro_root = "~/data/ePydRo",process_index = 79,overwrite=FALSE,vdatum_domain_override=NULL,is_quiet = FALSE)
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

scrape_and_transform = function(ehydro_root,process_index,overwrite=FALSE,vdatum_domain_override=NULL,is_quiet = FALSE) {
  # sinew::moga(file.path(getwd(),"R/scrape_and_transform.R"),overwrite = TRUE)
  #
  # devtools::load_all()
  # ehydro_root = "~/data/ePydRo"
  # is_quiet = FALSE
  # process_index <- 79

  ## -- Start --
  fn_time_start <- Sys.time()

  # Light input parsing
  if(is.null(vdatum_domain_override)) {
    vdatum_domain_override <- "contiguous"
  }
  if(vdatum_domain_override == "prvi") {
    t_v_frame_override <- "prvd02"
  }

  # Load sheet
  epydro_tracker <- sf::st_read(file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep))
  epydro_tracker_entry <- epydro_tracker[epydro_tracker$sheet_index==process_index,]
  epydro_tracker_sourcedata_name <- gsub('.{4}$', '', basename(epydro_tracker_entry$sourcedata))

  # Remove if overwrite
  if(dir.exists(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,fsep = .Platform$file.sep)) & overwrite) {
    unlink(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,fsep = .Platform$file.sep), recursive=TRUE)
  }

  # Scrape
  ##################################################
  if(!file.exists(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,basename(epydro_tracker_entry$sourcedata),fsep = .Platform$file.sep))) {
    dir.create(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,fsep = .Platform$file.sep))
    httr::GET(epydro_tracker_entry$sourcedata,
              httr::write_disk(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,basename(epydro_tracker_entry$sourcedata),fsep = .Platform$file.sep)), overwrite=TRUE)
    unzip(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,basename(epydro_tracker_entry$sourcedata),fsep = .Platform$file.sep),
          exdir = file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,fsep = .Platform$file.sep))
  }

  # And transform
  ##################################################
  if(!file.exists(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,"epydro","ept.json",fsep = .Platform$file.sep))) {

    # Loose priority list
    full_survey_path <- file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,paste0(epydro_tracker_sourcedata_name,"_FULL.XYZ"),fsep = .Platform$file.sep)
    gdb_path <-file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,paste0(epydro_tracker_sourcedata_name,".gdb"),fsep = .Platform$file.sep)
    tin_path <- file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,paste0(epydro_tracker_sourcedata_name,"_tin"),fsep = .Platform$file.sep)

    full_survey_present <- file.exists(full_survey_path)
    gdb_present <- file.exists(gdb_path)
    tin_present <- file.exists(tin_path)

    # Table
    xyz_pts <- NULL
    if(full_survey_present) {
      message('Using XYZ')
      full_pts <- data.table::fread(full_survey_path)
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

    # TIN
    tin_pts <- NULL
    if(tin_present) {
      message('Using TIN - Not yet implemented.')
      # tin_pts <- data.table::fread(tin_path)
    }

    # Final points frame

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
    # rm(resp)
    # resp <- httr::GET(datum_url)
    # jsonRespParsed <- httr::content(resp, as = "parsed")
    # jsonRespParsed$t_z

    resp <- httr::GET(datum_url)
    if (httr::http_error(resp)) {
      print_warning_block()
      message(paste('poorly formed url - Request URL:', datum_url))
    }
    jsonRespParsed <- httr::content(resp, as = "parsed")

    if(is.null(jsonRespParsed$t_z) || jsonRespParsed$t_z=="-999999") {
      print_warning_block()
      message(glue::glue("You broke someting on {epydro_tracker_sourcedata_name}, vdatum of 0 is used..."))
      jsonRespParsed$t_z <- 0
    } else {
      if(!is_quiet) { message(glue::glue("Ran {datum_url} and got {jsonRespParsed$t_z}")) }
    }

    out_z = gdb_pts$surveyPointElev * elev_unit_norm + as.numeric(jsonRespParsed$t_z)
    corrected_xy <- sf::st_transform(gdb_pts,sf::st_crs("EPSG:6349"))
    out_x <- sf::st_coordinates(corrected_xy)[, 1]
    out_y <- sf::st_coordinates(corrected_xy)[, 2]
    data <- data.frame(X = out_x, Y = out_y, Z = out_z)

    # Write out etp
    dir.create(file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,"epydro_las",fsep = .Platform$file.sep))
    cloud <- lidR::LAS(data)
    # lidR::las_check(cloud)
    lidR::writeLAS(cloud, file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,"epydro_las","epydro_las.laz",fsep = .Platform$file.sep))

    data.table::fwrite(data,file.path(ehydro_root,"survey",epydro_tracker_sourcedata_name,"epydro_xyz.csv",fsep = .Platform$file.sep))

    # Write out tracking sheet
    epydro_tracker[epydro_tracker$sheet_index==process_index,]$vertical_datum <- in_datum
    epydro_tracker[epydro_tracker$sheet_index==process_index,]$vertical_datum_offset <- as.numeric(jsonRespParsed$t_z)
    epydro_tracker[epydro_tracker$sheet_index==process_index,]$pos_z_convention <- 'up'
    epydro_tracker[epydro_tracker$sheet_index==process_index,]$vertical_datum_conversion <- glue::glue('{t_v_frame_override} = {in_datum} * {elev_unit_norm} + {as.numeric(jsonRespParsed$t_z)}')
    epydro_tracker[epydro_tracker$sheet_index==process_index,]$epydro_entwine <- T
    unlink(file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep))
    sf::st_write(epydro_tracker,file.path(ehydro_root,"epydro_SurveyJob.fgb",fsep = .Platform$file.sep),append=FALSE)
  }

  # wrap up
  return(TRUE)
}
