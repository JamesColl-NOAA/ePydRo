#' @title Hello world
#' @description A hello world tester
#' @return a hello world message (polo!)
#' @family helper
#' @details Used to verify that ePydRo is successfully loaded.  Copy-paste-d from Google (https://ascii.co.uk/art/scuba).
#' @examples
#' ePydRo::marco()
#' @rdname marco
#' @export
marco <- function() {
  # sinew::moga(file.path(getwd(),"R/marco.R"),overwrite = TRUE)
  # devtools::document()
  # pkgdown::build_site(new_process=TRUE)
  # devtools::load_all()
  # devtools::check()

  ## -- Start --
  message("          O   ")
  message("        O     ")
  message("   ~   0  ]   ")
  message("     ___  ]   ")
  message("    (o o) ]   ")
  message("      C===|   ")
  message("##############")
  message("##  POLO!   ##")
  message("##############")

  return(TRUE)
}

