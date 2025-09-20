#' Launch the richCluster Shiny app
#'
#' Starts the app shipped in inst/application.
#'
#' @return No return value.
#' @examplesIf interactive() && requireNamespace("shiny", quietly = TRUE)
#' # run_shiny_app()
#' @export
run_shiny_app <- function(...) {
  if (!interactive()) {
    stop("Shiny app can only be launched in an interactive session.")
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is not installed. Install it to launch the app.")
  }

  app_dir <- system.file("application", package = "richCluster")
  if (identical(app_dir, "")) {
    stop("App directory not found under inst/application in this package.")
  }

  # Optionally source modules from inst/application/modules
  mod_dir <- file.path(app_dir, "modules")
  if (dir.exists(mod_dir)) {
    rfiles <- list.files(mod_dir, pattern = "[.][Rr]$", full.names = TRUE)
    for (f in rfiles) sys.source(f, envir = parent.frame())
  }

  shiny::shinyAppDir(appDir = app_dir, ...)
}

# Backward compatible alias
#' @export
launch_shiny <- function(...) run_shiny_app(...)
