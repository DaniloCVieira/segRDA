#' Generate Report
#'
#' #' Generates a report for the analyses
#' @param som.graphics - the list of graphics resulted from the \code{\link{exploreSOM}} function
#' @param rf.graphics - the list of graphics resulted from the \code{\link{exploreRF}} function
#' @param output_dir - The output directory for the rendered output_file. If is missing the file will be saved in the working directory.
#' @return The function reders a HTML in the specified output directory
#' @export
ecoMLreport <- function( som.graphics, rf.graphics,output_file="ReportEcoML",  output_dir){
  file <- system.file("rmd", "template.Rmd", package = "ecoML")
  if (missing(output_dir)) {
    output_dir <- getwd()
  }
  rmarkdown::render(file, envir = list(som.graphics=som.graphics, rf.graphics=rf.graphics, firstletters=c(1,2)), output_dir = output_dir,output_file=output_file)
}

