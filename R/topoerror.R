#' Tophological error from SOM
#'
#' @description - In construction
#' @param m an objected of class \code{'kohonen'}, resulted from the functions \code{\link[kohonen:supersom]{som}} and \code{\link[kohonen]{supersom}}.
#' @return The tophological error
#' @import kohonen
#' @export
topoerror <- function(m, type = c("nodedist", "bmu")) {

  type = match.arg(type,c("nodedist", "bmu"))
  dmat <- switch(type,
                 nodedist = {
                   as.matrix(dist(m$codes[[1]]))
                 },
                 bmu = {

                   ## following code adapted from map.kohonen, may be put in a
                   ## separate auxiliary function lateron
                   nd <- nrow(m$data[[1]])
                   ncodes <- nrow(m$codes[[1]])
                   np <- ncol(m$codes[[1]])
                   distances <-  m$distances

                   matrix(distances, nd, ncodes, byrow = TRUE)
                 })

  ## here we assume that there are no cases where more than one zero
  ## exists, something that would be _very_ strange
  dmat.ordered <- t(apply(dmat, 1, order))

  dists.2D <- unit.distances(m$grid, m$grid$toroidal)
  mean(mapply(function(i, j) dists.2D[i,j],
              dmat.ordered[,1], dmat.ordered[,2]))
}
