#' Data ordering
#'
#' Ordinates both community and explanatory matrices based on the first RDA score.
#' @param x explanatory matrix;
#' @param y community matrix;
#' @param axis the RDA axis in which the ordering should be based. Defaults to \code{axis=1}.
#' @param ... parameters passed to vegan::\code{\link[vegan]{cca}};
#' @return An object of class \code{"ord"}, which is a list consisting of:
#' \enumerate{
#' \item{xo}: the ordered explanatory matrix
#' \item{yo}: the ordered community matrix
#' \item{x}: the original explanatory matrix
#' \item{y}: the original community matrix
#' }
#' @author Danilo Candido Vieira
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(x=sim1$envi, y=sim1$comm)
#' @export
#' @import vegan
OrdData<-function(x,y,axis=1,...)
{
  x<-as.matrix(x)
  y<-as.matrix(y)

  if(is.null(rownames(y))){ rownames(y)<-1:nrow(y)}
  if(is.null(colnames(y))){ colnames(y)<-1:ncol(y)}
  if(is.null(rownames(x))){ rownames(x)<-1:nrow(x)}
  if(is.null(colnames(x))){ colnames(x)<-1:ncol(x)}
  model<-rda(y~.,data.frame(x),...)
  site.order<-order(scores(model,axis,"sites"))
  sp.order<-order(scores(model,axis,"species"))
  yo<-y[site.order,sp.order]
  xo<-x[site.order,]
  output<-list(xo=xo,yo=yo, x=x, y=y)
  class(output)<-"ord"
  return(output)
}

