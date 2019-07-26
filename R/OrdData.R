#' Data ordering
#'
#' Ordinates both community and explanatory matrices based on the first RDA score.
#' @param x explanatory matrix;
#' @param y community matrix;
#' @param axis the RDA axis in which the ordering should be based. Defaults to \code{axis=1}
#' @param method standardization method (described in \code{\link[vegan]{decostand}}) to be applied on y before the RDA analysis. If \code{NA} (default), no transformation is performed.
#' @param ... furhter parameters passed to \code{\link[vegan]{decostand}} and \code{vegan::rda};
#' @return An object of class \code{"ord"}, which is a list consisting of:
#' \enumerate{
#' \item{xo}: the ordered explanatory matrix
#' \item{yo}: the ordered community matrix (non-transformed)
#' \item{x}: the original explanatory matrix
#' \item{y}: the original community matrix
#' }
#' @author Danilo Candido Vieira
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(x=sim1$envi, y=sim1$comm)
#' sim1.o<-OrdData(x=sim1$envi, y=sim1$comm, method="hellinger")
#' @export
#' @import vegan
#' @importFrom vegan decostand
#' @importFrom utils getS3method
OrdData<-function(x,y,axis=1,method=NA,...)
{ argg<-list(...)
  rda.args <- names(formals(getS3method("rda","default")))

  rda.args<-argg[names(argg)%in%rda.args]
  dec.args<-names(formals(decostand))
  dec.args<-argg[names(argg)%in%dec.args]
  x<-as.matrix(x)
  y<-as.matrix(y)
  if(!is.na(method))
  { yt<-do.call(decostand,  c(list(x=y, method=method),dec.args))} else {yt=y}
  if(is.null(rownames(y))){ rownames(y)<-1:nrow(y)}
  if(is.null(colnames(y))){ colnames(y)<-1:ncol(y)}
  if(is.null(rownames(x))){ rownames(x)<-1:nrow(x)}
  if(is.null(colnames(x))){ colnames(x)<-1:ncol(x)}
  model<-do.call(rda,  c(list(Y=data.frame(x), X=yt),rda.args))
  site.order<-order(scores(model,axis,"sites"))
  sp.order<-order(scores(model,axis,"species"))
  yo<-y[site.order,sp.order]
  xo<-x[site.order,]
  output<-list(xo=xo,yo=yo, x=x, y=y)
  class(output)<-"ord"
  return(output)
}


