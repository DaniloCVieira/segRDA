#' Plot the heatmap  throung a loop of variables

#' @description Plot the heatmap for multiple variables across the kohonen map.
#' @param  m an object of class \code{'kohonen'}, resulted from the functions \code{\link[kohonen:supersom]{som}} and \code{\link[kohonen]{supersom}}
#' @param pic  Integer or string indicating the variables (columns) to plot.
#' @return A graphical output
#' @import kohonen
#' @rdname mapCells
#' @export
#'
#' @examples
#' somC<-SOMwrap(nema_PCRBS,dist=c("BrayCurtis"),groups=NULL, seed=1,n_iterations = 500, xdim=5, ydim=5)
#' indicators<-rownames(somInd(somC$som))
#' mapCells(somC$som.model, pic=indicators)

mapCells<-function(m,pic)
{

  plotIndicators<-list()
  for (j in 1:length(pic)){
    plot(m,type="property",property=do.call(cbind,m$codes)[,pic][,j],palette.name=coolBlueHotRed,main=pic[j],cex=0.5,shape="straight" )
    plotIndicators[[j]]<-recordPlot()
    graphics.off()
  }
  names(plotIndicators)<-pic
  return(plotIndicators)
}

