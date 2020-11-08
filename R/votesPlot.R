#' Function to plot the best number of clusters
#'
#' @description  Barplot indicating the best number of clusters.
#' @param somC Object resulted from \code{SOMWrap} function.
#' @param col1 Color bars
#' @param col1 Color to highlight the number of clusters with the highest vote count.
#' @param main Title of the plot
#' @param xlab Label of the axis x
#' @param ylab Label of the axis y
#' @return A graphical output
#' @export
#'
#' @examples
#' somC<-SOMwrap(nema_PCRBS,dist=c("BrayCurtis"),groups=NULL, seed=1,n_iterations = 500, xdim=5, ydim=5)
#' votesPlot(somC)

votesPlot<-function(somC,col1="steelblue",col2="salmon",main='Votes on the best number of clusters',xlab='Best number of clusters',ylab="counting")

{
  colbars<-rep(col1,length(somC$groups.result[[2]]))
  colbars[which(names(somC$groups.result[[2]])==somC$groups)]<-col2
  barplot(somC$groups.result[[2]],las=1, main=main,xlab=xlab, ylab=ylab, col=colbars)
}

