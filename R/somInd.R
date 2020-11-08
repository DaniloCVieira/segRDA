#' Indicators from SOM
#'
#' @description  Plot the indicators variables along the kohonen cell map, with  the arrows indicating the direction and strength of the correlation between variables and cells.
#' @param m an object of class \code{'kohonen'}, resulted from the functions \code{\link[kohonen:supersom]{som}} and \code{\link[kohonen]{supersom}}.
#' @param indicate Calculation type. If \code{indicate=="cor"}, the function calculates the correlation of the variables with cells (weigthed by number of instances) and plot \code{npic} variables with the highest correlations, considering the different directions of the map. If \code{indicate=="var"}, the function calculates the variance of each variable weigthed by number of instances, and plot \code{npic} variables with the highest variances along the cell map. Variables with higher variances are those wich induced the strongest contrasts between the cells of the map.
#' @param somC object resulted from \code{\link{wrapSOM}} function.
#' @param npic Number of indicators variables to plot
#' @param pic  Integer or string indicating the variables (columns) to plot.
#' @param scaling Integer. Scaling factor for the arrows.
#' @return  The function invisibly returns the indicator variables, and their respective measure (settled by the argument \code{indicate})
#' @import kohonen
#' @importFrom graphics text rect strheight strwidth par
#' @export
#'
#' @examples
#' somC<-wrapSOM(nema_PCRBS,dist=c("BrayCurtis"),groups=NULL, seed=1,n_iterations = 500, xdim=5, ydim=5)
#' somInd(somC$som)
somInd<-function(m,npic=5, indicate=c("cor","var"), col.arrow="gray80", scaling=1.4)
{

  indicate=match.arg(indicate,c("cor","var"))
  colcodes<- coolBlueHotRed(nrow(somC$som.model$grid$pts))
  grid.size<-nrow(m$grid$pts)
  nb <- table(factor(m$unit.classif, levels=1:grid.size))
  CORMAP <- apply(do.call(cbind,m$codes),2,weighted.correlation,w=nb,grille=m)
  sigma2  <- sqrt(apply(do.call(cbind,m$codes),2,function(x,effectif){m2<-sum(effectif*(x- weighted.mean(x,effectif))^2)/(sum(effectif)-1)},effectif=nb))

  if(indicate=="cor")
  { indicadores<-names(sort(sigma2,decreasing=T))[1:npic]
  scores<-t(CORMAP)
  Xsp<-  scores[,1]
  Ysp<- scores[,2]
  A<- rowSums(abs(data.frame(S1=scores[names(which(Xsp<0&Ysp>0)),1],S2=scores[names(which(Xsp<0&Ysp>0)),2])))
  B<-rowSums(abs(data.frame(S1=scores[names(which(Xsp>0&Ysp>0)),1],S2=scores[names(which(Xsp>0&Ysp>0)),2])))
  C<-rowSums(abs(data.frame(S1=scores[names(which(Xsp<0&Ysp<0)),1],S2=scores[names(which(Xsp<0&Ysp<0)),2])))
  D<-rowSums(abs(data.frame(S1=scores[names(which(Xsp>0&Ysp<0)),1],S2=scores[names(which(Xsp>0&Ysp<0)),2])))
  Xsp1<-names(A)[order(A, decreasing = T)][1:npic]
  Xsp2<-names(B)[order(B, decreasing = T)][1:npic]
  Ysp1<-names(C)[order(C, decreasing = T)][1:npic]
  Ysp2<-names(D)[order(D, decreasing = T)][1:npic]
  indicadores<- unlist(lapply(data.frame(t(matrix(c(Xsp1,Xsp2,Ysp1,Ysp2),ncol=4))),cbind))
  indicadores= na.omit(indicadores[1:npic])
  result<-scores[indicadores,]
  colnames(result)<-c("cor.x", "cor.y")
  } else  { indicadores<-na.omit(names(sort(sigma2,decreasing=T))[1:npic])
    result<-data.frame(sigma2[indicadores])
    colnames(result)<-"Variance"


    }
  bp=data.frame(na.omit(as.matrix(t(CORMAP*scaling))))
  mul=1.15
  par(mar=c(7,6,2,1))
  plot(m,"codes",shape="straight", border="white", main="",keepMargins =T, bgcol=colcodes,codeRendering=F)
  text(m$grid$pts[,1],m$grid$pts[,2], labels=1:nrow(m$grid$pts))
  arrows(colMeans(m$grid$pts)[1],colMeans(m$grid$pts)[2],(bp[indicadores,1]*scaling)+colMeans(m$grid$pts)[1],(bp[indicadores,2]*scaling)+colMeans(m$grid$pts)[2], length=0.1, col=col.arrow, lwd=2)
  boxtext(x =(bp[indicadores,1]*scaling*mul)+colMeans(m$grid$pts)[1], y = (bp[indicadores,2]*mul*scaling)+colMeans(m$grid$pts)[2], labels = colnames(CORMAP[,indicadores]), col.bg = "white",  cex=.8, border.bg = col.arrow)
  #plotind(m,indicadores)

  return(invisible(result))

}



weighted.correlation <- function(v,w,grille)
{

  x <- grille$grid$pts[,"x"]
  y <- grille$grid$pts[,"y"]
  mx <- weighted.mean(x,w)
  my <- weighted.mean(y,w)
  mv <- weighted.mean(v,w)
  numx <- sum(w*(x-mx)*(v-mv))
  denomx <- sqrt(sum(w*(x-mx)^2))*sqrt(sum(w*(v-mv)^2))
  numy <- sum(w*(y-my)*(v-mv))
  denomy <- sqrt(sum(w*(y-my)^2))*sqrt(sum(w*(v-mv)^2)) #correlation for the two axes
  res <- c(numx/denomx,numy/denomy)
  return(res)
}


boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA,
                    border.bg = NA, adj = 1, pos = 4, offset = 0,
                    padding = c(0.5, 0.5), cex = 1, font = par('font')){

  ## The Character expansion factro to be used:
  theCex <- par('cex')*cex

  ## Is y provided:
  if (missing(y)) y <- x

  ## Recycle coords if necessary:
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }
  }
  strheight
  ## Width and height of text
  textHeight <- strheight(labels, cex = theCex, font = font)
  textWidth <- strwidth(labels, cex = theCex, font = font)

  ## Width of one character:
  charWidth <- strwidth("e", cex = theCex, font = font)

  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)
    }
  } else {
    adj <- c(0.5, 0.5)
  }

  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }
  } else {
    offsetVec <- c(0, 0)
  }

  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }

  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]

  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth
  rect(xleft = xMid - rectWidth/2,
                 ybottom = yMid - rectHeight/2,
                 xright = xMid + rectWidth/2,
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg)

  ## Place the text:
  text(xMid, yMid, labels, col = col.text, cex = theCex, font = font,
                 adj = c(0.5, 0.5))

  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }
}




