#' Confusion Matrix plot
#'
#' Plot the confusion matrix resulted from \code{RFwrap} ...
#' @param rf  Random forest model resulted from \code{\link{wrapRF}}
#' @param somC  object resulted from \code{\link{wrapSOM}} function
#' @param main title
#' @param cex.text size for text labels
#' @param line Integer for adjusting the title
#' @return A graphical ouptut
#' @importFrom caret confusionMatrix
#' @export
#'
#' @examples
#' rf=wrapRF(nema_PCRBS,somC)
#' plotCM(rf, somC)
plotCM<-function(rf,somC,main="",cex.text=1,line=2.8)
{

  ConfMat<-confusionMatrix(rf)
  CM<-CMerror<-ConfMat$table
  groups=somC$groups
  CM3<-as.data.frame(matrix(NA,nrow(CM),ncol(CM), dimnames=list(1:nrow(CM),1:ncol(CM))))
  CM3[rownames(CM),colnames(CM)]<-CM
  diag(CMerror)<-0

  CM3$class.error=(rowSums(CMerror)/rowSums(CM))
  CM<-CM3
  CM3<-as.data.frame(matrix(NA,groups,groups, dimnames=list(1:groups,1:groups)))
  CM3[rownames(CM),colnames(CM)]<-CM


  groups=nrow(CM)
  op <- par(no.readonly = TRUE)
  par(mar=c(5,5,7,7))
  CM3<-as.matrix(CM)
  CM3<-CM3[,nrow(CM3):1]
  bg.col <- somC$colhabs
  image(CM3, las=1,axes=F, col=rev(  gray.colors (100)))

  title(main, line=line+1.5)

  lim <- par("usr")

  x.bp<-seq(lim[1],lim[2], length.out = groups+1)
  bc <-x.bp
  y.bp<-rev(x.bp)

  bcy <- c(y.bp, lim[1])
  bg.adj <- adjustcolor(bg.col, alpha.f = 1)
  bc[groups+1]<-lim[4]
  for (i in 1:(length(bc) - 1)) {
    rect(bc[i], bcy[i], bc[i + 1], bcy[i + 1],
         border = NA, col = bg.adj[i])
  }

  factor<-abs(x.bp[1]-x.bp[2])/2
  seq<-x.bp+factor
  CM2<-  round(CM[nrow(CM):1,],2)
  CM2[,ncol(CM2)]<-paste0((CM2[,ncol(CM2)])*100,'%')
  par(xpd=T)

  for (x in 1:ncol(CM2)){
    for (y in nrow(CM2):1){
      points(seq[x], seq[y], pch=15, col= adjustcolor("white", alpha.f = .5), cex=4)
      text(seq[x], seq[y], CM2[y,x],cex=cex.text)
    }}


  abline(v=x.bp, xpd=F,col="white")
  abline(h=x.bp, xpd=F,col="white")
  i=1
  for (i in 1:(length(bc) - 1)) {
    rect(bc[i], lim[1], bc[i + 1], lim[1]-.05,
         border = NA, col = bg.adj[i], xpd=T)
  }

  for (i in 1:(length(bc) - 1)) {
    rect(lim[1], bcy[i], lim[1]-.05, bcy[i + 1],
         border = NA, col = bg.adj[i], xpd=T)
  }

  axis(2,at=seq(1,0,length.out = groups),(rownames(CM3)),las=1, tck=0, mgp=c(5,1.1,0))
  axis(3,at=seq,c(rownames(CM3),"class.error"),las=1,tck=0, mgp=c(5,1.1,0), xpd=T)
  accu<-(round(rf$results$Accuracy[which(rf$results$mtry%in%rf$bestTune)]*100,2))
  par(xpd=T)
  mtext(paste0("Acuracia: ",accu,"%"),3,line=line)

  legend("bottom", inset=-.2,legend=(ConfMat$text), bty="n", xpd=T,cex=.8,text.font=3 )

}

