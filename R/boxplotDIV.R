#' Diversity indexes by cluster
#'
#' Generates graphical outputs
#' @param data  data.frame containing the variable to be plotted.
#' @param somC  object resulted from \code{\link{wrapSOM}} function.
#' @return A graphical Output
#' @author Danilo Candido Vieira
#' @examples
#' boxplotDIV(nema_PCRBS, somC)
#' @importFrom vegan diversity specnumber
#' @export
boxplotDIV<-function(data,somC,labels=c("Ind/10cm2","Nº gêneros","Shannon-Winner","Simpson","Jaccard"))
{
  graphics.off()

  Ddata<-divI(data)
  IndNames<-c("Densidade","Riqueza","Diversidade","Diversidade","Equitabilidade")


  prevfac<-as.factor(somC$somC)
  colsfac<-somC$colhabs

  boxes.plot<-list()
  for( i in 1:ncol(Ddata))
  {
    par(mar=c(5,5,3,2))
    boxplot(Ddata[,i]~prevfac, las=1,col=colsfac, xlab="Habitats", ylab=labels[i],  cex.lab=1.6,cex.axis=1.4, main=IndNames[i], cex.main=2)
    boxes.plot[[i]]<-recordPlot()
    graphics.off()
  }
  return(boxes.plot)


}





### Calcular índices de diversidade
divI<-function(abund){ S<-vegan::specnumber(abund)
N<-rowSums(abund)
H<-vegan::diversity(abund)
D<-vegan::diversity(abund, index="simpson")
J <- H/log(S)
return(data.frame(N,S,H,D,J))}

