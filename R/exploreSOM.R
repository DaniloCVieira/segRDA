#' Exploring the SOMwrap results
#'
#' Creates a list of graphical outputs using a \code{\link{wrapSOM}} object
#' @param somC  object resulted from \code{\link{wrapSOM}} function.
#' @param coords data.frame containing the Spatial coordinates
#' @param spd An object of class \link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}.It will be used to interpolate results.
#' @param layer1 Optional, polygon. An object of class \code{\link[sp:SpatialPolygons]{SpatialPolygons}}
#' @param layer2 Optional, points. An object of class \code{\link[sp:SpatialPolygons]{SpatialPointsDataFrame}}
#' @param layer3 Optional, points. An object of class \code{\link[sp:SpatialPolygons]{SpatialPointsDataFrame}}
#' @param div Logical. Explore diversity indexes?
#' @param npic Integer. Number of indicator variables to plot
#' @param indicate Defaults to \code{indicate="cor"}. See \code{\link{somInd}}
#' @param scaling Integer. Scaling factor for the arrows when plotting \code{\link{somInd}}
#' @param res Integer. Interpolation resolution
#' @param main Optional, string. Overhall title of the plot
#' @param leg Optional, string.Title of the legend
#' @param crs.info string with CRS information. More details in \code{sp::\link[sp:CRS-class]{CRS}}
#' @return A graphical output
#' @return  Returns a list with graphical outputs from the \code{\link{wrapSOM}} analyses;
#' @aliases exploreSOM
#' @author Danilo Candido Vieira
#' @import kohonen
#' @rdname exploreSOM
#' @export
#'
#' @examples
#' exploreSOM(somC,coords_PCRBS,spd_PCRBS)
exploreSOM<-function(somC,coords,spd,layer1=NULL,layer2=NULL,layer3=NULL, div=T,npic=5,indicate="cor",scaling=1.4,res=2000,  crs.info="+proj=longlat",...)

{

  data=dataplot<-as.data.frame(do.call(cbind,somC$som.model$data))
  col_vector=somC$colhabs
  m=somC$som.model
  som_cluster_adj=somC$som.hc
  groups.result=somC$groups.result
  groups<-somC$groups

  # Fig. 1a - Train
  {

    plot(m$changes[,1], type="n", las=1,ylab="",xlab="", main="", mgp=c(3.2,1,0), cex.axis=1.2)
    for( i in 1:ncol(m$changes)) {lines(m$changes[,i], col= rainbow(ncol(m$changes))[i])}
    legend("topr", col=rainbow(ncol(m$changes)), lty=1, legend=colnames(m$changes), bty='n',cex=.7)
    legend("bottoml", legend=c(paste("topoerror:",round(topoerror(m),2)),paste("quant.error",round(mean(m$distances),2))), bty='n',cex=.7)
    mtext("Distancia media ate a unidade mais proxima",2, line=5,cex=1.2)
    mtext("Interacao",1, line=2.5,cex=1.2)
    changes.plot <- recordPlot()
    graphics.off()
  }

  # Fig. 1b - Count

  plot(m,"counts",shape="straight", border="white", main="",keepMargins =T)
  text(m$grid$pts[,1],m$grid$pts[,2], labels= table(factor(m$unit.classif, levels=1:(m$grid$xdim*m$grid$ydim))))
  counts.plot<-recordPlot()
  graphics.off()

  ###
  # Fig. 2a - Indicadores
  icell<-somInd(m, npic=npic, scaling=scaling, indicate=indicate)
  indicadores<-rownames(icell)
  icell<-recordPlot()
  graphics.off()


  # Fig. 2b - SOM interpolado

  bmu<-data.frame(SOM=somC$som.model$unit.classif)
  rownames(bmu)<- rownames(somC$som.model$data[[1]])
  colorkey=list(title="NeuID",   labels=list(labels=1:nrow( somC$som.model$grid$pts),at=(1:nrow( somC$som.model$grid$pts)),cex=.8), at=c(1:nrow( somC$som.model$grid$pts)-0.5,nrow( somC$som.model$grid$pts)+.5))
  SOM.map<-mapOne(data=bmu,coords=coords,spd=spd,layer1=layer1,layer2=layer2,layer3=layer3,leg="NeuID", crs.info = crs.info, res=res)
  graphics.off()


  # Fig. 3a

  Ind_som<-mapCells(somC$som.model,indicadores)
  Ind_map<-mapLoop(data,indicadores,coords=coords,spd=spd,layer1=layer1,layer2=layer2,layer3=layer3, indicadores, crs.info = crs.info, res=res,...)


  # Fig. 4 - Votos
  votesPlot(somC)
  votes.plot <- recordPlot()
  graphics.off()

  ###
  plot.new()
  par(mar=c(0,0,1,0))
  getResults(m,groups)
  info.plot<-recordPlot()
  graphics.off()
  ###

  plot(m, type="codes", main = "", bgcol = col_vector[som_cluster_adj], pchs = NA,border="white",shape="straight",codeRendering = "stars")
  add.cluster.boundaries(m, som_cluster_adj)
  codes.plot<- recordPlot()
  graphics.off()


  plot(m,"mapping",bgcol=col_vector[som_cluster_adj],col="black",pch=16, cex=1.2,shape="straight", border="white",main="")
  add.cluster.boundaries(m, som_cluster_adj)
  clust_cells<- recordPlot()
  graphics.off()


  ###
  habs.plot<-mapClust(somC,coords=coords, spd=spd,layer1=layer1,layer2=layer2,layer3=layer3,hlegend=1,res=res, crs.info=crs.info)

  graphics.off()
  ##
  if(isTRUE(div)){boxes.plot<-boxplotDIV(data,somC)} else{boxes.plot=NULL}
  return(c(list(info=info.plot,changes=changes.plot,counts=counts.plot,icell=icell),Ind_som=list(Ind_som),Ind_map=list(Ind_map),list(SOMmap=SOM.map,votes=votes.plot,codes=codes.plot,clust_cells=clust_cells,habs=habs.plot),list(boxes=boxes.plot)))
}

getResults<-function(m,groups)
{
  totnodes=nrow(do.call(cbind,m$codes))
  totclass<-length(table(m$unit.classif))
  nempty<-totnodes-totclass
  pempty<-round(nempty/totnodes*100,2)

  plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
  summ = c(paste("Grid:", paste0(m$grid$xdim,"x",m$grid$ydim)),
           paste("N camadas:",length(m$data)) ,
           paste("Distancia:",m$dist.fcts),
           paste("N amostras:",nrow(m$data[[1]])),
           paste("N agrupamentos:",groups),
           paste("Pesos:",m$distance.weights),
           paste("N Interacoes:",nrow(m$changes)),
           paste("Raio:",m$radius),
           paste0("Empty Nodes: ",nempty," (",pempty,"%)"))


  for (i in seq_along(summ)) {
    text(0, 100 - i*8, pos=4, summ[i], cex = 1, font=2)
  }
}
















