#' Functions to spatialize clusters resulted from  analysis
#'
#' Spatialize the clusters resulted from the \code{\link{wrapSOM}} analysis
#' @param somC  Object resulted from \code{\link{wrapSOM}} function.
#' @param coords Data.frame containing the Spatial coordinates
#' @param spd An object of class \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}}. It will be used to interpolate results.
#' @param layer1 Optional, polygon. An object of class \code{\link[sp:SpatialPolygons]{SpatialPolygons}}
#' @param layer2 Optional, points. An object of class \code{\link[sp:SpatialPolygons]{SpatialPointsDataFrame}}
#' @param layer3 Optional, points. An object of class \code{\link[sp:SpatialPolygons]{SpatialPointsDataFrame}}
#' @param res Integer. Interpolation resolution
#' @param main Optional, string. Overhall title of the plot
#' @param leg Optional, string. Legend title.
#' @param crs.info string with CRS information. More details in \code{sp::\link[sp:CRS-class]{CRS}}
#' @return A graphical output
#' @importFrom sp spTransform spsample coordinates gridded fullgrid proj4string sp.polygons sp.polygons
#' @importFrom gstat idw
#' @importFrom raster raster mask
#' @importFrom rasterVis GrTheme levelplot
#' @importFrom grid grid.rect grid.text
#' @rdname mapClust
#' @export
#'
#' @examples
#' somC<-wrapSOM(nema_PCRBS,dist=c("BrayCurtis"),groups=NULL, seed=1,n_iterations = 500, xdim=5, ydim=5)
#' mapClust(somC,coords=coords_PCRBS ,spd=spd_PCRBS,crs.info="+proj=longlat")
mapClust<-function(somC,coords,spd, layer1=NULL, layer2=NULL,layer3=NULL,res=20000,main="",hlegend=.5,  crs.info="+proj=longlat", leg="")
{
  suppressWarnings({
  prev=as.factor(somC[[1]])
  colhabs=somC$colhabs
  require(gstat)
  require(rasterVis)
  data_spat<-to_spatial(data.frame(coords,hab=as.numeric(prev)), crs.info=crs.info)
  shape=spTransform(spd,CRS=raster::crs(data_spat))
  if(is.null(layer1)==F){layer1=spTransform(layer1,CRS=raster::crs(data_spat))}
  if(is.null(layer2)==F){layer2=spTransform(layer2,CRS=raster::crs(data_spat))}
  if(is.null(layer3)==F){layer3=spTransform(layer3,CRS=raster::crs(data_spat))}
  lista<-list(shape,layer1,layer2, layer3)
  lista[which(unlist(lapply(lista,is.null)))]<-NULL
  limits<-do.call(rbind,  lapply(lapply(lista, bbox),t))


  xlimits<-range(limits[,1])
  ylimits<-range(limits[,2])
  colnames(limits)<-colnames(coords)
  newgrid<- as.data.frame(spsample(to_spatial(rbind(coords, limits), crs.info = crs.info), "regular", n=res))
  names(newgrid)       <- c("X", "Y")
  coordinates(newgrid) <- c("X", "Y")
  gridded(newgrid)     <- TRUE  # Create SpatialPixel object
  fullgrid(newgrid)    <- TRUE  # Create Spatialnewgrid object
  proj4string(newgrid) <- proj4string(to_spatial(coords, crs.info = crs.info))



  df_idw = idw(hab~1, data_spat, newgrid,nmax=1)
  df_idwraster <- raster::raster(df_idw)
  df_raster<-  raster::mask(df_idwraster, spd)
  my_rst<-  raster::ratify(df_raster)
  my_rst@data@attributes[[1]]$label <- levels((prev))
  names(colhabs)<-levels(prev)


  myTheme <- GrTheme ()
  myTheme$panel.background$col = 'white'
  colorkey=list(title=leg,space='right')

  levelplot(my_rst, par.settings = myTheme,col.regions=colhabs,margin=F,main=main, colorkey=list(labels = list(at = as.numeric(levels(prev)), labels = paste("Habitat",  levels(prev))),height=hlegend), xlim=xlimits, ylim=ylimits) +
    layer(sp.polygons(layer1, fill='gray80'), data=list(layer1=layer1)) +
    layer(sp.points(layer2,col="gray20",cex=.8), data=list(layer2=layer2))+
    layer(sp.points(layer3,pch=16,col="white", cex=1.5), data=list(layer3=layer3))+
    layer({
      xs <- seq(181000, 181400, by=100)
      grid.rect(x=xs, y=330500,
                width=100, height=30,
                gp=gpar(fill=rep(c('transparent', 'black'), 2)),
                default.units='native')
      grid.text(x= xs - 50, y=330560, seq(0, 400, by=100),
                gp=gpar(cex=0.5), rot=30,
                default.units='native')
    })
})
}


################################################
########    Iiternal functions
#' @export
to_spatial<-function(coords,  crs.info)
{
  colnames(coords)[1:2]<-c("Long","Lat")
  coordinates(coords)<-~Long+Lat
  proj4string(coords) <-CRS(crs.info)
  return(coords)
}

