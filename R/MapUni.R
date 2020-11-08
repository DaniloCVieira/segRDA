#' Functions to spatialize one variable from a dataset.
#'
#' @description  The function \code{mapOne} spatialize one  variable from a given dataset. The function `mapLoop` spatialize multiple variables.
#' @param data - data.frame containing the variable to be plotted.
#' @param coords data.frame containing the Spatial coordinates.
#' @param get - integer or a string containing the variable name to be plotted.
#' @param spd -  An object of clas \code{SpatialPolygonsDataFrame}. It will be used to interpolate results.
#' @param layer1 Optional, polygon. An object of clas \code{SpatialPolygons}
#' @param layer2 Optional, points.  An object of class \code{SpatialPointsDataFrame}
#' @param layer3 Optional, points.  An object of class \code{SpatialPointsDataFrame}
#' @param res Integer. Resolution of the interpolations
#' @param main Optinal, string. Title of the plot
#' @param leg Optinal, string. Legend title of the plot
#' @param crs.info string with CRS information. More details in \code{sp::\link[sp:CRS-class]{CRS}}
#' @param ... Arguments passed to \code{\link[ecoML:mapUni]{mapOne}} and \code{\link[gstat:krige]{idw}} functions.
#' @param pic  Integer or string indicating the variables (columns) to plot.
#' @return The function \code{mapOne} returns a list containing one single plot, whereas \code{mapLoop} returns a list containing multiple plots.
#' @author Danilo Candido Vieira
#' @import sp
#' @importFrom gstat idw
#' @importFrom raster raster mask
#' @importFrom rasterVis GrTheme levelplot
#' @importFrom grid grid.rect grid.text
#' @importFrom latticeExtra layer
#' @rdname mapUni
#' @aliases mapOne mapLoop
#' @export
#'
#'@examples
#' mapOne(nema_PCRBS,coords=coords_PCRBS,spd=spd_PCRBS)
#' mapLoop(nema_PCRBS, coords=coords_PCRBS, pic=1:5,spd=spd_PCRBS)
mapOne<-function(data,coords,get=1,spd,layer1=NULL, layer2=NULL, layer3=NULL, res=20000 ,main="",leg="",  crs.info="+proj=longlat",...)
{

  suppressWarnings({
    prev=unlist(data[get])
    colhabs=  colorRampPalette(c("blue", "green","yellow","orange","red"))
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
    data_idw = gstat::idw(hab~1, data_spat, newgrid,...)
    data_idwraster <- raster::raster(data_idw)
    my_rst<-  raster::mask(data_idwraster, spd)

    myTheme <- GrTheme ()
    myTheme$panel.background$col = 'white'
    colorkey=list(title=leg,space='right')

    levelplot(my_rst, par.settings = myTheme,main=main,margin=F,col.regions=colorRampPalette(c("blue", "green","yellow","orange","red"))(100),colorkey=colorkey, xlim=xlimits, ylim=ylimits) +
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

#' @rdname mapUni
#' @export
mapLoop<-function(data, coords, pic,...)
{
  Map_Inds<-list()
  newlist<-list()
  for(i in 1:length(pic)){
    Map_Inds[[i]]<-mapOne(data,coords=coords, get=pic[i], main=pic[i], leg=attr(data,"lab")[which(colnames(data)%in%pic[i])],...)
  }
  names(Map_Inds)<-pic
  for (i in 1:length(pic))
  { newlist[[i]]<- Map_Inds[which(names(Map_Inds)==pic[i])]}
  newlist<-unlist(newlist, recursive = F)
  return(newlist)
}




#### Internal Functions
to_spatial<-function(coords,  crs.info)
{
  colnames(coords)[1:2]<-c("Long","Lat")
  coordinates(coords)<-~Long+Lat
  proj4string(coords) <-CRS(crs.info)
  return(coords)
}
