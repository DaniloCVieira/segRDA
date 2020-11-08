#'  Wrap for Random Forest
#'
#' @description  Wrap of functions to perform Random Forest analysis with cross validation using the function \code{caret::\link[caret]{trainControl}}
#' @param data numerical data.frame
#' @param somC  object resulted from \code{\link{wrapSOM}} function
#' @param ntree Integer. Number of trees
#' @param seed Integer. The random number generator
#' @param prev.idw Logical. Should the supervisor use interpolated data?
#' @param trainControl.args List with arguments passed to \code{\link[caret:trainControl]{trainControl}}
#' @param crs.info string with CRS information. More details in \code{sp::\link[sp:CRS-class]{CRS}}
#' @return A list is returned of class \code{\link[caret]{train}}
#' @importFrom caret trainControl train
#' @importFrom gstat idw
#' @export
#'
#' @examples
#' rf=wrapRF(nema_PCRBS,somC)
wrapRF<-function(data,somC, ntree=500, seed=NULL,prev.idw=F,trainControl.args = list(method = 'repeatedcv', number = 10,repeats=20),coords=NULL,crs.info="+proj=longlat",...)
{
  suppressWarnings({
    controle_treinamento=do.call(trainControl,trainControl.args)

    if(is.null(seed==F)){set.seed(seed)}

    if(prev.idw==TRUE&is.null(coords)) stop("coords must be specified when prev.idw=T")
    if(prev.idw==TRUE){
      prev=somC$somC
      data_spat<-to_spatial(data.frame(coords[names(prev),],hab=as.numeric(prev)), crs.info=crs.info)
      data_idw = idw(hab~1, data_spat, to_spatial(coords[rownames(data),], crs.info=crs.info),nmax=1)
      data_idw$var1.pred<- round(data_idw$var1.pred)
      habpred<-data.frame(data_idw$var1.pred)

      if(isTRUE(prev.idw)){rownames(habpred)<-rownames(coords[rownames(data),])} else {rownames(habpred)<-rownames(data)}

      base<-na.omit(cbind(na.omit(round(habpred)),data[rownames(habpred),]))
      colnames(base)[1]<-'HAB'
    } else {
      base<-na.omit(cbind(HAB=as.factor(somC[[1]]),data[names(somC$somC),]))
    }

    lev<-names(table(base$HAB)>0)[which(table(base$HAB)>0)]
    base$HAB<-factor(base$HAB,levels=as.numeric(lev), labels = lev)
    modelo = train(HAB ~ ., data = base, trControl = controle_treinamento, method ="rf",ntree=ntree, localImp = TRUE,...)
    return(modelo)
  })


}


