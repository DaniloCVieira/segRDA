#' Plot ....
#'
#' Plot results from \code{RFwrap} ...
#' @param rf  Random forest model resulted from \code{\link{wrapRF}}
#' @param somC  object resulted from \code{\link{wrapSOM}} function
#' @param data numerical data.frame
#' @param k Number of variables to plot
#' @return  The function returns ...
#' @details ...
#' @seealso ...
#' @aliases plot.RFw
#' @author Danilo Candido Vieira
#' @importFrom randomForestExplainer measure_importance min_depth_distribution plot_min_depth_distribution plot_multi_way_importance
#' @examples
#' data(sim1)
#' @rdname plot.RFw
#' @export
#'
#' @examples
#' rf=wrapRF(nema_PCRBS,somC)
#' exploreRF(rf, somC, nema_PCRBS)
exploreRF<-function(rf,somC, data, k=10,...)
{
  suppressWarnings({


    plotCM(rf, somC)
    CM.plot<-recordPlot()
    graphics.off()

    forest=rf$finalModel
    multi_imps = measure_importance(forest)
    indicadores<-as.character(multi_imps[order(multi_imps$mean_min_depth),][1:k,"variable"])
    min_depth_frame <- min_depth_distribution(forest)
    multi_imps<-multi_imps[order(multi_imps$mean_min_depth),]
    BoldImportance<-multi_imps[which(multi_imps$p_value<0.05),]
    red.labels <- rev(ifelse(as.character( multi_imps$variable[1:k]) %in%as.character(BoldImportance$variable), yes = "red", no = "black"))
    MDD.plot<-plot_min_depth_distribution(min_depth_frame, k=k)+theme(axis.text.y = element_text(colour  = red.labels))
    MWI1.plot<-plot_multi_way_importance(multi_imps, no_of_labels = k)
    MWI2.plot<-plot_multi_way_importance(multi_imps, x_measure = "accuracy_decrease", y_measure = "gini_decrease", size_measure = "p_value", no_of_labels = k)
    graphics.off()

    indicators<-as.character(BoldImportance$variable)

    ind_rf<-mapLoop(data=data, coords=coords_PCRBS[rownames(data),],pic=indicators,...)


    return(c(CM=list(CM.plot),RFexplainer=list(list(MDD.plot,MWI1.plot,MWI2.plot)),ind_rf=list(ind_rf), Sigs=list(BoldImportance)))

  })

}
