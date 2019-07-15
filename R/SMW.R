#' Split moving window analysis
#'
#' Function \code{SMW} performs split moving window analysis (SMW) with randomizations tests. It may compute dissimilarities for a single window size or  for several windows sizes.
#'
#' @param yo  The ordered community matrix.
#' @param ws The window sizes to be analyzed. Either a single value or a vector of values.
#' @param dist The dissimilarity index used in vegan::\code{\link[vegan]{vegdist}}. Defaults to \code{'bray'}.
#' @param n.rand The number of randomizations.
#' @param rand The type of randomization for significance computation (Erdös et.al, 2014):
#' \itemize{
#'   \item \code{"shift"}: restricted randomization in which data belonging to the same species are randomly shifted along the data series ("Random shift");
#'   \item \code{"plot"}: unrestricted randomization: each sample is randomly repositioned along the data series ("Random plot").
#' }
#' @return  A two-level list object (\code{class smw}) describing the SMW results for each window \code{w} analyzed. The \code{smw} object is of length \code{ws}, and each of the  \code{w} slots is a list of SMW results:
#' \itemize{
#'   \item \code{..$dp}: The raw dissimilarity profile (DP). The DP is a data frame giving the positions, labels, values of dissimilarity and z-scores for each sample;
#'   \item \code{..$rdp}: data frame containing the randomized DP;
#'   \item \code{..$md}: mean dissimilarity of the randomized DP;
#'   \item \code{..$sd}: standard deviation for each sample position;
#'   \item \code{..$oem}: overall expected mean dissimilarity;
#'   \item \code{..$osd}: average standard deviation for the dissimilarities;
#'   \item \code{..$params}: list with input arguments
#' }
#' Available methods for class \code{"smw"} are \code{print}, \code{extract} and \code{plot}.
#' @seealso \code{\link{plot.smw}}, \code{\link{extract}}.
#' @author Danilo Candido Vieira
#' @references
#' \itemize{
#'   \item Erdos, L., Z. Bátori, C. S. Tölgyesi, and L. Körmöczi. 2014. The moving split window (MSW) analysis in vegetation science - An overview. Applied Ecology and Environmental Research 12:787–805.
#'   \item Cornelius, J. M., and J. F. Reynolds. 1991. On Determining the Statistical Significance of Discontinuities with Ordered Ecological Data. Ecology 72:2057–2070.
#' }
#' @examples
#' \dontrun{
#' data(sim1)
#' sim1o<-OrdData(sim1$envi,sim1$comm)
#' ws20<-SMW(yo=sim1o$yo,ws=20, n=10)
#' head(print(ws20))
#' pool<-SMW(yo=sim1o$yo,ws=c(10,20,30), n.rand=10)
#' head(print(pool))
#' }
#' @export
#' @import vegan
SMW<-function(yo,ws, dist="bray",rand=c("shift","plot"), n.rand=99)
{
  if(n.rand<2){stop ("number of randomizations not alowed")}
  if(any(ws%%2==1)){stop("all Window sizes must be enven")}
  rand<-match.arg(rand,c("shift","plot"))
  argg <- c(as.list(environment()), list())
  smw<-list()
  for( j in 1:length(ws))
  {
    w1<-ws[j]
    message("\n SMW analysis (",j,"/", length(ws), '); w =',w1,"\n")
    DPtable<-smw.root(yo,w1,dist)
    OB<-DPtable[,3]
    rdp<-matrix(NA,length(OB),n.rand)
    pb <- txtProgressBar(min = 0, max = n.rand, style = 3)
    for( b in 1:n.rand)
    { if(rand=='shift')
    { comm.rand<-apply(yo,2,function(sp)sp[sample(1:length(sp))])
    rdp[,b]<-smw.root(data.frame(comm.rand),w1,dist)[,3]
    rdp[,b]<-as.matrix(rdp[,b])}
      if(rand=='plot')
      { comm.rand<-t(apply(yo,1,function(sp)sp[sample(1:length(sp))]))
      rdp[,b]<-smw.root(data.frame(comm.rand),w1,dist)[,3]}
      rdp[,b]<-as.matrix(rdp[,b])
      setTxtProgressBar(pb, b)}
    rownames(rdp)<-DPtable[,1]
    Dmean<-apply(rdp,1,mean)
    SD<-apply(rdp,1,sd)
    oem<-sum(Dmean)/(nrow(yo)-w1)
    osd<-sum(SD)/(nrow(yo)-w1)
    Dz<-(OB-oem)/osd
    DPtable$zscore<-Dz
    #significance
    smw[[j]]<-list(dp= data.frame(DPtable),rdp=matrix(rdp),md=Dmean,
                   sd=SD,oem=oem,osd=osd, params=argg)

    class(smw[[j]])<-c('smw')
  }
  names(smw)<-paste("w",ws,sep="")

  class(smw)<-c("smw")
  return(smw)
}



smw.root<-function(yo,w, dist)
{
  if(w%%2==1){stop('window size should be even')}
  diss<-NULL
  yo<-data.frame(yo)
  for(i in 1:(nrow(yo)-w+1)){
    wy.ord<-yo[i:(i+(w-1)),]

    half.a<-apply(wy.ord[1:(nrow(wy.ord)/2),],2,sum)
    half.b<-apply(wy.ord[-c(1:(nrow(wy.ord)/2)),],2,sum)
    d<-vegdist(rbind(half.a,half.b),dist)
    diss[i]<-d
    k<-(w/2)
    for(i in 1:((nrow(yo)-w)))
    {k[i+1]<-(w/2)+i}
  }
  result<-data.frame(positions=k,sampleID= rownames(yo)[k],diss=diss)
  return(invisible(result))
}

