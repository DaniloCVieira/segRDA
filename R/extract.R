#' Extract results and breakpoints from a \code{smw} object
#'
#' Functions to extract results and breakpoints from a \code{smw} object
#'
#' @param smw An object of class \code{smw} resulted from the \code{\link{SMW}} analysis.
#' @param dp An object of class \code{dp} (see Details).
#' @param w  Numeric. A "target" window size from which results will be extracted (see Details). Only effective if the \code{smw} object contains results from multiple window sizes.
#' @param index The result to be extracted:
#' \itemize{
#'   \item \code{"dp"}: The dissimilarity profile (DP) table containing significant discontinuities and suggested breakpoints.  The returned DP has class \code{"dp"}, with an own generic \code{print} method. The \code{print} command prints and returns invisibly the DP.
#'   \item \code{"rdp"}: data frame containing the randomized DP;
#'   \item \code{"md"}: mean dissimilarity of the randomized DP;
#'   \item \code{"sd"}: standard deviation for each sample position;
#'   \item \code{"oem"}: overall expected mean dissimilarity;
#'   \item \code{"osd"}: average standard deviation for the dissimilarities;
#'   \item \code{"params"}: list with input arguments
#' }
#' @param sig Significance test for detecting dissimilarity values that differs significantly from those appearing in a random pattern. The following tests are considered with default to \code{sig.test="z"}:
#' \itemize{
#'   \item \code{"dp"} consider normalized dissimilarity (z-scores) discontinuities that exceed a "z" critical value;
#'   \item \code{"sd"} consider dissimilarity discontinuities that exceed mean plus one standard deviation;
#'   \item \code{"SD2"} consider dissimilarity discontinuities that exceed mean plus two standard deviation;
#'   \item \code{"tail1"} Consider dissimilarity discontinuities that exceed 95 percent confidence limits.
#' }
#' @param z The critical value for the significance of z-values. Defaults to \code{'z=1.85'} (Erdös et.al, 2014).
#' @param BPs Defines if the breakpoints should be chosen as those sample positions corresponding to the maximum dissimilarity in a sequence of significant values (\code{"max"}) or as those sample positions corresponding to the median position of the sequence (\code{"median"}). Defaults to \code{BPs="max"}. If \code{NULL} the breakpoints are not computed.
#' @param seq.sig The maximum length of consecutive, significant values of dissimilarity that will be considered in defining the community breakpoints. Defaults to \code{seq.sig=3};
#' @return
#' \itemize{
#'   \item \code{"extract"} returns a result called by the argument \code{index} (see Details);
#'   \item \code{"bp"} returns the locations of the breakpoints.
#'   }
#' @details If the \code{smw} object contains results from multiple window sizes, the DP table will be based on the average Z-score over the set of analysed window sizes.  Available methods for class \code{"dp"} are \code{print}, \code{bp} and \code{plot}.
#'
#' The argument \code{w} is optional. If the \code{smw} object is length 1, \code{w} is ignored. If \code{length(smw)>1} and \code{w} is \code{NULL}, the function will extract the dissimilarity profile averaged over the set of window sizes.
#' @seealso \code{\link{plot.smw}}.
#' @aliases extract print.dp print.smw bp
#' @examples
#'\dontrun{
#' data(sim1)
#' sim1o<-OrdData(sim1$envi,sim1$comm)
#' ws20<-SMW(yo=sim1o$yo,ws=20, n.rand=99)
#' ws20_dp<-extract(ws20)
#' head(ws20)
#'}

#' @importFrom stats anova cor median na.omit pf pnorm sd
#' @importFrom utils head capture.output setTxtProgressBar txtProgressBar

#' @usage NULL
#' @export
#' @rdname extract
extract<-function (smw,w, index,sig, z, BPs, seq.sig)(UseMethod("extract"))

#' @method extract smw
#' @export
#' @rdname extract
extract.smw<-function(smw,w=NULL, index="dp",sig="z", z=1.85, BPs="max", seq.sig=3) {

  index<-match.arg(index,c('dp','rdp','md','sd','oem','osd','params'))
  check(smw=smw,w=w, index=index)
  sig<-match.arg(sig,c("z","sd","sd2",'tail1'))
  if(!is.null(BPs)){BPs<-match.arg(BPs,c("max","median"))}
  if(length(smw)>1&is.null(w)&sig!="z") stop(paste("The object '",substitute(smw), "' cointains results from multiple window sizes.\n The statistical test '",substitute(sig) ,"'  does not apply when analyzing multiple window sizes.\n Please specify argument 'w' to implement the statistical tests 'sd', 'sd2' and 'tail1' '", sep=""))
  if(!is.null(BPs)&index=="dp"){
    extract01<-BP(smw=smw,w=w, sig=sig,z=z, index=index,seq.sig=seq.sig, peaks.choice=BPs)
    class(extract01)<-c("dp", "data.frame")}
  if(is.null(BPs)|index!="dp")
  {extract01<-extractROOT(smw=smw,w=w, index=index,sig=sig, z=z) }
  if(index=="dp"){ class(extract01)<-c("dp", "data.frame")} else{
    attr(extract01,"params")<-NULL}
  return((extract01))}


#' @export
print.dp<-function(x,...)
{attr(x,"params")<-NULL
attr(x,"class")<-NULL
class(x)<-'data.frame'
print(x)
return(invisible(x))}


#' @export
bp<-function(dp){UseMethod("bp")}

#' @rdname extract
#' @export
#' @method bp dp
bp.dp<-function(dp)
{return(dp$positions[dp$bp!="-"])}


##################################################################################
#### Internal function
####
dpr<-function(...){UseMethod("dpr")}
print.dpr<-function(dp){
  x<-dp[[1]]
  pa<-attr(dp,'params')
  attr(x,"params")<-pa
  return(x)}
check<-function(smw,w, index)
{
  message(paste("Selected index: '", index,"'", sep=""))

  if(length("w")>1){stop( "w must be a single value", call. = FALSE)}
  ws<-names(smw)[!is.na(suppressWarnings( as.numeric(gsub("w","",names(smw)))))]
  ws.names<-  paste(ws, collapse=", ")
  ws.call<-paste(gsub("w","",ws), collapse=", ")
  m0<- paste(length(ws)," window size(s) available in '", substitute(smw), "': ",ws.call, sep="")
  m1<-paste("Window size selected:", w)
  m2<-paste("Window size selected:", ws.call)
  m3<-paste("Window sizes averaged: ", ws.call, sep="")
  if(!is.null(w)&sum(names(smw)%in%paste('w',w,sep=""))==0){ stop(m0, call. = FALSE)}
  if(length(smw)>1&is.null(w)){return(message(m3))}
  message(m0)
  if(length(smw)==1&is.null(w)){return(message(m2))}
  if(length(smw)==1|length(smw)>1&!is.null(w)){ return(message(m1))}
}
sigSMW<-function (smw,  w,sig, z)
{

  yo<-smw[[1]]$params$yo
  if(length(smw)==1){wn=1}
  if(length(smw)>1){ if(!is.null(w)){wn<-which(names(smw)==paste("w",w,sep=""))}}
  if(length(smw)==1|length(smw)>1&!is.null(w))
  {
    smw.w<-smw[[wn]]
    DPtable<-smw.w$dp
    OB<-DPtable[,3]
    md<-smw.w$md
    SD<-smw.w$sd
    Sigif.Tests<-DPtable[,1:2]
    Sigif.Tests$sd<- OB>md+(SD)
    Sigif.Tests$sd2<-  OB>md+(2*SD)
    Dz<-(OB-smw.w$oem)/smw.w$osd
    Sigif.Tests$z<-Dz>z
    p.value<-NULL
    for(i in 1:length(OB))
    {p.value[i]<-pnorm(OB[i],mean=md[i], sd=SD[i], lower.tail = F)}
    Sigif.Tests$tail1<-p.value<0.05
    Sigif.Tests$p.value<-round(p.value,10)
    DPtable$sig<-  Sigif.Tests[,which(colnames(Sigif.Tests)==sig)]
    DPtable$sig[which(DPtable$sig)]<-"*"
    DPtable$sig[which(DPtable$sig==F)]<-"ns"
    colnames(DPtable)[ colnames(DPtable)=="sig"]<-paste("sig:",sig,sep="")
    if(sig=="tail1"){DPtable$p.value<-round(p.value,10)}

    return(invisible(DPtable))
  }
  if(length(smw)>1&is.null(w))
  {
    mess<-paste("'pooled' \n","The object '",substitute(smw), "' contains results from different window sizes. \n ","Significance tests were based on the mean Z-scores of ", length(smw), " window sizes: \n ", paste(gsub("w","",names(smw)),collapse=", "),sep="")
    smw.tab<-data.frame(positions=(smw[[1]]$params$ws[1]/2):(nrow(yo)-(smw[[1]]$params$ws[1]/2)))
    smw.tab$SiteID<-rownames(yo)[smw.tab[,1]]
    meansZ<-z.scores<- diss<-smw.tab
    for( i in 1:length(smw))
    { z.scores[which(smw.tab[,1]%in%smw[[i]]$dp[,1]),i+2]<- smw[[i]]$dp$zscore
    diss[which(diss[,1]%in%smw[[i]]$dp[,1]),i+2]<-smw[[i]]$dp$diss}
    meansZ$Average_diss<-apply(diss[,-(1:2)],1,mean,na.rm=T)
    meansZ$Average_z<-apply(z.scores[,-(1:2)],1,mean,na.rm=T)
    meansZ$sig<-meansZ$Average_z>z
    meansZ$sig[which(meansZ$sig)]<-"*"
    meansZ$sig[which(meansZ$sig==F)]<-"ns"
    colnames(meansZ)[ colnames(meansZ)=="sig"]<-paste("sig:","Z",sep="")
    return(invisible(meansZ))

  }
}
extractROOT<-function(...)(UseMethod("extractROOT"))
extractROOT.smw<-function(smw,w, index,sig, z)
{


  for(j in 1:length(smw)) {smw[[j]]$dp<-sigSMW(smw[j],w=NULL ,sig=sig, z=z)}
  if(length(smw)>1&is.null(w)&index!="dp") {stop(paste("w must be selected when index='",index,"'", sep=""))}
  if(length(smw)==1)
  {
    result<-list(smw[[1]][[index]])
    attr(result, "params")<-c(smw[[1]]$params,z=z)
    class(result)<-"dpr"
    x<-print(result)
    class(x)<-c("smw","data.frame")
    return(invisible(x))} else {
      smw$pooled_DP<-sigSMW(smw, w=NULL,sig=sig, z=z)
      if(!is.null(w))
      {
        result<-list(smw[[which(names(smw)==paste("w",w,sep=""))]][[index]])
        attr(result, "params")<-c(smw[[1]]$params,z=z)
        class(result)<-"dpr"
        x<-print(result)
        class(x)<-c("smw","data.frame")
        if(index!="dp") {class(x)<-NULL}

        return(invisible(x))}
      result<-list(smw$pooled_DP)
      attr(result, "params")<-c(smw[[1]]$params,z=z)
      class(result)<-"dpr"
      x<-print(result)
      class(x)<-c("smw","data.frame")
      return(invisible(x))}
}
BP<-function(smw,...) UseMethod("BP")
BP.smw<-function (smw,w, sig,z, index,seq.sig, peaks.choice)
{

  SMWs<-extractROOT(smw=smw,w=w, index=index,sig=sig, z=z)

  sigs<-SMWs[,5]
  if(sum(SMWs[,5]!="ns")==0)
  { message("DP without any significant dissimilarity value: no breakpoint could be determined.")
    SMWs$bp<-"-"
    return(invisible(SMWs))}
  pk<-peaks.choice
  SMWs[,5]<-SMWs[,5]!="ns"
  re<-split(SMWs, cumsum(c(1, diff(SMWs$sig) !=0)))
  groups<-re[unlist(lapply(re,function(x) sum(x$sig)!=0))]
  max.mess<-max(unlist(lapply(groups[ unlist(lapply(groups,nrow))>=1],nrow)))
  summ_leng<-summary(unlist(lapply(groups[ unlist(lapply(groups,nrow))>=1],nrow)))

  groups<-groups[ unlist(lapply(groups,nrow))>=seq.sig]
  mess.seq<-paste( "Summary stats of consecutive, significant dissimilarity values for detecting breakpoints in '", substitute(smw), "': \n",sep="" )
  if(length(groups)==0&sum(SMWs[,5])>0)
  {
    warning("DP shows '", sum(SMWs[,5]),"' significant dissimilarity values but no breakpoint could be determined for \n '",substitute(seq.sig),"=", seq.sig,"' \n \n", mess.seq,paste0(capture.output(summ_leng), collapse = "\n") ,sep="",  call. = FALSE)
    SMWs[,5]<-sigs
    SMWs$bp<-"-"
    return(invisible(SMWs))}
  point<-NULL
  if(pk=="max")
  {point<-unlist(lapply(groups,function(x) x[which(x[,4]==max(x[,4])),1]))}
  if(pk=="median"){point<-round(unlist(lapply(lapply(lapply(groups, "[",1),unlist),median)))}
  names(point)<-NULL
  SMWs$bp<-"-"
  ro<-which(SMWs[,1]%in%point)
  SMWs$bp[ro]<-seq(1:length(point))
  menss0<-paste("Window sizes averaged:",": ", paste(gsub("w","",names(smw)[-length(smw)]), collapse = ", "), sep="" )
  mess<-paste("Number of breakpoints detected:", length(point))
  mess2<-paste("Breakpoint positions: ", paste(point, collapse=", "))
  message(mess)
  message(mess2)
  SMWs[,5]<-sigs
  return((SMWs))

}
