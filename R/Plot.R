#' Plot the dissimilarity profiles
#'
#' Plot results from \code{smw} and \code{dp} objects. The command is a shortcut for extracting and plotting \code{\link{SMW}} resuts. Auxiliary arguments from \code{\link{extract}} (i.e.  \code{sig}, \code{z}, \code{BPs} and \code{seq.sig}) can be passed to \code{plot.smw}.  The auxilary method \code{bgDP} is available for the returned \code{dp} object when the argument \code{bg} is not NULL (see Details).
#' @param x  An object of class \code{"smw"} resulted from the function \code{\link{SMW}}.
#' @param dp An object of class \code{dp}.
#' @param w The window size from which results will be plotted. Only effective if \code{length(smw)>1}.
#' @param w.effect Logical, if \code{TRUE} draws a  dissimilarity profile using different windows sizes and returns an invisible data.frame with the breakpoint frequencies. Only effective if \code{length(smw)>1}. The function uses \code{\link{extract}} with defaults parameters to define the breakpoint positions for each of the evaluated window sizes.
#' @param values Character. \code{"zscore"} for plotting z-scores, \code{"diss"} for plotting dissimilarity values.
#' @param sig Significance test for detecting dissimilarity values that differs significantly from those appearing in a random pattern. If \code{NULL} the significance test is ommited from the plot.The following tests are considered with default to \code{sig.test="z"}:
#' \itemize{
#'   \item \code{'z'} consider normalized dissimilarity (z-scores) discontinuities that exceed a "z" critical value;
#'   \item \code{'sd'} consider dissimilarity discontinuities that exceed mean plus one standard deviation;
#'   \item \code{'SD2'} consider dissimilarity discontinuities that exceed mean plus two standard deviation;
#'   \item \code{'tail1'} Consider dissimilarity discontinuities that exceed 95 percent confidence limits.
#' }
#' @param z The critical value for the significance of z-values. Defaults to \code{'z=1.85'} (Erdös et.al, 2014).
#' @param BPs Defines if the breakpoints should be chosen as those sample positions corresponding to the maximum dissimilarity in a sequence of significant values (\code{"max"}) or as those sample positions corresponding to the median position of the sequence (\code{"median"}). Defaults to \code{BPs="max"}. If \code{NULL} the breakpoints are not computed.
#' @param seq.sig The maximum length of consecutive, significant values of dissimilarity that will be considered in defining the community breakpoints. Defaults to \code{seq.sig=3};
#' @param pchs A numerical vector of the form \code{c(d,s,b)} which modifies the default symbols of the plot. The default \code{pch = c(16,16,17)} describes respectively the dissimilarity values, significant dissimilarity values and breakpoints.
#' @param cols Vector of length 3 specifying the colors of the plot in the same way as the pch argument. Defaults to \code{colors=c("black","red","blue")}.
#' @param bg Optional. Sets background colors according to the breakpoints. It can be expressed either by a vector of colors or by the name of a pallet function (e.g. \code{"rainbow"}).
#' @param bg_alpha Factor modifying the opacity alpha of the backgroud [0,1].
#' @param wcols  Sets the colors for the window sizes (lines) when \code{w.effect=TRUE}. It can be expressed by a vector of colors or by the name of a pallet function. Defaults to \code{wcol="rainbow"} which uses the colour palette \code{rainbow} from R stats.
#' @return  The function returns invisibly an object of class \code{"dp"} (see Details).
#' @details If \code{bg} is not \code{NULL}, the attribute \code{params$bg} is added to the returned \code{dp}. This attribute contains the  sample colors used by the argument \code{bg}. The auxilary method \code{bgDP} can be used for accessing this color vector.
#' @seealso \code{\link{SMW}}, \code{\link{extract}}.
#' @param legend Logical. Should a default legend appear?
#' @param ... Further graphical parameters.
#' @aliases plot.smw bgDP
#' @author Danilo Candido Vieira
#' @importFrom grDevices adjustcolor
#' @importFrom graphics lines par plot plot.default points rect
#' @examples
#' data(sim1)
#' sim1o<-OrdData(sim1$envi,sim1$comm)
#' \dontshow{
#' pool<-SMW(yo=sim1o$yo,ws=c(40,50), n.rand=3)
#' plot(pool)
#' }
#' \donttest{
#' ws20<-SMW(yo=sim1o$yo,ws=20)
#' pool<-SMW(yo=sim1o$yo,ws=c(20,30,40))
#' plot(ws20)
#' plot(pool, w.effect=TRUE)
#' }
#' @rdname plot.smw
#' @export
plot.smw<-function(x, w=NULL, sig="z", z=1.85, BPs="max", seq.sig=3,w.effect=F,values=c("zscore","diss"), pchs=c(16,16,17),cols=c("black","red","blue"),bg=NULL,bg_alpha=0.1, wcols="rainbow", legend=TRUE,...){

   DP(input=x,w=w, sig=sig, z= z, BPs=BPs, seq.sig=seq.sig,w.effect=w.effect,values=values, pchs=pchs,cols=cols,bg=bg,bg_alpha=bg_alpha, wcols=wcols, legend=legend,... )}


####
DP<-function(input,w=NULL, sig="z", z=1.85, BPs="max", seq.sig=3,w.effect,values,pch=16, pchs,cols,bg,bg_alpha, wcols, legend,... )
{

  index="dp"
  sig.def<-sig
  BPs.def<-BPs
  peaks.choice<-match.arg(BPs,c("max","median"))
  argg<-list(...)
  class.in<-match.arg(class(input), c("smw", "dp"))
  values<-match.arg(values,c("zscore","diss"))

  dp<-dp.std<-suppressMessages(extract(smw=input,w=w, sig=sig,z=z, index="dp",seq.sig=seq.sig, BPs=peaks.choice))
  class(dp)<-'data.frame'
  if(isTRUE(w.effect)){
    if(length(input)==1){ stop("'smw' object must be length > 1 when 'w.effect=TRUE'. ")}
    ws<-as.numeric(gsub("w","",names(input)))
    breaks.list<-list(NULL)
    suppressWarnings(suppressMessages(for(i in 1:(length(input)))
    { breaks.list[[i]]<-suppressMessages(extract(smw=input,w=ws[i], sig=sig,z=z, index=index,seq.sig=seq.sig, BPs=peaks.choice))}))
    result.pooled<-input
    y.list<-lapply(result.pooled, function(x) x[[1]][,c("positions",values)])
    x<-breaks.list[[1]]
    for(i in 1:length(breaks.list)) {class(breaks.list[[i]])="data.frame"}
    bp.list<-lapply(breaks.list, function(x) x[,c("positions",values,"bp")])
    bp.list<-  lapply(bp.list, function(x) x[x$bp!="-",c(1,2)])
    y.limits<-range(na.omit(unlist(lapply(y.list, function(x) x[,2]))))}


  if(colnames(dp)[3]=="Average_diss"){colnames(dp)[c(3,4)]<-c("pooled diss", "pooled z-score")} else {
    colnames(dp)[c(3,4)]<-c("diss", "z-score")}

  {
    x<-dp[,1]
    if(values=="zscore"){y=dp[,4]} else {y=dp[,3]}
    y.sig<-dp[,5]
    x.null<-1:nrow(attr(dp, "params")$yo)
    y.null<-seq(min(y),max(y), length.out = length(x.null))
    x.sig<-x[dp[,5]!="ns"]
    y.sig<-y[dp[,5]!="ns"]
    if(w.effect==F){y.limits<-range(y.null)}
    plot.args <- c(names(formals(plot.default)), names(par()))
    plot.args<-argg[names(argg)%in%plot.args]
  }

  if(is.null(plot.args$ylab)){plot.args$ylab<-values}
  if(is.null(plot.args$xlab)){plot.args$xlab<-{xlab="sample positions"}}
  if(is.null(plot.args$las)){plot.args$las=1}
  if(!is.null(plot.args$type)){
    args.default<-plot.args
    args.default$type<-NULL
  } else{args.default=plot.args}
  x.limits<-range(x.null)
  do.call(plot,  c(list(x=x.null, y=y.null, xlim=x.limits,ylim=y.limits, type="n"),args.default))


  if(w.effect==T){
    if(length(wcols)==1){wcols<-get(get("wcols"))(length(ws))}
    for( i in 1:(length(ws)))
    {
      do.call(lines, c(list(y.list[[i]], col=wcols[i]),args.default))
      do.call(points, c(list(bp.list[[i]], pch=pch,col=wcols[i]),args.default))}
    if(legend==T){legend("topright" , col=wcols,lty=1,legend=names(input),y.intersp=.75, bty="n")}
    return(invisible(dp.std))}


  if(!is.numeric(BPs)) {
    x.bp<-x[dp$bp!="-"]
    y.bp<-y[dp$bp!="-"]
  }
  if(is.numeric(BPs))  {
    if(sum(BPs%in%x)!=length(BPs)) {
      stop("Mismatch between the settled BPs and the available sample positions \n check the dissimilarity profile for correct sample positions")}
    x.bp<-BPs
    y.bp<-y[which(x%in%x.bp)]}

  if(length(y.bp)>0) {
    if(!is.numeric(BPs)){
      if(!is.null(bg)) {
        if(length(bg)==1){bg.col<-get(get("bg"))(length(x.bp)+1)} else{bg.col=bg}
        cols.bp<-rep(bg.col, times=c(diff(c(0,x.bp,nrow(attr(dp.std, "params")$yo)))))
        lim <- par("usr")
        bc<-c(lim[1],x.bp,lim[2])

        bc<-c(lim[1],x.bp,lim[2])
        bg.adj<-adjustcolor(bg.col, alpha.f = bg_alpha)
        for(i in 1:(length(bc)-1)){rect(bc[i], bc[1]-100, bc[i+1], bc[2], border = "white", col = bg.adj[i])}
        rect(lim[1],lim[3],lim[2],lim[4])

      }
    } }

  opar <- par(new = par("new"))
  on.exit(par(opar))
  par(new=TRUE)
  do.call(plot,  c(list(x=x, y=y, xlim=x.limits, ylim=y.limits, pch=pchs[1], col=cols[1], ann=F, axes=F),plot.args))
  if(!is.null(sig.def)){do.call(points,  c(list(x=x.sig, y=y.sig, xlim=x.limits, ylim=y.limits, pch=pchs[2],col=cols[2]),args.default))}
  if(!is.null(BPs.def)){do.call(points,  c(list(x=x.bp, y=y.bp, xlim=x.limits, ylim=y.limits, pch=pchs[3],col=cols[3]),args.default))}
  if(legend==T)
  {sig.test<-colnames(dp)[5]
  if(sig.test=="sig:Z")
  {  sig.title<-paste("Significance test: \n",sig.test," >= ",z, sep="")} else {
    sig.title<-paste("Significance test:\n",sig.test)}
  legend("topleft",legend=sig.title,bty="n", cex=.8)}

  if(!is.null(bg)){

    attr(dp.std,"params")$bg<-cols.bp
    return(invisible(dp.std))}
  return(invisible(dp.std))


}

#' @export
bgDP<-function(dp){UseMethod("bgDP")}

#' @rdname plot.smw
#' @export
#' @method bgDP dp
bgDP.dp<-function(dp)
{return(attr(dp,'params')$bg)}
