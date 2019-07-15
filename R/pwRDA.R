#' Piecewise redundancy analysis (pwRDA)
#'
#' Perform a pwRDA using the specified breakpoints
#' @param x.ord ordered explanatory matrix
#' @param y.ord ordered community matrix
#' @param BPs community breakpoints
#' @param n.rand The number of randomizations for significance computation
#' @return  Returns an invisible list of length 4:
#'  \enumerate{
#'   \item \code{..$summ}: summary statistics of the pwRDA analysis;
#'   \item \code{..$rda.0}: full model cca object, which is described separately in vegan::\code{\link[vegan]{cca.object}}
#'   \item \code{..$rda.pw}: pw model cca object, which is described separately in vegan::\code{\link[vegan]{cca.object}}
#' }
#' @examples
#' \dontrun{
#' data(sim1)
#' sim1o<-OrdData(sim1$envi,sim1$comm)
#' pool<-SMW(sim1o$yo, ws=c(10,20,30), n.rand=10)
#' sim1.pw<-pwRDA(sim1o$xo,sim1o$yo, BPs=bp(extract(pool)))
#'}
#' @author Danilo Candido Vieira
#' @import vegan
#' @export
pwRDA<-function(x.ord,y.ord,BPs, n.rand=99)
{
  x.ord<-as.matrix(x.ord)
  y.ord<-as.matrix(y.ord)

  if(is.null(rownames(x.ord))){rownames(x.ord)<-1:nrow(x.ord)}
  if(is.null(rownames(y.ord))){rownames(y.ord)<-1:nrow(y.ord)}
  if(is.null(colnames(x.ord))){colnames(x.ord)<-1:ncol(x.ord)}
  if(is.null(colnames(y.ord))){colnames(y.ord)<-1:ncol(y.ord)}
  R.boot<-NULL
  pw.Models<-pwRDA.source(x.ord,y.ord,BPs)
  pw.obs<-pw.Models$summ
  obs<-pw.obs[2]
  rownames(x.ord)<-NULL
  rownames(y.ord)<-NULL
  pb <- txtProgressBar(min = 0, max = n.rand, style = 3)
  for( b in 1:n.rand)
  {
    sample<-sample(1:nrow(y.ord), replace=T)
    suppressWarnings(comm.rand<-y.ord[sample,])
    suppressWarnings(new.x<-x.ord[sample,])
    R.boot[b]<-pwRDA.source(new.x,comm.rand,BPs)$summ[2]
    setTxtProgressBar(pb, b) }
  p.value<-pnorm(obs,mean=mean(R.boot), sd=sd(R.boot), lower.tail = F)
  summ<-rbind(c(pw.obs[1],anova(pw.Models$rda.0)[1,4]),  c(pw.obs[2],p.value ),  c(pw.obs[3],pw.obs[4] ))
  summ<-round(summ,10)
  rownames(summ)<-c("FULL","PW","F")
  colnames(summ)<-c("Statistic","P.value")
  options(scipen=999)
  message("\n     pwRDA analysis       ",
          "\n ------------------------ ",
          "\n ** Summary statistics **",
          "\n ------------------------ \n",
          paste(capture.output(print(summ)), collapse = "\n"),
          "\n ------------------------\n ")
  pw.Models[[1]]<-summ
  class(pw.Models)<-"pw"
  return(invisible(pw.Models))
}

pwRDA.source <-function(x.ord,y.ord, BPs) {
  y.ord <- as.matrix(y.ord)
  x.ord <- as.matrix(x.ord)
  n <- nrow(x.ord)
  k <- ncol(x.ord)
  bks <- c(0, BPs,nrow(x.ord))
  nBPs <- length(bks)-1
  Xb <- matrix(0, ncol = nBPs * k, nrow = n)
  for (i in 1:(nBPs)){
    Xb[(bks[i]+1):bks[i+1], ((k*(i-1))+1):((k*(i-1))+k)] <- x.ord[(bks[i]+1):bks[i+1],]
  }
  Xb<-jitter(Xb)
  Xbc = scale(Xb, center = T, scale = F)
  rda.0<-vegan::rda ( data.frame( y.ord) ~ .,data.frame (x.ord))
  rda.pw<-vegan::rda ( data.frame( y.ord) ~ .,data.frame (Xb))

  Yc = scale(y.ord, center = T, scale = F)
  Y.avg = matrix(rep(apply(Yc,2,mean), times = nrow(Yc)), ncol = ncol(Yc), byrow = T)
  B.pw = solve(t(Xbc) %*% Xbc) %*% (t(Xbc) %*% Yc)
  # Species-environment correlations of PWRDA
  coord<-rda.pw$CCA$biplot
  bew.bp<-t(cor(coord,t(cor(x.ord,Xb))))
  rda.pw$CCA$biplot<-bew.bp
  Ypred.pw = Xbc%*%B.pw
  Yres <- Yc - Ypred.pw
  TSS.pw= sum((Yc)^2)
  RSS.pw= sum((Ypred.pw - Y.avg)^2)
  r2.pw<-RSS.pw/TSS.pw
  n.pw<-nrow(Xbc)
  k.pw<-ncol(Xbc)
  Radj.pw<-1-((1-r2.pw)*((n.pw-1)/(n.pw-k.pw-1)))
  Xc = scale(jitter(x.ord), center = T, scale = F)
  B.full = solve(t(Xc) %*% Xc) %*% (t(Xc) %*% Yc)
  Ypred.full = Xc%*%B.full
  Yres <- Yc - Ypred.full
  TSS.full= sum((Yc)^2)
  RSS.full= sum((Ypred.full - Y.avg)^2)
  r2.full<-RSS.full/TSS.full
  n.full<-nrow(Xc)
  k.full<-ncol(Xc)
  Radj.full<-1-((1-r2.full)*((n.full-1)/(n.full-k.full-1)))
  F.stat<-((RSS.full-RSS.pw)/(k.full-k.pw))/(RSS.full/(n.pw-k.full))
  dg1<- k.pw-k.full
  dg2<-n.pw-k.pw
  F.stat<-((RSS.pw-RSS.full)/(dg1))/(RSS.pw/(dg2))
  p.value<-1-pf(F.stat, dg1, dg2, lower.tail=T)
  summ<-c(Radj.full=Radj.full,Radj.pw=Radj.pw,F.stat=F.stat,p.value=p.value)
  pw<-list(summ=summ,rda.0=rda.0,rda.pw=rda.pw)
  class(pw)<-"pw"
  return(invisible(pw))
}

