## ----include = FALSE-----------------------------------------------------
library(segRDA)

## ----Install, eval=FALSE-------------------------------------------------
#  # to install the package via CRAN
#  install.packages("segRDA")
#  
#  # to install the package via GitHub
#  devtools::install_github("DaniloCVieira/segRDA")
#  

## ------------------------------------------------------------------------
data(sim3) ##Simulated data
x<-sim3$envi ## matrix of explanatory variables
y<-sim3$comm ## matrix of response variables


## ----Data ordering-------------------------------------------------------
sim1o<-OrdData(x=x,y=y, axis=1)
xo<-sim1o$xo ## ordered explanatory matrix
yo<-sim1o$yo ## ordered response matrix

## ---- fig.height=3, fig.width=6,fig.show='hold', tidy=TRUE, tidy.opts=list(width.cutoff=60)----

par(mfrow=c(1,2), mgp=c(1,1,0), cex=.9)
image(y, main="Original community data", col=topo.colors(100), axes=F, xlab="Sites", ylab="Species abundance")
image(yo, main="Ordered comunity data", col=topo.colors(100), axes=F, xlab="Sites")

## ----a_ws, results='hide',collapse = TRUE--------------------------------
ws20<-SMW(yo=yo,ws=20, n.rand=10)
pool<-SMW(y=yo,ws=c(10,20,30,40), n.rand=10)

## ----collapse = T--------------------------------------------------------
class(ws20)
length(ws20)
names(ws20)

class(pool)
length(pool)
names(pool)

## ----collapse = T--------------------------------------------------------
ws20_dp<-extract(ws20) 
ws20_dp[1:6,]

## ----collapse = T--------------------------------------------------------
pool_dp<-extract(pool) 
head(pool_dp)

## ----collapse = T--------------------------------------------------------
ws10_dp<-extract(pool, w=10) 

## ----collapse = T--------------------------------------------------------
ws20_dp<-extract(ws20, sig="tail1", seq.sig=20)

## ----collapse = T--------------------------------------------------------
bp(ws10_dp)
bp(pool_dp)

## ----collapse = T--------------------------------------------------------
extract(pool, w=10, index="osd")

## ----collapse=T,fig.height=3.5, fig.width=7,fig.show='hold'--------------
par(mfrow=c(1,2), cex=.9)
plot(pool,w=20, main="DP from a single window (10)", cex.main=.8) 
plot(pool, main="DP from pooled windows (10, 20, 30 and 40)",  bg=c("rainbow"),cex.main=.8)

## ----collapse=TRUE,fig.height=4, fig.width=4,fig.show='hold'-------------
plot(pool, w.effect = TRUE, main="Window size effect")

## ----include=F-----------------------------------------------------------
pw.sim<-pwRDA(x.ord=xo,y.ord=yo, BPs=bp(pool_dp))

## ----eval=F--------------------------------------------------------------
#  pw.sim<-pwRDA(x.ord=xo,y.ord=yo, BPs=bp(pool_dp))

## ----fig.height=3, fig.width=9, collapse=TRUE,fig.show='hold'------------
head(pw.sim$summ)
par(mfrow=c(1,3), cex=.75)
# plotting the full rda model:
plot(pw.sim$rda.0, main="full RDA model", las=1)
# plotting the DP profile and keeping the sample colors:
dp<-plot(pool, main="DP from pooled windows \n (10, 20, 30 and 40)",  bg='rainbow',cex.main=.8)

# plotting the pwRDA colored according to the breakpoints:
plot(pw.sim$rda.pw,type="n", scaling=3, main="pwRDA model")
points(pw.sim$rda.pw, pch=16, col=bgDP(dp), cex=1.2)
text(pw.sim$rda.pw,  display="bp",pch=16,col="steelblue4",lwd=2)

