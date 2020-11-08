
usethis::use_data(nema_PCRBS, overwrite = T)
usethis::use_data(envi_PCRBS, overwrite = T)
usethis::use_data(coords_PCRBS, overwrite = T)
usethis::use_data(spd_PCRBS, overwrite = T)
usethis::use_data(cont_PCRBS, overwrite = T)
usethis::use_data(points_PCRBS, overwrite = T)

usethis::use_data(envi_araca, overwrite = T)
usethis::use_data(nema_araca, overwrite = T)
usethis::use_data(coords_araca, overwrite = T)
usethis::use_data(spd_araca, overwrite = T)

#usethis::use_data(micro_PCRBS, overwrite = T)
#usethis::use_data(micro_PCRBS)

roxygen2::roxygenise()

rm(list=ls())
library(ecoML)
source(list.files(paste(path.package("ecoML"),"include/", sep = "/"), full.names = T))

nema_log<-vegan::decostand(nema_PCRBS,"log")
somC<-wrapSOM(nema_log,dist=c("BrayCurtis"),groups=NULL, seed=1,n_iterations = 1000, xdim=5, ydim=5)
som.graphics<-exploreSOM(somC,coords_PCRBS[rownames(nema_PCRBS),],spd_PCRBS)
save(som.graphics,file="som.graphics")
RF<-wrapRF(data=envi_PCRBS,somC=somC, prev.idw=T,seed=100,coords=coords_PCRBS, ntree=500)
rf.graphics<-exploreRF(RF,somC,envi_PCRBS,spd=spd_PCRBS, layer1=cont_PCRBS,layer2=points_PCRBS)
rf.graphics[1]
save(rf.graphics,file="rf.graphics")
ecoMLreport(som.graphics,rf.graphics)
library(randomForestExplainer)
require(usethis)
## 6. creates vignette folder
use_vignette("ecoML")

## 7. Build vignette in the doc folder and site to

devtools::build_vignettes()
pkgdown::build_reference()


## 8. Build manual to doc
library("qpdf")
devtools::build_manual()


## 9 Build readme - leave only the md file in the root directory
use_readme_md ('segRDA')
use_package_doc()


##10 Creates a Licence file
require(usethis)
use_mit_license("segRDA")
3
## 11 Release to the CRAN - checks
usethis::use_release_issue()
usethis::use_cran_comments()
devtools::check(remote=T)
devtools::check(run_dont_test=T)
devtools::check_win_devel()
rhub::check_for_cran()
rhub::check_on_ubuntu()
rhub::check_with_sanitizers()


usethis::use_version()
usethis::use_dev_version()
usethis::use_github_release()
devtools::submit_cran()
1


install.packages("segRDA")

# 12 After acceptance
usethis::use_github_release()
usethis::use_dev_version()
usethis::use_news_md()



# 11 - Push to git hub
devtools::install_github("hadley/devtools")
devtools::install_github("tinytex")

cat("GITHUB_PAT=97b305e603f549ef14a49a4ed355657df54672ab\n",
    file=file.path(normalizePath("~/"), ".Renviron"),
    append=TRUE)
Sys.getenv()
Sys.getenv("GITHUB_PAT")




git init
git add .
git config user.email "vieiradc@yahoo.com.br"
git commit -m "First commit package"
git remote add origin https://github.com/DaniloCVieira/segRDA.git
git push -f origin master
devtools::install_github("DaniloCVieira/segRDA")

git init
git add .
git config user.email "vieiradc@yahoo.com.br"
git commit -m "Development version"
git push -f origin master
