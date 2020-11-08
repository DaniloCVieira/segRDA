#' Number of clusters
#'
#' Compute the best number of clusters according to different indexes in \code{NbClust::\link[NbClust]{NbClust}} and defines the winning number as the one with the highest number of votes.
#' @param data - data.frame
#' @param seed Integer. The random number generator.
#' @param group_method the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans". More details in \code{\link[NbClust:NbClust]{Nbclust}} help page.
#' @return The function retunrs a list containing (1) the best number of cluster (the number with the highest number of votes) and (2) the counting table of votes.
#' @importFrom NbClust NbClust
#' @importFrom purrr safely
#' @rdname getGroups
#' @export
#'
#' @examples
#' getGroups(nema_PCRBS)
getGroups<-function(data, seed=NULL, group_method="kmeans"){
  suppressWarnings({
    set.seed(seed)
    indexes = c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "sdindex", "sdbw")
    results_nb = list()
    safe_nb = safely(NbClust)
    times = list()
    for(i in 1:length(indexes) ){
      t = lubridate::now()
      nb = safe_nb(as.dist(data)
                   , distance = "euclidean", method = group_method , min.nc = 2
                   , max.nc = 15
                   , index = indexes[i])
      results_nb[[i]] = nb$result$Best.nc[1]}
    votes<-table(unlist(results_nb))
    groups=as.numeric(names(votes)[which.max(votes)])
    return(list(groups=groups, votes=votes))
  })

}
