
#' Wrap for Self-organizing Maps
#'
#' @description  Wrap of functions to perform self-organizing maps analysis and determining the best number of clusters. The wrap uses the function \code{kohonen::\link[kohonen]{supersom}}
#' @param data numerical data.frame
#' @param dist distance function to be used for the data. Same as the argument \code{dist.fcts} from \code{kohonen::\link[kohonen:supersom]{som}} function, with the add that the wrap also implements the 'BrayCurtis' dissimilarity measure.
#' @param topo choose between a hexagonal or rectangular topology.
#' @param hmethod described in argument \code{method} from \code{\link[stats]{hclust}}
#' @param members described in argument \code{members} from \code{\link[stats]{hclust}}
#' @param groups Integer. The number of groups for clustering. If \code{groups=NULL} (default), the function will apply the function \code{\link{getGroups}} to determine the number of groups.
#' @param seed Integer. The random number generator
#' @param n_iterations Integer. The number of  interactions.
#' @param xdim number of cells in dimension X
#' @param ydim number of cells in dimension Y
#' @param scale Logical. Should the data be scaled?
#' @param ... Further arguments passed to \code{kohonen::\link[kohonen]{supersom}}
#' @return  Returns an list of length 7:
#'  \enumerate{
#'   \item \code{..$somC}: vector giving the cluster identification for each sample
#'   \item \code{..$colhabs}: the color attributed to each cluster
#'   \item \code{..$som.model}: object of class \code{kohonen}, resulted from \code{kohonen::\link[kohonen]{supersom}}
#'   \item \code{..$som.hc}: vector giving the cluster identification for each cell in the codebook
#'   \item \code{..$groups}: number of cluster defined/proposed
#'   \item \code{..$groups.result}: votes on the best number of clusters
#' }
#' @author Danilo Candido Vieira
#' @importFrom stats hclust cutree
#' @importFrom colorRamps matlab.like
#' @importFrom kohonen somgrid supersom classvec2classmat
#' @export
wrapSOM<-function(data,dist=c("BrayCurtis"),topo="hexagonal",hmethod='ward.D2',groups=NULL, seed=NULL,n_iterations = 1000,xdim=8,ydim=6,scale=F,members=NULL,...)
{

  set.seed(seed)
  numerics =  unlist(lapply(data,is.numeric))
  factors =  which(unlist(lapply(data,is.factor)))
  data_list = list()
  distances = vector()
  for (fac in factors)
    {
    data_list[[fac]] = classvec2classmat(data[[fac]] )
    rownames(data_list[[fac]])<-rownames(data)
    distances = c(distances, 'tanimoto')
    }

  if(scale==T){data_list[['numerics']] = scale(data[,numerics])}else{ data_list[['numerics']] = as.matrix(data[,numerics])}

  if( length(data_list[['numerics']])==0) {data_list[['numerics']]<-NULL} else{  distances = c( distances,dist)}

  { som_grid = somgrid(xdim = xdim, ydim=ydim, topo=topo, toroidal = F,neighbourhood.fct="gaussian")
    m = supersom( data_list, grid=som_grid, rlen= n_iterations, dist.fcts = distances, whatmap = c(factors, 'numerics')) }

  ## Clusters
  {
    dist_m = as.matrix(dist(do.call(cbind,m$codes)  ))
    dist_on_map = kohonen::unit.distances(som_grid)
    dist_adj =  dist_m ^ dist_on_map
    clust_adj = hclust(as.dist(dist_adj), method=hmethod,members = members)
    votesplot<-NULL
    groups.result<-getGroups(dist_adj)
    if(is.null(groups)){
      groups<-groups.result[[1]]
    } else { groups=groups}
    som_cluster_adj = cutree(clust_adj, groups)
    col_vector=colorRamps::matlab.like2(groups)
  }
  {
    newclass<-m$unit.classif
    for(i in 1:length(som_cluster_adj)) {newclass[newclass==i]<-rep(som_cluster_adj[i],sum(  newclass==i))}
    names(newclass)<-rownames(m$data[[1]])
  }
  somC=list(somC=newclass, colhabs=col_vector,som.model=m,som.hc=som_cluster_adj, groups=groups,groups.result=groups.result)
  class(somC)<-"somC"
  return(somC)
}

