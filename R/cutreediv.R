#' @title Cut the tree
#' @description 
#' This function cuts the tree into several cluters by specifying the desired number of clusters. 
#' @param tree the divclust object
#' @param K an integer with the desired number of clusters.
#' @return \item{clusters}{the list of observations in each  cluster}
#' @return \item{description}{the monothetic description of each cluster}
#' @return \item{which_cluster}{a vector of integers indicating the cluster of each observation} 
#' @return \item{B}{the proportion of inertia explained by the partition (between-cluster inertia/total inertia)}
#' @return \item{leaves}{an internal list of \strong{leaves}}
#' @export
#' @examples
#' data(protein) # pure quantitatives data
#' tree <- divclust(protein) # full clustering
#' p_5 <- cutreediv(tree,K=5) # partition in 5 clusters
#' p_5
#' 
#' data(dogs) # pure qualitative data
#' tree <- divclust(dogs) # full clustering
#' p_4 <- cutreediv(tree,K=4) # partition in 4 clusters

#' data(wine) # mixed data 
#' data <- wine[,1:29]
#' tree <- divclust(data) # full clustering
#' p_4 <- cutreediv(tree, 4) # 
#' 

cutreediv <- function(tree,K) {
  cl <- match.call()
  cluster <- tree
  tree_kmax <- cluster$kmax
  if (!inherits(cluster, "divclust")) 
    stop("use only with \"divclust\" objects")
  if (!(K %in% 1:(length(cluster$clusters)))) stop(paste("K must be an integer between 0 and",length(cluster$clusters)))
  if (K != tree_kmax){
    l <- list()
    f <- fifo_create()
    node_init <- cluster$tree
    inert_init <- node_init$v$inert
    fifo_add(f, list(node = node_init, path = NULL, inert = inert_init, stage = 1))
    
    while ( !fifo_is_empty(f)) {
      z <- fifo_remove(f)
      node <- z$node
      path <- z$path
      inertie <- z$inert
      stages <- z$stage
      
      sentence <- make_path(node, path)
      path_l <- sentence$l
      path_r <- sentence$r
      
      inert_l <- node$l$v$inert
      inert_r <- node$r$v$inert
      
      if (node$l$v$inert <= cluster$height[K]) {   
        l[[length(l) + 1]] <-  list(class = node$v$A_l, path = path_l, inert = inertie, stage = stages) 
      } else {
        fifo_add(f, list(node = node$l, path = path_l, inert = append(inertie, inert_l), stage = stages + 1))
      }
      if (node$r$v$inert <= cluster$height[K]) {
        l[[length(l) + 1]]  <- list(class = node$v$A_l_c, path = path_r, inert = inertie, stage = stages)
      } else {
        fifo_add(f, list(node = node$r, path = path_r, inert = append(inertie, inert_r), stage = stages + 1))
      }  
    }
    c <- rep(0, length(l)) 
    k <- 1 
    for (i in l) { 
      for (j in i$class) {
        c[j] <- k 
      } 
      k <- k + 1 
    }
    
    part <- list()
    part$description <- make_description(l,cluster)
    names(part$description) <- paste("C", 1:length(l), sep = "")
    part$clusters <- lapply(l, function(x) {cluster$rnames[x$class]})
    names(part$clusters) <- paste("C", 1:length(l), sep = "")
    part$height <- lapply(l, function(x) {x$inert})
    names(part$height) <- paste("C", 1:length(l), sep = "")
    part$stages <- lapply(l, function(x) {x$stage})
    names(part$stages) <- paste("C", 1:length(l), sep = "")
    #names(part$height) <- lapply(l, function(x) {cluster$rnames[x$class]})
    part$T <- cluster$T
    part$which_cluster <- c
    part$leaves <- l
    class(part) <- "cutreediv"
    return(part)}
  
  else{
    l <- list()
    f <- fifo_create()
    node_init <- cluster$tree
    inert_init <- node_init$v$inert
    fifo_add(f, list(node = node_init, path = NULL, inert = inert_init, stage = 1))
    
    while ( !fifo_is_empty(f)) {
      z <- fifo_remove(f)
      node <- z$node
      path <- z$path
      inertie <- z$inert
      stages <- z$stage
      
      sentence <- make_path(node, path)
      path_l <- sentence$l
      path_r <- sentence$r
      
      inert_l <- node$l$v$inert
      inert_r <- node$r$v$inert
      
      if (node$l$v$inert <= 0) {   
        l[[length(l) + 1]] <-  list(class = node$v$A_l, path = path_l, inert = inertie, stage = stages) 
      } else {
        fifo_add(f, list(node = node$l, path = path_l, inert = append(inertie, inert_l), stage = stages + 1))
      }
      if (node$r$v$inert <= 0) {
        l[[length(l) + 1]]  <- list(class = node$v$A_l_c, path = path_r, inert = inertie, stage = stages)
      } else {
        fifo_add(f, list(node = node$r, path = path_r, inert = append(inertie, inert_r), stage = stages + 1))
      }  
    }
    c <- rep(0, length(l)) 
    k <- 1 
    for (i in l) { 
      for (j in i$class) {
        c[j] <- k 
      } 
      k <- k + 1 
    }
    
    part <- list()
    part$description <- make_description(l,cluster)
    names(part$description) <- paste("C", 1:length(l), sep = "")
    part$clusters <- lapply(l, function(x) {cluster$rnames[x$class]})
    names(part$clusters) <- paste("C", 1:length(l), sep = "")
    part$height <- lapply(l, function(x) {x$inert})
    names(part$height) <- paste("C", 1:length(l), sep = "")
    part$stages <- lapply(l, function(x) {x$stage})
    names(part$stages) <- paste("C", 1:length(l), sep = "")
    #names(part$height) <- lapply(l, function(x) {cluster$rnames[x$class]})
    part$T <- cluster$T
    part$which_cluster <- c
    part$leaves <- l
    class(part) <- "cutreediv"
    return(part)}
}