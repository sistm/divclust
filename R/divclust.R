#' Monothetic divisive hierarchical clustering
#'
#' DIVCLUS-T is a divisive hierarchical clustering algorithm based on a monothetic bipartitional
#' approach allowing the dendrogram of the hierarchy to be read as a decision tree.
#' It is designed for numerical, categorical (ordered or not) or mixed data. Like the Ward agglomerative hierarchical
#' clustering algorithm and the k-means partitioning algorithm, it is based on the minimization of
#' the inertia criterion. However, it provides
#' a simple and natural monothetic interpretation of the clusters. Indeed, each cluster is decribed by set of
#' binary questions. The inertia criterion is calculated on all the principal components of PCAmix
#' (and then on standardized data in the numerical case).
#' @param data a data frame with numerical and/or categorical variables.
#' If the variable is ordinal, the column must be of class factor with the argument ordered=TRUE.
#' @param K the number of final clusters (leaves of the tree). By default, the complete dendrogram
#' is performed.
#' @param mtry the number of variables sampled at each tree node. By default, all the variables are selected.

#' @return \item{tree}{an internal \strong{tree}}
#' @return \item{clusters}{the list of observations in each final cluster (the leaves of the tree)}
#' @return \item{description}{the monothetic description of each final cluster (the leaves of the tree)}
#' @return \item{which_cluster}{a vector of integers indicating the final cluster of each observation}
#' @return \item{height}{the height of the clusters in the dendrogram of the tree}
#' @return \item{inertia}{the list of inertia nodes extracted from the dendrogram's height values, for each cluster}
#' @return \item{B}{the proportion of inertia explained by the final partition (between-cluster inertia/total inertia)}
#' @return \item{data_quanti}{the quantitative data set}
#' @return \item{data_quali}{the qualitative data set}
#' @return \item{mod_quali}{the list of categories of qualitative variables}
#' @return \item{vec_quali}{number of categories of each qualitative variable}
#' @return \item{kmax}{the number of different observations i.e. the maximal number of leaves}
#' @return \item{T}{The total inertia}
#' @return \item{mtry}{the number of variables sampled at each tree node}
#' @details The tree has K leaves corresponding to a partition in K clusters if K is specified
#' in input. Otherwise, each final cluster contains one observation and the tree is the
#' complete dendrogram. The between-cluster inertia of the final partition of the leaves is the sum of
#' the heights of the clusters in the tree. The total inertia for the quantitative dataset is equal
#' to p1 (the number of quantitative variables). The total inertia for the qualitative
#' dataset is m-p2 where m is the total number of categories and p2 is the number of qualitative
#' variables. For a mixture of quantitative and qualitative data, the total variance is
#' p1+m-p2. The quality of a partition is the proportion of inertia explained by the patition
#' which is the between-cluster inertia divided by the total inertia. The height of a cluster in the dendrogram  of divclust
#' is the inertia variation which is also the aggregation criterion of Ward used in ascendant
#' hierarchical clustering. This can be used, to help in the choice of the number of clusters
#' as for Ward hierarchical clustering.The inertia parameter provides, for each cluster, 
#' a list of inertia nodes extracted from the 'height' values of the dendrogram. 
#' This list captures the tree’s topology by marking the nodes that contribute to the cluster’s structure. 
#' To compute the distance between two clusters, you can sum the 'height' values from their last shared inertia 
#' node to the final node in the 'height' list, and then normalize this sum by dividing it by the total sum of all height values 
#'(representing the maximum inertia). For ordered qualitative variables (class factor with argument ordered
#' =TRUE), this order on the categories is used to reduce the number of possible binary questions.
#'
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @importFrom intervals Intervals interval_intersection
#' @useDynLib divclust, .registration = TRUE
#'
#' @seealso  \link{plot.divclust} \link{cutreediv}
#' @export
#' @examples
#' data(protein) # pure quantitatives data
#' tree <- divclust(protein) # full clustering
#' plot(tree)
#' plot(1:(tree$kmax-1),tree$height,xlab="number of cluster",ylab="height",main="Split levels")
#' c_5 <- divclust(protein, K=5) # stops clustering to 5 clusters
#' plot(c_5,nqbin=4)
#' c_5$B*100 #explained inertia
#' c_5$clusters  # retrieve the list of observations in each cluster
#' c_5$description # and their monothetic description
#'
#' data(dogs) # pure qualitative data
#' tree <- divclust(dogs) # full clustering
#' plot(tree)
#' plot(1:(tree$kmax-1),tree$height,xlab="number of cluster",ylab="height",main="Split levels")
#' c_4 <- divclust(dogs, K=4) # stops clustering to 4 clusters
#' plot(c_4)
#' c_4$clusters # retrieve the list of observations in each cluster
#' c_4$description # and their monothetic description
#' c_4$which_cluster # return a vector indicating to which cluster belongs each individual
#' c_4$B*100 #explained variance
#'
#' dogs2 <- dogs # take the order of categories into account (to reduce the complexity)
#' levels(dogs$Size)
#' size2 <- factor(dogs$Size,c("small","large","medium")) #changes the order of the levels
#' levels(size2)
#' dogs2$Size <- ordered(size2) #specify argument ordered=TRUE in the class factor
#' tree <- divclust(dogs2) # full clustering with variable Size considered as ordered.
#' plot(tree) #the constraint on the order changes the clustering


#' data(wine) # mixed data
#' data <- wine[,1:29]
#' c_tot <- divclust(data) # full clustering
#' plot(c_tot)
#' c_4 <- divclust(data, 4) # stops clustering to 4 clusters
#' plot(c_4)
#' p2 <- length(c_4$vec_quali)
#' p1 <- ncol(data)-p2
#' sum(c_4$height)/(p1+sum(c_4$vec_quali)-p2)*100 #explained variance
#' c_tot$tree$v #internal contain of the root node
#' c_tot$tree$r$v #internal contain of the right node of the root node

divclust <- function (data, K = NULL, mtry = ncol(data))
{
  cl <- match.call()
  if (is.null(data)) {
    stop("no problem to solve")
  }
  obj <- split_mix(data)
  data_quanti <- obj$data_quanti
  data_quali <- obj$data_quali

  # quanti preprocessing
  if (!is.null(data_quanti)) {
    mu_quanti <- colMeans(data_quanti)
    n <- nrow(data_quanti)
    nb_quanti <- ncol(data_quanti)
    rnames <- rownames(data_quanti)
    red <- sqrt((n-1)/n)
    sig_quanti <- apply(data_quanti,2, sd)*red
    #X_quanti <- as.matrix(data_quanti)
    X_quanti_init <- as.matrix(data_quanti)
    # BEWARE X_quanti is both centered and standardized
    X_quanti <- X_quanti_init - matrix(mu_quanti, n, nb_quanti, byrow=TRUE)
    X_quanti <- X_quanti /matrix(sig_quanti, n, nb_quanti, byrow=TRUE)
  } else {
    n <- nrow(data_quali)
    X_quanti <- X_quanti_init <-matrix(0,n,0)
    nb_quanti <- 0
    mu_quanti <- c()
    sig_quanti <- c()
  }
  # quali preprocessing
  if (!is.null(data_quali)) {
    if (sum(unlist(lapply(data_quali,is.character))) !=0)
      stop("All columns in data_quali must be of class factor")
    res <- format_cat_data(data_quali)
    rnames <- rownames(data_quali)
    mod_quali <- res$mod_quali
    vec_quali <- res$vec_quali
    cnames_quali <- res$cnames_quali
    X_quali <- as.matrix(res$Y)
    mu_quali <- colMeans(X_quali)
    nb_quali <- length(vec_quali)
    vec_order <- unlist(lapply(data_quali,function(x){class(x)[1]=="ordered"}))
  } else {
    X_quali <-matrix(0,n,0)
    mu_quali <- matrix(0,n,0)
    cnames_quali <- list()
    nb_quali <- 0
    mod_quali <- list()
    vec_quali <- c()
    vec_order <- c()
  }

  # preprocessing
  #X <- cbind(X_quanti, X_quali)
  X <- cbind(X_quanti_init, X_quali)
  p <- ncol(X)
  X_quali <- X_quali - matrix(mu_quali, n, length(mu_quali), byrow=TRUE)
  X_centre <- cbind(X_quanti, X_quali)
  maxdim <- min(n,ncol(X)-nb_quali)
  dec <- gsvd(X_centre, c(rep(1, nb_quanti), 1./mu_quali), rep(1./n, n))
  Z <- dec$u[,1:maxdim] %*% diag(dec$d[1:maxdim], nrow=maxdim, ncol=maxdim)
  #Z <- try(dec$u[,1:maxdim] %*% diag(dec$d)[1:maxdim, 1:maxdim, drop=FALSE])
  # if(inherits(Z, "try-error")){
  #   browser()
  # }
  w <- rep(1. / n, n)
  D <- rep(1, ncol(Z))

  # tree construction
  # algorithm's sketch :
  # start with one leaf consisting in best question chosen on all indices
  # leaves <- {Tree(choose_qestion(Z, all_indices))}
  # each iteration, compute "locally maximal partition"
  # leaves <- (leaves\\{Fmax}) U {Tree(A_l(F_max), Tree(Al_c(F_max))}
  
  #mtry à chaque noeud de l'arbre
  cut_var <- sample(1:ncol(data), size = mtry)
  #variable sample
  data_mtry <- data[, cut_var, drop=FALSE]
  obj_mtry <- split_mix(data_mtry)
  
  data_quanti_mtry <- obj_mtry$data_quanti
  if (!is.null(data_quanti_mtry)) {
    X_quanti_mtry <- as.matrix(data_quanti_mtry)
  } else {
    X_quanti_mtry <- matrix(0,n,0)
  }
  
  data_quali_mtry <- obj_mtry$data_quali
  if (!is.null(data_quali_mtry)) {
    res_mtry <- format_cat_data(data_quali_mtry)
    X_quali_mtry <- as.matrix(res_mtry$Y)
    vec_quali_mtry <- res_mtry$vec_quali
    vec_order_mtry <- unlist(lapply(data_quali_mtry,function(x){class(x)[1]=="ordered"}))
  } else {
    X_quali_mtry <-matrix(0,n,0)
    vec_quali_mtry <- c()
    vec_order_mtry <- c()
  }
  
  X_mtry <- cbind(X_quanti_mtry, X_quali_mtry)

  dendro <-new.env()
  dendro$v <- choose_question(X_mtry, Z, c(1:n), vec_quali_mtry, w, D,vec_order_mtry)
  leaves <- list(dendro)
  k <- 1
  kmax <- nrow(unique(X))
  if (is.null(K)) K <- kmax
  if (!(K %in% 2:kmax)) K <- kmax
  height <- rep(NA,K-1)

  while (k < K) {
    ind_max <- which.max(lapply(leaves, function(x){x$v$inert}))
    height[k] <- leaves[[ind_max]]$v$inert
    Fmin <- leaves[[ind_max]]
    leaves[[ind_max]] <- NULL
    
    cut_var_l <- sample(1:ncol(data), size = mtry)
    #variable sample
    data_mtry_l <- data[, cut_var_l, drop=FALSE]
    obj_mtry_l <- split_mix(data_mtry_l)
    
    data_quanti_mtry_l <- obj_mtry_l$data_quanti
    if (!is.null(data_quanti_mtry_l)) {
      X_quanti_mtry_l <- as.matrix(data_quanti_mtry_l)
    } else {
      X_quanti_mtry_l <- matrix(0,n,0)
    }
    
    data_quali_mtry_l <- obj_mtry_l$data_quali
    if (!is.null(data_quali_mtry_l)) {
      res_mtry_l <- format_cat_data(data_quali_mtry_l)
      X_quali_mtry_l <- as.matrix(res_mtry_l$Y)
      vec_quali_mtry_l <- res_mtry_l$vec_quali
      vec_order_mtry_l <- unlist(lapply(data_quali_mtry_l,function(x){class(x)[1]=="ordered"}))
    } else {
      X_quali_mtry_l <-matrix(0,n,0)
      vec_quali_mtry_l <- c()
      vec_order_mtry_l <- c()
    }
    
    X_mtry_l <- cbind(X_quanti_mtry_l, X_quali_mtry_l)

    cqg <- choose_question(X_mtry_l, Z, Fmin$v$A_l, vec_quali_mtry_l, w, D,vec_order_mtry_l)
    Fmin$l <- new.env()
    Fmin$l$v <- cqg
    leaves <- append(leaves, Fmin$l)
    
    cut_var_r <- sample(1:ncol(data), size = mtry)
    #variable sample
    data_mtry_r <- data[, cut_var_r, drop=FALSE]
    obj_mtry_r <- split_mix(data_mtry_r)
    
    data_quanti_mtry_r <- obj_mtry_r$data_quanti
    if (!is.null(data_quanti_mtry_r)) {
      X_quanti_mtry_r <- as.matrix(data_quanti_mtry_r)
    } else {
      X_quanti_mtry_r <- matrix(0,n,0)
    }
    
    data_quali_mtry_r <- obj_mtry_r$data_quali
    if (!is.null(data_quali_mtry_r)) {
      res_mtry_r <- format_cat_data(data_quali_mtry_r)
      X_quali_mtry_r <- as.matrix(res_mtry_r$Y)
      vec_quali_mtry_r <- res_mtry_r$vec_quali
      vec_order_mtry_r <- unlist(lapply(data_quali_mtry_r,function(x){class(x)[1]=="ordered"}))
    } else {
      X_quali_mtry_r <-matrix(0,n,0)
      vec_quali_mtry_r <- c()
      vec_order_mtry_r <- c()
    }
    
    X_mtry_r <- cbind(X_quanti_mtry_r, X_quali_mtry_r)

    cqd <- choose_question(X_mtry_r, Z, Fmin$v$A_l_c, vec_quali_mtry_r, w, D,vec_order_mtry_r)
    Fmin$r <- new.env()
    Fmin$r$v<- cqd
    leaves <- append(leaves, Fmin$r)

    k <- k+1
  }

  cluster <-list(tree=dendro)
  cluster$X <- X
  cluster$mod_quali <- mod_quali
  cluster$vec_quali <- vec_quali
  cluster$mu_quali <- mu_quali
  cluster$mu_quanti <- mu_quanti
  cluster$sig_quanti <- sig_quanti
  cluster$nb_quanti <- nb_quanti
  cluster$rnames <- rnames
  cluster$cnames_quanti <- colnames(data_quanti)
  cluster$cnames_quali <- cnames_quali
  cluster$kmax <- kmax
  cluster$height <- height
  cluster$T <- ncol(X)-nb_quali
  cluster$B <- sum(height)/cluster$T
  ret_leaves <- list_leaves(cluster)
  #cluster$description <-Reduce(append,lapply(ret_leaves$leaves, function(x) {x$path}))
  # suppress the first character of a string
  #cut_1 <- function(s) { Reduce(paste0,sapply(strsplit(s,""), identity)[-1])}
  #cluster$description <- lapply(cluster$description, cut_1)
  #cluster$description <- lapply(cluster$description, cut_1) # Yes we apply two times!

  cluster$mtry <- mtry
  cluster$description <- make_description(ret_leaves$leaves,cluster)
  names(cluster$description) <- paste("C", 1:K, sep = "")
  cluster$clusters <- lapply(ret_leaves$leaves, function(x) {rnames[x$class]})
  names(cluster$clusters) <- paste("C", 1:K, sep = "")
  cluster$inertia <- lapply(ret_leaves$leaves, function(x) {x$inert})
  names(cluster$inertia) <- paste("C", 1:K, sep = "")
  cluster$which_cluster <- ret_leaves$classes
  cluster$data_quanti <- data_quanti
  cluster$data_quali <- data_quali
  cluster$call <- cl
  class(cluster) <- "divclust"
  return(cluster)
}

#' Wines of Val de Loire data
#' @format A data frame with 21 rows (the number of wines) and
#' 31 columns: the first column corresponds to the label of origin,
#' the second column corresponds to the soil,
#' and the others correspond to sensory descriptors.
#' @source Centre de recherche INRA d'Angers
#' @description data refering to 21 wines of Val de Loire.
#' @name wine
NULL

#' Protein data
#' @format A data frame with 25 rows (the European countries) and 9 columns (the food groups)
#' @source Originated by A. Weber and cited in Hand et al., A Handbook of Small Data Sets, (1994, p. 297).
#' @description The data measure the amount of protein consumed for nine food groups in
#' 25 European countries. The nine food groups are red meat (RedMeat), white meat (WhiteMeat),
#' eggs (Eggs), milk (Milk), fish (Fish), cereal (Cereal), starch (Starch), nuts (Nuts), and fruits and vegetables (FruitVeg).
#' @name protein
NULL

#' Breeds of Dogs data
#' @format A data frame with 27 rows (the breeds of dogs) and  7 columns: their size, weight and speed with 3 categories
#' (small, medium, large), their intelligence (low, medium, high), their affectivity and aggressiveness
#' with 3 categories (low, high), their function (utility, compagny, hunting).
#' @source Originated by A. Brefort (1982) and cited in Saporta G. (2011).
#' @description Data refering to 27 breeds of dogs.
#' @name dogs
NULL

#' Equality case data
#' @format A numerial data frame with 20 rows and  2 columns simulated from two gaussien distributions.
#' @description  These data illstrate the case where two binary questions give the same bipartition and how
#' the binary question with the most discriminant variable (X1 here) is chosen.
#' @name equality_case
#' @examples
#' data(equality_case)
#' plot(equality_case) #X1 discriminates the bipartition better than X2
#' tree <- divclust(equality_case,K=3)
#' plot(tree,nqbin=1) # the binary question with X1 is chosen
NULL
