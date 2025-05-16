
#' Get_Clusters
#' Constructs clustering structures from both the correlation matrix and the
#' partial correlation matrix estimated using sparse GGM.
#' @param Zscores A numeric matrix of Z-scores (rows = SNPs, columns = traits).
#' or correlation matrix of Z scores or Correlation matrix of the Phenotypes
#'
#' @param penalty_function Penalty type for GGMncv (default: "lasso")
#'
#'@param n_rows This will be number of SNPs or number of Individuals if supplied
#'the correlation matrix from phenotypes directly

#' @returns A list with:
#' B_cor`: Cluster matrix from correlation matrix
#' B_Pcor`: Cluster matrix from partial correlation matrix

#' @examples
#' Z <- MASS::mvrnorm(100, mu = rep(0, 10), Sigma = diag(10))
#' res <- Get_Clusters(Z)
#' @export

Get_Clusters <- function(Zscores, penalty_function = "lasso",n_rows=NULL) {
  # Validate Zscores input
  if (is.null(dim(Zscores)) || is.vector(Zscores)) {
    stop("Error: Zscores must be a two-dimensional matrix or array.\n ")
  }
  # check if any entries are missing in correlation matrix or Z scores data
  if (any(is.na(Zscores))) {
    stop("Error:Input matrix contains NA values. Please handle them before proceeding.\n ")
  }

  # Check penalty function
  penalty_list<-c("atan","mcp","scad","exp","selo","log","lasso","sica","lq","adapt")

  if (!penalty_function %in% penalty_list) {
      stop(cat("\n penalty function  not found. current options are:", penalty_list,"\n"))
    }

   Zscores<-as.matrix(Zscores)

  # Check If correlation matrix supplied instead of full matrix
  is_correlation_matrix <- function(mat) {
    is.matrix(mat) &&
      nrow(mat) == ncol(mat) &&
      is.numeric(mat) &&
      all(abs(mat - t(mat)) < .Machine$double.eps^0.5) &&
      all(diag(mat) == 1)
  }

  if (is_correlation_matrix(Zscores)) {
    Cor <- Zscores
    if (is.null(n_rows)) {
      stop("\n The number of rows 'n' must be provided when supplying a correlation matrix.\n ")
    }
    num_row <- n_rows
  } else {
    Cor <- aSPU::estcov(Zscores)
    num_row <- nrow(Zscores)
  }

  # Clustering based on correlation matrix
  D1 <- as.dist(1 - Cor)
  hc1 <- hclust(D1)
  n_cluster_cor <- get_n_cluster(hc1, Cor)$n_cluster
  index_cor <- cutree(hc1, n_cluster_cor)
  B_cor <- sapply(1:n_cluster_cor, function(t) as.integer(index_cor == t))

  # Estimate sparse partial correlation using GGMncv
  Model <- GGMncv::ggmncv(
    R = Cor, n = num_row,
    penalty = penalty_function, n_lambda = 100,
    ic = "ebic", progress = FALSE
  )
  Pcor <- Model$P
  # Clustering based on partial correlation matrix
  D2 <- as.dist(1 - Pcor)
  hc2 <- hclust(D2)
  n_cluster_pcor <- get_n_cluster(hc2, Pcor)$n_cluster
  index_pcor <- cutree(hc2, n_cluster_pcor)
  B_Pcor <- sapply(1:n_cluster_pcor, function(t) as.integer(index_pcor == t))

  # Return list of results
  return(list(
    B_cor = B_cor,
    B_Pcor = B_Pcor
  ))
}

#' Determine Optimal Number of Clusters Using Modularity
#'
#' Computes the optimal number of clusters in hierarchical clustering
#' using a modularity-based criterion.Maximizing modularity
#'
#' @param hc A hierarchical clustering object (from `hclust()`).
#' @param Sigma A simialrity matrix used for clustering (e.g., correlation or partial correlation matrix).
#' @param m Maximum number of clusters to evaluate. Default is number of columns in Sigma.
#' @param between_cluster Threshold: if all entries in Sigma > threshold, returns 1 cluster
#'
#' @returns A list with:
#' -- `n_cluster`: optimal number of clusters
#' --`Qmodularity`: modularity scores for each tested cluster number
#' @export
#' @examples
#' set.seed(123)
#' K<-6
#' R=matrix(0.5,nrow=K,ncol=K)
#' diag(R)<-1
#' hc<-hclust(as.dist(1-R))
#' res<-get_n_cluster(hc,R)

get_n_cluster <- function(hc, Sigma, m=ncol(Sigma), between_cluster = 0.8){
  if (min(Sigma) > between_cluster){
    IND = 1
    Q = 1
  } else {
    Q <- c()
    if (ncol(Sigma) < 10){m = ncol(Sigma)}
    for(i in 1:m)
    {
      index=cutree(hc,i)
      B=sapply(1:i, function(t) as.numeric(index==t))
      Q[i] <- get_modularity(Sigma, B)
    }

    IND=which(Q==max(Q))
    L=length(IND)
    if (L>1) IND=IND[1]
  }
  return(list("n_cluster" = IND,
              "Qmodularity" = Q))
}


get_modularity <- function(Weight, B){
  # ------ Calculate modularity ----------
  if (dim(Weight)[1] == 1){
    Q <- 0
  } else {
    W_pos <- Weight * (Weight > 0)
    W_neg <- Weight * (Weight < 0)
    N <- dim(Weight)[1]
    K_pos <- colSums(W_pos)
    K_neg <- colSums(W_neg)
    m_pos <- sum(K_pos)
    m_neg <- sum(K_neg)
    m <- m_pos + m_neg
    cate <- B %*% t(B)
    if (m_pos == 0 & m_neg == 0){
      Q <- 0
    } else {
      if (m_pos == 0){
        Q_positive <- 0
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      } else if (m_neg == 0){
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- 0
      } else {
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      }
    }
    Q <- m_pos / m * Q_positive - m_neg / m * Q_negative
  }
}
