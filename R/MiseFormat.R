#' Title
#'
#' @param M
#' @param N
#' @param sigma
#'
#' @keywords internal
Image_transfo <- function(M,N, sigma){
  col = c(1:(M*N))%%N
  col[which(col==0)]=N
  ligne=rep(NA,M*N)
  for (i in 1:M){
    ligne[(N*(i-1)+1):(N*i)]=rep(i,N)
  }
  D= exp(-(as.matrix(dist( cbind(ligne,col), diag = TRUE))^2)/(2*sigma^2))/(sigma*sqrt(2*pi))
  vp = eigen(D)
  return( vp$vectors%*% diag(sqrt(vp$values)) %*% t(vp$vectors))
}

#' Title
#'
#' @param X
#' @param Q
#'
#' @keywords internal
transfo_lin <- function(X,Q){
  X_lin = matrix(0, nrow(X), ncol(X))
  for (i in 1:nrow(X_lin)){
    X_lin[i,]= Q%*%X[i,]
  }
  return(X_lin)
}


#' Factor partitions finder
#'
#' This function is used to find all the unique partitions of k factors into 2 groups
#'
#' @param Factor
#' @param id
#'
#' @keywords internal
Fact.partitions <- function(Factor, id){

  U <- unique(Factor)
  P <- Part.facts[[length(U)]]
  L <- list()
  for (k in 1:nrow(P)){
    w <- which(P[k,]==0)
    U_courant <- U[w]
    W <- NULL
    for (m in U_courant){
      W <- c(W,which(Factor==m))
    }
    L[[k]] <- id[W]
  }
  return(L)
}


#' Ordonne
#'
#' @param X
#' @param time
#' @param id
#'
#'
#' @keywords internal
ordonne <- function(X , time , id){
  mat  <- matrix(NA, length(unique(id)), length(unique(time)))
  for( i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    t_w <- time[w]
    w_time <- NULL
    for (j in 1:length(w)){
      w_time <- c(w_time, which(unique(time)==t_w[j]))
    }
    mat[i,w_time] <- X[w]
  }
  return(mat)
}
