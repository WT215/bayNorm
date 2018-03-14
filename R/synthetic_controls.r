

#' Generate synthetic control
#'
#' Input raw data and a vector of capture efficiencies of cells.
#' @param Data A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec A vector of capture efficiencies of cells.
#' @return  List containing 2D matrix of synthetic control, \code{BETA_vec} used and \code{lambda} used in \code{rpois}.
#'
#' @details Simulate control data (based on Poisson distribution).
#'
#'
#' @export
#'
SyntheticControl<-function(Data,BETA_vec){
  nGenes<-dim(Data)[1]
  nCells<-dim(Data)[2]
  beta_c<-BETA_vec
  Mean_depth<-mean(colSums(Data)/beta_c)
  Data_norm<-t(t(Data)/colSums(Data))*Mean_depth
  MU<-rowMeans(Data_norm)
  mu_mat<- MU%*%t(rep(1,length(MU)))

  N_c <- matrix(rpois(nGenes * nCells, lambda =mu_mat ),nrow = nGenes, ncol = nCells,byrow=F)

  N_c<-DownSampling(N_c ,beta_c)

  return(list(N_c=N_c,beta_c=beta_c,lambda=MU))
}
