#' @title Adjust MME size estimate
#'
#' @description  This function adjusts MME estimated size parameter
#' of prior, which is a negative binomial distribution,
#' using estimates from maximizing marginal distirbution
#' (\code{BB_SIZE}). Simulation studies has shown this hybrid
#' method of using adjusted MME size estimates is the most
#' robust (see bayNorm paper). Hence, this is the default
#' option for estimating size in bayNorm.
#'
#' @param BB_SIZE size estimated from \code{BB_Fun}.
#' @param  MME_MU mu estimated from EstPrior.
#' @param  MME_SIZE size estimated from EstPrior.
#' @return MME_SIZE_adjust: A vector of estimated size.
#' Adjusted MME_SIZE based on BB_SIZE
#' (size estimated by maximizing marginal distribution)
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' MME_MU<-rlnorm(100,meanlog=5,sdlog=1)
#' MME_SIZE<-rlnorm(100,meanlog=1,sdlog=1)
#' BB_SIZE<-rlnorm(100,meanlog=0.5,sdlog=1)
#' adjustt<-AdjustSIZE_fun(BB_SIZE, MME_MU, MME_SIZE)
#' @import stats
#'
#' @export
AdjustSIZE_fun <- function(BB_SIZE, MME_MU, MME_SIZE) {
  fitind <- which(BB_SIZE < max(BB_SIZE) & BB_SIZE > min(BB_SIZE))
  lmfit <- lm(log(BB_SIZE)[fitind] ~ log(MME_SIZE)[fitind])
  MME_SIZE_adjust <- coef(lmfit)[1] + coef(lmfit)[2] * log(MME_SIZE)
  MME_SIZE_adjust <- exp(MME_SIZE_adjust)
  return(MME_SIZE_adjust)
}

#' @title Rcpp version's as.matrix 
#'
#' @description  Transform sparse matrix to matrix.
#'
#' @param mat Sparse matrix.
#' @return Matrix.
#'
#' @examples
#' aa<-matrix(seq(1,6),nrow=2,ncol=3)
#' qq<-as(as.matrix(aa), "dgCMatrix")
#' all.equal(unname(as_matrix(qq)),unname(as.matrix(qq)))
#' @import stats
#'
#' @export
as_matrix <- function(mat){
  
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}






#' @title Check input 
#'
#' @description  Check input 
#'
#' @param Data Input data.
#' @return Matrix (of type matrix in R).
#'
#' @examples
#' aa<-matrix(seq(1,6),nrow=2,ncol=3)
#' Check_input(aa)
#'
#' @export
Check_input <- function(Data){
  
  if(!is(Data, 'matrix')){
    if(is(Data, 'sparseMatrix')){
      if(!is(Data, 'dgCMatrix')){
        Data <- as(as.matrix(Data), "dgCMatrix")
        Data<-as_matrix(Data)
      } else{
        Data<-as_matrix(Data)
      }
    } else{
      
      if (!(methods::is(Data, "SummarizedExperiment")) &
          !(methods::is(Data, "SingleCellExperiment"))) {
        Data <- as.matrix(Data)
      }
    }
  }
  return(Data)
}

