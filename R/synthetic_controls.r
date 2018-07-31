#' Generate synthetic control
#'
#' Input raw data and a vector of capture
#' efficiencies of cells.
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} or just matrix.
#' @param  BETA_vec A vector of capture efficiencies
#' (probabilities) of cells.
#' @return  List containing 2D matrix of synthetic control,
#' \code{BETA_vec} used and \code{lambda} used
#' in \code{rpois}.
#'
#' @details Simulate control data (based on Poisson
#' distribution).
#' @examples
#' data("EXAMPLE_DATA_list")
#' \dontrun{
#' SC_output<-SyntheticControl(Data=
#' EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta)
#' }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#' @export
#'
SyntheticControl<-function(Data,BETA_vec){

    if(methods::is(Data, "SummarizedExperiment")){

        if (
            is.null(
                SummarizedExperiment::assayNames(Data)) ||
            SummarizedExperiment::assayNames(Data)[1] != "Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"

            if (
                is.null(
                    colnames(SummarizedExperiment::assays(Data)[["Counts"]]))) {
                stop("Must supply sample/cell names!")
                }

        }
        Data<-SummarizedExperiment::assays(Data)[["Counts"]]
    }

    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
    }

    nGenes<-dim(Data)[1]
    nCells<-dim(Data)[2]
    beta_c<-BETA_vec
    Mean_depth<-mean(colSums(Data)/beta_c)
    Data_norm<-t(t(Data)/colSums(Data))*Mean_depth
    MU<-rowMeans(Data_norm)
    #mu_mat<- MU%*%t(rep(1,length(MU)))
    mu_mat<- MU%*%t(rep(1,nCells))

    N_c <- matrix(rpois(nGenes * nCells, lambda =mu_mat ),nrow = nGenes, ncol = nCells,byrow=FALSE)

    N_c<-DownSampling(N_c ,beta_c)

    return(list(N_c=N_c,beta_c=beta_c,lambda=MU))
}
