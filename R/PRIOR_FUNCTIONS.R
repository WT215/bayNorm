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
#' \dontrun{
#' MME_MU<-rlnorm(100,meanlog=5,sdlog=1)
#' MME_SIZE<-rlnorm(100,meanlog=1,sdlog=1)
#' BB_SIZE<-rlnorm(100,meanlog=0.5,sdlog=1)
#' adjustt<-AdjustSIZE_fun(BB_SIZE, MME_MU, MME_SIZE)
#' }
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


#' @title Estimate capture efficiency for cells
#'
#' @description  This function estimates cell specific
#' capture efficiencies (\code{BETA_vec}) using mean raw counts of
#' a subset of genes that is an input for bayNorm. A specific
#' method is used to exclude genes with high expression or high
#' drop-out are excluded.
#'
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix and
#' is named "Counts") or just matrix.
#' @param MeanBETA Mean capture efficiency of the scRNAseq data.
#'  This can be estimated via spike-ins or other methods.
#' @return List containing: \code{BETA}: a vector of capture efficiencies,
#'  which is of length number of cells;
#'  \code{Selected_genes}: a subset of
#'   genes that are used for estimating BETA.
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' \dontrun{
#' BETA_out<-BetaFun(Data=EXAMPLE_DATA_list$inputdata,
#' MeanBETA=0.06)
#' }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#' @export
BetaFun <- function(Data, MeanBETA) {

    if (methods::is(Data, "SummarizedExperiment")) {

        if (
            is.null(SummarizedExperiment::assayNames(Data)) ||
            SummarizedExperiment::assayNames(Data)[1] !="Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"

            if (
                is.null(
                    colnames(SummarizedExperiment::assays(Data)[["Counts"]]))) {
                stop("Must supply sample/cell names!")
            }

        }
        Data <- SummarizedExperiment::assays(Data)[["Counts"]]
    }

    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
    }


    Normcount <- t(t(Data)/colSums(Data)) * mean(colSums(Data))
    means <- rowMeans(Normcount)
    lmeans <- log(means)
    med <- apply(log(Normcount + 1), 1, function(x) {
        median(x)
    })
    mad <- apply(log(Normcount + 1), 1, function(x) {
        mad(x)
    })
    bound <- med + 3 * mad
    maxlogGene <- apply(log(Normcount + 1), 1, max)
    ind <- which(maxlogGene < bound)
    dropout = apply(Data, 1, function(x) {
        length(which(x == 0))/length(x)
    })
    Select_ind <- intersect(ind, which(dropout < 0.35))
    Selected_genes <- rownames(Data)[Select_ind]

    temppp <- colSums(Data[Select_ind, ])
    BETA <- temppp/mean(temppp) * MeanBETA
    names(BETA) <- colnames(Data)
    return(list(BETA = BETA, Selected_genes = Selected_genes))
}

#' @title Estimate size and mu for Negative Binomial distribution
#' for each gene using MME method
#'
#' @description  Input raw data and return
#' estimated size and mu for each gene using the MME method.
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix and
#' is named "Counts") or just matrix.
#' @param  parallel If TRUE, 5 cores will be used for parallelization.
#' Default is TRUE.
#' @param  NCores number of cores to use, default is 5.
#' This will be used to set up a parallel
#' environment using either MulticoreParam (Linux, Mac)
#' or SnowParam (Windows) with NCores
#' using the package BiocParallel.
#' @param verbose print out status messages. Default is TRUE.
#'
#' @details mu and size are two parameters of the prior that
#' need to be specified for each gene in bayNorm.
#' They are parameters of negative binomial distribution.
#' The variance is \eqn{mu + mu^2/size} in this parametrization.
#'
#' @return  List containing estimated mu and
#' size for each gene.
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' #Return estimated mu and size for each gene using MME method.
#' \dontrun{
#' MME_est<-EstPrior(Data=EXAMPLE_DATA_list$inputdata,
#' verbose=TRUE)
#' }
#' @import  fitdistrplus
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#' @export
EstPrior <- function(Data,parallel=FALSE,NCores=5, verbose = TRUE) {

    if (methods::is(Data, "SummarizedExperiment")) {

        if (
            is.null(
                SummarizedExperiment::assayNames(Data)) ||
            SummarizedExperiment::assayNames(Data)[1] !=
            "Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"

            if (
                is.null(
                    colnames(SummarizedExperiment::assays(Data)[["Counts"]]))) {
                stop("Must supply sample/cell names!")
            }

        }
        Data <- SummarizedExperiment::assays(Data)[["Counts"]]
    }

    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
    }

    i <- NULL

    if(parallel){
        cluster = makeCluster(NCores, type = "SOCK")
        registerDoSNOW(cluster)
        getDoParWorkers()

        iterations <- dim(Data)[1]
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        CoefDat <- foreach(
            i = seq_len(nrow(Data)),
            .combine = rbind,
            .options.snow = opts) %dopar% {

                suppressWarnings(
                    qq <- fitdistrplus::fitdist(
                        Data[i, ], "nbinom",
                        method = "mme",
                        keepdata = FALSE))
                return(coef(qq))
            }
        close(pb)
        stopCluster(cluster)
    } else{
        iterations <- dim(Data)[1]
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)

        CoefDat <- foreach(
            i = seq_len(nrow(Data)),
            .combine = rbind,
            .options.snow = opts) %do% {
                setTxtProgressBar(pb, i)

                suppressWarnings(
                    qq <- fitdistrplus::fitdist(
                        Data[i, ], "nbinom",
                        method = "mme",
                        keepdata = FALSE))
                return(coef(qq))
            }
        close(pb)
    }

    rownames(CoefDat) <- rownames(Data)
    M_ave_ori <- CoefDat[, 2]
    size_est <- CoefDat[, 1]
    names(M_ave_ori)<-rownames(Data)
    names(size_est)<-rownames(Data)

    if (verbose) {
        message("Priors estimation based on MME method has completed.")
    }

    return(list(MU = M_ave_ori, SIZE = size_est))
}
#' @title   A wrapper function of \code{EstPrior}
#' and \code{AdjustSIZE_fun}
#'
#' @description   A wrapper function for estimating the parameters
#' of prior using the hybrid method adjusted MME estimates based
#' on maximization of marginal likelihood. Input raw data and a
#' vector of capture efficiencies of cells.
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix and
#' is named "Counts") or just matrix.
#' @param  BETA_vec A vector of capture efficiencies of cells.
#' @param  parallel If TRUE, 5 cores will be used for
#' parallelization. Default is TRUE.
#' @param  NCores number of cores to use, default is 5.
#' This will be used to set up a parallel environment
#' using either MulticoreParam (Linux, Mac) or
#' SnowParam (Windows) with \code{NCores} using the
#' package \code{BiocParallel}.
#' @param  FIX_MU If TRUE, then 1D optimization, otherwise
#' 2D optimization (slow). Default is TRUE.
#' @param  GR If TRUE, the gradient function will be used
#' in optimization. However since the gradient function
#' itself is very complicated, it does not help too much
#' in speeding up. Default is FALSE.
#' @param  BB_SIZE If TRUE, estimate BB size, and then use
#' it for adjusting MME SIZE. Use the adjusted MME size
#' for bayNorm. Default is TRUE.
#' @param verbose Print out status messages. Default is TRUE.
#'
#' @details By Default, this function will estimate mu and
#' size for each gene using MME method. If \code{BB_size}
#' is enable, spectral projected gradient method from BB
#' package will be implemented to estimate 'BB size' by
#' maximizing marginal likelihood function. MME estimated
#' size will be adjusted according to BB size. BB size itself
#' will not be used in bayNorm this is because that in
#' our simulation we found that MME estimated mu and size
#' have more accurate relationship, but MME estimated
#' size deviates from the true value. BB size is overall
#' more close to the true size but it does not possess a
#' reasonable relationship with either MME estimated mu or
#' BB estimated mu.
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' \dontrun{
#' PRIOR_RESULT<-Prior_fun(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta,parallel=T,
#' NCores=5,FIX_MU=T,GR=F,BB_SIZE=T,verbose=T)
#' }
#'
#' @return  List of estimated parameters: mean
#' expression of genes
#' and size of each gene.
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' \dontrun{
#'PRIOR_RESULT<-Prior_fun(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta,parallel=TRUE,
#' NCores=5,FIX_MU=TRUE,GR=FALSE,BB_SIZE=TRUE,verbose=TRUE)
#'
#' }
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#' @export
#'
Prior_fun <- function(
    Data, BETA_vec, parallel = TRUE,
    NCores = 5, FIX_MU = TRUE,GR = FALSE,
    BB_SIZE = TRUE, verbose = TRUE) {

    if (methods::is(Data, "SummarizedExperiment")) {

        if (is.null(
            SummarizedExperiment::assayNames(Data))
            || SummarizedExperiment::assayNames(Data)[1] !="Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"

            if (
                is.null(
                    colnames(SummarizedExperiment::assays(Data)[["Counts"]]))) {
                stop("Must supply sample/cell names!")
            }

        }

        Data <- SummarizedExperiment::assays(Data)[["Counts"]]
    }

    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
    }

    normcount_N <- t(t(Data)/colSums(Data)) * mean(colSums(Data)/BETA_vec)
    Priors <- EstPrior(normcount_N, verbose = verbose,parallel=parallel,NCores=NCores)
    M_ave_ori <- Priors$MU
    size_est <- Priors$SIZE
    size_est[is.na(size_est)] <- min(size_est[!is.na(size_est)])
    MME_prior <- cbind(M_ave_ori, size_est)
    MME_prior <- as.data.frame(MME_prior)
    rownames(MME_prior) <- rownames(Data)
    colnames(MME_prior) <- c("MME_MU", "MME_SIZE")


    if (BB_SIZE) {
        if (verbose) {
            message("Start optimization using spg from BB package.
This part may be time-consuming.")
        }
        BB_size <- BB_Fun(
            Data, BETA_vec,
            INITIAL_MU_vec = MME_prior$MME_MU,
            INITIAL_SIZE_vec = MME_prior$MME_SIZE,
            MU_lower = min(MME_prior$MME_MU),
            MU_upper = max(MME_prior$MME_MU),
            SIZE_lower = min(MME_prior$MME_SIZE),
            SIZE_upper = ceiling(
                max(MME_prior$MME_SIZE)),
            parallel = parallel,
            NCores = NCores,
            FIX_MU = FIX_MU,
            GR = GR)

        if (FIX_MU) {

            BB_prior <- cbind(MME_prior$MME_MU, BB_size)
            rownames(BB_prior) <- rownames(Data)
            colnames(BB_prior) <- c("MME_MU", "BB_SIZE")

            MME_SIZE_adjust <- AdjustSIZE_fun(
                BB_prior[, 2],
                MME_prior$MME_MU,
                MME_prior$MME_SIZE)
        } else {
            BB_prior <- BB_size
            rownames(BB_prior) <- rownames(Data)
            colnames(BB_prior) <- c("BB_SIZE", "BB_MU")
            MME_SIZE_adjust <- AdjustSIZE_fun(
                BB_prior[, 1],
                MME_prior$MME_MU,
                MME_prior$MME_SIZE)
        }

        if (verbose) {
            message("Prior estimation has completed!")
        }
        names(MME_SIZE_adjust)<-rownames(Data)
        rownames(MME_prior)<-rownames(Data)
        rownames(BB_prior)<-rownames(Data)

        return(list(
            MME_prior = MME_prior, BB_prior = BB_prior,
            MME_SIZE_adjust = MME_SIZE_adjust,BETA_vec=BETA_vec))
    }

    return(list(MME_prior = MME_prior,BETA_vec=BETA_vec))
}


#' @title Estimating size for each gene by
#' either 1D or 2D maximization of marginal distribution
#'
#' @description  Estimating parameters of the prior distribution
#' for each gene by maximizing marginal distribution: 1D
#' (optimize with respect to size using MME estimate of mu,
#' 2D (optimize with respect to both mu and size)
#'
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix and
#' is named "Counts") or just matrix.
#' @param  BETA_vec A vector of capture efficiencies
#' (probabilities) of cells.
#' @param  INITIAL_MU_vec Mean expression of genes,
#' can be estimated from \code{EstPrior}.
#' @param  INITIAL_SIZE_vec size of genes (size is a parameter in
#' NB distribution), can come from EstPrior.
#' @param  MU_lower The lower bound for the mu.(Only need it when
#' you want to do 2D optimization). Default is 0.01.
#' @param  MU_upper The upper bound for the mu.(Only need it when
#' you want to do 2D optimization). Default is 500.
#' @param  SIZE_lower The lower bound for the size. Default is 0.01.
#' @param  SIZE_upper The upper bound for the size. Default is 30.
#' @param  parallel If TRUE, 5 cores will be used for parallelization.
#' Default is TRUE.
#' @param  NCores number of cores to use, default is 5. This will
#' be used to set up a parallel environment using either
#' MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores
#' using the package BiocParallel.
#' @param  FIX_MU If TRUE, then 1D optimization, otherwise 2D
#' optimization (slow).
#' @param  GR If TRUE, the gradient function will be used in
#' optimization. However since the gradient function itself is
#' very complicated, it does not help too much in speeding up.
#' Default is FALSE.
#' @return  BB estimated size (1D optimization) or size and mu (2D optimization).
#'
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' \dontrun{
#'
#'BB_RESULT<-BB_Fun(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta,INITIAL_MU_vec=
#' EXAMPLE_DATA_list$mu,
#' INITIAL_SIZE_vec=EXAMPLE_DATA_list$size,
#' MU_lower=0.01,MU_upper=500,SIZE_lower=0.01,
#' SIZE_upper=30,parallel=FALSE,NCores=5,FIX_MU=TRUE,GR=FALSE)
#' }
#'
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import utils
#' @import iterators
#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#' @useDynLib bayNorm
#' @importFrom Rcpp sourceCpp
#'
#' @export
#'
#'
BB_Fun <- function(
    Data, BETA_vec, INITIAL_MU_vec, INITIAL_SIZE_vec,
    MU_lower = 0.01, MU_upper = 500, SIZE_lower = 0.01,
    SIZE_upper = 30,parallel = FALSE, NCores = 5,
    FIX_MU = TRUE, GR = FALSE) {

    if (methods::is(Data, "SummarizedExperiment")) {

        if (
            is.null(SummarizedExperiment::assayNames(Data))
            || SummarizedExperiment::assayNames(Data)[1] !=
            "Counts") {
            message("Renaming the first element in assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"

            if (
                is.null(
                    colnames(SummarizedExperiment::assays(
                        Data)[["Counts"]]))) {
                stop("Must supply sample/cell names!")
            }

        }
        Data <- SummarizedExperiment::assays(Data)[["Counts"]]
    }

    if (!(methods::is(Data, "SummarizedExperiment"))) {
        Data <- data.matrix(Data)
    }


    Geneind <- NULL

    if(!GR){
        GR_pass<-NULL
    } else if(GR){
        if(FIX_MU){
            GR_pass<-GradientFun_NB_1D
        } else{
            GR_pass<-GradientFun_NB_2D
        }


    }


    if (FIX_MU) {
        # 1D

        lower_input = SIZE_lower
        upper_input = SIZE_upper

        if (parallel) {
            cluster = makeCluster(NCores, type = "SOCK")
            registerDoSNOW(cluster)
            getDoParWorkers()

            iterations <- dim(Data)[1]
            pb <- txtProgressBar(max = iterations, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)

            BB_parmat <- foreach(
                Geneind = seq_len(dim(Data)[1]),
                .combine = c,
                .options.snow = opts) %dopar% {
                    # print(Geneind)
                    mu <- INITIAL_MU_vec[Geneind]
                    size <- INITIAL_SIZE_vec[Geneind]
                    m_observed = Data[Geneind, ]

                    BB_opt <- BB::spg(
                            par = size, fn = MarginalF_NB_1D,
                            gr = GR_pass,MU = mu,
                            m_observed = m_observed,
                            BETA = BETA_vec,
                            control = list(
                                maximize = TRUE,
                                trace = FALSE,
                                maxfeval = 500),
                            lower = lower_input,
                            upper = upper_input)


                    optimal_par <- BB_opt$par
                    return(optimal_par)
                }

            close(pb)
            stopCluster(cluster)


        } else {

            iterations <- dim(Data)[1]
            pb <- txtProgressBar(max = iterations, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)


            BB_parmat <- foreach(
                Geneind = seq_len(dim(Data)[1]),
                .combine = c,.options.snow = opts) %do% {

                    # print(Geneind)
                    setTxtProgressBar(pb, Geneind)

                    mu <- INITIAL_MU_vec[Geneind]
                    size <- INITIAL_SIZE_vec[Geneind]
                    m_observed = Data[Geneind, ]

                    BB_opt <- BB::spg(
                            par = size, fn = MarginalF_NB_1D,
                            gr = GR_pass,MU = mu,
                            m_observed = m_observed,
                            BETA = BETA_vec,
                            control = list(
                                maximize = TRUE,
                                trace = FALSE,
                                maxfeval = 500),
                            lower = SIZE_lower,
                            upper = SIZE_upper)


                    #
                    optimal_par <- BB_opt$par
                    return(optimal_par)
                }
            close(pb)

        }
        names(BB_parmat) <- rownames(Data)
    } else {
        # 2D

        lower_input = c(SIZE_lower, MU_lower)
        upper_input = c(SIZE_upper, MU_upper)


        if (parallel) {
            cluster = makeCluster(NCores, type = "SOCK")
            registerDoSNOW(cluster)
            getDoParWorkers()

            iterations <- dim(Data)[1]
            pb <- txtProgressBar(max = iterations, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)


            BB_parmat <- foreach(
                Geneind = seq_len(dim(Data)[1]),
                .combine = rbind,
                .options.snow = opts) %dopar%
                {

                    mu <- INITIAL_MU_vec[Geneind]
                    size <- INITIAL_SIZE_vec[Geneind]
                    m_observed = Data[Geneind, ]

                    BB_opt <- BB::spg(
                            par = c(size, mu),
                            fn = MarginalF_NB_2D,
                            gr = GR_pass,
                            m_observed = m_observed,
                            BETA = BETA_vec,
                            control = list(
                                maximize = TRUE,
                                trace = FALSE,
                                maxfeval = 500),
                            lower = lower_input,
                            upper = upper_input)


                    optimal_par <- BB_opt$par
                    return(optimal_par)
                }
            close(pb)
            stopCluster(cluster)


        } else {


            iterations <- dim(Data)[1]
            pb <- txtProgressBar(max = iterations, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
            BB_parmat <- foreach(
                Geneind = seq_len(dim(Data)[1]),
                .combine = rbind) %do%

                {
                    setTxtProgressBar(pb, Geneind)
                    mu <- INITIAL_MU_vec[Geneind]
                    size <- INITIAL_SIZE_vec[Geneind]
                    m_observed = Data[Geneind, ]

                    BB_opt <- BB::spg(
                            par = c(size, mu), fn = MarginalF_NB_2D,
                            gr = GR_pass,
                            m_observed = m_observed,
                            BETA = BETA_vec,
                            control = list(
                                maximize = TRUE,
                                trace = FALSE,
                                maxfeval = 500),
                            lower = c(SIZE_lower, MU_lower),
                            upper = c(SIZE_upper, MU_upper))

                    optimal_par <- BB_opt$par
                    return(optimal_par)
                }

            close(pb)


        }


        rownames(BB_parmat) <- rownames(Data)

    }

    return(BB_parmat)
}

