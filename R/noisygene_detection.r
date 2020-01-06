#' A wrapper function for noisy gene detection from raw data.
#' his produces synthetic control, performs bayNorm on both real
#' cell data and synthetic controls and does noisy gene detection.
#'
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix,
#' is named "Counts"), just \code{matrix} or sparse matrix.
#' @param  BETA_vec A vector of capture efficiencies
#' of cells.
#' @param PRIORS (Need to be specified for efficiency
#' if \code{bayNorm} has already been applied)
#' A list of estimated prior
#' parameters obtained from bayNorm. Default is NULL.
#' @param  input_params (Need to be specified for efficiency
#' if \code{bayNorm} has already been applied)
#' A list of input parameters which have been used:
#' \code{BETA_vec}, \code{Conditions},
#' \code{UMI_sffl}, \code{Prior_type},
#' \code{FIX_MU}, \code{BB_SIZE} and \code{GR}.
#' @param  mode_version If TRUE, bayNorm return mode version
#' normalized data which is of 2D matrix instead of 3D array.
#' Default is FALSE.
#' @param  mean_version If TRUE, bayNorm return
#' mean version normalized data which
#' is of 2D matrix instead of 3D array.
#' Default is FALSE.
#' @param S The number of samples you would like
#' to generate from estimated posterior distribution
#' (The third dimension of 3D array).
#' Default is 20. S needs to be specified if
#' \code{mode_version}=FALSE.
#' @param  parallel If TRUE, \code{NCores} cores will be
#' used for parallelization. Default is TRUE.
#' @param  NCores number of cores to use, default is 5.
#' This will be used to set up a parallel environment
#' using either MulticoreParam
#' (Linux, Mac) or SnowParam (Windows) with
#' NCores using the package BiocParallel.
#' @param  FIX_MU Whether fix mu when estimating
#' parameters by
#' maximizing marginal distribution.
#' If TRUE, then 1D optimization,
#' otherwise 2D optimization (slow).
#' @param  GR If TRUE, the gradient function will
#' be used in optimization. However since
#' the gradient function itself is very complicated,
#' it does not help too much in speeding up.
#' Default is FALSE.
#' @param  BB_SIZE If TRUE, estimate BB size,
#' and then use it for
#' adjusting MME SIZE. Use the adjusted MME size
#' for bayNorm. Default is TRUE.
#' @param verbose Print out status messages.
#' Default is TRUE.
#' @param  plot.out If TRUE, show CV^2 vs Mean
#' expression plot.
#' Default is FALSE.
#' @return A list of objects.
#' @details A wrapper function for noisy gene detection
#' from raw scRNA-seq data.
#'
#' @examples
#' data("EXAMPLE_DATA_list")
#' noisy_out<-noisy_gene_detection(Data=
#' EXAMPLE_DATA_list$inputdata[,seq(1,30)],BETA_vec
#' =EXAMPLE_DATA_list$inputbeta[seq(1,30)], mode_version = FALSE,
#' mean_version=FALSE,
#' S = 20,parallel = FALSE, NCores = 5,
#' FIX_MU = TRUE, GR = FALSE,
#' PRIORS=NULL,
#' BB_SIZE = TRUE,
#' verbose = TRUE, plot.out = TRUE)
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @importFrom Matrix colSums
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#' @export
#'

noisy_gene_detection<-function(
    Data,BETA_vec=NULL,
    mode_version=FALSE,
    mean_version=FALSE,
    S=20,
    parallel=TRUE,NCores=5,
    FIX_MU=TRUE,GR=FALSE,
    BB_SIZE=TRUE,verbose=TRUE,
    plot.out=FALSE,
    PRIORS=NULL,
    input_params=NULL){

    if(mode_version & mean_version){
        stop("Only one of mode_version and mean_version
             should be specified to be TRUE, otherwise both
             should be set to FALSE so that 3D array
             normalized data will be returned.")
    }

    if(!mode_version & !mean_version){
        myFunc <- Main_NB_Bay
    }else if(mode_version & !mean_version){
        myFunc <- Main_mode_NB_Bay
    }else if(!mode_version & mean_version){
        myFunc <-  Main_mean_NB_Bay
    }



    #Set matrix as object for input data
    Data<-Check_input(Data)


    message("Apply bayNorm on the real cell datasets.")

    if(is.null(PRIORS)){
        bayNorm_N_out<-bayNorm(
            Data=Data,BETA_vec=BETA_vec,
            mode_version=mode_version,
            mean_version=mean_version,
            Conditions=NULL,UMI_sffl=NULL,
            Prior_type = NULL,S=S,
            parallel = parallel,NCores=NCores,
            FIX_MU = FIX_MU,GR=GR,
            BB_SIZE=BB_SIZE,verbose=verbose)


    } else if(!is.null(PRIORS) & !is.null(input_params)){

        BETA_vec<-PRIORS$BETA_vec

        bayNorm_N_out<-bayNorm_sup(
            Data=Data,
            PRIORS=PRIORS,
            input_params = input_params,
            S=S,
            mode_version=mode_version,
            mean_version=mean_version,
            parallel = parallel,NCores=NCores,
            BB_SIZE=BB_SIZE,verbose=verbose)
        }


    message("Generate synthetic control for the input data.")
    synthetic_out<-SyntheticControl(Data=Data,BETA_vec=BETA_vec)

    message("Apply bayNorm on the synthetic control.")

    mu_c<-bayNorm_N_out$PRIORS$MME_prior$MME_MU
    if(BB_SIZE){
        size_c<-bayNorm_N_out$PRIORS$MME_SIZE_adjust

    }else if(!BB_SIZE){
        size_c<-bayNorm_N_out$PRIORS$MME_prior$MME_SIZE
    }


    #tempppdat<-as(as.matrix(synthetic_out$N_c), "dgCMatrix")
    tempppdat<-Check_input(synthetic_out$N_c)
    bayNorm_C_array<-myFunc(Data=tempppdat,
                            BETA_vec=synthetic_out$beta_c,size=size_c,
                            mu=mu_c,S=S,
                            thres=max(synthetic_out$N_c)*2)

    NOISE_out<-NOISY_FUN(
        bayNorm_N_out[[1]],
        bayNorm_C_array,plot.out=plot.out)

    return(list(
        adjusted_Pvals=NOISE_out,
        synthetic_output=synthetic_out,
        bayNorm_N_out=bayNorm_N_out,
        bayNorm_C_array=bayNorm_C_array))

}





#' Noisy gene detection
#'
#' This function detects noisy genes using trends observed
#' in a set of synthetic controls. Input bayNorm normalized data
#' of real data (\code{bay_array_N}) and synthetic control
#' (\code{bay_array_C}) respectively.
#' @param bay_array_N A 2D matrix or 3D array of normalized
#' data(real cells).
#' @param  bay_array_C A 2D matrix or 3D array of normalized
#' data(synthetic control).
#' @param  plot.out If TRUE, show CV^2 vs Mean
#' expression plot. Default is FALSE.
#' @return A vector of adjusted P-values.
#' @details \code{bay_array_N} and \code{bay_array_C}
#' should be of the same dimension.
#'
#' @examples
#' bay_array_N<-array(rpois(1000*50*2,17),dim=c(1000,50,2))
#' bay_array_C<-array(rpois(1000*50*2,58),dim=c(1000,50,2))
#' noisy_output<-NOISY_FUN(bay_array_N,bay_array_C)
#'
#' @import foreach
#' @import MASS
#' @import locfit
#' @import grDevices
#' @import graphics
#' @export
#'



NOISY_FUN<-function(bay_array_N,bay_array_C,plot.out=FALSE){

    if(length(dim(bay_array_N))==3 &  length(dim(bay_array_C))==3){
        mean_1 = apply(bay_array_N, 1, mean)
        mean_1c = apply(bay_array_C, 1, mean)

        sd_1 = apply(apply(bay_array_N, c(1,3), sd), 1, mean)
        sd_1c = apply(apply(bay_array_C, c(1,3), sd), 1, mean)

    } else if(length(dim(bay_array_N))==2 &  length(dim(bay_array_C))==2){
        mean_1 = rowMeans(bay_array_N)
        mean_1c = rowMeans(bay_array_C)

        sd_1 = apply(bay_array_N, 1,sd)
        sd_1c = apply(bay_array_C, 1,sd)

    }

    DIM=dim(bay_array_N)
    Gene_names<-rownames(bay_array_N)

    noise_1 = (sd_1/mean_1)^2
    noise_1c = (sd_1c/mean_1c)^2

    Select<-!is.na(noise_1)
    x=log(mean_1)[Select]
    y=log(noise_1)[Select]

    Select_c<-!is.na(noise_1c)
    xc = log(mean_1c)[Select_c]
    yc = log(noise_1c)[Select_c]

    loessfit<-locfit.robust(yc ~ xc)
    h<-predict(loessfit, newdata=x, se.fit=TRUE)
    zval<-(y-h$fit)/h$residual.scale

    #####script###
    # use kernel density estimate to find the peak
    dor <- density(zval, kernel = "gaussian")
    distMid <-dor$x[which.max(dor$y)]
    dist2 <- zval - distMid
    tmpDist <- c(dist2[dist2 <= 0], abs(dist2[dist2 < 0])) + distMid
    distFit <- fitdistr(tmpDist, "normal")
    pRaw <- pnorm(
        zval, mean = distFit$estimate[1],
        sd = distFit$estimate[2], lower.tail = FALSE)
    pAdj <- p.adjust(pRaw, 'BH')
    noisy<-which(pAdj<0.05)
    noisy_ind<-match(names(noisy),Gene_names[Select])
    noisy_name<-names(noisy)

    col_N<-1
    col_C<-rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
    col_noisy<-rgb(red = 0, green = 1, blue = 0, alpha = 0.5)

    if(plot.out){
        plot(
            x, y,pch=16,main='Noisy genes detection',
            col=col_N,xlab='log: Mean expression',ylab='log: CV^2')
        lines(loessfit,col=6,lwd=3)
        points(xc, yc, col = col_C,pch=16)
        points(x[noisy_ind], y[noisy_ind], col = col_noisy,pch=16)
        legend(
            'topright',
            legend=c('Genes: real cell data',
                'Genes: synthetic control',
                'Noisy genes'),
            col=c(col_N,col_C,col_noisy),pch=16,cex=1)
    }


    pAdj2<-rep(1,length(Gene_names))
    pAdj2[Select]<-pAdj
    names(pAdj2)<-Gene_names
    return(pAdj2)
}





