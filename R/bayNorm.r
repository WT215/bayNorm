#' @title  A wrapper function of prior estimation and bayNorm function
#'
#' @description   This is the main wrapper function for bayNorm.
#' The input is a matrix of raw scRNA-seq data and a vector of
#' capture efficiencies of cells. You can also specify the
#' condition of cells for normalizing multiple groups
#' of cells separately.
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix and
#' is named "Counts") or just matrix.
#' @param  BETA_vec A vector of capture efficiencies
#' (probabilities) of cells.
#' If it is null, library size (total count) normalized to
#' 0.06 will be used
#' as the input \code{BETA_vec}. \code{BETA_vec} less and
#' equal to 0 or greater and equal to 1 will be replaced
#' by the minimum and maximum of the BETA_vec which
#' range between (0,1) respectively.
#' @param  Conditions vector of condition labels,
#' this should correspond to the columns of the Data. D
#' efault is NULL, which assumes that all cells
#' belong to the same group.
#' @param UMI_sffl Scaling factors are required only for
#' non-UMI based data for which \code{Data} is devided by
#' \code{UMI_sffl}. If non-null and \code{Conditions} is non-null,
#' then UMI_sffl should be a vector of length equal
#' to the number of groups. Default is \code{NULL}.
#' @param  Prior_type Determines what groups of cells is used
#' in estimating prior using \code{Conditions}.
#' Default is \code{NULL}.
#' If \code{Conditions} is \code{NULL},
#' priors are estimated based on all cells.
#' If \code{Conditions} is not \code{NULL} and
#' if \code{Prior_type} is LL,
#' priors are estimated within each group respectively.
#' If \code{Prior_type} is GG, priors are estimated based on cells
#' from all groups. LL is suitable for DE detection.
#' GG is preferred if reduction of batch effect between
#' samples are desired for example for technical replicates
#' (see bayNorm paper).
#' @param  mode_version If TRUE, bayNorm return modes of
#' posterior estimates as normalized data which is a 2D matrix
#' rather than samples from posterior which is a 3D array.
#' Default is FALSE.
#' @param  mean_version If TRUE, bayNorm return means of
#' posterior estimates as normalized data, which is a 2D matrix
#' rather than samples from posterior which is a 3D array.
#' Default is FALSE.
#' @param S The number of samples you would like to
#' generate from estimated posterior distribution
#' (The third dimension of 3D array). Default is 20.
#'  S needs to be specified if \code{mode_version}=FALSE.
#' @param  parallel If TRUE, multiple cores will be used
#' for parallelization.
#' @param  NCores number of cores to use, default is 5.
#' This will be used to set up a parallel environment
#' using either MulticoreParam (Linux, Mac) or
#' SnowParam (Windows) with NCores using
#' the package BiocParallel.
#' @param  FIX_MU Whether fix mu (the mean parameter of prior
#' distribution) to its MME estimate, when estimating prior
#' parameters by maximizing marginal distribution. If TRUE,
#' then 1D optimization is used, otherwise 2D optimization
#' for both mu and size is used (slow). Default is TRUE.
#' @param  GR If TRUE, the gradient function will be used
#' in optimization. However since the gradient function
#' itself is very complicated, it does not help too much
#' in speeding up. Default is FALSE.
#' @param  BB_SIZE If TRUE, estimate size parameter of
#' prior using maximization of marginal likelihood,
#' and then use it for adjusting MME estimate of SIZE Default is TRUE.
#' @param verbose print out status messages. Default is TRUE.
#' @return  List containing 3D arrays of normalized
#' expression (if \code{mode_version}=FALSE) or 2D matrix
#' of normalized expression (if \code{mode_version}=TRUE
#' or \code{mean_version}=TRUE),
#' a list contains estimated priors and a list contains
#' input parameters used: \code{BETA_vec},
#' \code{Conditions} (if specified),
#' \code{UMI_sffl} (if specified), \code{Prior_type},
#' \code{FIX_MU}, \code{BB_SIZE} and \code{GR}.
#'
#' @details A wrapper function of prior estimation
#' and bayNorm function.
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' #Return 3D array normalzied data:
#' bayNorm_3D<-bayNorm(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta,mode_version=FALSE)
#'
#' #Return 2D matrix normalized data:
#' bayNorm_2D<-bayNorm(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta
#' ,mode_version=TRUE)
#'
#' @references
#' Wenhao Tang, Francois Bertaux, Philipp Thomas,
#' Claire Stefanelli, Malika Saint, Samuel
#' Blaise Marguerat, Vahid Shahrezaei
#' bayNorm: Bayesian gene expression recovery,
#' imputation and normalisation for single cell RNA-sequencing data
#' bioRxiv 384586; doi: https://doi.org/10.1101/384586
#'
#'
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#'
#' @export
#'
bayNorm <- function(
    Data, BETA_vec, Conditions = NULL,
    UMI_sffl = NULL,Prior_type = NULL,
    mode_version = FALSE,
    mean_version=FALSE,
    S = 20,
    parallel = TRUE,NCores = 5,
    FIX_MU = TRUE, GR = FALSE,
    BB_SIZE = TRUE, verbose = TRUE) {

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


    input_params<-list(BETA_vec=BETA_vec,
                       Conditions=Conditions,
                       UMI_sffl=UMI_sffl,
                       Prior_type=Prior_type,
                       FIX_MU=FIX_MU,
                       BB_SIZE=BB_SIZE,
                       GR=GR)

    # Adapted from SCnorm
    if (methods::is(Data, "SummarizedExperiment") | methods::is(Data, "SingleCellExperiment")) {
        Data <- methods::as(Data, "SummarizedExperiment")
        if (is.null(SummarizedExperiment::assayNames(Data))
            || SummarizedExperiment::assayNames(Data)[1] !=
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
    if (!(methods::is(Data, "SummarizedExperiment")) &
        !(methods::is(Data, "SingleCellExperiment"))) {
        Data <- data.matrix(Data)
    }



    if (is.null(BETA_vec)) {
        BETA_vec <- colSums(Data)/mean(colSums(Data)) * 0.06
    }
    if (length(which(BETA_vec >= 1)) > 0) {
        BETA_vec[BETA_vec >= 1] = max(BETA_vec[BETA_vec < 1])
    }
    if (length(which(BETA_vec <= 0)) > 0) {
        BETA_vec[BETA_vec <= 0] = min(BETA_vec[BETA_vec > 0])
    }

    # Some pre-checkings:
    if (!is(Data,"matrix")) {
        stop("Input data should be of class matrix")
    }
    if (sum(duplicated(rownames(Data))) > 0) {
        warning("There are duplicated row names in Data")
    }
    if (sum(duplicated(colnames(Data))) > 0) {
        warning("There are duplicated column names in Data")
    }


    if (min(BETA_vec) <= 0 | max(BETA_vec) >= 1) {
        stop("The range of BETA must be within (0,1).")
    }
    if (ncol(Data) != length(BETA_vec)) {
        stop(
            "The number of cells(columns) in Data
        is not consistent with the number of elements
        in BETA_vec.")
    }


    if (is.null(Conditions)) {

        if (is.null(UMI_sffl)) {
            # Data_s<-Data
            Data_sr <- Data
        } else {
            Data_sr <- round(Data/UMI_sffl)
        }

        PRIORS = Prior_fun(
            Data = Data_sr,BETA_vec = BETA_vec,
            parallel = parallel,NCores = NCores,
            FIX_MU = FIX_MU,GR = GR,BB_SIZE = BB_SIZE,
            verbose = verbose)

        if (BB_SIZE) {
            MU_input = PRIORS$MME_prior$MME_MU
            SIZE_input = PRIORS$MME_SIZE_adjust
        } else {
            MU_input = PRIORS$MME_prior$MME_MU
            SIZE_input = PRIORS$MME_prior$MME_SIZE
        }


        Bay_out <- myFunc(
            Data = Data_sr,
            BETA_vec = BETA_vec,
            size = SIZE_input,
            mu = MU_input, S = S,
            thres = max(Data_sr)*2)
        rownames(Bay_out) <- rownames(Data)
        colnames(Bay_out) <- colnames(Data)

        return(list(
            Bay_out = Bay_out,
            PRIORS = PRIORS,
            input_params=input_params))




        if (verbose) {
            message("bayNorm has completed!")
        }


    } else {
        # multiple groups
        if (ncol(Data) != length(Conditions)) {
            stop("Number of columns in
                expression matrix must
                match length of conditions vector!")
        }
        if (is.null(Prior_type)) {
            warning("Prior_type needs to be specified
                    when Conditions are specified,
                    now Prior_type is set to be LL")
            Prior_type = "LL"
        }
        if (is.null(names(Conditions))) {
            names(Conditions) <- colnames(Data)
        }
        Levels <- unique(Conditions)



        DataList <- lapply(seq_along(Levels), function(x) {
            Data[, which(Conditions == Levels[x])]
        })
        BETAList <- lapply(seq_along(Levels), function(x) {
            BETA_vec[which(Conditions == Levels[x])]
        })


        if (is.null(UMI_sffl)) {
            # UMI
            DataList_sr <- lapply(seq_along(Levels), function(x) {
                Data[, which(Conditions == Levels[x])]
            })
        } else {
            # non-UMI
            DataList_sr <- lapply(seq_along(Levels), function(x) {
                round(Data[, which(Conditions == Levels[x])]/UMI_sffl[x])
            })
        }


        if (Prior_type == "LL") {
            PRIORS_LIST <- list()
            for (i in seq_len(Levels)) {
                PRIORS_LIST[[i]] <- Prior_fun(
                    Data = DataList_sr[[i]],
                    BETA_vec = BETAList[[i]],
                    parallel = parallel,
                    NCores = NCores,
                    FIX_MU = FIX_MU,
                    GR = GR, BB_SIZE = BB_SIZE,
                    verbose = verbose)
            }
        } else if (Prior_type == "GG") {
            PROPRS_TEMP <- Prior_fun(
                Data = do.call(cbind, DataList_sr),
                BETA_vec = do.call(c, BETAList),
                parallel = parallel,NCores = NCores,
                FIX_MU = FIX_MU, GR = GR,
                BB_SIZE = BB_SIZE,
                verbose = verbose)

            PRIORS_LIST <- list()
            for (i in seq_len(Levels)) {
                PRIORS_LIST[[i]] <- PROPRS_TEMP

            }

        }

        names(PRIORS_LIST) <- paste("Group", Levels)
        Bay_out_list <- list()
        for (i in seq_len(Levels)) {

            if (BB_SIZE) {
                MU_input = PRIORS_LIST[[i]]$MME_prior$MME_MU
                SIZE_input = PRIORS_LIST[[i]]$MME_SIZE_adjust
            } else {
                MU_input = PRIORS_LIST[[i]]$MME_prior$MME_MU
                SIZE_input = PRIORS_LIST[[i]]$MME_prior$MME_SIZE
            }


            Bay_out_list[[i]] <- myFunc(
                Data = DataList_sr[[i]],
                BETA_vec = BETAList[[i]],
                size = SIZE_input,
                mu = MU_input,
                S = S,
                thres = max(Data)*2)

            rownames(Bay_out_list[[i]]) <- rownames(DataList[[i]])
            colnames(Bay_out_list[[i]]) <- colnames(DataList[[i]])
        }
        names(Bay_out_list) <- paste("Group", Levels)

        return(list(Bay_out_list = Bay_out_list,
                    PRIORS_LIST = PRIORS_LIST,
                    input_params=input_params))


    }  # end for multiple groups

    if (verbose) {
        message("bayNorm has completed!")
    }

}

#' @title bayNorm with estimated parameters as input
#'
#' @description This is a supplementary wrapper function
#' for bayNorm. It is useful if one has already estimated
#' prior parameters and wants to simulate 2D or 3D
#' normalized output using the same prior estimates.
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells). \code{Data}
#' can be of class \code{SummarizedExperiment} (the
#' assays slot contains the expression matrix and
#' is named "Counts") or just matrix.
#' @param  PRIORS A list of estimated prior parameters
#' obtained from bayNorm.
#' @param  input_params A list of input parameters
#' which have been used: \code{BETA_vec}, \code{Conditions},
#' \code{UMI_sffl}, \code{Prior_type},
#' \code{FIX_MU}, \code{BB_SIZE} and \code{GR}.
#' @param  mode_version If TRUE, bayNorm return mode
#' version normalized data which is of 2D matrix
#' instead of 3D array. Default is FALSE.
#' @param  mean_version If TRUE, bayNorm return mean version
#' normalized data which is of 2D matrix instead of 3D array.
#' Default is FALSE.
#' @param S The number of samples you would like to
#' generate from estimated posterior distribution
#' (The third dimension of 3D array). Default is 20.
#' S needs to be specified if \code{mode_version}=FALSE.
#' @param  parallel If TRUE, 5 cores will be used for
#' parallelization.
#' @param  NCores number of cores to use, default is 5.
#' This will be used to set up a parallel environment
#' using either MulticoreParam (Linux, Mac) or
#' SnowParam (Windows) with NCores using the package
#' BiocParallel.
#' @param  BB_SIZE If TRUE (default), use adjusted size for
#' normalization. The adjusted size is obtained by adjusting
#' MME estimated size by a factor. The factor is
#' calculated based on both MME estimated size and BB
#' estimated size. If FALSE, use MME estimated SIZE.
#' @param verbose print out status messages. Default is TRUE.
#' @return  List containing 3D arrays of normalized
#' expression (if \code{mode_version}=FALSE) or 2D matrix
#' of normalized expression (if \code{mode_version}=TRUE
#' or \code{mean_version}=TRUE),
#' a list contains estimated priors and a list contains
#' input parameters used: \code{BETA_vec},
#' \code{Conditions} (if specified),
#' \code{UMI_sffl} (if specified), \code{Prior_type},
#' \code{FIX_MU}, \code{BB_SIZE} and \code{GR}.
#'
#' @details If you have run bayNorm before and obtained a
#' list of estimated prior parameters, then you may not want
#' to run parameter estimation again. You can just use
#' previous estimated parameters for obtaining 3D or
#' 2D normalized data.
#'
#' @examples
#' data('EXAMPLE_DATA_list')
#' #Return 3D array normalzied data:
#' bayNorm_3D<-bayNorm(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta
#' ,mode_version=FALSE)
#'
#' #Now if you want to generate 2D matrix using the same prior
#' #estimates as generated before:
#' bayNorm_2D<-bayNorm_sup(Data=EXAMPLE_DATA_list$inputdata
#' ,PRIORS=bayNorm_3D$PRIORS,
#' input_params = bayNorm_3D$input_params
#' ,mode_version=TRUE)
#'
#' @references
#' Wenhao Tang, Francois Bertaux, Philipp Thomas,
#' Claire Stefanelli, Malika Saint, Samuel
#' Blaise Marguerat, Vahid Shahrezaei
#' bayNorm: Bayesian gene expression recovery,
#' imputation and normalisation for single cell RNA-sequencing data
#' bioRxiv 384586; doi: https://doi.org/10.1101/384586
#'
#' @import parallel
#' @import foreach
#' @import SingleCellExperiment
#' @import doSNOW
#' @importFrom SummarizedExperiment SummarizedExperiment
#' assayNames assays colData
#'
#' @export
#'
bayNorm_sup <- function(
    Data,
    PRIORS = NULL,
    input_params=NULL,
    mode_version = FALSE,
    mean_version=FALSE,
    S = 20,
    parallel = TRUE, NCores = 5,
    BB_SIZE = TRUE, verbose = TRUE) {

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



    Conditions=input_params$Conditions
    UMI_sffl=input_params$UMI_sffl
    BETA_vec=input_params$BETA_vec

    if(!input_params$BB_SIZE & BB_SIZE){
        stop("Previous priors does not contain Adjusted MME size.
             Try to run bayNorm with BB_SIZE=TRUE.")
    }



    if (methods::is(Data, "SummarizedExperiment")
        | methods::is(Data, "SingleCellExperiment")) {
        Data <- methods::as(Data, "SummarizedExperiment")

        if (
            is.null(
                SummarizedExperiment::assayNames(Data)
            )
            || SummarizedExperiment::assayNames(Data)[1] !=
            "Counts") {
            message("Renaming the
                    firstelement in
                    assays(Data) to 'Counts'")
            SummarizedExperiment::assayNames(Data)[1] <- "Counts"

            if (is.null(colnames(
                SummarizedExperiment::assays(Data)[["Counts"]]))) {
                stop("Must supply sample/cell names!")
            }

        }
        Data <- SummarizedExperiment::assays(Data)[["Counts"]]
    }

    if (!(methods::is(Data, "SummarizedExperiment"))
        & !(methods::is(Data, "SingleCellExperiment"))) {
        Data <- data.matrix(Data)
    }




    if (is.null(Conditions)) {


        if (is.null(UMI_sffl)) {
            # Data_s<-Data
            Data_sr <- Data
        } else {
            Data_sr <- round(Data/UMI_sffl)
        }
        if (BB_SIZE) {
            MU_input = PRIORS$MME_prior$MME_MU
            SIZE_input = PRIORS$MME_SIZE_adjust
        } else {
            MU_input = PRIORS$MME_prior$MME_MU
            SIZE_input = PRIORS$MME_prior$MME_SIZE
        }

        Bay_out <- myFunc(
            Data = Data_sr,
            BETA_vec = BETA_vec,
            size = SIZE_input,
            mu = MU_input, S = S,
            thres = max(Data_sr)*2)
        rownames(Bay_out) <- rownames(Data)
        colnames(Bay_out) <- colnames(Data)

        return(list(
            Bay_out = Bay_out,
            PRIORS = PRIORS,
            input_params=input_params))



        if (verbose) {
            message("bayNorm has completed!")
        }


    } else {
        # multiple groups

        if (ncol(Data) != length(Conditions)) {
            stop("Number of columns in
                expression matrix must match
                length of conditions vector!")
        }
        if (is.null(names(Conditions))) {
            names(Conditions) <- colnames(Data)
        }
        Levels <- unique(Conditions)

        DataList <- lapply(seq_along(Levels), function(x) {
            Data[, which(Conditions == Levels[x])]
        })
        BETAList <- lapply(seq_along(Levels), function(x) {
            BETA_vec[which(Conditions == Levels[x])]
        })

        if (is.null(UMI_sffl)) {
            # UMI
            DataList_sr <- lapply(seq_along(Levels), function(x) {
                Data[, which(Conditions == Levels[x])]
            })
        } else {
            # non-UMI
            DataList_sr <- lapply(seq_along(Levels), function(x) {
                round(Data[, which(Conditions == Levels[x])]/UMI_sffl[x])
            })
        }

        ## use existing PRIORS
        PRIORS_LIST <- PRIORS
        Bay_out_list <- list()
        for (i in seq_len(Levels)) {
            if (BB_SIZE) {
                MU_input = PRIORS_LIST[[i]]$MME_prior$MME_MU
                SIZE_input = PRIORS_LIST[[i]]$MME_SIZE_adjust
            } else {
                MU_input = PRIORS_LIST[[i]]$MME_prior$MME_MU
                SIZE_input = PRIORS_LIST[[i]]$MME_prior$MME_SIZE

            }
            Bay_out_list[[i]] <- myFunc(
                Data = DataList_sr[[i]],
                BETA_vec = BETAList[[i]],
                size = SIZE_input, mu = MU_input,
                S = S, thres = max(Data) * 2)

            rownames(Bay_out_list[[i]]) <- rownames(DataList[[i]])
            colnames(Bay_out_list[[i]]) <- colnames(DataList[[i]])
        }
        names(Bay_out_list) <- paste("Group", Levels)
        names(BETAList)<-paste("Group", Levels)

        return(list(
            Bay_out_list = Bay_out_list,
            PRIORS_LIST = PRIORS_LIST,
            input_params=input_params))
    }  # end for multiple groups

    if (verbose) {
        message("bayNorm has completed!")
    }
}
