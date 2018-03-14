

#' @title  A wrapper function of prior estimation and bayNorm function
#'
#' @description   Input raw data and a vector of capture
#' efficiencies of cells. You can also specify the condition
#' of cells for normalizing multiple groups of cells separately.
#' @param Data A matrix of single-cell expression where rows
#' are genes and columns are samples (cells).
#' This object should be of class matrix rather than data.frame.
#' @param  BETA_vec A vector of capture efficiencies of cells.
#' If it is null, library size normalized to 0.06 will be used
#' as the input BETA_vec. BETA_vec less and equal to 0 or
#' greater and equal to 1 will be replaced by the minimum
#' and maximum of the BETA_vec which range between (0,1) respectively.
#' @param  Conditions vector of condition labels, this should correspond to the columns of the Data. Default is NULL, which assumes that all cells belong to the same group.
#' @param UMI_sffl (scaling factors for non UMI based data:
#' divide Data by UMI_sffl) Only needed when the input data
#' is non UMI based. If non-null and Conditions is non-null,
#' then UMI_sffl should be a vector of length equal to the
#' number of groups. Default is set to be NULL.
#' @param  Prior_type Default is NULL. If \code{Conditions}
#' is NULL, priors are estimated based on all cells. If
#' \code{Conditions} is not NULL: if \code{Prior_type} is \code{LL}, priors are estimated within each group respectively.
#' If \code{Prior_type} is \code{GG}, priors are estimated
#' based on cells from all groups. Basically, \code{LL} is suitable for DE detection. \code{GG} is prefered if there is
#' a prior knowledge about the data such that there should
#' not exist biological variation between groups.
#' @param  mode_version If TRUE, bayNorm return mode version normalized data which is of 2D matrix instead of 3D array. Default is FALSE.
#' @param S The number of samples you would like to generate from estimated posterior distribution (The third dimension of 3D array). Default is 20. S needs to be specified if \code{mode_version}=FALSE.
#' @param  parallel If TRUE, 5 cores will be used for parallelization.
#' @param  NCores number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  FIX_MU Whether fix mu when estimating parameters by maximizing marginal distribution. If TRUE, then 1D optimization, otherwise 2D optimization (slow).
#' @param  GR If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @param  BB_SIZE If TRUE, estimate BB size, and then use it for adjusting MME SIZE. Use the adjusted MME size for bayNorm. Default is TRUE.
#' @param verbose print out status messages. Default is TRUE.
#' @return  List containing 3D arrays of normalized expression (if \code{mode_version}=FALSE) or 2D matrix of normalized expression (if \code{mode_version}=TRUE), estimated parameters and input \code{BETA_vec}.
#'
#' @details A wrapper function of prior estimation and bayNorm function.
#'
#' @examples
#' data("EXAMPLE_DATA_list")
#' \dontrun{
#' #Return 3D array normalzied data:
#' bayNorm_3D<-bayNorm(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta,mode_version=F)
#'
#' #Return 2D matrix normalized data:
#' bayNorm_2D<-bayNorm(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta
#' ,mode_version=T)
#'
#' }
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
bayNorm<-function(Data,BETA_vec,Conditions=NULL,UMI_sffl=NULL,Prior_type=NULL,mode_version=FALSE,S=20,parallel=TRUE,NCores=5,FIX_MU=TRUE,GR=FALSE,BB_SIZE=TRUE,verbose=TRUE){

  if(is.null(BETA_vec)){
        BETA_vec<-colSums(Data)/mean(colSums(Data))*0.06
  }

  if(length(which(BETA_vec>=1))>0){
    BETA_vec[BETA_vec>=1]=max(BETA_vec[BETA_vec<1])
  }
  if(length(which(BETA_vec<=0))>0){
    BETA_vec[BETA_vec<=0]=min(BETA_vec[BETA_vec>0])
  }





#Some pre-checkings:
  if(class(Data)!='matrix'){stop("Input data should be of class matrix")}
  if(sum(duplicated(rownames(Data)))>0){warning("There are duplicated row names in Data")}
  if(sum(duplicated(colnames(Data)))>0){warning("There are duplicated column names in Data")}


  if(min(BETA_vec)<=0 | max(BETA_vec)>=1){stop("The range of BETA must be within (0,1).")}
  if(ncol(Data)!=length(BETA_vec)){stop("The number of cells (columns) in Data is not consistent with the number of elements in BETA_vec.")}


  if(is.null(Conditions)){

    if(is.null(UMI_sffl)){
      #Data_s<-Data
      Data_sr<-Data
    }else{
      Data_sr<-round(Data/UMI_sffl)
      }

    PRIORS=Prior_fun(Data=Data_sr,BETA_vec=BETA_vec,parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)
    if(BB_SIZE){
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_SIZE_adjust
    }else{
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_prior$MME_SIZE
    }

    if(!mode_version){
    Bay_array<-Main_Bay(Data=Data_sr,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data_sr)*2)
    rownames(Bay_array)<-rownames(Data)
    colnames(Bay_array)<-colnames(Data)
    return(list(Bay_array=Bay_array,PRIORS=PRIORS,BETA=BETA_vec))
    }else{ #mode
      Bay_mat<-Main_mode_Bay(Data=Data_sr,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data_sr)*2)
      rownames(Bay_mat)<-rownames(Data)
      colnames(Bay_mat)<-colnames(Data)
      return(list(Bay_mat=Bay_mat,PRIORS=PRIORS,BETA=BETA_vec))
    }

    if(verbose){
      message("bayNorm has completed!")
    }


  } else{# multiple groups




   if (ncol(Data) != length(Conditions)) {stop("Number of columns in
      expression matrix must match length of conditions vector!")}
   if(is.null(Prior_type)){warning("Prior_type needs to be specified when Conditions are specified, now Prior_type is set to be LL")
     Prior_type='LL'
     }
   if(is.null(names(Conditions))) {names(Conditions) <- colnames(Data)}
   Levels <- unique(Conditions)


   if(is.null(UMI_sffl)){#UMI
     DataList<- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})
     DataList_sr <- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})

     BETAList <- lapply(seq_along(Levels), function(x){BETA_vec[which(Conditions == Levels[x])]})
   }else{#non-UMI

     DataList<- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})
     DataList_sr <- lapply(seq_along(Levels), function(x){round(Data[,which(Conditions == Levels[x])]/UMI_sffl[x])})
     BETAList <- lapply(seq_along(Levels), function(x){BETA_vec[which(Conditions == Levels[x])]})
   }


   if(Prior_type=='LL'){
     PRIORS_LIST<-list()
     for(i in 1:length(Levels)){
       PRIORS_LIST[[i]]<-Prior_fun(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)

     }
   }else if (Prior_type=='GG'){
     PROPRS_TEMP<-Prior_fun(Data=do.call(cbind,DataList_sr),BETA_vec=do.call(c,BETAList),parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)

     PRIORS_LIST<-list()
     for(i in 1:length(Levels)){
       PRIORS_LIST[[i]]<-PROPRS_TEMP

     }

   }




   names(PRIORS_LIST)<-paste('Group',Levels)


   if(!mode_version){
   Bay_array_list<-list()
   for(i in 1:length(Levels)){

     if(BB_SIZE){
       MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
       SIZE_input=PRIORS_LIST[[i]]$MME_SIZE_adjust
     }else{
       MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
       SIZE_input=PRIORS_LIST[[i]]$MME_prior$MME_SIZE
     }


     Bay_array_list[[i]]<-Main_Bay(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2)

     rownames(Bay_array_list[[i]])<-rownames(DataList[[i]])
     colnames(Bay_array_list[[i]])<-colnames(DataList[[i]])
   }
   names(Bay_array_list)<-paste('Group',Levels)
   return(list(Bay_array_list=Bay_array_list,PRIORS_LIST=PRIORS_LIST,BETA=BETAList))
   }else{#mode

     Bay_mat_list<-list()
     for(i in 1:length(Levels)){

       if(BB_SIZE){
         MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
         SIZE_input=PRIORS_LIST[[i]]$MME_SIZE_adjust
       }else{
         MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
         SIZE_input=PRIORS_LIST[[i]]$MME_prior$MME_SIZE
       }

       Bay_mat_list[[i]]<-Main_mode_Bay(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2)

       rownames(Bay_mat_list[[i]])<-rownames(DataList[[i]])
       colnames(Bay_mat_list[[i]])<-colnames(DataList[[i]])
     }
     names(Bay_mat_list)<-paste('Group',Levels)
     return(list(Bay_mat_list=Bay_mat_list,PRIORS_LIST=PRIORS_LIST,BETA=BETAList))

   } #end of mode for multiple groups
 }# end for multiple groups

  if(verbose){
    message("bayNorm has completed!")
  }

}


#' @title bayNorm with estimated parameters as input
#'
#' @description This is a supplementary function for \code{bayNorm}. It is useful if you have already run \code{bayNorm} before and try to simulate 3D or 3D matrix using the same prior estimates.
#' @param Data A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec A vector of capture efficiencies of cells.
#' @param  PRIORS A list of estimated prior parameters obtained from bayNorm.
#' @param  Conditions vector of condition labels, this should correspond to the columns of the Data. Default is NULL, which assumes that all cells belong to the same group.
#' @param UMI_sffl (scaling factors for non UMI based data: divide Data by UMI_sffl) Only needed when the input data is non UMI based. If non-null and Conditions is non-null, then UMI_sffl should be a vector of length equal to the number of groups. Default is set to be NULL.
#' @param  mode_version If TRUE, bayNorm return mode version normalized data which is of 2D matrix instead of 3D array. Default is FALSE.
#' @param S The number of samples you would like to generate from estimated posterior distribution (The third dimension of 3D array). Default is 20. S needs to be specified if \code{mode_version}=FALSE.
#' @param  parallel If TRUE, 5 cores will be used for parallelization.
#' @param  NCores number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  BB_SIZE If TRUE, use adjusted size for normalization. The adjusted size is obtained by adjusting MME estimated size by a factor. The factor is calculated based on both MME estimated size and BB estimated size.
#' @param verbose print out status messages. Default is TRUE.
#' @return  List containing 3D arrays of normalized expression (if \code{mode_version}=FALSE) or 2D matrix of normalized expression (if \code{mode_version}=TRUE), estimated parameters and input \code{BETA_vec}.
#'
#' @details If you have run bayNorm before and obtained a list of estimated prior parameters, then you may not want to run parameter estimation again. You can just use previous estimated parameters for obtaining 3D or 2D normalized data.
#'
#' @examples
#' data("EXAMPLE_DATA_list")
#' \dontrun{
#' #Return 3D array normalzied data:
#' bayNorm_3D<-bayNorm(Data=EXAMPLE_DATA_list$inputdata,
#' BETA_vec = EXAMPLE_DATA_list$inputbeta
#' ,mode_version=F)
#'
#' #Now if you want to generate 2D matrix using the same prior
#' #estimates as generated before:
#' bayNorm_2D<-bayNorm_p(Data=EXAMPLE_DATA_list$inputdata
#' ,BETA_vec= bayNorm_3D$BETA,PRIORS=bayNorm_3D$PRIORS_LIST
#' ,mode_version=T)
#'
#' #If previous bayNorm was applied for normalizing multiple
#' #groups of cells (is.null(Origin_Conditions)=T), then:
#' inputbeta2<-unlist(bayNorm_3D$BETA)
#' bayNorm_2D<-bayNorm_p(Data=inputdata,BETA_vec = inputbeta2
#' ,PRIORS=bayNorm_3D$PRIORS_LIST,mode_version=T,Conditions
#' =Origin_Conditions)
#'
#' #You can also generate 3D array using the same prior estimates
#' #as generated before.
#' }
#'
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
bayNorm_sup<-function(Data,BETA_vec,PRIORS=NULL,Conditions=NULL,UMI_sffl=NULL,mode_version=FALSE,S=20,parallel=TRUE,NCores=5,BB_SIZE=TRUE,verbose=TRUE){


  if(is.null(Conditions)){
    if(is.null(UMI_sffl)){
      #Data_s<-Data
      Data_sr<-Data
    }else{
      Data_sr<-round(Data/UMI_sffl)
    }

    #PRIORS=Prior_fun(Data=Data_sr,BETA_vec=BETA_vec,parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)
    if(BB_SIZE){
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_SIZE_adjust
    }else{
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_prior$MME_SIZE
    }

    if(!mode_version){
      Bay_array<-Main_Bay(Data=Data_sr,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data_sr)*2)
      rownames(Bay_array)<-rownames(Data)
      colnames(Bay_array)<-colnames(Data)
      return(list(Bay_array=Bay_array,PRIORS=PRIORS,BETA=BETA_vec))
    }else{ #mode
      Bay_mat<-Main_mode_Bay(Data=Data_sr,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data_sr)*2)
      rownames(Bay_mat)<-rownames(Data)
      colnames(Bay_mat)<-colnames(Data)
      return(list(Bay_mat=Bay_mat,PRIORS=PRIORS,BETA=BETA_vec))
    }

    if(verbose){
      message("bayNorm has completed!")
    }


  } else{# multiple groups

    if (ncol(Data) != length(Conditions)) {stop("Number of columns in
                                                expression matrix must match length of conditions vector!")}
    # if(is.null(Prior_type)){warning("Prior_type needs to be specified when Conditions are specified, now Prior_type is set to be LL")
    #   Prior_type='LL'
    # }
    if(is.null(names(Conditions))) {names(Conditions) <- colnames(Data)}
    Levels <- unique(Conditions)

    if(is.null(UMI_sffl)){#UMI
      DataList<- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})
      DataList_sr <- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})

      BETAList <- lapply(seq_along(Levels), function(x){BETA_vec[which(Conditions == Levels[x])]})
    }else{#non-UMI

      DataList<- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})
      DataList_sr <- lapply(seq_along(Levels), function(x){round(Data[,which(Conditions == Levels[x])]/UMI_sffl[x])})
      BETAList <- lapply(seq_along(Levels), function(x){BETA_vec[which(Conditions == Levels[x])]})
    }


    # if(Prior_type=='LL'){
    #   PRIORS_LIST<-list()
    #   for(i in 1:length(Levels)){
    #     PRIORS_LIST[[i]]<-Prior_fun(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)
    #
    #   }
    # }else if (Prior_type=='GG'){
    #   PROPRS_TEMP<-Prior_fun(Data=do.call(cbind,DataList_sr),BETA_vec=do.call(c,BETAList),parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)
    #
    #   PRIORS_LIST<-list()
    #   for(i in 1:length(Levels)){
    #     PRIORS_LIST[[i]]<-PROPRS_TEMP
    #
    #   }
    #
    # }
    #names(PRIORS_LIST)<-paste('Group',Levels)
    PRIORS_LIST<-PRIORS

    if(!mode_version){
      Bay_array_list<-list()
      for(i in 1:length(Levels)){

        if(BB_SIZE){
          MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
          SIZE_input=PRIORS_LIST[[i]]$MME_SIZE_adjust
        }else{
          MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
          SIZE_input=PRIORS_LIST[[i]]$MME_prior$MME_SIZE
        }
        Bay_array_list[[i]]<-Main_Bay(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2)

        rownames(Bay_array_list[[i]])<-rownames(DataList[[i]])
        colnames(Bay_array_list[[i]])<-colnames(DataList[[i]])
      }
      names(Bay_array_list)<-paste('Group',Levels)
      return(list(Bay_array_list=Bay_array_list,PRIORS_LIST=PRIORS_LIST,BETA=BETAList))
    }else{#mode

      Bay_mat_list<-list()
      for(i in 1:length(Levels)){

        if(BB_SIZE){
          MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
          SIZE_input=PRIORS_LIST[[i]]$MME_SIZE_adjust
        }else{
          MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
          SIZE_input=PRIORS_LIST[[i]]$MME_prior$MME_SIZE
        }

        Bay_mat_list[[i]]<-Main_mode_Bay(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2)

        rownames(Bay_mat_list[[i]])<-rownames(DataList[[i]])
        colnames(Bay_mat_list[[i]])<-colnames(DataList[[i]])
      }
      names(Bay_mat_list)<-paste('Group',Levels)
      return(list(Bay_mat_list=Bay_mat_list,PRIORS_LIST=PRIORS_LIST,BETA=BETAList))

    } #end of mode for multiple groups
  }# end for multiple groups

  if(verbose){
    message("bayNorm has completed!")
  }
}
