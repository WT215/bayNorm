

#' A wrapper function of prior estimation and bayNorm function
#'
#' Input raw data and a vector of capture efficiencies of cells. You can
#' also need to specify the condition of cells.
#' @param Data A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec A vector of capture efficiencies of cells.
#' @param  mode_version If TRUE, bayNorm return mode version normalized data which is of 2D matrix instead of 3D array. Defaut is FALSE.
#' @param S The number of samples you would like to generate from estimated posterior distribution (The third dimension of 3D array). Default is 20. S needs to be specified if \code{mode_version}=FALSE.
#' @param  parallel If TRUE, 5 cores will be used for parallelization.
#' @param  NCores number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  FIX_MU If TRUE, then 1D optimization, otherwise 2D optimization (slow).
#' @param  GR If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @param  Conditions vector of condition labels, this should correspond to the columns of the Data. Default is NULL, which assumes that all cells belong to the same group.
#' @param  BB_SIZE If TRUE, estimate BB size, and then use it for adjusting MME SIZE. Use the adjusted MME size for bayNorm. Defaut is TRUE.
#' @param UMI_sffl (scaling factors for full-length data: divide Data by UMI_sffl) Only needed when the input data is full-length protocol based. If non-null and Conditions is non-null, then flsf should be a vector of length equal to the number of groups. Defaut is set to be NULL.
#' @param  Prior_type Default is NULL. If \code{Condition} is NULL, priors are estimated based on all cells. If \code{Condition} is not NULL: if \code{Prior_type} is \code{LL}, priors are estimated within each group respectively. If \code{Prior_type} is \code{GG}, priors are estimated based on cells from all groups. Basically, \code{LL} is suitable for DE detection. \code{GG} is prefered if there is a prior knowledge about the data such that there should not exist biological variation between groups.
#' @param verbose print out status messages. Default is TRUE.
#' @return  List containing 3D arrays of normalized expression (if \code{mode_version}=FALSE) or 2D matrix of normalized expression (if \code{mode_version}=TRUE) and estimated parameters.
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
bayNorm<-function(Data,BETA_vec,S=20,parallel=T,NCores=5,FIX_MU=T,GR=F,Conditions=NULL,BB_SIZE=T,mode_version=F,UMI_sffl=NULL,Prior_type=NULL,verbose=T){

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
    Bay_array<-Main_Bay(Data=Data_sr,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data_sr)*2,Mean_depth=1000000)
    rownames(Bay_array)<-rownames(Data)
    colnames(Bay_array)<-colnames(Data)
    return(list(Bay_array=Bay_array,PRIORS=PRIORS))
    }else{ #mode
      Bay_mat<-Main_mode_Bay(Data=Data_sr,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data_sr)*2,Mean_depth=1000000)
      rownames(Bay_mat)<-rownames(Data)
      colnames(Bay_mat)<-colnames(Data)
      return(list(Bay_mat=Bay_mat,PRIORS=PRIORS))
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


     Bay_array_list[[i]]<-Main_Bay(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000)

     rownames(Bay_array_list[[i]])<-rownames(DataList[[i]])
     colnames(Bay_array_list[[i]])<-colnames(DataList[[i]])
   }
   names(Bay_array_list)<-paste('Group',Levels)
   return(list(Bay_array_list=Bay_array_list,PRIORS_LIST=PRIORS_LIST))
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

       Bay_mat_list[[i]]<-Main_mode_Bay(Data=DataList_sr[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000)

       rownames(Bay_mat_list[[i]])<-rownames(DataList[[i]])
       colnames(Bay_mat_list[[i]])<-colnames(DataList[[i]])
     }
     names(Bay_mat_list)<-paste('Group',Levels)
     return(list(Bay_mat_list=Bay_mat_list,PRIORS_LIST=PRIORS_LIST))

   } #end of mode for multiple groups
 }# end for multiple groups

  if(verbose){
    message("bayNorm has completed!")
  }

}



