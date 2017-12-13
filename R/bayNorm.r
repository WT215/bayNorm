

#' A wrapper function of prior estimation and bayNorm function
#'
#' Input raw data and a vector of capture efficiencies of cells. You can
#' also need to specify the condition of cells.
#' @param Data: A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec: A vector of capture efficiencies of cells.
#' @param S: The number of samples you would like to genrate. Default is 20.
#' @param  parallel: If TRUE, 5 cores will be used for parallelization.
#' @param  NCores: number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  FIX_MU: If TRUE, then 1D optimization, otherwise 2D optimization (slow).
#' @param  GR: If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @param  Conditions: vector of condition labels, this should correspond to the columns of the Data. Default is NULL, which assumes that all cells belong to the same group.
#' @param  BB_SIZE: If TRUE, estimate BB size, and then use it for adjusting MME SIZE. Use the adjusted MME size for bayNorm. Defaut is TRUE.
#' @param  mode_version: If TRUE, bayNorm return mode version normalized data which is of 2D matrix instead of 3D array. Defaut is FALSE.
#' @return  A list of objects.
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
bayNorm<-function(Data,BETA_vec,S=20,parallel=T,NCores=5,FIX_MU=T,GR=F,Conditions=NULL,BB_SIZE=T,mode_version=F){

  if(is.null(Conditions)){
    PRIORS=Prior_fun(Data=Data,BETA_vec=BETA_vec,parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE)
    if(BB_SIZE){
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_SIZE_adjust
    }else{
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_prior$MME_SIZE
    }

    if(!mode_version){
    Bay_array<-Main_Bay(Data=Data,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000)
    rownames(Bay_array)<-rownames(Data)
    colnames(Bay_array)<-colnames(Data)
    return(list(Bay_array=Bay_array,PRIORS=PRIORS))
    }else{ #mode
      Bay_mat<-Main_mode_Bay(Data=Data,BETA_vec=BETA_vec,size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000)
      rownames(Bay_mat)<-rownames(Data)
      colnames(Bay_mat)<-colnames(Data)
      return(list(Bay_mat=Bay_mat,PRIORS=PRIORS))
    }


  } else{# multiple groups

   if (ncol(Data) != length(Conditions)) {stop("Number of columns in
      expression matrix must match length of conditions vector!")}
   if(is.null(names(Conditions))) {names(Conditions) <- colnames(Data)}
   Levels <- unique(Conditions)

   DataList <- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})
   BETAList <- lapply(seq_along(Levels), function(x){BETA_vec[which(Conditions == Levels[x])]})

   PRIORS_LIST<-list()
   for(i in 1:length(Levels)){
     PRIORS_LIST[[i]]<-Prior_fun(Data=DataList[[i]],BETA_vec=BETAList[[i]],parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE)
     if(BB_SIZE){
       MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
       SIZE_input=PRIORS_LIST[[i]]$MME_SIZE_adjust
     }else{
       MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
       SIZE_input=PRIORS_LIST[[i]]$MME_prior$MME_SIZE
     }

   }
   names(PRIORS_LIST)<-paste('Group',Levels)


   if(!mode_version){
   Bay_array_list<-list()
   for(i in 1:length(Levels)){
     Bay_array_list[[i]]<-Main_Bay(Data=DataList[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000,debug=F)

     rownames(Bay_array_list[[i]])<-rownames(DataList[[i]])
     colnames(Bay_array_list[[i]])<-colnames(DataList[[i]])
   }
   names(Bay_array_list)<-paste('Group',Levels)
   return(list(Bay_array_list=Bay_array_list,PRIORS_LIST=PRIORS_LIST))
   }else{#mode

     Bay_mat_list<-list()
     for(i in 1:length(Levels)){
       Bay_mat_list[[i]]<-Main_mode_Bay(Data=DataList[[i]],BETA_vec=BETAList[[i]],size=SIZE_input,mu=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000)

       rownames(Bay_mat_list[[i]])<-rownames(DataList[[i]])
       colnames(Bay_mat_list[[i]])<-colnames(DataList[[i]])
     }
     names(Bay_mat_list)<-paste('Group',Levels)
     return(list(Bay_mat_list=Bay_mat_list,PRIORS_LIST=PRIORS_LIST))

   }

 }

}



