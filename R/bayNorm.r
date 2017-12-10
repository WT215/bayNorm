

#' A wrapper function of prior estimation and bayNorm function
#'
#' Input raw data and a vector of capture efficiencies of cells. You may
#' also need to specify the condition of cells.
#' @param Data: A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec: A vector of capture efficiencies of cells.
#' @param S: The number of samples you would like to genrate. Default is 20.
#' @param  Para: If TRUE, 5 cores will be used for parallelization.
#' @param  NCores: number of cores to use, default is 5.
#' @param  FIX_MU: If TRUE, then 1D optimization, otherwise 2D optimization (slow).
#' @param  GR: If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @param  Conditions: Default is NULL.
#' @param  BB_SIZE:If TRUE, estimate BB size, and then use it for adjusting MME SIZE. Use the adjusted MME size for bayNorm. Defaut is TRUE.
#' @return  A list of objects.
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
bayNorm<-function(Data,BETA_vec,S=20,Para=T,NCores=5,FIX_MU=T,GR=F,Conditions=NULL,BB_SIZE=T){
  if(is.null(Conditions)){
    PRIORS=Prior_fun(Data=Data,BETA_vec=BETA_vec,Para=Para,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE)
    if(BB_SIZE){
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_SIZE_adjust
    }else{
      MU_input=PRIORS$MME_prior$MME_MU
      SIZE_input=PRIORS$MME_prior$MME_SIZE
    }
    Bay_array<-Main_Bay(Table=Data,Beta_origin=BETA_vec,size=SIZE_input,M_ave_ori=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000,debug=F)
    rownames(Bay_array)<-rownames(Data)
    colnames(Bay_array)<-colnames(Data)


    return(list(Bay_array=Bay_array,PRIORS=PRIORS))
  }

  #multiple groups
 else{

   if (ncol(Data) != length(Conditions)) {stop("Number of columns in
      expression matrix must match length of conditions vector!")}
   if(is.null(names(Conditions))) {names(Conditions) <- colnames(Data)}
   Levels <- unique(Conditions)

   DataList <- lapply(seq_along(Levels), function(x){Data[,which(Conditions == Levels[x])]})
   BETAList <- lapply(seq_along(Levels), function(x){BETA_vec[which(Conditions == Levels[x])]})

   PRIORS_LIST<-list()
   Bay_array_list<-list()
   for(i in 1:length(Levels)){
     PRIORS_LIST[[i]]<-Prior_fun(Data=DataList[[i]],BETA_vec=BETAList[[i]],Para=Para,NCores=NCores,FIX_MU=FIX_MU,GR=GR,BB_SIZE=BB_SIZE)
     if(BB_SIZE){
       MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
       SIZE_input=PRIORS_LIST[[i]]$MME_SIZE_adjust
     }else{
       MU_input=PRIORS_LIST[[i]]$MME_prior$MME_MU
       SIZE_input=PRIORS_LIST[[i]]$MME_prior$MME_SIZE
     }

     Bay_array_list[[i]]<-Main_Bay(Table=DataList[[i]],Beta_origin=BETAList[[i]],size=SIZE_input,M_ave_ori=MU_input,S=S,thres=max(Data)*2,Mean_depth=1000000,debug=F)

     rownames(Bay_array_list[[i]])<-rownames(DataList[[i]])
     colnames(Bay_array_list[[i]])<-colnames(DataList[[i]])


   }
   names(PRIORS_LIST)<-paste('Group',Levels)
   names(Bay_array_list)<-paste('Group',Levels)


   return(list(Bay_array_list=Bay_array_list,PRIORS_LIST=PRIORS_LIST))

 }





}

