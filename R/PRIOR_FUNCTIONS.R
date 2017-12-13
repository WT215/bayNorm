#' Adjust MME size
#'
#' @param BB_SIZE: size estimated from BB_Fun or BB_Fun_1D
#' @param  MME_MU: mu estimated from EstPrior.
#' @param  MME_SIZE: size estimated from EstPrior.
#' @return MME_SIZE_adjust: Adjusted MME_SIZE based on BB_size
#'
#'
#' @export
AdjustSIZE_fun<-function(BB_SIZE,MME_MU,MME_SIZE){
  fitind<-which(BB_SIZE< max(BB_SIZE) & BB_SIZE> min(BB_SIZE))
  lmfit<-lm(log(BB_SIZE)[fitind]~log(MME_SIZE)[fitind])
  MME_SIZE_adjust<-coef(lmfit)[1]+coef(lmfit)[2]*log(MME_SIZE)
  MME_SIZE_adjust<-exp(MME_SIZE_adjust)
  return(MME_SIZE_adjust)
}




#' Estimate capture efficiency for cells
#'
#' This function aims to select of subset of genes for estimating capture efficiency: BETA for bayNorm.
#'
#' @param Data: A matrix of single-cell expression where rows are
#' genes and columns are samples (cells). This object should be of
#' class matrix rather than data.frame.
#' @param MeanBETA: Mean capture efficiency of the scRNAseq data. This can be estimated via spike-ins.
#'
#' @return BETA: a vector of capture efficiencies, which is of length number of cells.
#' @return Selected_genes: a subset of genes that are used for estimating BETA.
#'
#'
#' @export
BetaFun<-function(Data,MeanBETA){
  Normcount<-t(t(Data)/colSums(Data))*mean(colSums(Data))
  means <- rowMeans(Normcount)
  lmeans <- log(means)
  med <- apply(log(Normcount+1),1,function(x){median(x)})
  mad <- apply(log(Normcount+1),1,function(x){mad(x)})
  bound <- med + 3 * mad
  maxlogGene<-apply(log(Normcount+1),1,max)
  ind<-which(maxlogGene<bound)
  dropout=apply(Data,1,function(x){length(which(x==0))/length(x)})
  Select_ind<-intersect(ind,which(dropout<0.35))
  Selected_genes<-rownames(Data)[Select_ind]

  temppp<-colSums(Data[Select_ind,])
  BETA<-temppp/mean(temppp)*MeanBETA
  names(BETA)<-colnames(Data)
  return(list(BETA=BETA,Selected_genes=Selected_genes))
}

#' Estimate size and mu for NB distribution for each gene using MME method
#'
#' Input raw data and return estimated size and mu for each gene.
#' @param Data: A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param verbose: print out status messages. Default is TRUE.
#' @return  Estimated size and mu for each gene.
#'
#' @import  fitdistrplus
#'
#' @export
EstPrior<-function(Data,verbose=T){
  CoefDat<-foreach(i=1:nrow(Data),.combine=rbind)%do%{
    qq<-fitdistrplus::fitdist(Data[i,],'nbinom',method='mme',keepdata=F)
    return(coef(qq))
  }
  rownames(CoefDat)<-rownames(Data)
  M_ave_ori<-CoefDat[,2]
  size_est<-CoefDat[,1]

  if(verbose){
    message("Priors estimation based on MME method has completed.")
  }

  return(list(MU=M_ave_ori,SIZE=size_est))
}
#' A wrapper function of EstPrior and AdjustSIZE_fun
#'
#' Input raw data and a vector of capture efficiencies of cells.
#' @param Data: A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec: A vector of capture efficiencies of cells.
#' @param  parallel: If TRUE, 5 cores will be used for parallelization. Defaut is TRUE.
#' @param  NCores: number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  FIX_MU: If TRUE, then 1D optimization, otherwise 2D optimization (slow). Defaut is TRUE.
#' @param  GR: If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @param  BB_SIZE: If TRUE, estimate BB size, and then use it for adjusting MME SIZE. Use the adjusted MME size for bayNorm. Defaut is TRUE.
#' @param verbose: Print out status messages. Default is TRUE.
#'
#' @details By defaut, this function will estimate mu and size for each gene using MME method. If \code{BB_size} is enable, spectral projected gradient method from BB package will be implemented to estimate "BB size" by maximizing marginal likelihood function. MME estimated size will be adjusted according to BB size. BB size itself will not be used in bayNorm this is because that in our simulation we found that MME estimated mu and size have more accurate relationship, but MME estimated size deviates from the true value. BB size is overall more close to the true size but it does not possess a reasonable relationship with either MME estimated mu or BB estimated mu.
#'
#' @examples
#' \dontrun{
#' Prior_fun(Data,BETA_vec,parallel=T,NCores=5,FIX_MU=T,GR=F,BB_SIZE=T,verbose=T)
#' }
#'
#' @return  A list of objects.
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'
Prior_fun<-function(Data,BETA_vec,parallel=T,NCores=5,FIX_MU=T,GR=F,BB_SIZE=T,verbose=T){

  normcount_N<-t(t(Data)/colSums(Data))*mean(colSums(Data)/BETA_vec)
  Priors<-EstPrior(normcount_N,verbose=verbose)
  M_ave_ori<-Priors$MU
  size_est<-Priors$SIZE
  size_est[is.na(size_est)]<-min(size_est[!is.na(size_est)])
  MME_prior<-cbind(M_ave_ori,size_est)
  MME_prior<-as.data.frame(MME_prior)
  rownames(MME_prior)<-rownames(Data)
  colnames(MME_prior)<-c("MME_MU","MME_SIZE")

  #BB_size<-BB_Fun_1D(Dat_mat=Data,BETA_vec=BETA_vec,INITIAL_MU_vec=MME_prior$MME_MU,INITIAL_SIZE_vec=MME_prior$MME_SIZE,SIZE_lower=min(MME_prior$MME_SIZE),SIZE_upper=ceiling(max(MME_prior$MME_SIZE)),parallel=parallel,NCores = NCores)

if(BB_SIZE){
  if(verbose){
    message("Begin to estimate size for each gene by maximizing the marginal distribution. The optimization method is spectral projected gradient implemented in BB package. The MME estimated size will be adjusted based on this.")
  }

  if(FIX_MU){
    BB_size<-BB_Fun(Data,BETA_vec,INITIAL_MU_vec=MME_prior$MME_MU,INITIAL_SIZE_vec=MME_prior$MME_SIZE,MU_lower=min(MME_prior$MME_MU),MU_upper=max(MME_prior$MME_MU),SIZE_lower=min(MME_prior$MME_SIZE),SIZE_upper=ceiling(max(MME_prior$MME_SIZE)),parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR)
    BB_prior<-cbind(MME_prior$MME_MU,BB_size)
    rownames(BB_prior)<-rownames(Data)
    colnames(BB_prior)<-c("MME_MU","BB_SIZE")

    MME_SIZE_adjust<-AdjustSIZE_fun(BB_prior[,2],MME_prior$MME_MU,MME_prior$MME_SIZE)
  }else{
    BB_size<-BB_Fun(Data,BETA_vec,INITIAL_MU_vec=MME_prior$MME_MU,INITIAL_SIZE_vec=MME_prior$MME_SIZE,MU_lower=min(MME_prior$MME_MU),MU_upper=max(MME_prior$MME_MU),SIZE_lower=min(MME_prior$MME_SIZE),SIZE_upper=ceiling(max(MME_prior$MME_SIZE)),parallel=parallel,NCores=NCores,FIX_MU=FIX_MU,GR=GR)
    BB_prior<-BB_size
    rownames(BB_prior)<-rownames(Data)
    colnames(BB_prior)<-c("BB_SIZE","BB_MU")
    MME_SIZE_adjust<-AdjustSIZE_fun(BB_prior[,1],MME_prior$MME_MU,MME_prior$MME_SIZE)
  }

  if(verbose){
    message("Prior estimation has completed!")
  }

  return(list(MME_prior=MME_prior,BB_prior=BB_prior,MME_SIZE_adjust=MME_SIZE_adjust))
}

  return(list(MME_prior=MME_prior))
}


#' Estimating size for each gene by either 1D or 2D optimization of marginal distribution
#'
#' @param Data: A matrix of single-cell expression where rows are genes
#' and columns are samples (cells). This object should be of class
#' matrix rather than data.frame.
#' @param  BETA_vec: A vector of capture efficiencies of cells.
#' @param  INITIAL_MU_vec: Mean expression of genes, can come from EstPrior.
#' @param  INITIAL_SIZE_vec: size of genes (size is a parameter in NB distribution), can come from EstPrior.
#' @param  MU_lower: The lower bound for the mu.(Only need it when you want to do 2D optimization). Defaut is 0.01.
#' @param  MU_upper: The upper bound for the mu.(Only need it when you want to do 2D optimization). Defaut is 500.
#' @param  SIZE_lower: The lower bound for the size. Defaut is 0.01.
#' @param  SIZE_upper: The upper bound for the size. Defaut is 30.
#' @param  parallel: If TRUE, 5 cores will be used for parallelization. Default is TRUE.
#' @param  NCores: number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  FIX_MU: If TRUE, then 1D optimization, otherwise 2D optimization (slow).
#' @param  GR: If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @return  A vector of estimated size.
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @useDynLib bayNorm
#' @importFrom Rcpp sourceCpp
#'
#' @export
#'
#'
BB_Fun<-function(Data,BETA_vec,INITIAL_MU_vec,INITIAL_SIZE_vec,MU_lower=0.01,MU_upper=500,SIZE_lower=0.01,SIZE_upper=30,parallel=F,NCores=5,FIX_MU=T,GR=F)
{
  if(FIX_MU){#1D

    lower_input=SIZE_lower
    upper_input=SIZE_upper

    if(parallel){
      cluster = makeCluster(NCores, type = "SOCK")
      registerDoSNOW(cluster)
      getDoParWorkers()

      BB_parmat<-foreach(Geneind=1:dim(Data)[1],.combine=c )%dopar%{
        print(Geneind)
        mu<-INITIAL_MU_vec[Geneind]
        size<-INITIAL_SIZE_vec[Geneind]
        m_observed=Data[Geneind,]

        if(!GR){
          BB_opt<-BB::spg(par=size,fn=MarginalF_1D,MU=mu,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,maxfeval=500),lower=lower_input,upper=upper_input)
        }else{
          BB_opt<-BB::spg(par=size,fn=MarginalF_1D,gr=GradientFun_1D,MU=mu,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,maxfeval=500),lower=lower_input,upper=upper_input)
        }

        optimal_par<-BB_opt$par
        return(optimal_par)
      }

      stopCluster(cluster)


    } else{

      BB_parmat<-foreach(Geneind = 1:dim(Data)[1],.combine=c)%do%{

        print(Geneind)
        mu<-INITIAL_MU_vec[Geneind]
        size<-INITIAL_SIZE_vec[Geneind]
        m_observed=Data[Geneind,]

        if(!GR){
          BB_opt<-BB::spg(par=size,fn=MarginalF_1D,MU=mu,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,maxfeval=500),lower=SIZE_lower,upper=SIZE_upper)
          }else{
            BB_opt<-BB::spg(par=size,fn=MarginalF_1D,gr=GradientFun_1D,MU=mu,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,maxfeval=500),lower=SIZE_lower,upper=SIZE_upper)

          }
        #
        optimal_par<-BB_opt$par
        return(optimal_par)
      }
    }
    names(BB_parmat)<-rownames(Data)
  }else{# 2D

    lower_input=c(SIZE_lower,MU_lower)
    upper_input=c(SIZE_upper,MU_upper)


    if(parallel){
      cluster = makeCluster(NCores, type = "SOCK")
      registerDoSNOW(cluster)
      getDoParWorkers()

      BB_parmat<-foreach(Geneind=1:dim(Data)[1],.combine=rbind )%dopar%{
        print(Geneind)
        mu<-INITIAL_MU_vec[Geneind]
        size<-INITIAL_SIZE_vec[Geneind]
        m_observed=Data[Geneind,]
        if(!GR){
          BB_opt<-BB::spg(par=c(size,mu),fn=MarginalF_2D,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,ftol=0.001,maxfeval=500),lower=lower_input,upper=upper_input)

        }else{
          BB_opt<-BB::spg(par=c(size,mu),fn=MarginalF_2D,gr=GradientFun_2D,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,maxfeval=500),lower=lower_input,upper=upper_input)
        }

        optimal_par<-BB_opt$par
        return(optimal_par)
      }

      stopCluster(cluster)


    } else{


      BB_parmat<-foreach(Geneind = 1:dim(Data)[1],.combine=rbind)%do%{

        print(Geneind)
        mu<-INITIAL_MU_vec[Geneind]
        size<-INITIAL_SIZE_vec[Geneind]
        m_observed=Data[Geneind,]
        if(!GR){
          BB_opt<-BB::spg(par=c(size,mu),fn=MarginalF_2D,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,ftol=0.001,maxfeval=500),lower=c(SIZE_lower,MU_lower),upper=c(SIZE_upper,MU_upper))
        }else{
          BB_opt<-BB::spg(par=c(size,mu),fn=MarginalF_2D,gr=GradientFun_2D,m_observed=m_observed,BETA=BETA_vec,control=list(maximize=T,trace=FALSE,maxfeval=500),lower=c(SIZE_lower,MU_lower),upper=c(SIZE_upper,MU_upper))
        }

        optimal_par<-BB_opt$par
        return(optimal_par)
      }
    }


    rownames(BB_parmat)<-rownames(Data)

  }

  return(BB_parmat)
}

