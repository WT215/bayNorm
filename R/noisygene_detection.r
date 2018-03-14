

#' A wrapper function of synthetic control generation, bayNorm on both real cell data and synthetic controls and noisy gene detection.
#'
#' @param Data A matrix of single-cell expression where rows are genes and columns are samples (cells). This object should be of class matrix rather than data.frame.
#' @param  BETA_vec A vector of capture efficiencies of cells.
#' @param  Conditions No need to specify.
#' @param UMI_sffl No need to specify. Currently this function is for UMI based data.
#' @param  Prior_type No need to specify. Default is NULL. If \code{Conditions} is NULL, priors are estimated based on all cells.
#' @param  mode_version If TRUE, bayNorm return mode version normalized data which is of 2D matrix instead of 3D array. Default is FALSE.
#' @param S The number of samples you would like to generate from estimated posterior distribution (The third dimension of 3D array). Default is 20. S needs to be specified if \code{mode_version}=FALSE.
#' @param  parallel If TRUE, 5 cores will be used for parallelization.
#' @param  NCores number of cores to use, default is 5. This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) or SnowParam (Windows) with NCores using the package BiocParallel.
#' @param  FIX_MU Whether fix mu when estimating parameters by maximizing marginal distribution. If TRUE, then 1D optimization, otherwise 2D optimization (slow).
#' @param  GR If TRUE, the gradient function will be used in optimization. However since the gradient function itself is very complicated, it does not help too much in speeding up. Default is FALSE.
#' @param  BB_SIZE If TRUE, estimate BB size, and then use it for adjusting MME SIZE. Use the adjusted MME size for bayNorm. Default is TRUE.
#' @param verbose print out status messages. Default is TRUE.
#' @return  List of objects. The first element in the list is the adjusted P-values (noisy genes detection).
#'
#' @details A wrapper function of synthetic control generation, bayNorm on both real cell data and synthetic controls and noisy gene detection.
#'
#' @examples
#' \dontrun{
#' data("noisy_gene_check")
#' noisy_out<-noisy_gene_detection(Data=inputdata,BETA_vec
#' =inputbeta, mode_version = F, S = 20,parallel = T, NCores = 5,
#' FIX_MU = T, GR = F, BB_SIZE = T, verbose = T, plot.out = T)
#'}
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
#'

noisy_gene_detection<-function(Data,BETA_vec,mode_version=F,S=20,parallel=T,NCores=5,FIX_MU=T,GR=F,BB_SIZE=T,verbose=T,plot.out=F){



message("Apply bayNorm on the real cell datasets.")
bayNorm_N_out<-bayNorm(Data=Data,BETA_vec=BETA_vec,Conditions=NULL,UMI_sffl=NULL,Prior_type = NULL,S=S,parallel = parallel,NCores=NCores,FIX_MU = FIX_MU,GR=GR,BB_SIZE=BB_SIZE,verbose=verbose)

message("Generate synthetic control for the input data.")
synthetic_out<-SyntheticControl(Data=Data,BETA_vec=BETA_vec)

message("Apply bayNorm on the synthetic control.")

mu_c<-bayNorm_N_out$PRIORS$MME_prior$MME_MU
size_c<-bayNorm_N_out$PRIORS$MME_SIZE_adjust

if(!mode_version){
  bayNorm_C_array<-Main_Bay(Data=synthetic_out$N_c,BETA_vec=synthetic_out$beta_c,size=size_c,mu=mu_c,S=S,thres=max(synthetic_out$N_c)*2)

}else{
  bayNorm_C_array<-Main_mode_Bay(Data=synthetic_out$N_c,BETA_vec=synthetic_out$beta_c,size=size_c,mu=mu_c,S=S,thres=max(synthetic_out$N_c)*2)
}


NOISE_out<-NOISY_FUN(bayNorm_N_out[[1]],bayNorm_C_array,plot.out=plot.out)

return(list(adjusted_Pvals=NOISE_out,synthetic_output=synthetic_out,bayNorm_N_out=bayNorm_N_out,bayNorm_C_array=bayNorm_C_array))

}





#' Noisy gene detection
#'
#' Input raw data and a vector of capture efficiencies of cells. You can
#' also need to specify the condition of cells.
#' @param bay_array_N A 2D matrix or 3D array of normalized data(real cells).
#' @param  bay_array_C A 2D matrix or 3D array of normalized data(synthetic control).
#' @details Noisy gene detection
#'
#' @import foreach
#' @import MASS
#' @import locfit
#'
#' @export
#'



NOISY_FUN<-function(bay_array_N,bay_array_C,plot.out=F){

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
  h<-predict(loessfit, newdata=x, se.fit=T)
  zval<-(y-h$fit)/h$residual.scale

  #####script###
  # use kernel density estimate to find the peak
  dor <- density(zval, kernel = "gaussian")
  distMid <-dor$x[which.max(dor$y)]
  dist2 <- zval - distMid
  tmpDist <- c(dist2[dist2 <= 0], abs(dist2[dist2 < 0])) + distMid
  distFit <- fitdistr(tmpDist, "normal")
  pRaw <- pnorm(zval, mean = distFit$estimate[1], sd = distFit$estimate[2], lower.tail = FALSE)
  pAdj <- p.adjust(pRaw, 'BH')
  noisy<-which(pAdj<0.05)
  noisy_ind<-match(names(noisy),Gene_names[Select])
  noisy_name<-names(noisy)

  col_N<-1
  col_C<-rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
  col_noisy<-rgb(red = 0, green = 1, blue = 0, alpha = 0.5)

  if(plot.out){
    plot(x, y,pch=16,main='Noisy genes detection',col=col_N,xlab='log: Mean expression',ylab='log: CV^2')
    lines(loessfit,col=6,lwd=3)
    points(xc, yc, col = col_C,pch=16)
    points(x[noisy_ind], y[noisy_ind], col = col_noisy,pch=16)
    legend('topright',legend=c('Genes: real cell data','Genes: sysnthetic control','Noisy genes'),col=c(col_N,col_C,col_noisy),pch=16,cex=1)
  }


  pAdj2<-rep(1,length(Gene_names))
  pAdj2[Select]<-pAdj
  names(pAdj2)<-Gene_names
  return(pAdj2)
}




