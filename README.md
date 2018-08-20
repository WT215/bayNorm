# bayNorm

bayNorm is an R package which is used to normalize single-cell RNA-seq data. 

## Code for bayNorm paper
The code for producing figures in bayNorm paper [1] can be found [here](https://github.com/WT215/bayNorm_papercode)

## Installation

Make sure to use the latest version of bayNorm by installing it from GitHub. 


```R
library(devtools)
devtools::install_github("WT215/bayNorm")
```


# Quick start: for either single or groups of cells
The main function is `bayNorm`, which is a wrapper function for  gene specific prior parameter estimation and normalization. The input is a matrix of scRNA-seq data with rows different genes and coloums different cells. The output is either point estimates from posterior (2D array) or samples from posterior (3D array). 

Essential input and parameters for running `bayNorm` are: 

* `Data`: a `SummarizedExperiment` object or matrix (rows: genes, columns: cells). 
* `BETA_vec`: a vector of probabilities which is of length equal to the number of cells. 
* `Conditions`: If `Conditions` is provided, prior parameters will be estimated within each group of cells (we name this kind of procedure as "LL" procedure where "LL" stands for estimating both $\mu$ and $\phi$ locally). Otherwise, bayNorm applied "GG" procedure for estimating prior parameters (estimating both $\mu$ and $\phi$ globally).
* `Prior_type`: Even if you have specified the `Conditions`, you can still choose to estimate prior parameters across all the cells by setting `Prior_type="GG"`.

***

   


```R
data('EXAMPLE_DATA_list')
#Return 3D array normalzied data:
bayNorm_3D<-bayNorm(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec = EXAMPLE_DATA_list$inputbeta,
    mode_version=F,
    mean_version = F)

#Return 2D matrix normalized data (MAP of posterior):
bayNorm_2D<-bayNorm(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec = EXAMPLE_DATA_list$inputbeta,
    mode_version=T,
    mean_version = F)

#Return 2D matrix normalized data (mean of posterior):
bayNorm_2D<-bayNorm(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec = EXAMPLE_DATA_list$inputbeta,
    mode_version=F,
    mean_version = T)
```

## Non-UMI scRNAseq dataset
bayNorm's mathematical model is suitable for UMI dataset. However it can be also applied on non-UMI dataset. In `bayNorm`, you need to specify the following parameter:
* `UMI_sffl`: bayNorm can also be applied on the non-UMI dataset. However, user need to provide a scaling factor. Raw data will be divided by the scaled number and bayNorm will be applied on the rounded scaled data. This scaling factor can be interpreted as the average number of times original mRNA molecules were sequenced after PCR amplification. It is chosen so that the Dropout vs Mean expression plots to be close to assymptotic expression ($e^{-mean}$).


## Output 3D array or 2D array with existing estimated prior parameters.
If you have run bayNorm on a dataset to estimate prior parameters, but want to output new posterior estimates (3D or 2D array), you can use the function `bayNorm_sup`. It is important to input the existing estimated parameters by specifying the following parameter in `bayNorm_sup`:
* `BETA_vec`: If `Conditions` has been specified previously, then input `unlist(bayNorm_output$BETA)`
* `PRIORS`: `input bayNorm_output$PRIORS_LIST`
* `Conditions`: make sure to specify the same `Conditions` as before.
You can find these two objects from the previous output of bayNorm function, which is a list.

```R
#Return 3D array normalzied data:
bayNorm_3D<-bayNorm(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec = EXAMPLE_DATA_list$inputbeta,
    mode_version=F,
    mean_version = F)

#Now if you want to generate 2D matrix (MAP) using the same prior
#estimates as generated before:
bayNorm_2D<-bayNorm_p(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec= bayNorm_3D$BETA,
    PRIORS=bayNorm_3D$PRIORS_LIST,
    mode_version=T,
    mean_version = F)

#Or you may want to generate 2D matrix 
#(mean of posterior) using the same prior
#estimates as generated before:
bayNorm_2D<-bayNorm_p(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec= bayNorm_3D$BETA,
    PRIORS=bayNorm_3D$PRIORS_LIST,
    mode_version=F,
    mean_version = T)
```

## References

- [1] <a href="https://www.biorxiv.org/content/early/2018/08/03/384586">Tang <em>et al.</em> (2018). bioRxiv. </a>
