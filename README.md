# bayNorm

bayNorm is an R package which is used to normalize single-cell RNA-seq data. 

The bayNorm in Julia language is now available: [https://github.com/WT215/bayNormJL.jl](https://github.com/WT215/bayNormJL.jl).

## Code for bayNorm paper
The code for producing figures in bayNorm paper [1] can be found [here](https://github.com/WT215/bayNorm_papercode)

## Installation

Make sure to use the latest version of bayNorm by installing it from GitHub. 


```R
library(devtools)
devtools::install_github("WT215/bayNorm")
```

bayNorm has been submitted to Bioconductor, once it is accepted, it can be installed via:
```R
library(BiocManager)
BiocManager::install("bayNorm")
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
rse <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=EXAMPLE_DATA_list$inputdata[,seq(1,30)]))
#SingleCellExperiment object can also be input in bayNorm:
#rse <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=EXAMPLE_DATA_list$inputdata))

#Return 3D array normalzied data, draw 20 samples from posterior distribution:
bayNorm_3D<-bayNorm(
    Data=rse,
    BETA_vec = NULL,
    mode_version=FALSE,
    mean_version = FALSE,S=20
    ,verbose =FALSE,
    parallel = TRUE)

#Return 2D matrix normalized data (MAP of posterior):
#Simply set mode_version=TRUE, but keep mean_version=FALSE


#Return 2D matrix normalized data (mean of posterior):
#Simply set mean_version=TRUE, but keep mode_version=FALSE
```

## Non-UMI scRNAseq dataset
bayNorm's mathematical model is suitable for UMI dataset. However it can be also applied on non-UMI dataset. In `bayNorm`, you need to specify the following parameter:
* `UMI_sffl`: bayNorm can also be applied on the non-UMI dataset. However, user need to provide a scaling factor. Raw data will be divided by the scaled number and bayNorm will be applied on the rounded scaled data. This scaling factor can be interpreted as the average number of times original mRNA molecules were sequenced after PCR amplification. It is chosen so that the Dropout vs Mean expression plots to be close to assymptotic expression ($e^{-mean}$).


## Output 3D array or 2D array with existing estimated prior parameters.
If you have run bayNorm on a dataset to estimate prior parameters, but want to output new posterior estimates (3D or 2D array), you can use the function `bayNorm_sup`. It is important to input the existing estimated parameters by specifying the following parameter in `bayNorm_sup`:
* `BETA_vec`: If `Conditions` has been specified previously, then input `unlist(bayNorm_output$BETA)`
* `PRIORS`: input `bayNorm_output$PRIORS`
* `input_params`: input `bayNorm_output$input_params`

```R
data('EXAMPLE_DATA_list')
#Return 3D array normalzied data:
bayNorm_3D<-bayNorm(
    Data=EXAMPLE_DATA_list$inputdata,
    BETA_vec = EXAMPLE_DATA_list$inputbeta,
    mode_version=FALSE,
    mean_version = FALSE)

#Now if you want to generate 2D matrix (MAP) using the same prior
#estimates as generated before:
bayNorm_2D<-bayNorm_sup(
    Data=EXAMPLE_DATA_list$inputdata,
    PRIORS=bayNorm_3D$PRIORS,
    input_params=bayNorm_3D$input_params,
    mode_version=TRUE,
    mean_version = FALSE)

#Or you may want to generate 2D matrix 
#(mean of posterior) using the same prior
#estimates as generated before:
bayNorm_2D<-bayNorm_sup(
    Data=EXAMPLE_DATA_list$inputdata,
    PRIORS=bayNorm_3D$PRIORS,
    input_params=bayNorm_3D$input_params,
    mode_version=FALSE,
    mean_version = TRUE)
```

## Work around Seurat object: clusters detection
```R
library(bayNorm)
library(Seurat)

data('EXAMPLE_DATA_list')

library(Seurat)
bay_out<-bayNorm(EXAMPLE_DATA_list$inputdata,mean_version = TRUE)
x.seurat <- CreateSeuratObject(counts =bay_out$Bay_out,assay = 'bayNorm')
# x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat)
#Specifying: assay='bayNorm'
x.seurat <- RunPCA(x.seurat, features = x.seurat@assays$bayNorm@var.features,
                   pcs.compute = 20,assay='bayNorm')

x.seurat <- RunUMAP(x.seurat, dims = 1:20,assay='bayNorm')
#x.seurat <- JackStraw(x.seurat, prop.freq = 0.06)
x.seurat <- FindNeighbors(x.seurat, dims = 1:20)
x.seurat <- FindClusters(x.seurat, resolution = 0.5)
head(Idents(x.seurat), 5)
#It is a toy example. Here only one cluster was found
plot(x.seurat@reductions$umap@cell.embeddings,pch=16,col=as.factor(Idents(x.seurat)))
#Double check that the assay we used comes from bayNorm
x.seurat@reductions$pca@assay.used
```

## References

- [1] <a href="https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz726/5581401">Tang <em>et al.</em> (2019). Bioinformatics. </a>
