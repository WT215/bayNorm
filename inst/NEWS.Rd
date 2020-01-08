\name{bayNorm: updating news}
\title{bayNorm News}
\encoding{UTF-8}

\section{Version 0.99.0}{
\itemize{
\item Pass all checks of both R CMD build and R CMD BiocCheck.
\item Future work:
\itemize{
   \item Vignette needs to be further improved.
   \item Need to Improve man pages.
   \item Submit to Bioconductor.}}
}

\section{Version 0.99.1}{
\itemize{
\item Fix two warnings reported by Bioconductor
\itemize{
   \item WARNING: Use TRUE/FALSE instead of T/F, Found in files: tests/testthat/test-bayNorm.r
   \item WARNING: Use is() or !is()
        }
    }
}



\section{Version 0.99.4}{
\itemize{
\item Apply corrections proposed by the reviewer, for more details please see \url{https://github.com/Bioconductor/Contributions/issues/878}
    }
}

\section{Version 0.99.8}{
\itemize{
\item Fix some errors.
    }
}

\section{Version 0.99.9}{
\itemize{
\item Let more examples to be runnable; Updated vignette; Including SingleCellExperiment Class.
    }
}

\section{Version 0.99.19}{
\itemize{
\item Replacing foreach with BiocParallel
\item Fix some errors in BB_fun
    }
}

\section{Version 1.1.4}{
\itemize{
\item For a smoother progress bar, use foreach instead of BiocParallel.
    }
}


\section{Version 1.1.6}{
\itemize{
\item Suppress warning message came from BiocParallel.
    }
}

\section{Version 1.5.0}{
\itemize{
\item Set default setting of BETA_vec=NULL.
\item Update citation.
    }
}

\section{Version 1.5.2}{
\itemize{
\item Vectorize computation in EstPrior function.
    }
}

\section{Version 1.5.3}{
\itemize{
\item Add Rcpp EstPrior function.
    }
}

\section{Version 1.5.4}{
\itemize{
\item Support sparse matrix: dgTMatrix and dgCMatrix.
    }
}

\section{Version 1.5.5}{
\itemize{
\item Make BB_fun faster: speed up the access of row of dgCMatrix object.
\item Modify EstPrior, using rowSums from Matrix package to calculated variance. Hence no need for Rcpp EstPrior function.
    }
}

\section{Version 1.5.6}{
\itemize{
\item as_matrix for transforming sparse matrix to matrix in R efficiently than as.matrix.
    }
}

\section{Version 1.5.7}{
\itemize{
\item Add useful function file.
\item Use Rcpp EstPrior function.
    }
}

\section{Version 1.5.8}{
\itemize{
\item Change back to rowSums instead of using Rcpp EstPrior function.
    }
}

\section{Version 1.5.9}{
\itemize{
\item Change back to Rcpp EstPrior function (sp_mat version).
\item Use dgCMatrix instead of matrix in R.
\item #define ARMA_64BIT_WORD 1 and so on allows for big sparse matrix (learn from the liger R package).
    }
}

\section{Version 1.5.10}{
\itemize{
\item Update RELEASE_3_10 1.4.2 to 1.4.3.
    }
}

\section{Version 1.5.11}{
\itemize{
\item Allow sparse matrix as output (valid only for mean version procedure).
\item Let version of RELEASE_3_10 to be 1.4.11.
    }
}

\section{Version 1.5.12}{
\itemize{
\item Deprecated: Allow sparse matrix as output (valid only for mean version procedure). Since sparse matrix is useful only when the matrix it is sparse.
\item Calibrate documentations.
    }
}

\section{Version 1.5.13}{
\itemize{
\item Activate "Allow sparse matrix as output (valid only for mean version procedure)". Hence you can work on dgCMatrix throughout the process. Input which is not of type dgCMatrix will be converted to dgCMatrix at first.
    }
}