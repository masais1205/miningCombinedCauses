\name{mhpc-package}
\alias{mhpc-package}
\alias{mhpc}
\docType{package}
\title{
Mining combined causes
}
\description{
Mining combined causes (based on bnlearn - hiton-pc)
}
\details{
\tabular{ll}{
Package: \tab mhpc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-10-02\cr
License: \tab GPLv3 \cr
}
This package implements to find single and combined causes
}
\author{
Saisai Ma \cr
University of South Australia

Maintainer: Saisai Ma \email{maysy020@mymail.unisa.edu.au}

}
\references{
Ma, Saisai, Jiuyong Li, Lin Liu, and Thuc Duy Le. "Mining combined causes in large data sets." Knowledge-Based Systems 92 (2016): 104-111.
}
\examples{
library(mhpc)
mhpc(data, variant = "B", alpha=0.05, supp=0.05, cor.pvalue=0.01) 
# binary data; variant: B for MH-PC-B or F for MH-PC-F; 
# alpha: threshold of the conditional independece test; 
# supp & cor.pvalue: threshold of the pruning scheme.
}
