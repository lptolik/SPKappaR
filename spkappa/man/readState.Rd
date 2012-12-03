\name{readState}
\alias{readState}
\title{readState}
\description{main function of the package. Reads content of SPKappa snapshot file and defines color of identified voxels.}
\usage{readState(file, anchor = c("GluR1_m:membrane", "PSD95:cytosol"), 
    colorDef = NA)}
\arguments{
  \item{file}{snapshot file to read}
  \item{anchor}{within content of the voxel. 
list of regular expressions to be identified in}
  \item{colorDef}{snapshot just read.
predefined color set. If NA, new color definition will be generated from content of}
}

\value{Returns data.frame sutable for updateStruct amend argument. Attribute 'palette' of the
result contains colorDef structure.}

\author{Anatoly Sorokin, Oksana Sorokina}
\note{Complex definitions in anchor list processed sequentially in exclusive manner, so second complex define color}




