\name{updateStruct}
\alias{updateStruct}
\title{updateStruct}
\description{Function incorporate data from SPkapp snapshot into the 3D model of the domain   }
\usage{updateStruct(struct = spine, amend, show = TRUE, alpha = 0.03)}
\arguments{
  \item{struct}{base model of the domain}
  \item{amend}{simulation results to be incorporated into the model}
  \item{show}{logical value defining should be the new structure sent to RGL package for visualilation}
  \item{alpha}{transparancy value for the replaced voxels}
}

\value{data.frame with columns X,Y,Z,Col,Alph that contains coordinates (X,Y,Z) color (Col) and transparancy value (Alph) for voxels}

\author{Anatoly Sorokin, Oksana Sorokina}





