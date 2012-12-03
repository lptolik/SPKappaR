\name{showStruct}
\alias{showStruct}
\title{showStruct}
\description{Modifyes and shows 3D structure of domain defined by logical function}
\usage{showStruct(struct = spine, color = "black", alpha = 1, def, show = FALSE, 
    ...)}
\arguments{
  \item{struct}{basic structure of the domain that have to be modified.}
  \item{color}{color of spheres, matched by the function}
  \item{alpha}{transparency of spheres, matched by the function}
  \item{def}{logical function that defines coordinates of the voxels to be modified}
  \item{show}{logical flag, if show=TRUE new domain will be visualised by the 'rgl' library}
  \item{\dots}{additional parameters required by the logical function}
}

\value{data.frame with columns X,Y,Z,Col,Alph that contains coordinates (X,Y,Z) color (Col) and transparancy value (Alph) for voxels}

\author{Anatoly Sorokin, Oksana Sorokina}
\note{Logica def function could loop like def=function(.x,...){return(.x$X>20&.x$Col=='green')}}




