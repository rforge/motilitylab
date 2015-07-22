\name{motilityLab-package}
\alias{motilityLab-package}
\alias{motilityLab}
\docType{package}
\title{
MotilityLab: Quantitative Motion Analysis
}
\description{
Quantitative analysis of tracks of moving
objects like, for instance, cells, bacteria or animals.
Compute various flavours of common measures like mean
square displacement, confinement ratio, autocorrelation,
asphericity and so on.
}
\details{
\tabular{ll}{
Package: \tab motilityLab\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-04\cr
License: \tab GPL-2\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Katharina Dannenberg, Johannes Textor

Maintainer: Johannes Textor <johannes.textor@gmx.de>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{cluster
spatial}
\examples{
## a simple hierarchical clustering based on 2D asphericity

## tag track IDs so we can identify them later
names(TCells) <- paste0("T",names(TCells))
names(BCells) <- paste0("B",names(BCells))
names(Neutrophils) <- paste0("N",names(Neutrophils))

## project all tracks down to 2D
cells <-  projectDimensions(
		c(TCells,BCells,Neutrophils), c("x","y") )

## compute asphericity
asph <- lapply( cells, asphericity )

## plot clustering
plot(hclust(dist(asph)))
}
