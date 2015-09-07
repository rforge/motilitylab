\name{motilityLab-package}
\alias{motilityLab-package}
\alias{motilityLab}
\docType{package}
\title{
MotilityLab: Quantitative Motion Analysis
}
\description{
Quantitative analysis of tracks of moving objects like
cells, bacteria or animals. Compute various 
common measures like mean square displacement, confinement ratio, 
autocorrelation, and asphericity. 
}
\details{
\tabular{ll}{
Package: \tab motilityLab\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-04\cr
License: \tab GPL-2\cr
}
The MotilityLab package is designed for analyzing cell tracks acquired by
time-lapse microscopy (like those provided in the included datasets 
\code{\link{TCells}}, \code{\link{BCells}} and \code{\link{Neutrophils}}). 
However, it can process any kind of x-y-(z)-t data, and
we hope that it may be useful for other purposes as well.

The basic data structure that most functions in this package operate on is a set of 
\emph{tracks}. A track is a list of spatial coordinates that are recorded at 
\emph{fixed} time intervals. Thus, we assume that each track is stored in a matrix whose
first column denotes a time interval (for example, seconds elapsed since the beginning
of the experiment), and whose remaining columns denote a spatial coordinate. A set of
tracks is stored as a \code{\link{list}} with S3 class \code{tracks}. MotilityLab provides
some S3 functions for this class, e.g. \code{\link{c.tracks}} to merge track sets and
\code{\link{plot.tracks}} to plot them.

A wide range of common track measures are included in the package. These are all functions
that take a single track as an input, and output one or several numbers that describe
the input track quantitatively. For instance, the function \code{\link{speed}} estimates
the average instantaneous speed of the track by linear interpolation, and 
\code{\link{straightness}} computes the start-to-end distance divided by the trajectory
length (this is a number between 0 and 1, where 1 indicates a perfectly straight track).

MotilityLab is designed to support various flavors of track analysis that have been
suggested in the literature. The simplest kind is a \emph{track-based} analysis, where
we compute a single statistic for each track in a dataset (Beltman et al, 2009). Because
track sets are lists, this is achieved simply by using \code{\link{lapply}} or
\code{\link{sapply}} together with the track measure (see Examples).
 
A second kind of track analysis is called \emph{step-based} (Beltman et al, 2009). Here, we
chop each track up into segments of the same length and then apply our measures to those
segments instead of the whole track. This can help to avoid biases that arise when 
the tracks have variable length (as is nearly always the case in cell tracking experiments).
In MotilityLab, step-based analyses are performed by using the \code{\link{subtracks}}
function. Often we want to perform such step-based analyses for all possible subtrack
lengths simultaneously, and plot the result as a function of the subtrack length. A 
famous example for such an analysis is the mean square displacement plot. This can be 
achieved by using the \code{\link{aggregate.tracks}} function, which has several options
to control which subtrack lengths are considered, whether overlapping subtracks are 
considered etc.

A third kind of analysis is where we analyze all subtracks of a single track. The result
is often plotted as a matrix. This can reveal dynamic patterns along a single track,
e.g. turning behaviour or local slowdowns. This so-called \emph{staggered} analysis 
(Mokhtari et al, 2013) is performed by the \code{\link{computeStaggered}} function.

Lastly, in addition to data analysis, the package contains some function to generate 
cell tracks by simulation. This is useful to develop and benchmark track analysis 
methodology (Textor et al, 2011), and for computational biology studies that try to 
extrapolate the long-term consequences of observed cell migration behaviour. Alongside
a simple uncorrelated random walk (\code{\link{brownianTrack}}), this package implements
a simulation model proposed by Beauchemin et al (2007) 
in the function \code{\link{beaucheminTrack}}. That model can also simulate 
directionally biased motion.
}
\author{
Katharina Dannenberg, Jeffrey Berry, Johannes Textor

Maintainer: Johannes Textor <johannes.textor@gmx.de>
}
\references{

Joost B. Beltman, Athanasius F.M. Maree and Rob. J. de Boer (2009).
Analysing immune cell migration. \emph{Nature Reviews Immunology} \bold{9},
789--798. doi:10.1038/nri2638

Zeinab Mokhtari, Franziska Mech, Carolin Zitzmann, Mike Hasenberg, Matthias Gunzer
and Marc Thilo Figge (2013), Automated Characterization and 
Parameter--Free Classification of Cell Tracks Based on Local Migration 
Behavior. \emph{PLoS ONE} \bold{8}(12), e80808. doi:10.1371/journal.pone.0080808

Johannes Textor, Antonio Peixoto, Sarah E. Henrickson, Mathieu
Sinn, Ulrich H. von Andrian and Juergen Westermann (2011),
Defining the Quantitative Limits of Intravital Two-Photon Lymphocyte Tracking.
\emph{PNAS} \bold{108}(30):12401--12406. doi:10.1073/pnas.1102288108

Catherine Beauchemin, Narendra M. Dixit and Alan S. Perelson (2007), Characterizing 
T cell movement within lymph nodes in the absence of antigen. \emph{Journal of Immunology}
\bold{178}(9), 5505-5512. doi:10.4049/jimmunol.178.9.5505
}
\keyword{cluster
spatial}
\examples{
## track-based speed comparison
boxplot( sapply( Neutrophils, straightness ), sapply( BCells, straightness ) )

## step-based turning angle comparison
boxplot( sapply( subtracks( Neutrophils, 2 ), overallAngle  ), 
	sapply( subtracks( BCells, 2 ), overallAngle ) )

## mean square displacement plot; a step-based displacement analysis for all step lengths
plot( aggregate( TCells, squareDisplacement )[,"value"] ) 

## 'staggered' analysis of displacement over whole track. Reveals that this track
## slows down near its beginning and near its end. 
filled.contour( computeStaggered(TCells[[4]], displacement, matrix=TRUE) )

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
