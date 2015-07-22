\name{Neutrophils}
\alias{Neutrophils}
\docType{data}
\title{
Two-Photon Data: Neutrophils homing towards an infection in the lung
}
\description{
Labelled neutrophils  were adoptively transfered and intravitally
imaged (using two-photon microscopy) inside the lung 
of the recipient mouse. These cells show a fairly directed kind of
motion, as they move towards infection foci.
}
\usage{data("Neutrophils")}
\format{
  An S3 object of class "tracks"; a list with 10 elements. Each 
  element name identifies a cell track. Each element is a matrix
  containing the following four columns.
  \describe{
    \item{\code{t}}{the time (in seconds)}
    \item{\code{x}}{The X coordinate (in micrometers)}
    \item{\code{y}}{The Y coordinate (in micrometers)}
    \item{\code{z}}{The Z coordinate (in micrometers)}
  }
}
\source{
Data were generated in 2012 in the Mark J. Miller Lab, 
Department of Medicine,
Washington University in St Louis, USA.
}
\examples{
## load the tracks
data(Neutrophils)
## visualize the tracks (calls function plot.tracks)
plot(Neutrophils)
}
\keyword{datasets}
