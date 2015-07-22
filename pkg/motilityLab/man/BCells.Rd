\name{BCells}
\alias{BCells}
\docType{data}
\title{
Two-Photon Data: B Cells in a Lymph Node
}
\description{
Labelled B cells were adoptively transfered and intravitally
imaged (using two-photon microscopy) inside a peripheral lymph
node of the recipient mouse. These data illustrate the characteristic
"random-walk-like" motion pattern of B cells in lymph nodes.
}
\usage{data("BCells")}
\format{
  An S3 object of class "tracks"; a list with 24 elements. Each 
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
data(BCells)
## visualize the tracks (calls function plot.tracks)
plot(BCells)
}
\keyword{datasets}
