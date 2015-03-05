\name{TCells}
\alias{TCells}
\docType{data}
\title{
Two-Photon Data: T Cells in a Lymph Node
}
\description{
Labelled T cells were adoptively transfered and intravitally
imaged (using two-photon microscopy) inside a peripheral lymph
node of the recipient mouse. These data represent the characteristic
"random-walk-like" motion pattern of T cells in lymph nodes.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("TCells")}
\format{
  An S3 object of class "tracks"; a list with 22 elements. Each 
  element name identifies a cell track. Each element is a data 
  frame on the following four variables.
  \describe{
    \item{\code{t}}{the time (in seconds)}
    \item{\code{x}}{The X coordinate (in micrometers)}
    \item{\code{y}}{The Y coordinate (in micrometers)}
    \item{\code{z}}{The Z coordinate (in micrometers)}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
## load the tracks
data(TCells)
## visualize the tracks (calls function plot.tracks)
plot(TCells)
}
\keyword{datasets}
