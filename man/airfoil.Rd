\name{airfoil}
\alias{airfoil}
\docType{data}
\title{
Airfoil self-noise
}
\description{
The Airfoil self-noise data set
}
\usage{data("airfoil")}
\format{
  A data frame with 1503 observations on the following 6 variables.
  \describe{
    \item{\code{V1}}{a numeric vector}
    \item{\code{V2}}{a numeric vector}
    \item{\code{V3}}{a numeric vector}
    \item{\code{V4}}{a numeric vector}
    \item{\code{V5}}{a numeric vector}
    \item{\code{V6}}{a numeric vector}
  }
}
\details{
The data set contains 1503 data points, including the 6 variables. Among them, the scaled sound pressure level is the dependent variable and the other five are independent variables.
}
\source{
The Airfoil Self-Noise data set is from the NASA data set in UCI database. 
}
\references{
T.F. Brooks, D.S. Pope, and A.M. Marcolini. Airfoil self-noise and prediction. Technical report, NASA RP-1218, July 1989.
}
\examples{
data(airfoil)
## maybe str(airfoil) ; plot(airfoil) ...
}
\keyword{datasets}