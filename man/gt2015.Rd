\name{gt2015}
\alias{gt2015}
\docType{data}
\title{
Gas turbine NOx emission
}
\description{
The gas turbine NOx emission data set.
}
\usage{data("gt2015")}
\format{
  A data frame with 7384 observations on the following 11 variables.
  \describe{
    \item{\code{AT}}{a numeric vector}
    \item{\code{AP}}{a numeric vector}
    \item{\code{AH}}{a numeric vector}
    \item{\code{AFDP}}{a numeric vector}
    \item{\code{GTEP}}{a numeric vector}
    \item{\code{TIT}}{a numeric vector}
    \item{\code{TAT}}{a numeric vector}
    \item{\code{TEY}}{a numeric vector}
    \item{\code{CDP}}{a numeric vector}
    \item{\code{CO}}{a numeric vector}
    \item{\code{NOX}}{a numeric vector}
  }
}
\details{
To predict nitrogen oxide emissions, we use the gas turbine NOx emission data set in UCI database, which contains 36,733 instances of 11,733 sensor measurements. The pollutant emission factors of gas turbines include 9 variables. We select 7,200 data points in 2015.
}
\source{
The gas turbine NOx emission data set is from UCI database.
}
\references{
NA
}
\examples{
data(gt2015)
## maybe str(gt2015) ; plot(gt2015) ...
}
\keyword{datasets}