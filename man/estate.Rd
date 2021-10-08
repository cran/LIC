\name{estate}
\alias{estate}
\docType{data}
\title{
Real estate valuation
}
\description{
The real estate valuation data set.
}
\usage{data("estate")}
\format{
  A data frame with 414 observations on the following 8 variables.
  \describe{
    \item{\code{No}}{a numeric vector}
    \item{\code{X1.transaction.date}}{a numeric vector}
    \item{\code{X2.house.age}}{a numeric vector}
    \item{\code{X3.distance.to.the.nearest.MRT.station}}{a numeric vector}
    \item{\code{X4.number.of.convenience.stores}}{a numeric vector}
    \item{\code{X5.latitude}}{a numeric vector}
    \item{\code{X6.longitude}}{a numeric vector}
    \item{\code{Y.house.price.of.unit.area}}{a numeric vector}
  }
}
\details{
Real estate valuation data set contains information about 414 real estate prices of 5 independent variables. The dependent variable is the price per unit area.
}
\source{
The data set is from Xindian District, New Taipei City, Taiwan.
}
\references{
Yeh, I. C., & Hsu, T. K. (2018). Building real estate valuation models with comparative approach through case-based reasoning. Applied Soft Computing, 65, 260-271.
}
\examples{
data(estate)
## maybe str(estate) ; plot(estate) ...
}
\keyword{datasets}