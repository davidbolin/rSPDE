\name{matern.covariance}
\alias{matern.covariance}
\title{
  The Matern covariance function
}
\description{
  \code{matern.covariance} evaluates the Matern covariance function
  \deqn{C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h)}{C(h) = (\sigma^2/(2^(\nu-1)\Gamma(\nu))(\kappa h)^\nu K_\nu(\kappa h)}
}
\usage{
  matern.covariance(h,
                    kappa,
                    sigma,
                    nu)
  }
\arguments{
  \item{h}{Distances to evaluate the covariance function at.}
  \item{kappa}{Range parameter.}
  \item{sigma}{Standard deviation.}
  \item{nu}{Shape parameter.}
}
\value{
\code{matern.covariance} return a vector with the values C(h).
}
\author{
  David Bolin \email{davidbolin@gmail.com}
}

\examples{

x = seq(from = 0, to = 1, length.out = 101)
plot(h,matern.covariance(abs(x-0.5),kappa=10,nu=1/5,sigma=1),type="l",ylab = "C(h)",xlab="h")

}
