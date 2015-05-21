\name{predict.speedlm}
\alias{predict.speedlm}


\title{Predict method for a speedlm object}
\description{
\code{summary} The method is currently under construction but some functionalities are available.
}
\usage{

\method{predict}{speedlm}(object, newdata, na.action = na.pass, ...)
}

\arguments{
  \item{object}{an object of class 'speedlm'.}
  \item{newdata}{An optional data frame with new data or the original data.}
  \item{na.action}{function determining what should be done with missing values in               
        \code{newdata}.}
  \item{\dots}{further optional arguments}
}
\details{If \code{newdata} is omitted prediction are based on the data used for the fit only if argument \code{fitted} was previously set to \code{TRUE} in the speedlm object.}
\value{
  \item{predictor}{a vector of predictions.}
}


\author{ Tomer Kalimi and Marco Enea }

\seealso{ \link{speedlm} }
\examples{
set.seed(10)
x <- round( matrix(rnorm(20*3),20,3),digits=3)
colnames(x) <-c("y","s1","s2")
da <- as.data.frame(x)
mod <- speedlm(y~s1+s2, data=da, fitted=TRUE)
predict(mod)

}

\keyword{models}

