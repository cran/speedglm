\name{predict.speedglm}
\alias{predict.speedglm}


\title{Predict method for a speedglm object}
\description{
\code{summary} The method is currently under construction but some functionalities are available.
}
\usage{

\method{predict}{speedglm}(object, newdata, type = c("link", "response"), 
                           na.action = na.pass, ...)
}

\arguments{
  \item{object}{an object of class 'speedlgm'.}
  \item{newdata}{An optional data frame with new data or the original data.}
  \item{type}{Type of prediction.}
  \item{na.action}{function determining what should be done with missing values in \code{newdata}.}
  \item{\dots}{further optional arguments}
}
\details{If \code{newdata} is omitted prediction are based on the data used for the fit only if argument \code{fitted} was previously set to \code{TRUE} in the speedglm object.
Currently the method does not work for function \code{shglm}.}
\value{
  \item{pred}{a vector of predictions.}
}

\author{ Tomer Kalimi and Marco Enea }

\seealso{ \link{speedglm} }
\examples{

data(data1)
mod <- speedglm(y~x1+x2+factor(fat1), data=data1, family=Gamma(log), fitted=TRUE)
predict(mod)

}

\keyword{models}

