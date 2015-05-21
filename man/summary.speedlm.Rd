\name{summary.speedlm}
\alias{summary.speedlm}
\alias{coef.speedlm}
\alias{vcov.speedlm}
\alias{logLik.speedlm}
\alias{AIC.speedlm}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods to summarize Linear Models fits}
\description{
\code{summary} method for class 'speedlm'.   
}
\usage{
\method{summary}{speedlm}(object, correlation = FALSE,...)
\method{coef}{speedlm}(object,...)
\method{vcov}{speedlm}(object,...)
\method{logLik}{speedlm}(object,...)
\method{AIC}{speedlm}(object,...,k = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class 'speedlm'.}
  \item{correlation}{logical. Do you want to print the correlation matrix? By default it is false.}
  \item{k}{numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
  \item{\dots}{further optional arguments}
}


\value{
  \item{coefficients}{the matrix of coefficients, standard errors, t-statistics and two-side p-values.}
  \item{rdf}{degrees of freedom of the fitted model. It is a component from \code{object}.}
  \item{call}{the component from \code{object}.}
  \item{r.squared}{R^2, the fraction of variance explained by the model.}
  \item{adj.r.squared}{the "adjusted" R^2 statistic, penalizing for higher p.}
  \item{fstatistic}{(for models including non-intercept terms) a 3-vector with the value of 
  the F-statistic with its numerator and denominator degrees of freedom.}
  \item{f.pvalue}{p-value of the F-statistic.}
  \item{RSS}{Residual sum of squares.}  
  \item{var.res}{estimated variance of residuals.}  
  \item{rank}{the component from \code{object}.}
  \item{correlation}{(only if \code{correlation} is true) the correlations of the estimated parameters.}
  \item{\dots}{the results from the functions \code{logLik}, \code{AIC} and \code{vcov}.}
}
%\references{ put references to the literature/web site here  }
\author{ Marco ENEA }
%\note{ further notes 

 %Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{ \link{speedlm} }
\examples{

y <- rnorm(100,1.5,1)
x <- round(matrix(rnorm(200), 100, 2), digits = 3)
colnames(x) <- c("s1","s2")
da <- data.frame(y, x)

m <- speedlm(y ~ s1 + s2,da)
summary(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

