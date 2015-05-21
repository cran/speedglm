\name{add1.speedlm}
\alias{add1.speedlm}
\alias{drop1.speedlm}
\alias{extractAIC.speedlm}
\alias{nobs.speedlm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrappers to the add1 and drop1 methods for speedlm objects}

\description{
These are adviced to be used for \code{speedlm} models fitted on moderately large data sets. It is also possible to use \link[MASS]{stepAIC} function from package \code{MASS}.
}
\usage{
\method{add1}{speedlm}(object, scope, scale = 0, test = c("none", "Chisq","F"), 
                       x = NULL, k = 2, data, ...)
\method{drop1}{speedlm}(object, scope, scale = 0, all.cols = TRUE,
                        test = c("none","Chisq", "F"), k = 2, data, ...)
\method{extractAIC}{speedlm}(fit, scale = 0, k=2,...)
\method{nobs}{speedlm}(object, use.fallback = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{speedlm} object for which \code{model=TRUE} was previously set to.}
    \item{fit}{a \code{speedlm} object}
  \item{scope}{see \link[stats]{add1} from package \code{stats}.}
  \item{scale}{see \link[stats]{add1} from package \code{stats}.}
  \item{all.cols}{see \link[stats]{drop1} from package \code{stats}.}
  \item{test}{see \link[stats]{add1} from package \code{stats}.}
  \item{x}{see \link[stats]{add1} from package \code{stats}.}
  \item{k}{see \link[stats]{add1} from package \code{stats}.}
  \item{data}{the data that the model was previously fitted to. 
              If not provided, these will be searched in the parent environment.}
  \item{use.fallback}{logical. Should fallback methods be used to try to guess the value?}
  \item{\dots}{further optional arguments.}
}


\author{
Ronen Meiri
}


\examples{
\dontrun{
library(MASS)
set.seed(10)
n <- 1000
k <- 3
x <- round(matrix(rnorm(n * k), n, k), digits = 3)
beta <- c(0.05,0.5,0.8,1.3,1.8)
y <- drop(tcrossprod(cbind(1,x,x[,2]*x[,3]),t(beta))) + rnorm(n,,0.2)
colnames(x) <- c("s1", "s2", "s3") 
da <- data.frame(y, x)

m0 <- speedlm(y ~ 1, data = da)
m0.1 <- add1(m0,scope=~(s1+s2+s3)^2, data = da)
m1 <- stepAIC(m0,scope=~(s1+s2+s3)^3)
m1
}
}

