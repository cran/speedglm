\name{add1.speedlm}
\alias{add1.speedlm}
\alias{drop1.speedlm}
\alias{extractAIC.speedlm}
\alias{nobs.speedlm}
\alias{add1.speedglm}
\alias{drop1.speedglm}
\alias{extractAIC.speedglm}
\alias{nobs.speedglm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrappers to the add1 and drop1 methods for speedlm and speedglm objects}

\description{
These are adviced to be used for \code{speedlm} and \code{speedglm} models fitted on moderately large data sets. It is also possible to use \link[MASS]{stepAIC} function from package \code{MASS}.
}
\usage{
## S3 method for class 'speedlm'
\method{add1}{speedlm}(object, scope, scale = 0, test = c("none", "Chisq","F"), 
                       x = NULL, k = 2, data, ...)
\method{drop1}{speedlm}(object, scope, scale = 0, all.cols = TRUE,
                        test = c("none","Chisq", "F"), k = 2, data, ...)
\method{extractAIC}{speedlm}(fit, scale = 0, k=2,...)
\method{nobs}{speedlm}(object, use.fallback = FALSE, ...)

## S3 method for class 'speedglm'
\method{add1}{speedglm}(object, scope, scale = 0, test = c("none", "Rao", "LRT", 
                         "Chisq", "F"), x = NULL, k = 2, weights=rep(1,object$n), ...) 
\method{drop1}{speedglm}(object, scope, scale = 0, test = c("none", "Rao", "LRT", 
                         "Chisq", "F"), k = 2, weights=rep(1,object$n), ...)
\method{extractAIC}{speedglm}(fit, scale = 0, k=2,...)
\method{nobs}{speedglm}(object, use.fallback = FALSE, ...) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{speedlm} or \code{speedglm} object for which \code{model=TRUE} was previously set to.}
    \item{fit}{a \code{speedlm} or \code{speedglm} object}
  \item{scope}{see \link[stats]{add1} from package \code{stats}.}
  \item{scale}{see \link[stats]{add1} from package \code{stats}.}
  \item{all.cols}{see \link[stats]{drop1} from package \code{stats}.}
  \item{test}{see \link[stats]{add1} from package \code{stats}.}
  \item{x}{see \link[stats]{add1} from package \code{stats}.}
  \item{k}{see \link[stats]{add1} from package \code{stats}.}
  \item{data}{the data that the model was previously fitted to. 
              If not provided, these will be searched in the parent environment.}
  \item{weights}{the model weights, if provided in the speedglm object}
  \item{use.fallback}{logical. Should fallback methods be used to try to guess the value?}
  \item{\dots}{further optional arguments.}
}

\details{It is possible to use functions step() and stepAIC() for both speedlm and speedglm objects but objects fitted using updateWithMoreData()}
\author{
Ronen Meiri and Marco Enea
}
\section{Warnings}{Note that these functions have been poorly tested and need to be checked out more carefully.}


\examples{
\dontrun{
set.seed(10)
n <- 1000
k <- 3
x <- round(matrix(rnorm(n * k), n, k), digits = 3)
beta <- c(0.05,0.5,0.8,1.3,1.8)
y <- drop(tcrossprod(cbind(1,x,x[,2]*x[,3]),t(beta))) + rnorm(n,,0.2)
colnames(x) <- c("s1", "s2", "s3") 
da <- data.frame(y, x)

m0 <- speedlm(y ~ 1, data = da,model=TRUE, y=TRUE)
m0.1 <- add1(m0,scope=~(s1+s2+s3)^2, data = da)
m1 <- step(m0,scope=~(s1+s2+s3)^3)
m1

m2 <- speedglm(y ~ 1, data = da, model=TRUE, y=TRUE)
m2.1 <- add1(m2,scope=~(s1+s2+s3)^2, data = da)
m3 <- step(m2,scope=~(s1+s2+s3)^3)
m3
}
}

