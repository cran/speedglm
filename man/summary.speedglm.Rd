\name{summary.speedglm}
\alias{summary.speedglm}
\alias{coef.speedglm}
\alias{vcov.speedglm}
\alias{logLik.speedglm}
\alias{AIC.speedglm}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods to summarize Generalized Linear Models fits}
\description{
\code{summary} method for the class 'speedglm'.
}
\usage{

\method{summary}{speedglm}(object,correlation=FALSE,...)
\method{coef}{speedglm}(object,...)
\method{vcov}{speedglm}(object,...)
\method{logLik}{speedglm}(object,...)
\method{AIC}{speedglm}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class 'speedglm'.}
  \item{correlation}{logical. Do you want to print the correlation matrix? By default it is false.}
  \item{\dots}{further optional arguments}
}
\value{
  \item{coefficients}{the matrix of coefficients, standard errors, z-statistics and two-side p-values.}
  \item{df.residual}{the component from object.}
  \item{df.null}{the component from object.}
  \item{null.deviance}{the component from object.}
  \item{deviance}{the component from object.}
  \item{family}{the component from object.}
  \item{call}{the component from object.}
  \item{AIC}{the Akaike Information Criterion.}
  \item{RSS}{Residuals sums of squares.}
  \item{correlation}{(only if \code{correlation} is true.) The correlations of the estimated coefficients.}
  \item{logLik}{the log-likelihood value.}
  \item{rank}{the component from object.}
  \item{dispersion}{the estimated dispersion parameter of the fitted model.}
  \item{convergence}{the component from object.}    
  \item{iter}{the component from object.}
  \item{tol}{the component from object.}

}
%\references{ put references to the literature/web site here  }
\author{ Marco ENEA }
%\note{ further notes 

 %Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{ \link{speedglm} }
\examples{
n<-1000
k<-5
y <- rgamma(n,1.5,1)
x <-round( matrix(rnorm(n*k),n,k),digits=3)
colnames(x) <-paste("s",1:k,sep = "")
da<- data.frame(y,x)
fo <- as.formula(paste("y~",paste(paste("s",1:k,sep=""),collapse="+")))

m4 <- speedglm(fo,data=da,family=Gamma(log))
summary(m4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

