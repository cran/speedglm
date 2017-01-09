\name{speedlm}
\alias{speedlm}
\alias{speedlm.fit}
\alias{speedlm.wfit}
\alias{update.speedlm}
\alias{updateWithMoreData}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Linear Models to Large Data Sets}
\description{
The functions of class 'speedlm' may speed up the fitting of LMs to large data
sets. High performances can be obtained especially if R is linked against an
optimized BLAS, such as ATLAS.
}
\usage{
# S3 method of class 'data.frame'
speedlm(formula, data, weights = NULL, offset = NULL, sparse = NULL, 
        set.default = list(), method=c('eigen','Cholesky','qr'), 
        model = FALSE, y = FALSE, fitted = FALSE, subset=NULL, ...)

# S3 method of class 'matrix'
speedlm.fit(y, X, intercept = FALSE, offset = NULL, row.chunk = NULL, 
            sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
            tol.solve = .Machine$double.eps, sparse = NULL, tol.values = 1e-07, 
            tol.vectors = 1e-07, method=c('eigen','Cholesky','qr'), ...)

speedlm.wfit(y, X, w, intercept = FALSE, offset = NULL, row.chunk = NULL, 
             sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
             tol.solve = .Machine$double.eps, sparse = NULL, tol.values = 1e-07, 
             tol.vectors = 1e-07, method=c('eigen','Cholesky','qr'), ...)
                      
# S3 method of class 'speedlm' (object) and 'data.frame' (data)                    
\method{update}{speedlm}(object, formula, data, add=TRUE, evaluate=TRUE, 
                         subset=NULL,...)

# S3 method of class 'speedlm' (object) and 'data.frame' (data)                    
updateWithMoreData(object, data, weights = NULL, offset = NULL, sparse = NULL, 
                   all.levels = FALSE, set.default = list(), subset=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
Most of arguments are the same of functions \link[stats]{lm} but with 
some difference.  
  \item{formula}{the same of function \code{lm}.}
  \item{data}{the same of function \code{lm}, but it must always specified.}
  \item{weights}{the same of function \code{lm}, but it must be specified as
       \code{data$weights}.}
  \item{w}{the same of \code{weights}.}  
  \item{intercept}{a logical value which indicates if an intercept is used.} 
  \item{offset}{the same of function \code{lm}.}
  \item{X}{the same of \code{x} in function \code{lm}.}
  \item{y}{the same of \code{lm},\code{lm.wfit} and \code{lm.fit} .}
  \item{sparse}{logical. Is the model matrix sparse? By default is NULL, so a
        quickly sample survey will be made.}
  \item{set.default}{a list in which to specify the parameters to pass to the
        functions \link{cp}, \link{control} and \link{is.sparse}.}
  \item{sparselim}{a value in the interval [0, 1]. It indicates the minimal 
        proportion of zeroes, in the model matrix X, in order to consider X as
        sparse.}
  \item{camp}{see function \code{is.sparse}.}
  \item{eigendec}{logical. Do you want to investigate on rank of X? You may set
        it to false if you are sure that X is full rank.}
  \item{row.chunk}{an integer, see the function \code{cp} for details.}
  \item{tol.solve}{see function \link{solve}.}
  \item{tol.values}{see function \link{control}.}
  \item{tol.vectors}{see function \link{control}.}
  \item{method}{see function \link{control}.}
  \item{object}{an object of class 'speedlm'.}
  \item{all.levels}{are all levels of eventual factors present in each data
        chunk? If so, set \code{all.levels} to true to speed up the fitting.}
  \item{model}{logical. Should the model frame be returned?}  
  \item{fitted}{logical. Should the fitted values be returned? }    
  \item{subset}{the same of function \code{lm}} 
  \item{add}{logical. Are additional data coming from a new chunk provided?}
  \item{evaluate}{logical. If true evaluate the new call else return the call.}
  \item{\dots}{further optional arguments.}
  
}
\details{Unlikely from \link[stats]{lm} or \link[biglm]{biglm}, the functions of class 'speedlm' do not use 
 the QR decomposition but directly solve the normal equations. Further, the most recent version of the package include method 'qr'. However such qr decomposition is not applied directly on matrix X, but on X'WX.   
 In some extreme case, this might have some problem of numerical stability but may take advantage from the use of 
 an optimized BLAS. The memory size of an object of class 'speedlm' is \eqn{O(p^2)}, where \eqn{p} is the number of covariates.
 If an optimized BLAS library is not installed, an attempt to speed up calculations may be done by setting \code{row.chunk} 
 to some value, usually less than 1000, in \code{set.default}. See the function \link{cp} for details. Factors are permitted 
 without limitations. \cr
In the most recent versions, function \code{update.speedlm} is now a wrapper to call either \code{updateWithMoreData}  (the new name of the old \code{update.speedlm}, for additional data chunks), or \link[stats]{update} from package \code{stats}. 
}
\value{
  \item{coefficients}{the estimated coefficients.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{XTX}{the product X'X (weighted, if the case).}
  \item{A}{the product X'X (weighted, if the case) not checked for singularity.}
  \item{Xy}{the product X'y (weighted, if the case).}
  \item{ok}{the set of column indeces of the model matrix where the model has
            been fitted.}
  \item{rank}{the numeric rank of the fitted linear model.}
  \item{pivot}{see the function \link{control}.}
  \item{RSS}{the estimated residual sums of squares of the fitted model.}
  \item{sparse}{a logical value indicating if the model matrix is sparse.}
  \item{deviance}{the estimated deviance of the fitted model.}
  \item{weigths}{the weights used in the last updating.}
  \item{zero.w}{the number of non-zero weighted observations.}
  \item{nobs}{the number of observations.}
  \item{nvar}{the number of independent variables.}
  \item{terms}{the \code{terms} object used.}
  \item{intercept}{a logical value which indicates if an intercept has been used.} 
  \item{call}{the matched call.}
  \item{model}{Either NULL or the model frame, if \code{model} was previously set to TRUE. }  
  \item{y}{Either NULL or the response variable, if \code{y} was previously set to TRUE. }  
  \item{fitted.values}{Either NULL or the fitted values, if \code{fitted} was previously set to TRUE. }    
  \item{\dots}{others values necessary to update the estimation.}
  
}
\references{Enea, M. (2009) Fitting Linear Models and Generalized Linear Models With Large Data Sets in R.
            In \emph{book of short papers, conference on ``Statistical Methods for the analysis of large data-sets'',
            Italian Statistical Society}, Chieti-Pescara, 23-25 September 2009, 411-414.\cr
            
            Klotz, J.H. (1995) Updating Simple Linear Regression. \emph{Statistica Sinica}, \bold{5}, 399-403.\cr
            
            Bates, D. (2009) Comparing Least Square Calculations. Technical report. Available at 
            \url{http://cran.rakanu.com/web/packages/Matrix/vignettes/Comparisons.pdf}\cr
            
            Lumley, T. (2009) biglm: bounded memory linear and generalized linear models.\emph{ R package version 0.7}
            \url{https://CRAN.R-project.org/package=biglm}.
             }
\author{ Marco Enea, with contribution from Ronen Meiri.}
\note{All the above functions make an object of class 'speedlm'.}
 %Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{ \link{summary.speedlm},\link{speedglm}, \link[stats]{lm}, and \link[biglm]{biglm} }
\examples{
\dontrun{
n <- 1000
k <- 3
y <- rnorm(n)
x <- round(matrix(rnorm(n * k), n, k), digits = 3)
colnames(x) <- c("s1", "s2", "s3") 
da <- data.frame(y, x)
do1 <- da[1:300,]
do2 <- da[301:700,]
do3 <- da[701:1000,]

m1 <- speedlm(y ~ s1 + s2 + s3, data = do1)
m1 <- update(m1, data = do2)
m1 <- update(m1, data = do3)

m2 <- lm(y ~ s1 + s2 + s3, data = da)
summary(m1)
summary(m2)
}

\dontrun{
  # as before but recursively
  make.data <- function(filename, chunksize,...){       
    conn <- NULL
    function(reset=FALSE, header=TRUE){
      if(reset){
        if(!is.null(conn)) close(conn)
         conn<<-file(filename,open="r") 
      } else{
        rval <- read.table(conn, nrows=chunksize,header=header,...)
        if (nrow(rval)==0) {
          close(conn)
          conn<<-NULL
          rval<-NULL
        }
        return(rval)
      }
    }
  }

  write.table(da,"da.txt",col.names=TRUE,row.names=FALSE,quote=FALSE) 
  x.names <- c("s1", "s2", "s3")
  dat <- make.data("da.txt",chunksize=300,col.names=c("y",x.names))
  dat(reset=TRUE) 
  da2 <- dat(reset=FALSE)
  
  # the first model runs on the first 300 rows.
  m3 <- speedlm(y ~ s1 + s2 + s3, data=da2)

  # the last three models run on the subsequent 300, 300 and 100 rows, respectively
  for (i in 1:3){
    da2 <- dat(reset=FALSE, header=FALSE)
    m3 <- update(m3, data=da2, add=TRUE)
  }  
  all.equal(coef(m1),coef(m3))  
  file.remove("da.txt")
}

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models 	}

