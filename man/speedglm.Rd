\name{speedglm}
\alias{speedglm.wfit}
\alias{speedglm}
\alias{shglm}
%\alias{shglm_data}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Generalized Linear Models for Large Data Sets}
\description{
\code{speedglm} and \code{speedglm.wfit} fit GLMs to medium-large 
data sets, that is those storable into the R memory. The highest performances, in terms of computation time,
are obtained when R is linked against an optimized BLAS, such as ATLAS. The function \code{shglm} 
is for a data set stored into a file of size greater than the available memory, and takes as 
argument a function to manipulate connections. 
}
\usage{
## S3 method for class 'data.frame':
speedglm(formula,data,family=gaussian(),weights=NULL,start=NULL,
         etastart=NULL,mustart=NULL,offset=NULL,maxit=25, k=2, 
         sparse=NULL,set.default=list(), trace=FALSE,
         method=c('eigen','Cholesky','qr'), model=FALSE, y=FALSE, 
         fitted=FALSE,...)

## S3 method for class 'matrix':       
speedglm.wfit(y, X, intercept=TRUE, weights=NULL,row.chunk=NULL,
              family=gaussian(), start=NULL, etastart=NULL,
              mustart=NULL, offset=NULL, acc=1e-08, maxit=25, k=2,
              sparselim=.9,camp=.01, eigendec=TRUE, tol.values=1e-7,
              tol.vectors=1e-7, tol.solve=.Machine$double.eps,
              sparse=NULL,method = c('eigen','Cholesky','qr'), 
              trace=FALSE,...)
              
## S3 method for class 'function':              
shglm(formula, datafun, family = gaussian(), weights.fo = NULL, start = NULL, 
      etastart = NULL, mustart = NULL, offset = NULL, maxit = 25, k = 2, 
      chunksize = 5000, sparse = NULL, trace = FALSE, all.levels = FALSE,
      set.default = list(),...)
     
}
%- maybe also 'usage' for other objects documented here.
\arguments{
Most of arguments are the same of \link[stats]{glm} or \link[biglm]{bigglm} 
but with some difference.     
  \item{formula}{the same of \code{glm}.}
  \item{data}{a data frame.}
  \item{datafun}{a function which uses connections. See the example below.}     
  \item{family}{the same of \code{glm}, but it must be specified with brackets.}
  \item{start}{the same of \code{glm}.}
  \item{weights}{the same of \code{glm}, but it must be specified as \code{data$weights}.}
  \item{weights.fo}{weights for the response. It must be specified as a formula 
  (see the example below).}
  \item{etastart}{the same of \code{glm}.}
  \item{mustart}{the same of \code{glm}.}
  \item{offset}{the same of \code{glm}.}
  \item{intercept}{the same of \code{glm}.}
  \item{X}{the same of \code{x} in \code{glm.fit}.}
  \item{y}{the same of \code{glm} and \code{glm.fit}.}
  \item{maxit}{the same of \code{glm}.}
  \item{k}{numeric, the penalty per parameter to be used; the default k = 2 is
  the classical AIC.}
  \item{trace}{logical. Do you want to be informed about the model estimation 
  progress?}
  \item{sparse}{logical. Is the model matrix sparse? By default is NULL, so a 
  quickly sample survey will be made.}
  \item{chunksize}{an integer indicates the number of rows of the data file to 
  read at time.}
  \item{all.levels}{logical, are all factor's levels present in each data chunk?}
  \item{set.default}{a list in which to specify the below parameters.}
  \item{sparselim}{a real in the interval [0, 1]. It indicates the minimal 
  proportion of zeroes in the data matrix X in order to consider X as sparse.}   
  \item{camp}{see the function \link{is.sparse}.}
  \item{eigendec}{logical. Do you want to check the rank of X? You may set it to 
  false if you are sure that X is full rank.}  
  \item{row.chunk}{an integer, see the function \link{cp} for details.}
  \item{acc}{tolerance to be used for the estimation.} 
  \item{tol.solve}{see the function \link{solve}.}
  \item{tol.values}{see the function \link{control}.}
  \item{tol.vectors}{see the function \link{control}.}
  \item{method}{the chosen method to detect for singulatity.}  
  \item{model}{logical. If TRUE the model frame will be returned. }  
  \item{fitted}{logical. If TRUE the fitted values will be returned.}  
  \item{\dots}{further optional arguments.}
  
}
\details{The function \code{shglm} works like
\code{biglm}, but it checks for singularity and does not impose restrictions on 
factors. Since during the IWLS estimation \code{shglm} uses repeated accesses to data file stored, for example, 
into the hard disk, the estimation time could be very long. 
Unlike from \code{glm} or \code{biglm}, the functions of class 'speedglm' do not use 
 the QR decomposition, but directly solve the equations in the form of Iterative(-ly) (Re-)Weighted Least Squares (IWLS). 
The memory size of an object of class 'speedglm' is \eqn{O(p^2)}, where \eqn{p} is the number of covariates, unless one or more of argument \code{model}, \code{y} and \code{fitted} are set to TRUE. If an optimized BLAS 
 is not installed, an attempt to speed up calculations might be done by setting \code{row.chunk} 
 to some value, usually less than 1000, in \code{set.default}. See the function \link{cp} for details.\cr
 If the model matrix is (very) sparse, the package \code{Matrix} could be used.
 Note that if method 'qr' is chosen, then the qr decomposition will not be applied on matrix X, as in \code{lm}, 
 but on X'WX.   
}
\value{
  \item{coefficients}{the estimated coefficients.}
  \item{logLik}{the log likelihood of the fitted model.}
  \item{iter}{the number of iterations of IWLS used.}
  \item{tol}{the maximal value of tolerance reached.}
  \item{convergence}{a logical value which indicates if convergence was reached.}
  \item{family}{the family object used.}
  \item{link}{the link function used.}
  \item{df}{the degrees of freedom of the model.}
  \item{XTX}{the product X'X (weighted, if the case).}
  \item{dispersion}{the estimated dispersion parameter of the model.}
  \item{ok}{the set of column indeces of the model matrix where the model has
           been fitted.}
  \item{rank}{the rank of the model matrix.}  
  \item{RSS}{the estimated residual sum of squares of the fitted model.}
  \item{aic}{the estimated Akaike Information Criterion.}
  \item{sparse}{a logical value which indicates if the model matrix is sparse.}
  \item{deviance}{the estimated deviance of the fitted model.}
  \item{nulldf}{the degrees of freedom of the null model.}
  \item{nulldev}{the estimated deviance of the null model.}
  \item{ngoodobs}{the number of non-zero weighted observations.}
  \item{n}{the number of observations.}
  \item{intercept}{a logical value which indicates if an intercept has been used.}
  \item{terms}{the terms object used.}
  \item{call}{the matched call.}
  \item{model}{Either NULL or, if \code{model} was previously set to TRUE, the model frame. }  
  \item{y}{Either NULL or, if \code{y} was previously set to TRUE, the response variable. }  
  \item{linear.predictors}{Either NULL or, if \code{fitted} was previously set to TRUE, the fitted values. }    
  \item{offset}{the model offset.}
}

\references{Enea, M. (2009) Fitting Linear Models and Generalized Linear Models with large data sets in R.
            In \emph{book of short papers, conference on ``Statistical Methods for the analysis of large data-sets'',
            Italian Statistical Society}, Chieti-Pescara, 23-25 September 2009, 411-414.\cr
            
            Bates, D. (2009) Comparing Least Square Calculations. Technical report. \cr
            
            Lumley, T. (2009) biglm: bounded memory linear and generalized linear models. \emph{R package version 0.7}.
            \url{https://CRAN.R-project.org/package=biglm}.
            }
            
\author{ Marco Enea. Ronen Meiri contributed with method 'qr'}
\note{All the above functions make an object of class 'speedglm'.\cr
In the current package version, arguments \code{start}, \code{mustart} and \code{etastart} of function \code{shglm} have been disabled. 
These will be restored in future.}


 %Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{ \link{speedlm}, \link[biglm]{bigglm}, \link[stats]{glm} }
\examples{

\dontrun{
# The following comparison among glm(), bigglm() and speedglm() cannot be considered rigorous 
# and exhaustive, but it is only to give an idea of the computation time. 
# It may take a long time.
require(biglm)
n<-50000
k<-80
y <- rgamma(n,1.5,1)
x <-round( matrix(rnorm(n*k),n,k),digits=3)
colnames(x) <-paste("s",1:k,sep = "")
da<- data.frame(y,x)
fo <- as.formula(paste("y~",paste(paste("s",1:k,sep=""),collapse="+")))   

system.time(m1 <- glm(fo,data=da,family=Gamma(log)))
system.time(m2 <- bigglm(fo,data=da,family=Gamma(log)))
system.time(m3 <- speedglm(fo,data=da,family=Gamma(log)))

# You may also try speedglm when R is linked against an optimized BLAS,
# otherwise try to run the following function. In some computers, it is
# faster for large data sets.
system.time(m4 <- speedglm(fo,data=da,family=Gamma(log),set.default=list(row.chunk=1000)))
}
##################
\dontrun{
## An example of function using a connection to an out-memory file  
## This is a slightly modified version of the function from the bigglm's help page 
 make.data<-function(filename, chunksize,...){       
     conn<-NULL
     function(reset=FALSE){
     if(reset){
       if(!is.null(conn)) close(conn)
       conn<<-file(filename,open="r")
     } else{
       rval<-read.table(conn, nrows=chunksize,...)
       if ((nrow(rval)==0)) {
            close(conn)
            conn<<-NULL
            rval<-NULL
       }
       return(rval)
     }
  }
}


# data1 is a small toy dataset
data(data1)
write.table(data1,"data1.txt",row.names=FALSE,col.names=FALSE)
rm(data1) 

da<-make.data("data1.txt",chunksize=50,col.names=c("y","fat1","x1","x2"))

# Caution! make sure to close the connection once you have run command #1
da(reset=T) #1: opens the connection to "data1.txt"
da(reset=F) #2: reads the first 50 rows (out of 100) of the dataset
da(reset=F) #3: reads the second 50 rows (out of 100) of the dataset
da(reset=F) #4: is NULL: this latter command closes the connection


require(biglm)        
# fat1 is a factor with four levels                                    
b1<-shglm(y~factor(fat1)+x1,weights=~I(x2^2),datafun=da,family=Gamma(log))
b2<-bigglm(y~factor(fat1)+x1,weights=~I(x2^2),data=da,family=Gamma(log))
summary(b1) 
summary(b2) 

file.remove("data1.txt")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models 	}

