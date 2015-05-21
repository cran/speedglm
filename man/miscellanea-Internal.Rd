\name{control}
\alias{control}
\alias{cp}
\alias{is.sparse}




%- Also NEED an '\alias' for EACH other topic documented here.
\title{Miscellanea of functions}
\description{
Utility functions for least squares estimation in large data sets.
}
\usage{
control(B, symmetric = TRUE, tol.values = 1e-7, tol.vectors = 1e-7,
       out.B = TRUE, method = c("eigen", "Cholesky"))
cp(X, w = NULL, row.chunk = NULL, sparse = FALSE)
is.sparse(X, sparselim = .9, camp = .05)               
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{a squared matrix.}
  \item{symmetric}{logical, is \code{B} symmetric?}
  \item{tol.values}{tolerance to be consider eigenvalues equals to zero.}
  \item{tol.vectors}{tolerance to be consider eigenvectors equals to zero.}
  \item{out.B}{Have the matrix B to be returned?}
  \item{method}{the method to check for singularity. By default is "eigen", and
                an eigendecomposition of X'X is made. The "Cholesky" method is 
                faster than "eigen" and does not use tolerance, but the former 
                seems to be more stable for opportune tolerance values.}
  \item{X}{the model matrix.}
  \item{w}{a weights vector.}
  \item{sparse}{logical, is \code{X} sparse?}  
  \item{sparselim}{a real in the interval [0; 1]. It indicates 
                the minimal proportion of zeroes in the data matrix X in order 
                to consider X as sparse } 
                eigendec Logical. Do you want to investigate on rank of X? You may set to
  \item{row.chunk}{an integer which indicates the total rows number 
       compounding each of the first g-1 blocks. If \code{row.chunk} is not a divisor 
       of nrow(X), the g-th block will be formed by the remaining data.}
  \item{camp}{the sample proportion of elements of X on which the survey will be based.}       
}
\details{
   Function \code{control} makes an eigendecomposition of B according established values of tolerance. 
   Function \code{cp} makes the cross-product X'X by partitioning X in row-blocks. 
   When an optimized BLAS, such as ATLAS, is not installed, the function represents an attempt 
   to speed up the calculation and avoid overflows with medium-large data sets loaded in R memory.
   The results depending on processor type. Good results are obtained, for example, with an AMD Athlon 
   dual core 1.5 Gb RAM by setting \code{row.chunk} to some value less than 1000. Try the example below 
   by changing the matrix size and the value of \code{row.chunk}. If the matrix X is sparse, it will have 
   class "dgCMatrix" (the package Matrix is required) and the cross-product will be made without 
   partitioning. However, good performances are usually obtained with a very 
   high zeroes proportion. 
   Function \code{is.sparse} makes a quick sample survey on sample proportion of zeroes in X.
   }
\value{
for the function \code{control}, a list with the following elements:
  \item{XTX}{the matrix product B without singularities (if there are).}
  \item{rank}{the rank of B}
  \item{pivot}{an ordered set of column indeces of B with, if the case, the last \eqn{rank+1,...,p}  
               columns which indicate possible linear combinations.}
for the function \code{cp}:  
  \item{new.B}{the matrix product X'X (weighted, if \code{w} is given).}
for the function \code{is.sparse}:  
  \item{sparse}{a logical value which indicates if the sample proportion of zeroes is 
                greater than \code{sparselim}, with the sample proportion as attribute.}
}
%\references{ put references to the literature/web site here  }
\author{Marco ENEA}
\seealso{ \link{eigen}, \link{chol}, \link{qr}, \link{crossprod}}
\examples{
#### example 1.

n <- 100000
k <- 100
x <- round(matrix(rnorm(n*k),n,k),digits=4)
y <- rnorm(n)

# if an optimized BLAS is not installed, depending on processor type, cp() may be 
# faster than crossprod() for large matrices.

system.time(a1 <- crossprod(x))
system.time(a2 <- cp(x,,row.chunk = 500))
all.equal(a1, a2)  

#### example 2.1.
n <- 100000
k <- 10
x <- matrix(rnorm(n*k),n,k)
x[,2] <- x[,1] + 2*x[,3]  # x has rank 9
y <- rnorm(n)

# estimation by least squares 
A <- function(){
  A1 <- control(crossprod(x))
  ok <- A1$pivot[1:A1$rank]
  as.vector(solve(A1$XTX,crossprod(x[,ok],y)))
}
# estimation by QR decomposition
B <- function(){
  B1 <- qr(x)
  qr.solve(x[,B1$pivot[1:B1$rank]],y)    
}  
system.time(a <- A())
system.time(b <- B())

all.equal(a,b)

###  example 2.2
x <- matrix(c(1:5, (1:5)^2), 5, 2)
x <- cbind(x, x[, 1] + 3*x[, 2])
m <- crossprod(x)
qr(m)$rank # is 2, as it should be
control(m,method="eigen")$rank # is 2, as it should be
control(m,method="Cholesky")$rank # is wrong


### example 3. 
n <- 10000
fat1 <- gl(20,500)
y <- rnorm(n)
da <- data.frame(y,fat1)
m <- model.matrix(y ~ factor(fat1),data = da)
is.sparse(m)





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

