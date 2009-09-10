
`summary.speedlm` <-
function(object,correlation=FALSE,...)
{
  if (!inherits(object, "speedlm")) stop("object must be an object of class speedlm")
  z <- object 
  n <- if (is.null(z$weights)) z$n.obs else z$n.obs - z$zero.w
  nvar <- z$nvar
  rdf <- z$df.residual
  var_res <-  as.numeric(z$RSS)/rdf
  se_coef <- rep(NA,z$nvar) 
  inv <- solve(z$XTX)
  se_coef[z$ok] <- sqrt(var_res * diag(inv))  
  t1 <- z$coefficients / se_coef
  p  <- 2 * pt(abs(t1),df=z$df.residual,lower.tail=FALSE) 
  ip <-  !is.na(p)  
  if (is.null(z$weights)) {  
    X1X <- z$X1X[z$ok]
  #  if (attr(z$terms, "intercept")){
    if (z$intercept){
      X1X <- matrix(kronecker(X1X,X1X),z$rank,z$rank)  
      mss <- crossprod(z$coef,z$XTX-X1X/n)%*%z$coef
    } else  mss <-crossprod(z$coef,z$XTX)%*%z$coef      
  } else {
    XW1 <- z$XW1[z$ok]
#    mss <- if (attr(z$terms, "intercept")){  
   mss <- if (z$intercept){  
             XWX <- matrix(kronecker(XW1,XW1),length(XW1),length(XW1)) 
             XW1 <- matrix(kronecker(XW1*z$SW,XW1),z$rank,z$rank,byrow=TRUE)     
             crossprod(z$coef,z$XTX - 2*XWX/z$SW + XW1/(z$SW^2))%*%z$coef
           } else crossprod(z$coef,z$XTX)%*%z$coef 
  }   
  rss <- z$RSS
  if (nvar != (z$intercept)) {
    df.int <- if (z$intercept) 1L else 0L
    r.squared <- as.numeric(mss/(mss + rss))
    adj.r.squared <- 1 - (1 - r.squared) * ((n - df.int)/rdf)
    fstatistic <- c(value = (as(mss,"numeric")/(z$rank - df.int))/var_res, 
                    numdf = z$rank - df.int, dendf = rdf)
  } else  r.squared <- adj.r.squared <- 0
  f.pvalue <- 1-pf(fstatistic[1],fstatistic[2],fstatistic[3]) 
  param <- data.frame("coef"=z$coefficients,"se"=se_coef,"t"=t1,"p-value"=p)                                            
  keep <- match(c("call", "terms", "frame", "ok","RSS","rank"), names(object),0)
  ans <- c(object[keep], list(coefficients =param,"var.res"=var_res, 
           "df.residuals"=z$df.residual,"n.obs"=z$n.obs,"r.squared"=r.squared,
           "adj.r.squared"=adj.r.squared,correlation=correlation,
           fstatistic=fstatistic,f.pvalue=f.pvalue, rdf=rdf, cov.scaled = inv))
  if (correlation) 
    ans$correl <- inv * as.numeric(var_res) / outer(se_coef[z$ok],se_coef[z$ok])          
  class(ans) <- "summary.speedlm"
  return(ans)
}

`coef.speedlm` <-
function(object,...)  object$coefficients


`vcov.speedlm` <-
function(object,...)
{
  z <- object
  var_res <- z$RSS/z$df.residual
  rval<- var_res * solve(z$XTX)
  rval
}

logLik.speedlm <- function(object,...){
    p <- object$rank
    N <- object$n.obs
    if (is.null(object$pw)) pw <- 1 else  {
      N <- object$n.obs - object$zero.w
      pw <- object$pw
    }  
    N0 <- N
    val <- 0.5 * (log(pw) - N * (log(2 * pi) + 1 - log(N) +  log(object$RSS)))
    attr(val, "nall") <- N0
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val)<-"logLik.speedlm"
    val
}

print.logLik.speedlm <- function(x, digits = getOption("digits"), ...) {
    cat("'log Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
        " (df=", format(attr(x, "df")), ")\n", sep = "")
    invisible(x)
}

AIC.speedlm<-function (object,...,k = 2){ 
    p <- object$rank
    N <- object$n.obs
    if (is.null(object$pw)) pw <- 1 else  {
      N <- object$n.obs - object$zero.w
      pw <- object$pw
    }  
    val <- -(log(pw) - N * (log(2 * pi) + 1 - log(N) + log(object$RSS)))+k*(p+1)
    val
}


`print.summary.speedlm` <-
function(x,digits=max(3,getOption("digits")-3),...){

  x$coefficients$coef <- if (any(abs(na.omit(x$coefficients$coef))<0.0001))
                           format(x$coefficients$coef,scientific=TRUE,
                           digits=4) else  round(x$coefficients$coef,digits=6)
  x$coefficients$se <- if (any(na.omit(x$coefficients$se)<0.0001))
                         format(x$coefficients$se,scientific=TRUE,digits=4) else
                         round(x$coefficients$se,digits=6)
  x$coefficients$t <- round(x$coefficients$t,digit=3)
  x$coefficients$p.value <- if (any(na.omit(x$coefficients$p.value)<0.0001))
                              format(x$coefficients$p.value,scientific=TRUE,
                              digits=4) else round(x$coefficients$p.value,
                              digits=6)
  s<-sum(Vectorize(is.na(x$coefficients$coef)))
  cat("Linear Regression Model of class 'speedlm':\n")  
  if (!is.null(x$call)) cat("\nCall: ", deparse(x$call), "\n\n")  
  if (length(x$coef)) {
        cat("Coefficients:\n")
        print(x$coefficients)
    } else cat("No coefficients\n")
  cat("---\n")
  cat("Residual standard error: ",round(sqrt(x$var.res),6)," on ",x$rdf,
      " degrees of freedom;\n",
      "observations: ",x$n.obs,";  R^2: ",format(x$r.squared,digits=3),
      ";  adjusted R^2: ",format(x$adj.r.squared,digits=3),";\n",
      "F-statistic: ",format(x$fstatistic[1],digits=4)," on ", x$fstatistic[2],
      " and ",x$fstatistic[3]," df;  p-value: ",format(x$f.pvalue,digits=6),
      ".\n",sep="")
  if (s==1) cat("One coefficient not defined because of singularities. \n")
  if (s>1) cat(s," coefficients not defined because of singularities. \n")
  invisible(x)
  if (x$correlation) {
    cat("---\n")
    cat("Correlation of Coefficients:\n")
    x$correl[upper.tri(x$correl,diag=TRUE)]<-NA
    print(x$correl[-1,-nrow(x$correl)],na.print="",digits=2)
  }
}

`print.speedlm` <-
function(x,digits = max(3, getOption("digits") - 3),...)
{
  cat("Linear Regression Model of class 'speedlm':\n")
  if (!is.null(x$call)) cat("\nCall: ", deparse(x$call), "\n\n")
  if (length(x$coef)) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2,
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}


