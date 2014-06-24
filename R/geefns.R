#' GEE estimating functions
#' 
#' Internal functions for computing the GEE. Should generally not be called by user.
#' @export
#' @param Y Vector of (correlated) outcomes
#' @param X Matrix of predictors
#' @param b Vector of coefficients
#' @param mu.Y Mean function
#' @param v.Y Variance function
#' @param aux Auxiliary function for coputing (co)variance parameters
#' @param id Vector of cluster IDs
#' @param uid Vector of unique subject IDs 
#' @param rows.indivs List of rows of \code{X} corresponding to each subject ID
#' @param corstr Working correlation structure
ee.GEE <- function(Y,X,b,mu.Y=function(eta){eta},v.Y=function(eta){rep(1,length(eta))},aux=NULL,id=1:length(Y),uid=sort(unique(id)),
                   rows.indivs=lapply(uid,function(j) { which(id==j)}),corstr="ind") {
  a <- aux[1]
  phi <- aux[2]
  
  eta <- t(X%*%b)
  
  if(corstr=="exch") {
    ##print("Setting up exchangeable working correlation matrix")    
    indivMats <- lapply(rows.indivs,function(rs) { 
      clust.size <- length(rs)
      sqA.i <- diag(sqrt(v.Y(eta[rs])))
      R <- stats::toeplitz(c(1,rep(a,clust.size-1)))
      phi*solve(sqA.i%*%R%*%sqA.i)
    })    
    Vinv <- as.matrix(Matrix::bdiag(indivMats))
  }
  else {
    Vinv <- phi*diag(as.vector(1/v.Y(eta)))
  }
  
  contrib.indiv <- function(ind,b,cols=1:length(b)) {
    
    if(length(ind)==1) {
      A <- v.Y(X[ind,cols]%*%b[cols])
      cont <- A*X[ind,cols]*Vinv[ind,ind]*(Y[ind] - mu.Y(X[ind,cols]%*%b[cols]))
    }
    else {
      A <- diag(as.vector(v.Y(X[ind,cols]%*%b[cols])))
      cont <- t(A %*% X[ind,cols]) %*% Vinv[ind,ind] %*% (Y[ind] - mu.Y(X[ind,cols]%*%b[cols]))
    }
    cont
  }
  
  L <- lapply(rows.indivs,contrib.indiv,b=b)
  
  apply(do.call("cbind",L),1,sum)
}

#' @describeIn ee.GEE
#' @export
ee.GEE.aux <- function(Y,X,b,mu.Y=function(eta){eta},v.Y=function(eta){rep(1,length(eta))},id=1:length(Y),uid=sort(unique(id)),
                       rows.indivs=lapply(uid,function(j) { which(id==j)})) {
  
  eta <- t(X%*%b)
  pearson.resids <- (Y - mu.Y(eta)) / sqrt(v.Y(eta))
  
  p <- sum(b!=0) ## Really, this p should be number of covariates in model
  
  num <- sum(unlist(lapply(rows.indivs,function(rs) {
    combs <- combn(rs,2)
    sum(pearson.resids[combs[1,]]*pearson.resids[combs[2,]])    
  })))
  
  ## This function could be moved "outside" because it doesn't need to be recomputed every time.
  denom <- sum(unlist(lapply(rows.indivs,function(rs) {
    n.i <- length(rs)
    (n.i*(n.i-1))/2
  })))- p
  
  phi <- 1/(sum(pearson.resids^2)/(length(Y)-p))
  
  alpha <- phi*num/denom ## This gives an estimate of pairwise correlation alpha
  return(c(alpha,phi)) ## We return alpha/phi since this is what we need in the original equation
  
}

expit <- function(x) { exp(x)/(1+exp(x))}

ee.GEELin <- function(...) { ee.GEE(...,mu.Y=function(eta){eta},v.Y=function(eta){rep(1,length(eta))}) }
ee.GEEBin <- function(...) { ee.GEE(...,mu.Y=function(eta) {expit(eta)},v.Y=function(eta) { expit(eta)/(1+exp(eta))}) }
ee.GEEPois <- function(...) {ee.GEE(...,mu.Y=function(eta){exp(eta)},v.Y=function(eta){exp(eta)})}

ee.GEELin.aux <- function(...) { ee.GEE.aux(...,mu.Y=function(eta){eta},v.Y=function(eta){rep(1,length(eta))}) }
ee.GEEBin.aux <- function(...) { ee.GEE.aux(...,mu.Y=function(eta) { expit(eta)},v.Y=function(eta) { expit(eta)/(1+exp(eta))}) }
ee.GEEPois.aux <- function(...) {ee.GEE.aux(...,mu.Y=function(eta){exp(eta)},v.Y=function(eta){exp(eta)})}
