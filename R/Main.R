#' Computes quasiconvex/quasiconcave least squares estimate with or without additional monotonicity constraint via CPLEX  or Gurobi.
#'
#' This function will take in the multivariate regressors and the response, and return the least squares estimates at the regressors, based on the response. X,y and Monotone must be supplied, though all other parameters have pre-set values the user can proceed with unless they wish to change the prior specification.
#' @param X              An n by d matrix of regressors.
#' @param y              An n by 1 vector of responses.
#' @param Shape       A categorical variable indicating the type of regression.
#'                       The user must input one of the following three types:
#'                       "QuasiConvex": means QuasiConvex regression -- with possible monotonicty constraint enforced by 'Monontone'
#'                       "QuasiConcave": means QuasiConcave regression -- with possible monotonicty constraint enforced by 'Monontone'
#' @param Monotone       A categorical variable indicating the type of  monotonicitty constraint on top of the shape constraint enforced by 'Shape' parameter
#'                       The user must input one of the following three types:
#'                       "non.inc": means additional nonincreasing.
#'                       "non.dec": means nondecreasing.
#'                       "non": means no additional monotonicity constraint
#' @param Max.b          Bound on the absolute values of the variables.
#' @param MM             The 'M' parameter in the big-M method applied by us.
#' @param ep             A small positive quantity to convert open constraints to close constraints.
#' @param time.limit     Time limit in seconds of call to optimizer. Can be any nonnegative number.
#' @param tol            Relative MIP optimality gap tolerance. Can be any nonnegative number. Default 1e-4.
#' @param Init.start     A matrix for warm start. Default is NULL. If provided should be matrix with n rows.
#' @param optimizer      Which optimizer to use when solving the MIQP.  The user must input one of the following two types: 'gurobi' or 'cplex'. Make sure they installed locally and accesbible frome the R envoirment.
#' @return               An object of class 'Quasi.reg' that contain among other things, the
#'                       least squares estimates at the regressors (within f.hat) and
#'                        the optimal objective function (within minvalue), list of solutions based on each of the intial starts (within FhatList), the initial start that lead to the the optimum (within K), and the corresponding convergence status (within status)
#' @export
#' @examples
#' n=20
#' d=4
#' library(MASS)
#' ## We demonstrate an example, where the regressor matrix X is Gaussian, and the response vector y
#' ## is a noisy version of a real-valued quasiconvex, decreasing function applied to the rows of X.
#' X = matrix(runif(n*d, 0, 3), nrow=n)
#' y = exp(-rowSums(X)) + rnorm(n)
#' ret = Quasi.reg(X, y,Shape ="QuasiConvex", Monotone = "non.inc",tol =  1e-06, optimizer= "gurobi")
#'
Quasi.reg <- function(X, y, Shape = c("QuasiConvex", "QuasiConcave"),  Monotone = c("non.inc", "non.dec","none"), Max.b = 10^4, MM =10^4, ep =0.001,  time.limit= 200, tol =  1e-06, Init.start = NULL, optimizer = c("cplex", "gurobi")){

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.vector(y)) y <- as.vector(y)
  dd <- ncol(X)
  if( !(Monotone %in% c("non.inc", "non.dec","none")) ){
    stop("Monotonicity constraint not correctly specified!")
  }
  if( !(Shape %in%  c("QuasiConvex", "QuasiConcave"))){
    stop("Shape must be either 'QuasiConvex' or 'QuasiConcave'.")
  }
  if( !(optimizer %in% c("cplex", "gurobi"))){
    stop("An optimizer (either 'cplex' or 'gurobi') must be chosen, installed, and loaded.")
  }

  if( optimizer == "cplex" && !is.null(Init.start)){
    warning("Currently 'Rcplex' does not support initial value. 'Init.start' is ignored.")
  }
  if(!is.matrix(Init.start) && !is.null(Init.start)) Init.start <- as.matrix(Init.start, nrow=nrow(X))
  if(Shape=="QuasiConvex"){
    Cons <-  QuasiConv.control(X, y, Max.b = Max.b, MM = MM, ep = ep, Monotone= Monotone)
  } else if (Shape == "QuasiConcave"){
    if( Monotone == "non.inc"){
      Monotone.inv <- "non.dec"
    } else if (Monotone == "non.dec"){
      Monotone.inv <- "non.inc"
    } else {
      Monotone.inv <- Monotone
    }
    Cons <-  QuasiConv.control(X, -y, Max.b = Max.b, MM = MM, ep = ep, Monotone= Monotone.inv )
  }
  out.list <- fhat.list <- vector("list", 1)
  obs.vec <-rep(0, 1)
  if (optimizer == "cplex"){
    t <- Sys.time()
    out.list[[1]]<- out <- Rcplex::Rcplex(cvec = Cons$cvec, Amat = Cons$Amat, bvec = Cons$bvec, Qmat = Cons$Qmat, lb = Cons$lb, ub = Cons$ub, objsense = "min",sense = "L", vtype = Cons$vtype, control = list(tilim = time.limit ,epgap = tol))
    fhat.list[[1]] <- opt.fhat <- out$xopt[(1:nrow(X))]
    least.squares.loss<- obs.vec[1] <- out.list[[1]]$obj+sum(y^2)
    # opt.zeta <- out$xopt[(nrow(X)+1) : (nrow(X)+ nrow(X)*dd)]
    out$xopt <- out.list[[1]]$xopt <- out.list[[1]]$xopt[(1:nrow(X))]## Removing the zetas for memory purposes (If you need to access the zeta's then comment this line)
    t1 <- Sys.time()
  } else if (optimizer== "gurobi"){
    t <- Sys.time()
    model<- list(obj = Cons$cvec, A = Cons$Amat, lb = Cons$lb, ub = Cons$ub, rhs = Cons$bvec, vtype = Cons$vtype, Qmat = 1/2*Cons$Qmat) ## QMAT had to be halved
    params <- list()
    params$TimeLimit <- time.limit
    ## uses Initial starts to compute all estimates and chooses one of them
    if(!is.null(Init.start)){
      Init.start<- as.matrix(Init.start, nrow= nrow(X))
      out.list  <- fhat.list <- vector("list", ncol(Init.start))
      obs.vec <-rep(0, ncol(Init.start))
      for (ii in 1:ncol(Init.start)){
        cat("\rMultistart ", ii, " of ", ncol(Init.start), " done!")
        cat("\n")
        model$start = c(Init.start[,ii], rep(NA, nrow(X)*dd + nrow(X)^2 ))
        out.list[[ii]] <- gurobi::gurobi(model, params)
        fhat.list[[ii]] <- out.list[[ii]]$x[(1:nrow(X))]
        out.list[[ii]]$x <- out.list[[ii]]$x[(1:nrow(X))] ## Removing the zetas for memory purposes (If you need to access the zeta's then comment this line)
        obs.vec[ii]<- out.list[[ii]]$objval+sum(y^2)
      }
      out <- out.list[[which.min(obs.vec)]]
      opt.fhat <- fhat.list[[which.min(obs.vec)]]
      least.squares.loss <- min(obs.vec)
    } else{
        out.list[[1]]<- out <- gurobi::gurobi(model, params)
        least.squares.loss<- obs.vec[1] <- out.list[[1]]$objval+sum(y^2)
        fhat.list[[1]] <- opt.fhat <- out.list[[1]]$x[(1:nrow(X))]
        out.list[[1]]$x <- out.list[[1]]$x[(1:nrow(X))] ## Removing the zetas for memory purposes (If you need to access the zeta's then comment this line)
    }

    t1 <- Sys.time()
  }
  if(Shape=="QuasiConvex"){
    f.hat <-  opt.fhat
  } else if (Shape=="QuasiConcave") {
    f.hat <-  -opt.fhat
    fhat.list <- lapply(fhat.list, function(x) -x)
  }
  ret <- NULL
  ret$out <-  out
  ret$fhat  <-   f.hat
  ret$minvalue  <-  least.squares.loss
  ret$time  <-   t1-t
  ret$tol <-   tol
  ret$method <-  optimizer
  ret$monotone <-  Monotone
  ret$Shape <-  Shape
  ret$status <-  out$status
  ret$K <- which.min(obs.vec)
  # Initial starts stored
  ret$ObjValVec <- obs.vec
  ret$OutList <- out.list
  ret$FhatList <- fhat.list
  ret$InitStart <- Init.start
  ret$call <- match.call()

  class(ret) <- "Quasi.reg"
  return(ret)
}

print.Quasi.reg <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("Estimate of f is:\n")
  print(x$fhat)
  cat("Number of Starting Vectors:\n")
  print(ncol(x$InitStart))
  cat("Initial vector leading to the minimum:\n")
  print(x$InitStart[,x$K])
  cat("Minimum Criterion Value Obtained:\n")
  print(x$minvalue)
  cat("Convergence Status:\n")
  print(x$out$status)
}





#'Checks whether the points in the vector y are actually the outputs of a single quasiconvex function applied to the corresponding rows of the matrix X.
#'
#'This function will take the matrix X and the vector y as inputs, and return a string,
#'indicating whether the vector y is or not a quasiconvex realization of X.
#'
#' @param X              An n by d matrix.
#' @param y              An n by 1 vector.
#' @param tol            A small positive quantity.
#'
#' @return               A string, indicating whether the vector y is, or not
#'                       a quasiconvex realization of X.
#' @export
#' @examples
#' \dontrun{
#' library(MASS)
#' ## We demonstrate an example, where the regressor matrix X is Gaussian, and the response vector y
#' ## is a real-valued quasiconvex function applied to the rows of X.
#' X = mvrnorm(20,rep(0,8), diag(8), tol=1e-06,empirical=FALSE)
#' y = exp(-rowSums(X))
#'
#' out = qconvcheck(X,y, tol = 1e-3)
#' }
qconvcheck<-function(X,y, tol = 1e-3){
  X<-as.matrix(X)
  y<-as.matrix(y)
  n<-nrow(X)
  d<-ncol(X)
  if(length(y)!=n)
  {
    print("ERROR in qconvcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  }
  else
  {
    output.qconv<-"Quasiconvex"
    for(i in 1:n)
    {
      if(sum(y<(y[i]-tol))>0)
      {
        Xchull<-X[(y<(y[i]-tol)),]
        if(sum(y<(y[i]-tol))==1)
        {
          Xchull<-t(as.matrix(Xchull))
        }
        flag<-chullcheck(Xchull,X[i,])
        if(flag==1)
        {
          output.qconv<-"Not Quasiconvex"
          break
        }
      }
    }
    return(output.qconv)
  }
}





#'Checks whether the point test is in the convex hull of the rows of X.
#'
#'This function will take the matrix X and the vector test as inputs,
#'and return a binary variable, which is 1 if the point test is in the
#'convex hull of the rows of X, and 0 otherwise.
#'
#' @param X              An n by d matrix.
#' @param test           An n by 1 vector.
#'
#' @return               A binary variable, which is 1 if the point test is in the
#'                       convex hull of the rows of X, and 0 otherwise.
#' @export
#' @importFrom lpSolve lp
#' @examples
#' \dontrun{
#' ## We demonstrate an example, where the point test is in the convex hull of the rows of X.
#' X<-matrix(c(0,0,1,1,0,1,0,1),4,2)
#' test<-c(0.3,0.4)
#'
#' out = chullcheck(X,test)
#' }
chullcheck<-function(X,test){
  X<-as.matrix(X)
  test<-as.matrix(test)
  n<-nrow(X)
  d<-ncol(X)
  if(length(test)!=d)
  {
    print("ERROR in chullcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  }
  else
  {
    objective.in<-rep(0,n)
    const.mat<-rbind(t(X),rep(1,n),diag(n))
    const.dir<-c(rep("=",(d+1)),rep(">=",n))
    const.rhs<-rbind(test,1,as.matrix(rep(0,n)))
    out<-lp(direction="min", objective.in,const.mat,const.dir,const.rhs)
    if(out$status==2)
    {
      return(0)
    }
    else
    {
      return(1)
    }
  }
}





#'Checks whether the point test is in the upper orthant of the convex hull of the rows of X.
#'
#'This function will take the matrix X and the
#'vector test as inputs, and return a binary variable, which is 1 if
#'the point test is in the upper orthant of the convex hull of the
#'rows of X, and 0 otherwise.
#'
#' @param X              An n by d matrix.
#' @param test           An n by 1 vector.
#'
#' @return               A binary variable, which is 1 if the point test
#'                       is in the upper orthant of the convex hull of the
#'                       rows of X, and 0 otherwise.
#' @export
#' @importFrom lpSolve lp
#' @examples
#' \dontrun{
#' ## We demonstrate an example, where the point test is in the upper orthant of the convex hull
#' ## of the rows of X.
#' X<-matrix(c(0,0,1,1,0,1,0,1),4,2)
#' test<-c(2,7)
#'
#' out = chullupcheck(X,test)
#' }
chullupcheck<-function(X,test){
  X<-as.matrix(X)
  test<-as.matrix(test)
  n<-nrow(X)
  d<-ncol(X)
  if(length(test)!=d){
    print("ERROR in chullcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  } else{
    objective.in<-rep(0,(n+d))
    const.mat1<-cbind(t(X),diag(d))
    const.mat2<-c(rep(1,n),rep(0,d))
    const.mat<-rbind(const.mat1,const.mat2,diag(n+d))
    const.dir<-c(rep("=",(d+1)),rep(">=",(n+d)))
    const.rhs<-rbind(test,1,as.matrix(rep(0,(n+d))))
    out<-lp(direction="min", objective.in,const.mat,const.dir,const.rhs)
    if(out$status==2){
      return(0)
    } else{
      return(1)
    }
  }
}





#'Checks whether the point test is in the lower orthant of the convex
#'hull of the rows of X.
#'
#'
#'
#'This function will take the matrix X and the
#'vector test as inputs, and return a binary variable, which is 1 if
#'the point test is in the lower orthant of the convex hull of the
#'rows of X, and 0 otherwise.
#'
#' @param X              An n by d matrix.
#' @param test           An n by 1 vector.
#'
#' @return               A binary variable, which is 1 if the point test
#'                       is in the lower orthant of the convex hull of the
#'                       rows of X, and 0 otherwise.
#' @export
#' @importFrom lpSolve lp
#' @examples
#' \dontrun{
#' ## We demonstrate an example, where the point test is in the lower orthant of the convex hull
#' ## of the rows of X.
#' X<-matrix(c(0,0,1,1,0,1,0,1),4,2)
#' test<-c(-1,-2)
#'
#' out = chulldowncheck(X,test)
#' }

chulldowncheck<-function(X,test){
  X<-as.matrix(X)
  test<-as.matrix(test)
  n<-nrow(X)
  d<-ncol(X)
  if(length(test)!=d){
    print("ERROR in chullcheck: Dimension of the testing point does not match the dimension of the rows of the base matrix!")
  } else{
    objective.in<-rep(0,(n+d))
    const.mat1<-cbind(t(X),-diag(d))
    const.mat2<-c(rep(1,n),rep(0,d))
    const.mat<-rbind(const.mat1,const.mat2,diag(n+d))
    const.dir<-c(rep("=",(d+1)),rep(">=",(n+d)))
    const.rhs<-rbind(test,1,as.matrix(rep(0,(n+d))))
    out<-lp(direction="min", objective.in,const.mat,const.dir,const.rhs)
    if(out$status==2){
      return(0)
    } else{
      return(1)
    }
  }
}





#'Interpolates function under shape and/or monotocity constraints
#'
#'Assumes that 'fin' satisfy some monotonicity and shape constraint
#'fit obtained from the matrix X of regressors and some response vector.
#'This function will take as inputs, the matrix X, the fitted vector fin,
#'a matrix nondatapoint_matrix whose rows may not be among those of X, and a
#'categorical variable Shape which denotes the type of regression in terms of
#'shape, and Monotone which denotes the type of regression in terms of
#'monotonicity. It will return a quasiconvex/quasiconcave/none and monotone (type of
#'monotonicity determined by the user entered option for the variable Monotone) function
#'evaluated at the rows of the matrix nondatapoint_matrix, which is an
#'extrapolation of the fitted values in the vector fin. This function is user friendly wrapper for 'interpol.func'.
#'
#' @param X                         An n by d matrix of regressors.
#' @param fin                An n by 1 vector of fitted values.
#' @param nondatapoint_matrix       An N by d matrix whose rows may not be among the rows of X.
#'@param Shape       A categorical variable indicating the type of regression.
#'                       The user must input one of the following three types:
#'                       "QuasiConvex": means QuasiConvex regression -- with possible monotonicty constraint enforced by 'Monontone'
#'                       "QuasiConcave": means QuasiConcave regression -- with possible monotonicty constraint enforced by 'Monontone'
#'                      "none": means there is neither
#' @param Monotone                  A categorical variable denoting the type of regression in terms of monotonicity, in addition to possible shape constraint imposed by the 'Shape' parameter
#'                                  The user must input one of the following three types:
#'                                  "non.inc": means nonincreasing
#'                                  "non.dec": means nondecreasing
#'                                  "non": means quasiconvex regression.
#'
#' @return                          A vector, representing a quasiconvex and monotone (type of monotonicity
#'                                  determined by the user entered option for the input variable Monotone)
#'                                  function evaluated at the rows of the matrix nondatapoint_matrix,
#'                                  which is an extrapolation of the fitted values in the vector fin.
#' @export
# 1
# # ' @examples
# # ' library(MASS)
# # ' X =  matrix(runif(40*3, 0, 3), nrow=40)
# # ' y = exp(-rowSums(X)) + rnorm(20)
# # '
# # ' ret = Quasi.reg(X, y,Shape ="QuasiConvex", Monotone = "non.inc", tol =  1e-06, optimizer= "gurobi")
# # ' Xnondata = matrix(runif(20*3, 0, 2), nrow=20)
# # ' out=interpolation.matrix.func(X,ret$f.hat,Xnondata,Monotone="non.inc")

interpolation.matrix.func<-function(X,fin,nondatapoint_matrix,Shape = c("QuasiConvex", "QuasiConcave", "none"), Monotone = c("non.inc", "non.dec","none"))
{
  if(Shape=="QuasiConvex")
  {
  int.f<-function(ndp)
  {
    return(interpol.func(X,fin,ndp,Shape="QuasiConvex", Monotone))
  }
  return(apply(nondatapoint_matrix,1,int.f))
  }
  else if(Shape=="QuasiConcave")
  {
    Mon<-"none"
    if(Monotone=="non.inc")
    {
      Mon<-"non.dec"
    }
    else if (Monotone=="non.dec")
    {
      Mon<-"non.inc"
    }
    int.f<-function(ndp)
    {
      return(-interpol.func(X,-fin,ndp, Shape= "QuasiConvex", Mon))
    }
    return(apply(nondatapoint_matrix,1,int.f))
  } else if (Shape == "none"){
     int.f<-function(ndp)
    {
      return(interpol.func(X,fin,ndp, Shape= "none", Monotone))
    }
    return(apply(nondatapoint_matrix,1,int.f))
  }
}

# ##Monotone Example
# monotonegrid<-function(n)
# {
# m<-matrix(0,n^2,2)
# for(i in 0:(n-1)){for(j in 1:n){m[(n*i)+j,2]<-i}}
# for(i in 0:(n-1)){for(j in 1:n){m[(n*i)+j,1]<-j-1}}
# return(m)
# }

# monotone_example<-function(m)
# {
#   u<-as.matrix(rowSums(m))
#   return(as.matrix(cbind(m,u)))
# }

# m<-monotonegrid(10)
# mfull<-monotone_example(m)


