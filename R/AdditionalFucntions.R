#Additional functions

#'Interpolates function under shape and/or monotocity constraints
#'
#'Interpolates/extrapolates a quasiconvex function with/without additional
#'monotonicity constraints (already defined on some data points), to a non-data point.
#'This function takes a matrix of data points, a vector of the corresponding functional
#'values at the data points, a non-data point, and shape and monotonicity specifications
#'as inputs, and returns the interpolated/extrapolated value of the function (subject to
#'the shape and monotonicity constraints specified) at the non-data point.
#'
#' @param X              An n by d matrix of data points.
#' @param cplex.fhat     An n by 1 vector of functional values evaluated at the rows of X.
#' @param nondatapoint   A d by 1 non-data point where the function needs to be interpolated/extrapolated.
#' @param Shape          A categorical variable indicating the shape of the function.
#'                       The user must input one of the following three types:
#'                       "QuasiConvex": means quasiconvex function,
#'                       "none": means no shape restriction is imposed.
#' @param Monotone       A categorical variable indicating the monotonicity pattern of the function.
#'                       The user must input one of the following three types:
#'                       "non.inc": means nonincreasing function,
#'                       "non.dec": means nondecreasing function,
#'                       "none": means no monotonicity restriction is imposed.
#' @return               The interpolated/extrapolated value of the function (subject to
#'                       the shape and monotonicity constraints specified) at the non-data point.
#' @export
interpol.func<-function(X,cplex.fhat,nondatapoint, Shape = c("QuasiConvex","none"), Monotone = c("non.inc", "non.dec","none"))
{
  if(Monotone=="non.inc")
  {
    chcheck<-function(A,x)
    {
      return(chullupcheck(A,x))
    }
  }
  else if(Monotone=="non.dec")
  {
    chcheck<-function(A,x)
    {
    return(chulldowncheck(A,x))
    }
  }
  else if (Monotone=="none")
  {
    chcheck<-function(A,x)
    {
    return(chullcheck(A,x))
    }
  }
  if (Shape == "QuasiConvex"){
    X<-as.matrix(X)
    z<-as.matrix(cplex.fhat)
    mat<-cbind(X,z)
    mat<-mat[order(mat[,ncol(mat)]),]
    ####Use Binary Search to arrive at the minimal set####################
    out<-mat[nrow(mat),ncol(mat)]
    Y<-mat
    counter<-1
    if(nrow(Y)>1){
      counter<-chcheck(Y[,1:(ncol(Y)-1)],nondatapoint)
    }
    else{
      counter<-chcheck(t(as.matrix(Y[,1:(ncol(Y)-1)])),nondatapoint)
    }
    #print(counter)
    #############Now will begin binary search##################
    l<-0
    u<-0
    if(counter)
    {
      l<-1
      u<-nrow(Y)
    }
    while(u-l)
    {
      mid<-floor((l+u)/2)
      if(mid>1)
      {
        counter<-chcheck(Y[1:mid,1:(ncol(Y)-1)],nondatapoint)
      }
      else
      {
        counter<-chcheck(t(as.matrix(Y[1,1:(ncol(Y)-1)])),nondatapoint)
      }
      if(counter)
      {
        u<-mid
      }
      else
      {
        l<-mid+1
      }
      out<-Y[u,ncol(Y)]
    }
    return(out)
  }else if (Shape == "none"){
    if(Monotone=="non.dec"){
      X<-as.matrix(X)
      z<-as.matrix(cplex.fhat)
      mat<-cbind(X,z)
      mat<-mat[order(mat[,ncol(mat)]),]
      matsmall<-mat[,1:(ncol(mat)-1)]
      prelog<-apply(t(t(matsmall)<=nondatapoint),1,prod)
      if(sum(prelog)){
        logicvec<-which(prelog==1)
        return(mat[logicvec[length(logicvec)],ncol(mat)])
      }
      else{
        return(mat[1,ncol(mat)])
      }
    }else if(Monotone=="non.inc"){
      Mon<-"non.dec"
      out<-interpol.func(X,-cplex.fhat,nondatapoint,Shape,Mon)
      return(-out)
    }else {
      return(0)
    }
  }
}




# ###################

#  else
#   {
#     if(Monotone=="non.dec")
#     {
#       X<-as.matrix(X)
#       z<-as.matrix(cplex.fhat)
#       mat<-cbind(X,z)
#       mat<-mat[order(mat[,ncol(mat)]),]
#       matsmall<-mat[,1:(ncol(mat)-1)]
#       prelog<-apply(t(t(matsmall)<=nondatapoint),1,prod)
#       if(sum(prelog))
#       {
#         logicvec<-which(prelog==1)
#         return(mat[logicvec[length(logicvec)],ncol(mat)])
#       }
#       else
#       {
#         return(mat[1,ncol(mat)])
#       }
#     }
#     else if(Monotone=="non.inc")
#     {
#       Mon<-"non.dec"
#       out<-interpol.func(X,-cplex.fhat,nondatapoint,Shape,Mon)
#       return(-out)
#     }
#     else
#     {
#       return(0)
#     }
#   }
# }


########################



#' Creates matrices, vectors and variable bounds corresponding to the objective and the constraints
#' of the mixed integer quadratic optimization.
#'
#' This function takes the matrix of regressors, the response
#' variables, the monotonicity specification and some other bounds and tuning parameters as inputs, and returns
#' a list of matrices, vectors and variable bounds corresponding to the objective function and the constraints
#' of the mixed integer quadratic optimization. The output of this function can be used by CPLEX/GUROBI for
#' implementing the mixed integer optimization.
#'
#' @param X              An n by d matrix of data points.
#' @param y              An n by 1 vector of responses.
#' @param Max.b          Bound on the absolute values of the variables.
#' @param MM             The 'M' parameter in the big-M method applied by us.
#' @param ep             A small positive quantity to convert open constraints to close constraints.
#' @param Monotone       A categorical variable indicating the monotonicity pattern of the function.
#'                       The user must input one of the following three types:
#'                       "non.inc": means nonincreasing function,
#'                       "non.dec": means nondecreasing function,
#'                       "none": means no monotonicity restriction is imposed.
#' @return               List of matrices, vectors and variable bounds corresponding to the objective function
#'                       and the constraints of the mixed integer quadratic optimization. The output of this
#'                       function can be used by CPLEX/GUROBI for implementing the mixed integer optimization.
#' @importFrom Matrix Matrix
#' @export
QuasiConv.control <- function(X,y, Max.b , MM , ep, Monotone = c("non.inc", "non.dec","none")) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.vector(y)) y <- as.vector(y)
  nn <- nrow(X)
  dd <- ncol(X)
  cvec<- c( -2*y, rep(0,  nn*dd+ nn^2))
  Qmat <- Matrix(0, nrow= nn+ nn*dd+ nn^2, ncol =nn+ nn*dd+ nn^2, sparse= TRUE)
  Qmat[1: nn, 1:nn] <- 2* diag(nn)
  Qmat <- methods::as(Qmat, "dgTMatrix") ### Qmat is now sparse
  Amat <- Matrix(0,(2*nn*(nn-1)),(nn+(nn*dd)+nn^2), sparse =TRUE)
  negid<- -diag(nn-1)
  jdmat<-cbind(matrix(1,nn-1,1),negid)
  for(i in 2:(nn-1))
  {
    jd<-cbind(negid[,1:(i-1)],matrix(1,nn-1,1),negid[,i:(nn-1)])
    jdmat<-rbind(jdmat,jd)
  }
  jd<-cbind(negid,matrix(1,nn-1,1))
  jdmat<-rbind(jdmat,jd)
  Amat[1:(nn*(nn-1)),1:nn]<-jdmat
  for(alp in 1:nn)
  {
    for(bet in setdiff(c(1:nn),alp))
    {
      if(bet<alp)
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
      else
      {
        Amat[(nn*(nn-1))+((nn-1)*(alp-1))+bet-1,(nn+(alp-1)*dd+1):(nn+alp*dd)]<-X[alp,]-X[bet,]
      }
    }
  }
  Amat[1:(nn-1),(nn+nn*dd+1):(nn+nn*dd+nn)]<-cbind(matrix(0,nn-1,1),MM*negid)
  for(i in 2:(nn-1))
  {
    Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<-cbind(MM*negid[,1:(i-1)],matrix(0,nn-1,1),MM*negid[,i:(nn-1)])
  }
  i<-nn
  Amat[((i-1)*(nn-1)+1):(i*(nn-1)),((nn*dd)+(i*nn)+1):(nn*dd + i*nn + nn)]<- cbind(MM*negid,matrix(0,nn-1,1))
  Amat[(nn*(nn-1)+1):(2*nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]<- - Amat[1:(nn*(nn-1)),(nn+nn*dd+1):(nn+nn*dd+nn^2)]
  Amat <- methods::as(Amat, "dgTMatrix") ### Qmat is now sparse

  ########################
  bvec <-c(rep(0, nn*(nn-1)), rep(MM-ep, nn*(nn-1)))
  if(Monotone=="none")
  {
  lb <- c(rep(-Max.b, nn), rep(-Max.b, nn*dd), rep(0, nn^2))
  ub <-  c(rep(Max.b, nn+ nn*dd), rep(1, nn^2))
  }
  else if(Monotone=="non.inc")
  {
    lb <- c(rep(-Max.b, nn), rep(0, nn*dd), rep(0, nn^2))
    ub <-  c(rep(Max.b, nn+ nn*dd), rep(1, nn^2))
  }
  else if(Monotone=="non.dec")
  {
    lb <- c(rep(-Max.b, nn+ nn*dd), rep(0, nn^2))
    ub <-  c(rep(Max.b, nn), rep(0, nn*dd), rep(1, nn^2))
  }
  else
  {
    stop("Monotonicity constraint not specified correctly!")
  }
  vtype <- c(rep("C", nn+ nn*dd), rep("B", nn^2))
  ###########################

  list(cvec = cvec,
       Amat = Amat,
       lb = lb,
       ub = ub,
       bvec = bvec,
       vtype = vtype,
       Qmat = Qmat,
       dd = dd,
       varlength = length(cvec) )
}




