#FSFA

#libraries
library(Rfast)
source("solveV.R")
library(doParallel)
library(purrr)


# Helpers -----------------------------------------------------------------
cfun = function(...) NULL

tr <- function(A=NULL) { sum(  diag(A)  ) }

pred = function(gpars=NULL, newdata=NULL, probs=FALSE)
{
  G = length(gpars$pi)
  dens = sapply(1:G, function(g){ logden( newdata, gpars$dpar[[g]] ) })
  
  vals = t( exp( log(gpars$pi) + t(dens) ) )
  vals = vals/rowsums(vals)
  
  if(!probs){  vals = rowMaxs(vals)  }
  
  
  return(vals)
}

norm_eig2 = function(del=NULL, eta=NULL, p=NULL)
{
  q = length(del)
  logdet = sum(  c(log(del),(p-q)*log(eta))  )
  logvals = c(  log(del)-logdet/p, log(eta)-logdet/p  )
  exp(logvals)
}

wbinv = function(B=NULL, D=NULL)
{
  dinv = D^(-1)
  dB = dinv*B
  
  val = diag(dinv) - dB%*% solve( diag(dim(B)[2]) + crossprod(dB, B) ) %*%t(dB)
  
  return(val)
}

dirdraw = function(alpha=NULL)
{
  y = rgamma(length(alpha), shape=alpha)
  val = y/sum(y)
  return(val)
}

gparbar = function(gpars=NULL, parbar=NULL, denom=NULL)
{
  G = length(gpars$dpar)
  for( i in seq_len(G))
  {
    parbar$dpar[[i]]$mu = parbar$dpar[[i]]$mu + 
      (gpars$dpar[[i]]$mu - parbar$dpar[[i]]$mu)/(denom+1)
    
    parbar$dpar[[i]]$lam1 = parbar$dpar[[i]]$lam1 + 
      (gpars$dpar[[i]]$lam1-parbar$dpar[[i]]$lam1)/(denom+1)
    
    parbar$dpar[[i]]$xi1 = parbar$dpar[[i]]$xi1 + 
      (gpars$dpar[[i]]$xi1-parbar$dpar[[i]]$xi1)/(denom+1)
    
    parbar$dpar[[i]]$lam2 = parbar$dpar[[i]]$lam2 + 
      (gpars$dpar[[i]]$lam2-parbar$dpar[[i]]$lam2)/(denom+1)
    
    parbar$dpar[[i]]$omg2 = parbar$dpar[[i]]$omg2 + 
      (gpars$dpar[[i]]$omg2-parbar$dpar[[i]]$omg2)/(denom+1)
    
    parbar$dpar[[i]]$eta2 = parbar$dpar[[i]]$eta2 + 
      (gpars$dpar[[i]]$eta2-parbar$dpar[[i]]$eta2)/(denom+1)
  }
  parbar$pi = parbar$pi + (gpars$pi - gpars$pi)/(denom+1)
  
  return(parbar)
}

relDiff = function(x=NULL)
{
  abs(x[-length(x)] - x[-1])/x[-length(x)]
}


# Density  ----------------------------------------------------------------

logden <- function(datax=NULL, par=NULL) 
{
  ## datax is an array with dimension (p1 x p2 x n) 
  logc = -prod(par$pb)*log(2*pi)/2
  chols1 = chol(par$lam1%*%t(par$lam1) + diag(par$xi1) )
  
  logDet1 = -( par$pb[2]/2 ) * ( 2 * sum(log(diag(chols1))) ) #verbose
  logDet2 = -( par$pb[1]/2 ) * 
  ( sum(log(par$omg2)) + (par$pb[2]-par$qd[2])*log(par$eta2) ) 
  
  trterm = -0.5*trv2( dt=datax, par=par )
  
  lval <- logc + logDet1 + logDet2 + trterm 
  lval
}

logdets <- function(datax=NULL, par=NULL) 
{
  ## datax is an array with dimension (p1 x p2 x n) 
  logc = -prod(par$pb)*log(2*pi)/2
  chols1 = determinant( 
    par$lam1%*%t(par$lam1) + diag(par$xi1), 
    logarithm = TRUE 
  )$modulus
  
  logDet1 = -(par$pb[2]/2) * chols1
  logDet2 = -(par$pb[1]/2) * 
  ( sum(log(par$omg2)) + (par$pb[2]-par$qd[2])*log(par$eta2) )  #not necessary
  
  
  lval <- logc + logDet1 + logDet2
  lval
}

trwz = function(dt=NULL, par=NULL, zs=NULL)
{
  nobs = dim(dt)[3]
  # r = (x - m)
  r = sweep(dt, c(1,2), par$mu, "-",  check.margin = FALSE) #correct
  r = sweep(r, 3, sqrt(zs), "*",  check.margin = FALSE) #correct
  
  di2 = diag( 1/par$omg2 - 1/par$eta2, nrow=length(par$omg2) )
  
  invS2 = par$lam2%*% di2 %*%t(par$lam2) + 1/par$eta2*diag(nrow=par$pb[2]) #correct
  invS1 = wbinv(par$lam1, par$xi1) #correct
  
  W2  = matrix( 
    rowmeans(apply( r, 3, function(z, A=NULL){ z %*% tcrossprod(A,z) }, A=invS2)),
    nrow=par$pb[1]
  ) #correct
  
  val = nobs*sum(invS1 * W2) #correct
  
  
  return( val )
  
}

trv2 = function(dt=NULL, par=NULL)
{
  nr = dim(dt)[3]
  r = sweep(dt, c(1,2), par$mu, "-",  check.margin = FALSE) #correct
  
  #invS1 = wbinv(par$lam1, par$xi1) #correct
  invS1 = solve(par$lam1 %*% t(par$lam1) + diag(par$xi1) )
  
  di2   = diag(1/par$omg2 - 1/par$eta2, nrow=length(par$omg2))
  invS2 = par$lam2%*% di2 %*%t(par$lam2) +  1/par$eta2*diag(nrow=par$pb[2]) #correct
  
  trtm = apply( r, 3, function(z, A=NULL, B=NULL){ sum((z%*%tcrossprod(A,z)) * B) }, A=invS2, B=invS1) #correct
  
  return(trtm)
}

mxll2 = function(dat=NULL, gpar=NULL)
{
  G = length(gpar$pi)
  ldens = sapply(1:G, FUN=function(g){ logden(datax=dat, gpar$dpar[[g]]) })
  # n times g matrix
  # rows are observations, cols are groups
  
  val = ldens + rep( log(gpar$pi), rep(nrow(ldens),G) ) #correct
  
  rmx = rowMaxs( val, value = TRUE )
  
  valr = log( rowsums( exp(val - rmx) ) ) + rmx
  valr = sum(valr) 
  
  return( valr )
}

dmvn = function(datx=NULL, par=NULL) #sigma needs to be fixed. not pd
{
  nr = dim(datx)[3]
  
  s1 = par$lam1%*%t(par$lam1) + diag(par$xi1, nrow=length(par$xi1)) 
  s2 = par$lam2 %*% diag(par$omg2, nrow=length(par$omg2)) %*% t(par$lam2) + 
    par$eta2 * (diag(nrow=par$pb[2]) - par$lam2%*%t(par$lam2))
  
  s = s2 %x% s1
  m = as.numeric(par$mu)
  
  val = apply(datx, 3, function(z){ dmvnorm(as.numeric(z), m, s, logged=TRUE) })
  
}

mxllfast = function(dat=NULL, gpar=NULL)
{
  G = length(gpar$pi)
  ldens = sapply(1:G, FUN=function(g){ dmvn(datx=dat, gpar$dpar[[g]]) })
  val = ldens + rep( log(gpar$pi), rep(nrow(ldens),G) ) #correct
  vrmx = rowMaxs(val, value = TRUE)
  valr = sum(log( rowsums( exp(val - vrmx) ) ) + vrmx)
  return(valr)
}

tr1 = function(x=NULL, par=NULL)
{
  invs1 = wbinv(par$lam1, par$xi1) #correct
  invs1 = (invs1 + t(invs1))/2
  di2   = diag(1/par$omg2 - 1/par$eta2, nrow=length(par$omg2))
  invs2 = par$lam2%*% di2 %*%t(par$lam2) +  1/par$eta2*diag(nrow=par$pb[2]) #correct
  
  #Got lazy. make this regular multiplication later
  tr = -0.5*sum(diag((invs1 %*% x %*% invs2 %*% t(x))))
  
  
}

emll = function(dat=NULL, gpar=NULL)
{
  mapz = rowMaxs(gpar$Z)
  
  nobs = dim(dat)[3]
  G = length(gpar$pi)
  
  #logdets each group
  lds = sapply(1:G, function(i){ logdets(datax=dat, par=gpar$dpar[[i]]) }) #assumes logdets is correct
  lgp = log(gpar$pi)
  lgs = lgp + lds
  
  val= numeric()
  for(i in seq_len(nobs))
  {
    val[i] = lgs[ mapz[i] ] + tr1( dat[,,i], gpar$dpar[[ mapz[i] ]] )
  }
  
  return( sum(val) )
}

ezll = function(dat=NULL, gpar=NULL)
{
  #mxll = sumsum zig*logpi + zig*logdets + zig*trace
  #     = sum_g ng*logpig + ng*logdets -0.5*trace{ W2 %*% S1 }
  nr = dim(dat)[3]
  G = length(gpar$pi)
  
  #logdets each group
  lds = sapply(1:G, function(i){ logdets(datax=dat, par=gpar$dpar[[i]]) }) #assumes logdets is correct
  lds = sum((log(gpar$pi) + lds) * gpar$pi * nr) #correct

  trs = sapply(1:G, function(i)
  {
    -0.5*trwz(dt=dat, par=gpar$dpar[[i]], zs=gpar$Z[,i])
  }
  )/G
  
  val = sum(lds) + G*sum(trs)
  return( val )
}

ent = function(dat=NULL, gpar=NULL, map=FALSE)
{
  #  map option gives Biernacki ICL, otherwise
  #  it is calculated using Zs at MLE 
  
  G = length(gpar$pi)
  ldens = sapply(1:G, FUN=function(g){ logden(datax=dat, gpar$dpar[[g]]) }) #correct if logden is
  # n times g matrix
  # rows are observations, cols are groups
  
  val = ldens + rep( log(gpar$pi), rep(nrow(ldens),G) ) #correct
  v2 = val - rowMaxs(val,value=TRUE)
  v2s = rowsums(exp(v2))
  v3 = v2 - log(v2s)
  Zs = exp(v3)
  
  #Biernacki ICL
  if( map )
  {
    ind = rowMaxs(Zs)
    ind = cbind(seq_along(ind), ind)
    
    Zs = matrix(0, nrow=nrow(Zs), ncol=G)
    Zs[ind] = 1
  }
  
  val = sum(Zs*v3) #This is too many computations in the case of map.
  
  return(val)
} #interpretation: uncertainty (0 implies low uncertainty in the model)


# Initialization ----------------------------------------------------------

#Gen one param set
# d = (row, col) = (p, b)
# q = (latent gaussian row dim, intrinsic column dim) = (q, dg)

#Random params
rgpar <- function(
  g=2, ddim=c(4,5), idim=NULL, nobs=NULL, hW=NULL, qs=NULL, ...
) {
  #q converted to a 2xG matrix
  #stop()
  useqs=TRUE
  
  if( is.null(qs)[1] ) #set idim if qs are not chosen
  {
    useqs = FALSE
    if( is.null(idim)[1] )  
    {
      idim = matrix(rep(c(2,2), times=g), nrow=g, byrow=TRUE)
    }
    if( length(idim) == length(ddim) ) 
    {  
      idim = matrix(rep(idim, times=g), nrow=g, byrow=TRUE)  
    }  
    if(!(dim(idim)[1] == g) ) stop("qd's must be in the form of a 2xg matrix")
  }
  
  if( is.null(hW) ){ hW = diag(pb[2]) }
  
  val=list()
  
  if( useqs )
  {
    val$dpar = lapply(1:g, FUN=function(gp)
    { 
      rpar(pb=ddim, auto=TRUE, qi=qs[gp], ...) 
    })
  }else
  {
    val$dpar = lapply(1:g, FUN=function(gp)
    { 
      rpar(pb=ddim, qd=as.numeric(idim[gp,]), auto=FALSE, ...) 
    })
  }
  #Diriclet parameter and draws
  a = matrix( runif(nobs*g, min=0.1, max=5), nrow=nobs, byrow=TRUE )
  Zm = t( apply(a, 1, dirdraw) )
  
  val$hW  = hW
  val$Z   = Zm
  val$pi  = colmeans(Zm)
  val$lls = -Inf
  
  return(val)
}	

# when auto is on, d is chosen large, so that it can be automated
# later in the cyc1 steps
rpar <- function(pb=c(4,5), qd=c(2,2), auto=TRUE, qi=NULL)
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  if(auto)
  { 
    if(is.null(qi)[1]){ qi = round(runif( 1, 1, (pb[1]-1) )) }
    qd = c(qi, pb[2]-1) 
  }
  
  val = list()
  s1  = cov(matrix(  rnorm( prod(2*pb[1]^2), sd=1   ), 2*pb[1], pb[1]  ))
  s2  = cov(matrix(  rnorm( prod(2*pb[2]^2), sd=0.1 ), 2*pb[2], pb[2]  )) 
  
  tempmu = matrix( rnorm(pb[1]*qd[2], mean=10, sd=2), pb[1], qd[2] ) #pxd
  
  temp1 = eigen.sym(s1, k=qd[1])
  val$lam1  = temp1$vectors * rep( sqrt(temp1$values), 
                                   rep(nrow(temp1$vectors), 
                                       length(temp1$values))
  ) #works
  val$xi1 =  diag(s1 - tcrossprod(val$lam1,val$lam1))
  
  temp2 = eigen(s2)
  val$lam2  = temp2$vectors[,1:qd[2],drop=FALSE]
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=temp2$values[1:qd[2]], eta=temp2$values[qd[2]+1], p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$mu = tcrossprod(tempmu, val$lam2)
  
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}


#Random Zs
genZ = function(datx=NULL, g=NULL)
{
  nobs = dim(datx)[3]
  
  #Diriclet parameter and draws
  a = matrix( runif(nobs*g, min=0.01, max=10), nrow=nobs, byrow=TRUE)
  Zm = t(apply(a, 1, dirdraw ))
  return(Zm)
}

# When auto is on, automatically initializes d based on the
# MLE for the row covariance matrix. q must be set. 
rzpar <- function(
  datx=NULL, zs=NULL, pb=NULL, qd=c(2,2), iter=10, auto=TRUE, thres=NULL, 
  qinit=NULL
) {
  #
  # pb = (row, col) = (p, b)
  # qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  #
  
  #
  # auto = TRUE: choose d automatically 
  #              (no need to set thresh, qi should be set)
  #
  
  val = list()
  
  ng=sum(zs)
  val$mu = apply( sweep(datx, 3, zs/ng, "*"), 
                  2, 
                  rowsums) #correct
  
  # r = (x - m)
  r = sweep(datx, c(1,2), val$mu, "-",  check.margin = FALSE) #correct
  r = sweep(r, 3, sqrt(zs), "*",  check.margin = FALSE)
  
  
  #ML estimation of covariances
  if(!is.null(iter))
  {
    s1  = diag(nrow=pb[1])
    s2  = diag(nrow=pb[2]) 
    
    for(i in seq_len(iter))
    {
      invs1 = solve(s1)
      
      s2 = matrix(
        rowmeans( 
          apply( r, 3, function(z, A=NULL){ t(z) %*% A %*% z }, A=invs1)
        ), 
        nrow=pb[2] 
      )/pb[1]
      s2 = (s2 + t(s2))/2 #+ diag(10, nrow=pb[2])

      invs2 = solve(s2)
      
      s1 = matrix(
        rowmeans( 
          apply( r, 3, function(z, A=NULL){ z %*% A %*% t(z) }, A=invs2)
        ), 
        nrow=pb[1] 
      )/pb[2]
      s1 = (s1 + t(s1))/2
      
    }
  }
  else
  {
    s1  = cov(matrix(  rnorm( prod(2*pb[1]^2), sd=1 ), 2*pb[1], pb[1]  ))
    s2  = cov(matrix(  rnorm( prod(2*pb[2]^2), sd=1 ), 2*pb[2], pb[2]  )) 
    
    invs1 = solve(s1)
    invs2 = solve(s2)
    
    s1 = matrix(
      rowmeans( 
        apply( r, 3, function(z, A=NULL){ z %*% A %*% t(z) }, A=invs2 )
      ), 
      nrow=pb[1] 
    )/pb[2]
    
    s2 = matrix(
      rowmeans( 
        apply( r, 3, function(z, A=NULL){ t(z) %*% A %*% z }, A=diag(nrow=pb[1]))
      ), 
      nrow=pb[2] 
    )/pb[1]
  }
  
  if(auto) #automatically choose d (not q)
  {
    if(is.null(thres)){ thres=0.68 }
    qd = c(0,0)
    eigv2 = eigen(s2)$values/sum(eigen(s2)$values)
    
    #rd = which.max(relDiff(eigv2)/eigv2[-1])
    
    #if(rd >= ((pb[2]-1)/2) ){ 
      props = cumsum(eigv2)
      rd = ifelse( which(props >= thres)[1] < pb[2], #statement
                   which(props >= thres)[1],         #TRUE
                   round((pb[2]-1)/2))               #FALSE
    #}
    
    qd[2] = rd
    qd[1] = ifelse(!is.null(qinit)[1], qinit, min(round(pb[1]/4), qd[2]) )
  }
  
  
  temp1 = eigen.sym(s1, k=qd[1])
  val$lam1  = temp1$vectors * rep( sqrt(temp1$values), 
                                   rep(nrow(temp1$vectors),length(temp1$values))
  ) #correct
  val$xi1 =  diag(s1 - tcrossprod(val$lam1,val$lam1))
  
  temp2 = eigen(s2)
  val$lam2  = temp2$vectors[,1:qd[2], drop=FALSE]
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=temp2$values[1:qd[2]], eta=temp2$values[qd[2]+1], p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  #val$mu = tcrossprod(tempmu, val$lam2)
  
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

rgzpar <- function(
  dt=NULL, g=2, ddim=c(4,5), idim=c(2,2), hW=NULL, qs=NULL, ...
) {
  #q converted to a 2xG matrix
  #stop()
  nobs=dim(dt)[3]
  pb = dim(dt)[1:2]
  useqs=TRUE
  
  if( is.null(qs)[1] ) #set idim if qs are not chosen
  {
    useqs = FALSE
    if( is.null(idim)[1] )  
    {
      idim = matrix(rep(c(2,2), times=g), nrow=g, byrow=TRUE)
    }
    if( length(idim) == length(ddim) ) 
    {  
      idim = matrix(rep(idim, times=g), nrow=g, byrow=TRUE)  
    }  
    if(!(dim(idim)[1] == g) ) stop("qd's must be in the form of a 2xg matrix")
  }
  
  if(is.null(hW)){ hW = diag(pb[2]) }
  
  
  nobs=dim(dt)[3]
  pb = dim(dt)[1:2]
  
  
  val=list()
  
  #Diriclet parameter and draws
  a = matrix( runif(nobs*g, min=0.01, max=10), nrow=nobs, byrow=TRUE)
  Zm = t(apply(a, 1, dirdraw ))
  
  #Generate parameters
  if(useqs)
  {
    val$dpar = lapply(1:g, FUN=function(gp)
    { 
      rzpar(datx=dt, zs=Zm[,gp], pb=ddim, qinit=qs[gp], ...) 
    })
  }
  else
  {
    val$dpar = lapply(1:g, FUN=function(gp)
    { 
      rzpar(datx=dt, zs=Zm[,gp], pb=ddim, qd=as.numeric(idim[gp,]), ...) 
    })
  }
  val$hW  = hW
  val$Z   = Zm
  val$pi  = colmeans(Zm)
  val$lls = -Inf
  
  return(val)
}	



#SEM init
typeSEM = function(z=TRUE) #specifies SEM to use two or one zigup
{
  if(z)
  { #two z updates
    f = function(datx=NULL, gpar=NULL, cutoff=NULL)
    {
      #if(length(cutoff)==1){ cutoff = c(cutoff, cutoff) }
      #cycle1
      gpar$Z = zig.up(datx, gpar)
      gpar$Z = zdraw(gpar$Z)
      gpar$dpar = lapply(seq_len(ncol(gpar$Z)), FUN=function(g)
      { 
        cyc1(datx, par=gpar$dpar[[g]], zs=gpar$Z[,g], thres=cutoff, hW=gpar$hW) 
      })
      
      #cycle2
      #gpar$Z = zig.up(datx, gpar)
      #gpar$Z = zdraw(gpar$Z)
      gpar$dpar = lapply(seq_len(ncol(gpar$Z)),FUN=function(g)
      { 
        cyc22(datx, gpar$dpar[[g]], gpar$Z[,g]) 
      })
      
      gpar$pi = colmeans(gpar$Z)
      gpar$lls = c(gpar$lls, mxll2(datx, gpar))
      
      return(gpar)
    }
  }
  else
  { #one z update
    f = function(datx=NULL, gpar=NULL, cutoff=NULL)
    {
      gpar$Z = zig.up(datx, gpar)
      gpar$Z = zdraw(gpar$Z)
      gpar = Mstep(datx, gpar, thresh=cutoff)
      gpar$pi = colmeans(gpar$Z)
      gpar$lls = c(gpar$lls, mxll2(datx, gpar))
    }
  }
  
  return(f)
}

SEMz = function(
  data=NULL, gpar0=NULL, niter=NULL, burnin=NULL, fullEM=TRUE, ...
) {
  #
  # SEMmean currently disabled, because an iteratively updated d potentially
  # makes dimensions of cyc2 parameters differ across iterations
  #
  itSEM = typeSEM(fullEM)
  try.SEM = possibly(itSEM, otherwise=NA)
  burnits = 0
  
  val = list("max"=NA, 
             "BIC"=NA,
             "fail"=FALSE,
             "burnlls"=NA)
  
  if( !is.null(burnin) )
  {
    
    for( i in seq_len(burnin) )
    {
      t=0
      gpna=NA
      
      while( is.na(gpna)[1] && t<10 )
      {
        gpna = try.SEM(data, gpar0, ...) #... == cutoff basically
        t=t+1
      }
      
      if( !is.na(gpna)[1] )
      { 
        gpar0 = gpna 
        burnits = burnits + 1
      }
      else
      { #Report the failure
        val$fail=TRUE
        val$message = "SEM stuck during burn in. 
                       Try different Z's. 
                       Returning last parameter update"
        break;
      }
    }
    
  }
  if( val$fail )
  {
    #parbar = gpar0
    val$max=gpar0
    val$BIC=BICfc(data, gpar0)
    val$burnits=burnits
  }
  else
  {
    burnlls   = gpar0$lls
    gpar0$lls = gpar0$lls[length(gpar0$lls)]
    val$max   = gpar0
    #parbar   = gpar0
    bll       = gpar0$lls
    
    for( i in seq_len(niter) )
    {
      t=0 
      gpna=NA
      
      while( is.na(gpna)[1] && t<10 ) #try to update the parameters with SEM
      {
        gpna = try.SEM(data, gpar0, ...)
        t=t+1
      }
      if( !is.na(gpna)[1] ) #update best param values
      {
        gpar0 = gpna
        if( gpar0$lls[i+1] > bll ){ val$max = gpar0; bll = gpar0$lls[i+1] }
        #parbar = gparbar(gpar0, parbar, denom=i)
        #parbar$lls = mxll2(data, parbar)
      }
      else
      {#If iteration doesnt succeed, report the failure and quit
        val$fail=TRUE
        val$message = "Burnin successful, but SEM stuck during iteration. 
                       Try different Z's. 
                       Returning last parameter update"
        break;
      }
      
    }
    
    #parbar$Z = zig.up(data, parbar)
    val$BIC = BICfc(data, val$max)
    val$fulls = c(burnlls, gpar0$lls)
  }
  
  return(val)
  
} #Make this faster in the future if possible



# Update Functions --------------------------------------------------------

zig.up = function(dat=NULL, gpar=NULL)
{
  # n times g matrix
  # rows are observations, cols are groups
  G = length(gpar$pi)
  ldens = sapply(
    seq_len(G), 
    FUN=function(g){ logden(datax=dat, gpar$dpar[[g]]) }
  ) 
  
  val = ldens + rep( log(gpar$pi), rep(nrow(ldens),G) ) #correct
  val = exp(val -  rowMaxs(val, value=TRUE)) 
  val = val / rowsums(val)
  
  return(val)
}

cyc1 <- function(datax=NULL, par=NULL, zs=NULL, auto=TRUE, thres=NULL, hW=NULL) 
{
  # datax is an array with dimension (p x b x n) 
  # zs the column of Z corresponding to the current group
  #if (is.null(zs)) zs = rep(1, dim(datax)[3])
  
  ng = sum(zs)

  #mu
  par$mu = apply(  sweep( datax, 3, as.numeric(zs)/ng, "*" ), 2, rowsums  )
  
  # r = (x - m)
  r = sweep( datax, c(1,2), par$mu, "-",  check.margin = FALSE )
  r = sweep( r, 3, sqrt(zs), "*",  check.margin = FALSE )
  
  #update column covariance parameters
  #invS1 = wbinv(par$lam1, par$xi1) #correct
  invS1 = solve( par$lam1%*%t(par$lam1) + diag(par$xi1) ) #correct
  invS1 = ( invS1 + t(invS1) )/2
  
  W1  = matrix(
    rowmeans( 
      apply( r, 3, function(z, A=NULL){ t(z) %*% A %*% z }, A=invS1 ) 
    ), 
    nrow=par$pb[2] 
  )
  W1 = ( W1 + t(W1) )/2
  W1 = hW %*% W1 %*% t(hW)
  
  if(auto)
  {
    if( is.null(thres) ){ thres=0.68  } 
    
    eigv = eigen(W1)$values/sum(eigen(W1)$values)
    
    props = cumsum(eigv)
    rd = ifelse( 
      which(props >= thres)[1] < par$pb[2], 
      which(props >= thres)[1], 
      par$qd[2]
    )

    par$qd[2] = rd
  }
  
  par$lam2 = eigen.sym( W1, par$qd[2] )$vectors
  
  com2 = colsums(  t( crossprod(par$lam2,W1) )*par$lam2  ) #correct
  cae2 = ( sum(diag(W1)) - sum(com2) )/( par$pb[2]-par$qd[2] ) #correct 
  c2g = ( prod(com2)*(cae2)^(par$pb[2]-par$qd[2]) )^(-1/ par$pb[2] ) 
  
  #eigenvalues
  par$omg2 = c2g * com2  
  par$eta2 = c2g * cae2 
  
  return(par)
}

cyc2 <- function( datax=NULL, par=NULL, zs=NULL ) #deprecated
{
  # datax is an array with dimension (p x b x n) 
  # zs the group column of Z matrix
  # if (is.null(zs)) zs = rep(1, dim(datax)[3])
  #stop()
  
  pig = mean(zs)
  
  # r = (x - m)
  r = sweep(datax, c(1,2), par$mu, "-",  check.margin = FALSE)
  r = sweep(r, 3, sqrt(zs), "*",  check.margin = FALSE)
  
  #update row covariance parameters
  di2   = diag(1/par$omg2 - 1/par$eta2, nrow=length(par$omg2))
  invS2 = par$lam2%*% di2 %*%t(par$lam2) +  1/par$eta2*diag(nrow=par$pb[2]) #correct
  
  invS1 = wbinv(par$lam1, par$xi1) #correct
  invS1 = (invS1 + t(invS1))/2
  
  W2  = matrix(
    rowmeans( 
      apply( r, 
             3, 
             function(z, A=NULL){ z%*%tcrossprod(A,z) }, 
             A=invS2 
      )
    ), 
    nrow=par$pb[1]
  ) #correct 

  W2 = (W2 + t(W2))/2
  
  #intermediate parameters
  beta1 = t(par$lam1) %*% invS1 #correct
  vee   = diag( par$pb[2]*pig, nrow=par$qd[1] ) - (par$pb[2]*pig)*(beta1 %*% par$lam1)  +  
    beta1%*%tcrossprod(W2, beta1) 
  
  invV = solve(vee) #Better way? 
  
  WxB  = tcrossprod(W2, beta1)
  
  #updates
  par$lam1 =  WxB %*% invV
  par$xi1  =  ( 1/(par$pb[2]*pig) ) * diag( W2 - par$lam1%*%beta1%*%W2 )

  return(par)
}

cyc22 <- function( datax=NULL, par=NULL, zs=NULL ) 
{
  # datax is an array with dimension (p x b x n) 
  # zs the the column of Z corresponding to current group
  # if (is.null(zs)) zs = rep(1, dim(datax)[3])
  
  ng = sum(zs)
  
  # r = (x - m)
  r = sweep( datax, c(1,2), par$mu, "-",  check.margin = FALSE )
  r = sweep( r, 3, sqrt(zs/ng), "*",  check.margin = FALSE ) #note ng
  
  #update row covariance parameters
  di2   = diag(  1/par$omg2 - 1/par$eta2, nrow=length(par$omg2)  )
  invS2 = par$lam2%*% di2 %*%t(par$lam2) +  
  1/par$eta2*diag(nrow=par$pb[2]) #correct
  
  invS1 = wbinv(par$lam1, par$xi1) #correct. 
  invS1 = ( invS1 + t(invS1) )/2
  
  
  
  W2  = matrix(
    rowsums( 
      apply( r, 3, function(z, A=NULL){ z %*% tcrossprod(A,z) }, A=invS2 )
    ), 
    nrow=par$pb[1]
  ) #correct 
  
  W2 = ( W2 + t(W2) )/2
  
  #intermediate parameters
  beta1 = t(par$lam1) %*% invS1 #correct
  vee   = diag( par$pb[2], nrow=par$qd[1] ) - 
  (par$pb[2])*(beta1 %*% par$lam1) +
  beta1 %*% W2 %*% t(beta1) 
  
  vee = ( vee + t(vee) )/2
  
  invV = solve(vee) #Better way? 
  
  #updates
  lam1 =  W2 %*% t(beta1) %*% invV #ngs in denominator cancel
  
  xi1  =diag( W2 - lam1%*%beta1%*%W2 )/par$pb[2]
  
  par$lam1 = lam1
  par$xi1 = xi1
  
  return(par)
}

par.up = function(dtx=NULL, gpar=NULL, ...)
{
  #Cycle 1
  gpar$dpar = lapply(1:ncol(gpar$Z), FUN=function(g)
  { 
    cyc1(dtx, gpar$dpar[[g]], gpar$Z[,g], hW=gpar$hW, ...) 
  })
  
  #Cycle 2
  #gpar$Z = zig.up(datx, gpar)
  gpar$dpar = lapply(1:ncol(gpar$Z), FUN=function(g)
  { 
    cyc22(dtx, gpar$dpar[[g]], gpar$Z[,g]) 
  })
  
  gpar$pi = colmeans(gpar$Z)
  
  return(gpar)
}

zdraw = function(Zs=NULL, rerolls=5)
{
  val = apply(Zs, 1, function(r){ rmultinom(1,1, r) })
  
  rerolls = NULL #not implemented yet
  
  return(t(val))
}


#Parsimonious specifications of the updates




# Em function -------------------------------------------------------------

Mstep = function(datx=NULL, gpars=NULL, thresh=NULL)
{
  #cycle1
  gpars$dpar = lapply(1:ncol(gpars$Z), FUN=function(g)
  { 
    cat("group",g)
    cyc1(datx, gpars$dpar[[g]], gpars$Z[,g], thres=thresh[1], hW=gpars$hW) 
  }
  )
  #cycle2
  gpars$dpar = lapply(1:ncol(gpars$Z),FUN=function(g)
  { 
    cyc2(datx, gpars$dpar[[g]], gpars$Z[,g], thres=thresh[2]) 
  }
  )
  
  gpars$lls = c(gpars$lls, mxll2(datx, gpars))
  
  return(gpars)
}

itAECM = function(datx=NULL, gpar=NULL, cutoff=NULL)
{
  #Cycle 1
  gpar$Z = zig.up( datx, gpar )
  gpar$dpar = lapply( 1:ncol(gpar$Z), FUN=function( g )
  { 
    cyc1( datx, gpar$dpar[[g]], gpar$Z[,g], thres=cutoff, hW=gpar$hW ) 
  })
  
  #Cycle 2
  #gpar$Z = zig.up(datx, gpar) (optional)
  gpar$dpar = lapply(1:ncol(gpar$Z), FUN=function( g )
  { 
    cyc22( datx, gpar$dpar[[g]], gpar$Z[,g] ) 
  })
  
  gpar$pi = colmeans( gpar$Z )
  gpar$lls = c( gpar$lls, mxll2(datx, gpar) )
  
  return(gpar)
}

AECMn = function(data=NULL, gpar0=NULL, niter=NULL, ...)
{
  for( i in seq_len(niter) )
  {
    gpar0 = itAECM(data, gpar0, ...)
  }
  
  return(gpar0)
}




# Model Selection ---------------------------------------------------------


BICfc = function(dt=NULL, fgpar=NULL, icl=FALSE, mapz=FALSE)
{
  G = length(fgpar$pi)
  nobs = nrow(fgpar$Z)
  pb = fgpar$dpar[[1]]$pb
  
  qds = sapply(1:G, function(g)
  {
    val = fgpar$dpar[[g]]$qd
    return(val)
  })
  
  k1 = (2*pb[1] + pb[2] +0.5)*sum(qds[2,]) + (pb[1]+1)*G
  k2 = -0.5*sum(qds[2,]^2)
  k3 = 0.5*sum(qds[1,]*(qds[1,]-1))
  
  #Calc entropy
  if( icl )
  { 
    zent=ent(dt, fgpar, map=mapz) 
    bic = -max(fgpar$lls) - zent + 0.5*(k1+k2+k3)*log(nobs)
  }
  else
  {
    bic = -2*max(fgpar$lls) + (k1+k2+k3)*log(nobs)
  }
  
  return(bic)
}

cplx = function(fgpar=NULL)
{
  G = length(fgpar$pi)
  nobs = nrow(fgpar$Z)
  pb = fgpar$dpar[[1]]$pb
  
  qds = sapply(1:G, function(g)
  {
    val = fgpar$dpar[[g]]$qd
    return(val)
  })
  
  k1 = (2*pb[1] + pb[2] +0.5)*sum(qds[2,]) + (pb[1]+1)*G
  k2 = -0.5*sum(qds[2,]^2)
  k3 = 0.5*sum(qds[1,]*(qds[1,]-1))
  
  k1+k2+k3
}


#Deprecated
ICLmfc = function(dat=NULL, fgpar=NULL)
{
  G = length(fgpar$pi)
  nobs = nrow(fgpar$Z)
  pb = fgpar$dpar[[1]]$pb
  
  qds = sapply(1:G, function(g)
  {
    val = fgpar$dpar[[g]]$qd
    return(val)
  })
  
  k1 = (2*pb[1] + pb[2] +0.5)*sum(qds[2,]) + (pb[1]+1)*G
  k2 = -0.5*sum(qds[2,]^2)
  k3 = sum(qds[1,]*(qds[1,]-1))
  
  #Requires calculated expected log likelihood instead
  icl = -2*emll(dat, fgpar) + (k1+k2+k3)*log(nobs)
  
  return(icl)
}

ICLzfc = function(dat=NULL, fgpar=NULL)
{
  G = length(fgpar$pi)
  nobs = nrow(fgpar$Z)
  pb = fgpar$dpar[[1]]$pb
  
  qds = sapply(1:G, function(g)
  {
    val = fgpar$dpar[[g]]$qd
    return(val)
  })
  
  k1 = (2*pb[1] + pb[2] +0.5)*sum(qds[2,]) + (pb[1]+1)*G
  k2 = -0.5*sum(qds[2,]^2)
  k3 = sum(qds[1,]*(qds[1,]-1))
  
  #Requires calculated expected log likelihood instead
  icl = -2*ezll(dat, fgpar) + (k1+k2+k3)*log(nobs)
  
  return(icl)
}

















# Garbagio ----------------------------------------------------------------






## Garbage shute
# #mean update
# mu.up = function(dat=NULL, par=NULL, zs=NULL)
# {
#   ng =sum(zs)
#   # zs the group column of Z matrix
#   mg = apply(
#     sweep(sweep(datax, c(1,2), par$mu, "-"), 
#           3, zs, "*", 
#           check.margin=FALSE
#     ),
#     2, rowsums
#   )
#   # %*%par$lam2 to find mu_gs -- rows are the mu_gs
#   
#   return(mg/ng)
# }
# #single group z update
# zg.up = function(dat=NULL, par=NULL, pig=NULL)
# {
#   ldens =  logden(par, datax=dat)
#   # n times g matrix
#   # rows are observations, cols are groups
#   
#   val = ldens + rep( log(pig), length(ldens) )
#   val = exp(val -  max(val)) 
#   val = val / sum(val)
#   
#   return(val)
# }


# rzpar = function(dat=NULL, zs=NULL, G=NULL, qd=NULL)
# {
#   pb = dim(dat)[1:2]
#   nobs = dim(dat)[3]
#   
#   #q converted to a 2xG matrix
#   if(is.null(qd))  qd = matrix(rep(c(2,2), times=G), nrow=G, byrow=TRUE)  
#   if(length(qd) == length(pb)) 
#   {  idim = matrix(rep(qd, times=G), nrow=G, byrow=TRUE)  }  
#   if(!(dim(idim)[1] == G) ) stop("qd's must be in the form of a 2xg matrix")
#   
#   val = list()
#   
#   #Diriclet parameter and draws
#   a = matrix( runif(nobs*G, min=0.1, max=5), nrow=nobs, byrow=TRUE)
#   Zm = t(apply(a, 1, dirdraw ))
#   
#   #Generate the parameters
#   val$dpar = lapply(seq_len(G), function(g)
#   {
#     cyc1(datax=dat, par)
#   })
#   
#   
# }
# 
# ranz = function(dat=NULL, g=2, qd=c(2,2))
# {
#   
#   nobs = dim(datax)[3]
#   val=list()
#   
#   #Diriclet parameter and draws
#   a = matrix( runif(nobs*g, min=0.1, max=5), nrow=nobs, byrow=TRUE)
#   Zinit = t(apply(a, 1, dirdraw ))
#   
#   #calc params from Z's
#   
#   
#   
# }
# 
# s2frmZ = function(dt=NULL, zs=NULL, q=NULL){
#   par = list()
#   
#   ng = sum(zs)
#   #mu
#   par$mu = apply( 
#     sweep(dt, 3, as.numeric(zs)/ng, "*"), 
#     2, 
#     rowsums
#   )
#   
#   #
#   
# }


#lol whoops
# mxll = function(dat=NULL, gpar=NULL)
# {
#   #mxll = sumsum zig*logpi + zig*logdets + zig*trace
#   #     = sum_g ng*logpig + ng*logdets -0.5*trace{ W2 %*% S1 }
#   nr = dim(dat)[3]
#   G = length(gpar$pi)
#   
#   #logdets each group
#   lds = sapply(1:G, function(i){ logdets(datax=dat, par=gpar$dpar[[i]]) }) #assumes logdets is correct
#   lds = sum((log(gpar$pi) + lds) * gpar$pi * nr) #correct
#   #lds = sum(  (log(gpar$pi) + lds)*colsums(Z)  )
#   
#   trs = sapply(1:G, function(i)
#   { 
#     -0.5*trwz(dt=dat, par=gpar$dpar[[i]], zs=gpar$Z[,i]) 
#   }
#   )/G
#   
#   val = sum(lds) + G*sum(trs) 
#   return( val )
# }



#
# Old SEM function. Replaced by typeSEM
#
# itSECM = function(datx=NULL, gpar=NULL, cutoff=NULL)
# {
#   gpar$Z = zig.up(datx, gpar)
#   gpar$Z = zdraw(gpar$Z)
#   print(colsums(gpar$Z))
#   gpar = Mstep(datx, gpar, thresh=cutoff)
#   
#   return(gpar)
# }

#
# Old SEM functions.
#
# SEMmm = function(data=NULL, gpar0=NULL, niter=NULL, burnin=NULL, ...)
# {
#   
#   if(!is.null(burnin))
#   {
#     for(i in seq_len(burnin))
#     {
#       gpar0 = itSECM(data, gpar0, ...) #... == cutoff basically
#     }
#   }
#   
#   gpar0$lls = gpar0$lls[length(gpar0$lls)]
#   best = gpar0
#   parbar = gpar0
#   bll=gpar0$lls
#   
#   for(i in seq_len(niter))
#   {
#     gpar0 = itSECM(data, gpar0, ...)
#     
#     if(gpar0$lls[i+1] > bll){ best = gpar0; bll = gpar0$lls[i+1]}
#     parbar = gparbar(gpar0, parbar, denom=i)
#     parbar$lls = mxll2(data, parbar)
#   }
#   
#   parbar$Z = zig.up(data, parbar)
#   
#   return(list("max" =best, 
#               "mean"=parbar))
# }
# SEMmean = function(data=NULL, gpar0=NULL, niter=NULL, burnin=NULL, ...)
# {
#   
#   if(!is.null(burnin))
#   {
#     for(i in seq_len(burnin))
#     {  
#       gpar0 = itSECM(data, gpar0, ...)  
#     }
#   }
#   
#   parbar = gpar0
#   
#   for(i in seq_len(niter))
#   {
#     gpar0 = itSECM(data, gpar0, ...) 
#     parbar = gparbar(gpar0, parbar, denom=i)
#   }
#   parbar$Z = zig.up(data, parbar)
#   parbar$lls = mxll2(data, parbar)
#   
#   return(parbar)
# }
# 
# SEMmmz = function(data=NULL, gpar0=NULL, niter=NULL, burnin=NULL, full=TRUE, ...)
# {
#   itSEM = typeSEM(full)
#   
#   if(!is.null(burnin))
#   {
#     for(i in seq_len(burnin))
#     {
#       gpar0 = itSEM(data, gpar0, ...) #... == cutoff basically
#     }
#   }
#   
#   gpar0$lls = gpar0$lls[length(gpar0$lls)]
#   best = gpar0
#   parbar = gpar0
#   bll=gpar0$lls
#   
#   for(i in seq_len(niter))
#   {
#     gpar0 = itSEM(data, gpar0, ...)
#     
#     if(gpar0$lls[i+1] > bll){ best = gpar0; bll = gpar0$lls[i+1]}
#     parbar = gparbar(gpar0, parbar, denom=i)
#     parbar$lls = mxll2(data, parbar)
#   }
#   
#   parbar$Z = zig.up(data, parbar)
#   
#   return(list("max" =best, 
#               "mean"=parbar))
# }


#
# zig update with alternative calculation strategies
#
# zig.up = function(dat=NULL, gpar=NULL)
# {
#   G = length(gpar$pi)
#   ldens = sapply(1:G, FUN=function(g){ logden(datax=dat, gpar$dpar[[g]]) }) #correct if logden is
#   # n times g matrix
#   # rows are observations, cols are groups
#   
#   val = ldens + rep( log(gpar$pi), rep(nrow(ldens),G) ) #correct
#   # logpi = rep( log(gpar$pi), rep(nrow(ldens),G) )
#   # val2 = ldens + matrix(logpi, ncol=G)
#   # val2 = val - max(val)
#   # val2 = exp(val2)
#   
#   val = exp(val -  rowMaxs(val, value=TRUE)) 
#   
#   
#   #val3 = exp(ldens + rep( log(gpar$pi), rep(nrow(ldens),G) ))
#   
#   val = val / rowsums(val)
#   # val2 = val2 /rowsums(val2)
#   # val3 = val3 / rowsums(val3)
#   
#   return(val)
# }

#
# Old BIC. ICL and BIC to close when entropy is negigible.
#
# BICfc = function(dt=NULL, fgpar=NULL, icl=FALSE, mapz=FALSE)
# {
#   G = length(fgpar$pi)
#   nobs = nrow(fgpar$Z)
#   pb = fgpar$dpar[[1]]$pb
#   
#   qds = sapply(1:G, function(g)
#   {
#     val = fgpar$dpar[[g]]$qd
#     return(val)
#   })
#   
#   k1 = (2*pb[1] + pb[2] +0.5)*sum(qds[2,]) + (pb[1]+1)*G
#   k2 = -0.5*sum(qds[2,]^2)
#   k3 = sum(qds[1,]*(qds[1,]-1))
#   
#   zent = 0
#   
#   #Calc entropy
#   if( icl ){ zent=ent(dt, fgpar, map=mapz) }
#   
#   bic = -2*max(fgpar$lls) + zent +  (k1+k2+k3)*log(nobs)
#   
#   return(bic)
# }