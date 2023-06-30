# Comparison Study Functions

#Simulation of choosing q and d using ICL and BIC
norm_eig2 = function(del=NULL, eta=NULL, p=NULL)
{
  q = length(del)
  logdet = sum(  c(log(del),(p-q)*log(eta))  )
  logvals = c(  log(del)-logdet/p, log(eta)-logdet/p  )
  exp(logvals)
}


# Generate Random Parameters--MFSF ----------------------------------------


rfpar <- function(pb=c(4,5), qd=c(2,2))
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
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
  val$gam2 = temp2$vectors
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=temp2$values[1:qd[2]], eta=temp2$values[qd[2]+1], p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$mu = tempmu %*% diag(sqrt(val$omg2)) %*% t(val$lam2)
  
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

#Generate random parameters
rfpar2 <- function(pb=c(4,5), qd=c(2,2), noise.scale=1)
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  s1  = cov(matrix(  rnorm( prod(2*pb[1]^2), sd=1 ), 2*pb[1], pb[1]  ))
  tmpo2 = sort(runif(qd[2], 50, 100), decreasing=TRUE)
  tmpe2 = runif(1, 0.5, 5)
  
  
  
  tempmu = matrix( rnorm(pb[1]*qd[2], mean=4, sd=2), pb[1], qd[2] ) #pxd
  
  temp1 = eigen.sym(s1, k=qd[1])
  val$lam1  = temp1$vectors * rep( sqrt(temp1$values), 
                                   rep(nrow(temp1$vectors), 
                                       length(temp1$values))
  ) #works
  val$xi1 =  noise.scale*diag(s1 - tcrossprod(val$lam1,val$lam1))
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2 , p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  val$mu = tempmu %*% 
  diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% 
  t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag( sqrt(val$omg2), nrow=length(val$omg2) ) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

#Generate parameters close to subspace clustering
rfparfin <- function(pb=c(4,5), qd=c(2,2), a=0, 
                     unif1=c(100,200), unif2=c(0.5,5) )
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  tmpo1 = sort(runif(qd[1], unif1[1], unif1[2]), decreasing=TRUE)
  tmpe1 = runif(1, unif2[1], unif2[2])
  
  tmpo2 = sort(runif(qd[2], 50, 100), decreasing=TRUE)
  tmpe2 = runif(1, 0.5, 5)
  
  
  tempmu = matrix( rnorm(pb[1]*qd[2], mean=1, sd=1), pb[1], qd[2] ) #pxd
  
  
  #Makes determinant equal to 1
  tv1 = norm_eig2( del=tmpo1, eta=tmpe1 , p=pb[1] )
  omg1 = tv1[1:qd[1]]
  eta1 = tv1[qd[1]+1]
  
  #val$gam1 = pracma::randortho(pb[1], type="orthonormal")
  val$gam1 = diag(pb[1])
  templam = val$gam1[, 1:qd[1], drop=FALSE] #p x d
  
  
  
  if(length(omg1)>=2)
  {
    somg1 = diag(sqrt(omg1 - eta1))
    val$lam1 = templam %*% somg1
  }else
  {
    val$lam1 = omg1*templam
  }
  
  #Draw bridge
  beta = a * (omg1[qd[1]] - eta1)/(pb[1] - qd[1]-1)
  xs = seq_len(pb[1] - qd[1]) + qd[1]
  val$xi1 = c( rep(eta1, qd[1]), eta1 + beta*(pb[1]-xs)) 
  
  
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2 , p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  
  
  val$mu = tempmu %*% diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

rfparfinTest <- function(pb=c(4,5), qd=c(2,2), a=0, 
                     unif1=c(100,200), unif2=c(0.5,5) )
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  tmpo1 = sort(runif(qd[1], unif1[1], unif1[2]), decreasing=TRUE)
  tmpe1 = runif(1, unif2[1], unif2[2])
  
  tmpo2 = sort(runif(qd[2], unif1[1], unif1[2]), decreasing=TRUE)
  tmpe2 = runif(1, unif2[1], unif2[2])
  
  
  tempmu = matrix( rnorm(pb[1]*qd[2], mean=0, sd=0.75), pb[1], qd[2] ) #pxd
  
  
  #Makes determinant equal to 1
  tv1 = norm_eig2( del=tmpo1, eta=tmpe1 , p=pb[1] )
  omg1 = tv1[1:qd[1]]
  eta1 = tv1[qd[1]+1]
  
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2 , p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  #val$gam1 = pracma::randortho(pb[1], type="orthonormal")
  val$gam1 = diag(pb[1])
  templam = val$gam1[, 1:qd[1], drop=FALSE] #p x d
  
  
  
  if(length(omg1)>=2)
  {
    somg1 = diag(sqrt(omg1 - eta1))
    val$lam1 = templam %*% somg1
  }else
  {
    val$lam1 = omg1*templam
  }
  
  #Draw bridge
  beta = a * (omg1[qd[1]] - eta1)/(pb[1] - qd[1]-1)
  xs = seq_len(pb[1] - qd[1]) + qd[1]
  val$xi1 = c( rep(eta1, qd[1]), eta1 + beta*(pb[1]-xs)) 
  
  
  val$gam2 = diag(pb[2])
  val$lam2 = val$gam2[, seq_len(qd[2]), drop=FALSE]
  # val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  # val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  
  
  val$mu = tempmu %*% diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}


#Generate parameters close to subspace clustering
rfparfinII <- function(pb=c(4,5), qd=c(2,2), a=0, 
                       unif1=c(100,200), unif2=c(0.5,5) )
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  tmpo1 = sort(runif(qd[1], unif1[1], unif1[2]), decreasing=TRUE)
  tmpe1 = runif(1, unif2[1], unif2[2])
  
  tmpo2 = sort(runif(qd[2], 50, 100), decreasing=TRUE)
  tmpe2 = runif(1, 0.5, 5)
  
  
  tempmu = matrix( rnorm(pb[1]*qd[2], mean=1, sd=1), pb[1], qd[2] ) #pxd
  
  
  #Makes determinant equal to 1
  #tv1 = norm_eig2(del=tmpo1, eta=tmpe1 , p=pb[1])
  omg1 = tmpo1
  eta1 = tmpe1
  
  #val$gam1 = pracma::randortho(pb[1], type="orthonormal")
  val$gam1 = diag(pb[1])
  templam = val$gam1[, 1:qd[1], drop=FALSE] #p x d
  
  
  
  if(length(omg1)>=2)
  {
    somg1 = diag(sqrt(omg1 - eta1))
    val$lam1 = templam %*% somg1
  }else
  {
    val$lam1 = omg1*templam
  }
  
  #Draw bridge
  beta = a * (omg1[qd[1]] - eta1)/(pb[1] - qd[1]-1)
  xs = seq_len(pb[1] - qd[1]) + qd[1]
  val$xi1 = c( rep(eta1, qd[1]), eta1 + beta*(pb[1]-xs)) 
  
  
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2 , p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  
  
  val$mu = tempmu %*% diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

rfparfinIII <- function( pb=c(4,5), qd=c(2,2) )
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  tempmu = matrix( rnorm(pb[1]*qd[2], mean=1, sd=1), pb[1], qd[2] ) #pxd
  
  s1  = cov(matrix(  rnorm( prod(2*pb[1]^2), sd=1   ), 2*pb[1], pb[1]  ))
 
  temp1 = eigen.sym(s1, k=qd[1])
  val$lam1  = temp1$vectors * rep( sqrt(temp1$values), 
                                   rep(nrow(temp1$vectors), 
                                       length(temp1$values))
  ) #works
  val$xi1 =  diag(runif(pb[1], min=min(temp1), max=max(temp1)))
  
  tmpo2 = sort(runif(qd[2], 50, 100), decreasing=TRUE)
  tmpe2 = runif(1, 0.5, 5)
  
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2 , p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  
  
  val$mu = tempmu %*% 
  diag( sqrt(val$omg2), nrow=length(val$omg2) ) %*% 
  t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

genpars = function(g=2, ddim=c(4,5), idim=c(2,2), nobs=NULL,...) 
{
  #q converted to a 2xG matrix
  #nobs a vector of group sizes
  if(is.null(idim))  idim = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)  
  if(length(idim) == length(ddim)) 
  {  idim = matrix(rep(idim, times=g), nrow=g, byrow=TRUE)  }  
  if(!(dim(idim)[1] == g) ) stop("qd's must be in the form of a 2xg matrix")
  
  val=list()
  val$dpar = lapply(
    1:g, 
    FUN=function(gp){  rfpar2(pb=ddim, qd=as.numeric(idim[gp,], ...))  })
  
  #Diriclet parameter and draws
  a = matrix( runif(sum(nobs)*g, min=0.1, max=5), nrow=sum(nobs), byrow=TRUE)
  Zm = t(apply(a, 1, dirdraw ))
  
  val$Z   = Zm
  val$pi  = colmeans(Zm)
  val$lls = -Inf
  
  return(val)
}	

genparsfin = function(g=2, ddim=c(4,5), idim=c(2,2), acept=0, ...)
{
  #q converted to a 2xG matrix
  #nobs a vector of group sizes
  if( is.null(idim) )  idim = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)
  if( length(idim) == length(ddim) )
  {  idim = matrix(rep(idim, times=g), nrow=g, byrow=TRUE)  }
  if(!(dim(idim)[1] == g) ) stop("qd's must be in the form of a 2xg matrix")


  val=list()
  val$dpar = lapply(seq_len(g), FUN=function(gp)
  {
    rfparfinTest(pb=ddim, qd=idim[gp,], a=acept, ... )
  })

  return(val)
}



# Generate Mus according to specified distance between them
rfpar_nomu <- function(pb=c(4,5), qd=c(2,2), a=0, 
                       unif1=c(100,200), unif2=c(0.5,5))
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  cond=FALSE
  
  while(!cond)
  {
    
    tmpo1 = sort(runif(qd[1], unif1[1], unif1[2]), decreasing=TRUE)
    tmpe1 = runif(1, unif2[1], unif2[2])
    
    tmpo2 = sort(runif(qd[2], unif1[1], unif1[2]), decreasing=TRUE)
    tmpe2 = runif(1, unif2[1], unif2[2])
    
    #Eigenvalues -------------------------------------------
    
    # Makes determinant equal to 1
    tv1 = norm_eig2( del=tmpo1, eta=tmpe1, p=pb[1] )
    omg1 = tv1[1:qd[1]]
    eta1 = tv1[qd[1]+1]
    
    tv2 = norm_eig2(del=tmpo2, eta=tmpe2, p=pb[2])
    omg2 = tv2[1:qd[2]]
    eta2 = tv2[qd[2]+1]
    
    
    omg12 = sort( omg2%x%omg1, decreasing=TRUE )
    om = omg12[length(omg12)] # ???
    
    vtest = sort( c(eta2*omg1, eta1*omg2), decreasing=TRUE, index.return=TRUE ) 
    
    minomg = min(omg12)
    maxnom = max(c(eta2*omg1, eta1*omg2))
    
    if( minomg >= maxnom ) cond = TRUE
  }
  
  val$omg2 = omg2
  val$eta2 = eta2
  
  
  #val$gam1 = pracma::randortho(pb[1], type="orthonormal")
  val$gam1 = diag(pb[1])
  templam = val$gam1[, 1:qd[1], drop=FALSE] #p x d
  
  
  
  if(length(omg1)>=2)
  {
    somg1 = diag(sqrt(omg1 - eta1))
    val$lam1 = templam %*% somg1
  }else
  {
    val$lam1 = omg1*templam
  }
  
  #Draw bridge
  beta = a * (omg1[qd[1]] - eta1)/(pb[1] - qd[1]-1)
  xs = seq_len(pb[1] - qd[1]) + qd[1]
  val$xi1 = c( rep(eta1, qd[1]), eta1 + beta*(pb[1]-xs)) 
  
  
  val$gam2 = diag(pb[2])
  val$lam2 = val$gam2[, seq_len(qd[2]), drop=FALSE]
  # val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  # val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}

genmu <- function(pb=c(4,5), G=2, d=0.75)
{
  #this function currently does not work for G!=2
  
  m1 = matrix( rnorm( prod(pb)), nrow=pb[1], ncol=pb[2] )
  
  rnorm(prod(pb)) |> 
    (\(u){ d * u/sqrt(sum(u^2))  })() |> 
    (\(u){  m1 + u   })() -> 
    m2
  
  list("mu1"=m1, 
       "mu2"=m2
  )
  
  
}

genpars_co <- function(g=2, ddim=c(4,5), idim=c(2,2), accept=0, di=1, ...)
{
  #q converted to a 2xG matrix
  #nobs a vector of group sizes
  #stop()
  if( is.null(idim) )  idim = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)
  if( length(idim) == length(ddim) )
  {  idim = matrix(rep(idim, times=g), nrow=g, byrow=TRUE)  }
  if(!(dim(idim)[1] == g) ) stop("qd's must be in the form of a 2xg matrix")
  
  val = list()
  
  val$dpar = lapply(seq_len(g), FUN=function(grp)
  {
    rfpar_nomu(pb=ddim, qd=idim[grp,], a=accept, ...)
  })
  
  mulist = genmu(pb=ddim, G=g, d=di)
  
  for( i in seq_along(val$dpar))
  {
    val$dpar[[i]]$mu = mulist[[i]]
  }
  
  return( val )
  
}



#Generate Z matrix for each group
genz4par = function(gpars=NULL, nobs=NULL)
{
  #Diriclet parameter and draws
  g = length(gpars$dpar)

  stopifnot("nobs must specify sample size for each component"=length(nobs)==g)
  
  
  a = matrix( runif(sum(nobs)*g, min=0.1, max=5), nrow=sum(nobs), byrow=TRUE)
  Zm = t(apply(a, 1, dirdraw ))
  
  gpars$Z   = Zm
  gpars$pi  = colmeans(Zm)
  gpars$lls = -Inf
  
  return(gpars)
}

#Generate the data
pardat = function(par=NULL, ng=NULL)
{
  udat = rnorm(par$qd[1]*par$pb[2]*ng)
  eps  = rnorm(prod(par$pb)*ng)
  uray = array(udat, dim=c(par$qd[1], par$pb[2], ng))
  eray = array(eps, dim=c(par$pb, ng))
  
  hs2 = diag( 
    sqrt( 
      c(   par$omg2, rep( par$eta2, times=(par$pb[2]-length(par$omg2)) )   )
    ) ) %*% 
    t(par$gam2)
  
  udat = apply( uray, 3, function(u, A=NULL, B=NULL)
  {
    A %*% u %*% B
  }, 
  A=par$lam1, B=hs2
  ) #each column an observation
  
  eps = apply(eray, 3, function(e, A=NULL, B=NULL)
  {
    A %*% e %*% B
  }, 
  A=diag(sqrt(par$xi1)), B=hs2
  )
  
  
  dat = c(par$mu) + udat + eps #correct
  #dat = array(dat, c(par$pb, ng))
  
  #returns data as ng columns, each corresponding to a vectorized obs.
  return(dat) 
  
}

gpardat = function(gpar=NULL, nobs=NULL)
{
  #nobs a vector of group sizes
  sapply(1:length(nobs), FUN=function(g)
  {
    pardat(gpar$dpar[[g]], nobs[g])
  }) |> 
    array( dim=c(gpar$dpar[[1]]$pb, sum(nobs)) )
}




# Generate Random Parameters--Their Model ---------------------------------

#Generate parameters close to subspace clustering
rspar <- function( pb=c(4,5), qd=c(2,2), unif1=c(50,100), unif2=c(0.5,5), a=0 )
{
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  tmpo1 = sort(runif(qd[1], unif1[1], unif1[2]), decreasing=TRUE)
  tmpe1 = runif(1, unif2[1], unif2[2])
  
  tmpo2 = sort(runif(qd[2], unif1[1], unif1[2]), decreasing=TRUE)
  tmpe2 = runif(1, unif2[1], unif2[2])
  
  
  tempmu = matrix(rnorm(pb[1]*qd[2], mean=0, sd=0.75), pb[1], qd[2] ) #pxd
  
  #Eigenvalues -------------------------------------------
  
  # Makes determinant equal to 1
  tv1 = norm_eig2( del=tmpo1, eta=tmpe1, p=pb[1] )
  omg1 = tv1[1:qd[1]]
  eta1 = tv1[qd[1]+1]
  
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2, p=pb[2])
  omg2 = tv2[1:qd[2]]
  eta2 = tv2[qd[2]+1]
  
  
  omg12 = sort( omg2%x%omg1, decreasing=TRUE )
  om = omg12[length(omg12)]
  
  vtest =  sort( c(eta2*omg1, eta1*omg2), decreasing=TRUE, index.return=TRUE ) 
  
  indies = c( rep(1, times=qd[1]), rep(2, times=qd[2])  )[vtest$ix]
  tts = list( seq_len(pb[2]-qd[2]), seq_len(pb[1]-qd[1]) )
  
  omgval = vector( "list", length=length(vtest$x) )
  
  for( i in seq_along(vtest$x) )
  {
    t = tts[[ indies[i] ]]
    beta = a * (om - vtest$x[i])/t[length(t)]
    
    omgval[[i]] = vtest$x[i] + beta*(t[length(t)]-t)
    om = omgval[[i]][length(t)]
  }
  
  
  val$omg = c(omg12, unlist(omgval))
  val$eta = eta1*eta2
  
  val$omg1 = omg1
  val$omg2 = omg2
  val$eta1 = eta1
  val$eta2 = eta2
  
  
  #Eigenvectors -------------------------------------------
  
  gam1 = diag(pb[1])
  #rotation matrix
  val$r1 = gam1[,seq_len(qd[1]), drop=FALSE]
  
  #Factor loadings
  if(length(omg1)>=2)
  {
    val$lam1 = val$r1 %*% diag(sqrt(omg1 - eta1))
  }else
  {
    val$lam1 = omg1*val$r1
  }
  
  gam2 = diag(pb[2])
  val$lam2 = gam2[, seq_len(qd[2]), drop=FALSE]
  
  val$lam = val$lam2 %x% val$lam1
  val$gam = gam2 %x% gam1
  
  
  
  val$mu = tempmu %*% 
  diag( sqrt(val$omg2), nrow=length(val$omg2) ) %*% 
  t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  
  
  val$pb  = pb
  val$qd  = qd
  
  val$p = prod(val$pb)
  val$d = length(val$omg)
  
  return( val )
}

grspar = function(g=2, ddim=c(4,5), idim=c(2,2), accept=0, ...)
{
  if( is.null(idim) ) idim = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)
  else if( length(idim) == length(ddim) )
  {
    idim =  matrix( rep(idim, times=g), nrow=g, byrow=TRUE )
  }
  if( !(dim(idim)[1] == g) ) stop("qds must be in the form of a 2xg matrix")
  
  val = list()
  val$dpar = lapply( seq_len(g), FUN=function(gp)
  {
    rspar(pb=ddim, qd=idim[gp,], a=accept, ...)  
  })
  
  return( val )
}

#Generate mus according to specified distance
rspar_nomu <- function( 
  pb=c(4,5), qd=c(2,2), unif1=c(50,100), unif2=c(0.5,5), a=0 
) {
  ## pb = (row, col) = (p, b)
  ## qd = (latent gaussian dim, intrinsic basis dim) = (q, d)
  
  val = list()
  
  cond=FALSE
  
  while(!cond)
  {
    
    tmpo1 = sort(runif(qd[1], unif1[1], unif1[2]), decreasing=TRUE)
    tmpe1 = runif(1, unif2[1], unif2[2])
    
    tmpo2 = sort(runif(qd[2], unif1[1], unif1[2]), decreasing=TRUE)
    tmpe2 = runif(1, unif2[1], unif2[2])
    
    #Eigenvalues -------------------------------------------
    
    # Makes determinant equal to 1
    tv1 = norm_eig2( del=tmpo1, eta=tmpe1, p=pb[1] )
    omg1 = tv1[1:qd[1]]
    eta1 = tv1[qd[1]+1]
    
    tv2 = norm_eig2(del=tmpo2, eta=tmpe2, p=pb[2])
    omg2 = tv2[1:qd[2]]
    eta2 = tv2[qd[2]+1]
    
    
    omg12 = sort( omg2%x%omg1, decreasing=TRUE )
    om = omg12[length(omg12)] # ???
    
    vtest =  sort( c(eta2*omg1, eta1*omg2), decreasing=TRUE, index.return=TRUE )
    
    minomg = min(omg12)
    maxnom = max(c(eta2*omg1, eta1*omg2))
    
    if( minomg >= maxnom ) cond = TRUE
  }
  
  #Create intermediate values
  indies = c( rep(1, times=qd[1]), rep(2, times=qd[2])  )[vtest$ix]
  tts = list( seq_len(pb[2]-qd[2]), seq_len(pb[1]-qd[1]) )
  
  omgval = vector( "list", length=length(vtest$x) )
  
  for( i in seq_along(vtest$x) )
  {
    t = tts[[ indies[i] ]]
    beta = a * (om - vtest$x[i])/t[length(t)]
    
    omgval[[i]] = vtest$x[i] + beta*(t[length(t)]-t)
    om = omgval[[i]][length(t)]
  }
  
  
  val$omg = c(omg12, unlist(omgval))
  val$eta = eta1*eta2
  
  val$omg1 = omg1
  val$omg2 = omg2
  val$eta1 = eta1
  val$eta2 = eta2
  
  
  #Eigenvectors -------------------------------------------
  
  gam1 = diag(pb[1])
  #rotation matrix
  val$r1 = gam1[,seq_len(qd[1]), drop=FALSE]
  
  #Factor loadings
  if(length(omg1)>=2)
  {
    val$lam1 = val$r1 %*% diag(sqrt(omg1 - eta1))
  }else
  {
    val$lam1 = omg1*val$r1
  }
  
  gam2 = diag(pb[2])
  val$lam2 = gam2[, seq_len(qd[2]), drop=FALSE]
  
  val$lam = val$lam2 %x% val$lam1
  val$gam = gam2 %x% gam1
  
  
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  
  
  val$pb  = pb
  val$qd  = qd
  
  val$p = prod(val$pb)
  val$d = length(val$omg)
  
  return( val )
}

grspar_co = function(g=2, ddim=c(4,5), idim=c(2,2), accept=0, di=1, ...)
{
  if( is.null(idim) ) idim = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)
  else if( length(idim) == length(ddim) )
  {
    idim =  matrix( rep(idim, times=g), nrow=g, byrow=TRUE )
  }
  if( !(dim(idim)[1] == g) ) stop("qds must be in the form of a 2xg matrix")
  
  val = list()
  
  val$dpar = lapply( seq_len(g), FUN=function(gp)
  {
    rspar_nomu(pb=ddim, qd=idim[gp,], a=accept, ...)  
  })
  
  mulist = genmu(pb=ddim, G=g, d=di)
  
  
  for( i in seq_along(val$dpar))
  {
    val$dpar[[i]]$mu = mulist[[i]]
  }
  
  return( val )
}

p2dat = function(par=NULL, ng=NULL, as_array=FALSE)
{
  vm = as.numeric(par$mu)
  par |> 
  ( \(u)
    {
      u$gam %*% 
      diag(   sqrt(  c( u$omg, rep(u$eta, times=(u$p-u$d)) )  )   ) %*% 
      t(u$gam)  
    }
  )() ->
    sig
  
  rnorm( par$p*ng ) |> 
    matrix( nrow=par$p )->
    zdat

  val = sig %*% zdat + rep(vm, times=ng)
  
  
  if(as_array){  val = array(val, dim=c(par$pb, ng))  }
  
  val
}

gp2dat = function(gpar=NULL, ngs=NULL)
{
  sapply(seq_along(ngs), FUN=function(g)
  {
    p2dat(gpar$dpar[[g]], ng=ngs[g], as_array=FALSE)  
  }) |> 
    array(  dim=c( gpar$dpar[[1]]$pb, sum(ngs) )  ) ->
    dat
}









# Generate Functions ------------------------------------------------------

h1 = Vectorize(function(t)
{
  if(abs(t-15) < 6) v = 6-abs(t-15)
  else v=0
  
  v
})

h2 = Vectorize(function(t)
{
  if(abs(t-7) < 6) v = 6-abs(t-7)
  else v=0
  
  v 
})


gentridat = function(scale=1, seed=90053)
{
  g=4
  n=250
  tt=101
  mn=1
  mx=21
  
  set.seed(seed)
  
  u = runif( g*n, 0 , 0.1 )
  ones = rep(1, times=tt)
  t = seq(mn, mx, length.out=tt)
  h1t = scale*h1(t)
  h2t = scale*h2(t)
  
  
  g1 = cbind(u[1:n], 1-u[1:n]) %*% rbind(ones, h1t) + rnorm(n*tt, 0, sqrt(0.25))
  g2 = cbind( u[ (n+1):(2*n) ], 1-u[ (n+1):(2*n) ] ) %*% rbind( ones, h2t )+ 
  rnorm( n*tt, 0, sqrt(0.25) )
  g3 = cbind( u[ (2*n+1):(3*n) ], 0.5-u[ (2*n+1):(3*n) ] ) %*% rbind(ones, h1t)+
  rnorm(n*tt, 0, sqrt(0.25))
  g4 = cbind( u[ (3*n+1):(4*n) ], 0.5-u[ (3*n+1):(4*n) ] ) %*% rbind(ones, h2t)+
  rnorm( n*tt, 0, sqrt(0.25) )
  
 val1 = t(  do.call(rbind, list(g1,g2,g3,g4))  )
 
 
 #dim 2
 u = runif( g*n, 0 , 0.1 )
 ones = rep(1, times=tt)
 t = seq(mn,mx, length.out=tt)
 h1t = scale*h1(t)
 h2t = scale*h2(t)
 
 
 g1 = cbind(u[1:n], 0.5-u[1:n]) %*% rbind(ones, h1t)+ rnorm(n*tt, 0, sqrt(0.25))
 g2 = cbind(u[ (n+1):(2*n) ], 0.5-u[ (n+1):(2*n) ]) %*% rbind( ones, h2t )+ 
 rnorm(n*tt, 0, sqrt(0.25))
 g3 = cbind( u[ (2*n+1):(3*n) ], 1-u[(2*n+1):(3*n) ]) %*% rbind( ones, h2t )+ 
 rnorm(n*tt, 0, sqrt(0.25))
 g4 = cbind(u[(3*n+1):(4*n) ], 1-u[(3*n+1):(4*n) ]) %*% rbind( ones, h1t )+ 
 rnorm(n*tt, 0, sqrt(0.25))
 
 val2 = t(  do.call(rbind, list(g1,g2,g3,g4))  )
 
 val = list( "fundat" = list(val1, val2),
             "spars"  = list("g"=g, "n"=n, "tt"=tt, "range"=c(mn,mx), "time"=t)
 )
}


plot_mat = function(fdat=NULL, ti=NULL, nobs=NULL, strname="rplot" )
{
  
  png(file=paste0(strname, ".png"),
      width=600, height=400)
  matplot(ti, 
          fdat, 
          type="l",
          col=rep(c("red", "black", "blue", "green"), each=nobs )
  )
  dev.off()
  
  return("done plotting")
}


fit_trifuns = function(fdat=NULL, ts=NULL, nbase=NULL)
{
  
  pseq = seq_along(fdat)
  #fnames = paste0("value", pseq)
  
  bbasis = create.bspline.basis(
    rangeval = c(ts[1], ts[length(ts)]), 
    nbasis   = nbase
  )
  
  funs = lapply(pseq, FUN=function(i, ...)
  {
    smooth.basis( y=fdat[[i]], ... )$fd
  }, 
    argvals=ts, fdParobj=bbasis
  )
  
  
  funs
}


plot_fd = function(fdat=NULL, interval=NULL, nobs=NULL, pname=NULL,cl=NULL)
{
  if( is.null(cl) )
  {
    cl = rep( c("red", "black", "blue", "green"), each=nobs )
  }
  pseq = seq_along(fdat)
  pnames = paste0(pname, pseq)
  
  for(i in pseq)
  {
    sdat = eval.fd( fdat[[i]], interval )
    
    png( file=paste0(pnames[i], ".png"),
         width=600, 
         height=400
    )
    
    graphics::matplot( interval, sdat, type="l", col=cl )
    dev.off()
  }
  
  return("done plotting")
}


fit_funhddc = function( fdat=NULL, seed=NULL, pll=TRUE, ...)
{
  set.seed(seed)
  if(pll) ncls = parallel::detectCores()-2
  else    ncls = 1
  
  fdat |> 
    funHDDC(mc.cores=ncls,...)
  
}

fd_2_coef = function(fdat=NULL)
{
  #I already know this is slow, but it works for now.
  pseq = seq_along(fdat)
  ary = array(  dim=c( length(fdat), dim(fdat[[1]]$coefs) )  )
  
    for( i in pseq ) 
    {
      ary[i,,] = fdat[[i]]$coefs
    }
  
  ary
} 

fd2coef_v2 = function(fdat=NULL)
{
  for(i in seq_along(fdat))
  {
    fdat[[i]] = fdat[[i]]$coefs
  }
  
  bn = dim(fdat[[1]])
  p = length(fdat)
  
  array( unlist(fdat), dim=c(bn,p) ) |> 
    aperm(c(3,1,2))->
    val
  val
}

fd2inprod = function(fdlist=NULL)
{
  abasis = fdlist[[1]]$basis
  val    = pracma::sqrtm( inprod(abasis,abasis) )
}

coef2fd = function( coefs=NULL, abasis=NULL )
{
  fdlist = lapply( seq_len(dim(coefs)[1]), FUN=function(u)
  {
    fd(coef=coefs[u,,], basisobj=abasis)  
  })
  
  fdlist
}

fdlist2vals = function(fdl=NULL, ts=NULL)
{
  lapply( seq_along(fdl), FUN=function(i)
  {
    eval.fd(ts, fdl[[i]]) |> t()
  })
}

fdlist2gmfd = function(fdl=NULL, ts=NULL)
{
  lapply( seq_along(fdl), FUN=function(i)
  {
    eval.fd(ts, fdl[[i]]) |> t()
  }) |> 
    ( \(u){ gmfd::funData(grid=ts, data=u) } )() ->
    gmfdat
  
  gmfdat
}



# Model Functions ---------------------------------------------------------

dq_2_k = function(qd=NULL, pb=NULL)
{
  return( qd[2]*pb[1] + qd[1]*pb[2] - prod(qd) )
}

# //