# Parameter Recovery Data Functions

#Simulation of choosing q and d using ICL and BIC
norm_eig2 = function(del=NULL, eta=NULL, p=NULL)
{
  q = length(del)
  logdet = sum(  c(log(del),(p-q)*log(eta))  )
  # if(q > (p-q) & logdet <0 ){
  #   logvals = c(  log(del)-logdet/(2*(p-q)), log(eta)-logdet/(2*(q))  )
  # }
  logvals = c(  log(del)-logdet/p, log(eta)-logdet/p  )
  exp(logvals)
}

#generator for single group
rfpar <- function(pb=c(4,5), qd=c(2,2), noise.scale=1)
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
                                   rep(nrow(temp1$vectors),length(temp1$values))
  ) #works
  val$xi1 =  noise.scale*rep( 0.25, dim(val$lam1)[1] )
  
  #Makes determinant equal to 1
  tv2 = norm_eig2(del=tmpo2, eta=tmpe2 , p=pb[2])
  val$omg2 = tv2[1:qd[2]]
  val$eta2 = tv2[qd[2]+1]
  
  val$gam2 = pracma::randortho(pb[2], type="orthonormal")
  val$lam2 = val$gam2[, 1:qd[2], drop=FALSE]
  
  val$mu = tempmu %*% 
    diag(
      sqrt(val$omg2), nrow=length(val$omg2)
    ) %*% 
    t(val$lam2)
  val$oldmu = tempmu
  val$extra = diag(sqrt(val$omg2), nrow=length(val$omg2)) %*% t(val$lam2)
  val$pb  = pb
  val$qd  = qd
  
  
  return( val )
}


#Use this as a general looping function for the generators above.
genpars = function(g=2, ddim=c(4,5), idim=c(2,2), nobs=NULL,...) 
{
  #q converted to a 2xG matrix
  #nobs a vector of group sizes
  if(is.null(idim))  idim = matrix(rep(c(1,2), times=g), nrow=g, byrow=TRUE)  
  if(length(idim) == length(ddim)) 
  {  idim = matrix(rep(idim, times=g), nrow=g, byrow=TRUE)  }  
  if(!(dim(idim)[1] == g) ) stop("qd's must be in the form of a 2xg matrix")
  
  val=list()
  val$dpar = lapply( seq_len(g), FUN=function(gp)
    { 
      rfpar(pb=ddim, qd=as.numeric(idim[gp,]), ...) 
    }
  )
  
  #Diriclet parameter and draws
  a = matrix( runif(sum(nobs)*g, min=0.1, max=5), nrow=sum(nobs), byrow=TRUE)
  Zm = t(apply(a, 1, dirdraw ))
  
  val$Z   = Zm
  val$pi  = colmeans(Zm)
  val$lls = -Inf
  
  return(val)
}	



#Generate the data
pardat = function(par=NULL, ng=NULL)
{
  #returns data as ng columns, each corresponding to a vectorized obs.
  udat = rnorm( par$qd[1]*par$pb[2]*ng )
  eps  = rnorm( prod(par$pb)*ng )
  uray = array( udat, dim=c(par$qd[1], par$pb[2], ng) )
  eray = array( eps, dim=c(par$pb, ng) )
  
  hs2 = diag( 
      sqrt( 
        c( par$omg2, 
          rep(  par$eta2, times=( par$pb[2]-length(par$omg2) )  )  
        )
      ) 
    ) %*%
    t(par$gam2)
  
  udat = apply(uray, 3, function(u, A=NULL, B=NULL)
  {
    A %*% u %*% B
  }, 
  A=par$lam1, B=hs2 ) #each column an observation
  
  eps = apply(eray, 3, function(e, A=NULL, B=NULL)
  {
    A %*% e %*% B
  }, 
  A=diag(sqrt(par$xi1)), B=hs2 )
  
  dat = c(par$mu) + udat + eps #correct
  
  return(dat) 
  
}

gpardat = function(gpar=NULL, nobs=NULL)
{
  #nobs is a vector of group sample sizes
  dat = sapply( 1:length(nobs), function(g)
  {
    pardat(gpar$dpar[[g]], nobs[g])
  },
  simplify=FALSE )
}


#
