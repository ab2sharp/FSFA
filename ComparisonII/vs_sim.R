#
# filename: vs_sim.R
#

# Output: Classification results of our method and funHDDC on the parameters 
#         found in the file `params_d#.R`
#

source("R/frankenclustAuto.R")
source("R/comps_funs.R")
library(purrr)
library(gmfd)
library(fda)
library(funHDDC)
library(pbapply)
library(doParallel)
load("params_d7.RData") 



# Helpers -----------------------------------------------------------------

scale.array <- function(x=NULL)
{
  xx = apply(x, c(1,2), scale)
  xx = aperm(xx, c(2,3,1))  
  return(xx)
}

try.AECM   = possibly(AECMn, otherwise=NA)
try.init   = possibly(rgzpar, otherwise=NA)
try.fkmeans = possibly(gmfd_kmeans, otherwise=NA)


aecm_gOG = function(cdat=NULL,Bmat=NULL, Gee=NULL, qds=NULL, autoq=FALSE)
{
  bll=-Inf;  best=NA;  its=0;  fail=FALSE; #loop parameters
  
  pbn = dim(cdat)
  
  while(is.na(best) && its <20)
  {
    best = try.init(
      dt   = cdat, 
      g    = Gee, 
      ddim = c(pbn[1], pbn[2]), 
      idim = qds, 
      hW   = Bmat, 
      auto = autoq
    )
    its  = its+1
  }
  if(is.na(best)[1])
  {
    fail=TRUE
    res = NA;
  }else
  {
    #mini em
    for(k in seq_len(30))
    {
      gpars = try.init(cdat, g=Gee, ddim=c(pbn[1], pbn[2]), 
                       idim=qds, hW=Bmat, auto=autoq
      )
      res   = try.AECM(cdat, gpars, niter=15)
      
      if( !is.na(res)[1] )
      {
        if(max(res$lls) > bll) { bll=max(res$lls); best=res }
      }
      
    }
    
    #full EM
    res = try.AECM(cdat, best, niter=100)
    
  }
  
  if(is.na(res)[1])
  {
    BIC = Inf
    pclass = rep("nay", times=pbn[3])
    pG = 1
    npars = NA
  }else
  {
    BIC = BICfc(cdat, res)
    pclass = rowMaxs(res$Z)
    pG = length(unique(pclass))
    npars = cplx(res)
  }
  
  val = list( "BIC"        = BIC, 
              "class"      = pclass, 
              "G"          = pG, 
              "pars"       = res, 
              "complexity" = npars )
}
aecm_g = function(cdat=NULL,Bmat=NULL, Gee=NULL, qds=NULL, autoq=FALSE)
{
  #qds should now be a list
  
  #loop parameters
  bll=-Inf;  best=NA;  its=0;  fail=FALSE; bBIC=Inf; bestmod=NA;
  
  pbn = dim(cdat)
  
  if(pbn[1] < 10) thold=0.68
  else thold=0.40
  
  for( i in seq_along(qds) )
  {
    
    while(is.na(best) && its <20)
    {
      best = try.init(
        dt   = cdat, 
        g    = Gee, 
        ddim = c(pbn[1], pbn[2]), 
        idim = qds[[i]], 
        hW   = Bmat, 
        auto = autoq
      )
      its  = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = NA;
    }else
    {
      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(cdat, g=Gee, ddim=c(pbn[1], pbn[2]), 
                         idim=qds[[i]], hW=Bmat, auto=autoq
        )
        res = try.AECM(cdat, gpars, niter=15, cutoff=thold)
        
        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }
        
      }
    }
    
    if(is.na(res)[1])
    {
      BIC = Inf
    }else
    {
      BIC = BICfc(cdat, res)
    }
    
    if( BIC <= bBIC ){ bestmod=res; bBIC=BIC }
    
    bll = -Inf; best=NA; its=0; fail=FALSE;
    
  }
    
    #full EM
    res = try.AECM(cdat, bestmod, niter=100, cutoff=thold)
    
  
  
  if(is.na(res)[1])
  {
    BIC = Inf
    pclass = rep("nay", times=pbn[3])
    pG = 1
    npars = NA
    qdl = rep(NA, times=4)
  }else
  {
    BIC = BICfc(cdat, res)
    pclass = rowMaxs(res$Z)
    pG = length(unique(pclass))
    npars = cplx(res)
    qdl = c(res$dpar[[1]]$qd, res$dpar[[2]]$qd)
  }
  
  val = list( "BIC"        = BIC, 
              "class"      = pclass, 
              "G"          = pG, 
              "pars"       = res, 
              "complexity" = npars,
              "qdl"        = qdl)
}
sim_fun = function( dstr=NULL, nstr=NULL, pset=NULL, nset=NULL, 
                    ngrp=NULL, qdset=NULL, sb=FALSE, lT=25 )
{
  gpar = pset[[dstr]]
  nobs = nset[[nstr]]
  ddim = dim(gpar$dpar[[1]]$mu)
  
  if(!sb){ mdat = gpardat( pset[[dstr]], nset[[nstr]] ) }
  else   { mdat = gp2dat(  pset[[dstr]], nset[[nstr]] ) }
  
  
  fb = create.fourier.basis( nbasis=dim(mdat)[2] )
  fdlist = coef2fd( coefs=mdat, abasis=fb )
  
  bipmat = diag(dim(mdat)[2])
  
  fb$rangeval |> (\(u){ seq(u[1], u[2], length.out=lT) })() -> tts
  gmfdat = fdlist2gmfd(fdl=fdlist, ts=tts)
  
  s1 = Sys.time()
  v1 = aecm_g( cdat=mdat, Bmat=bipmat, Gee=ngrp, qds=qdset )
  f1 = Sys.time()
  t1 = f1-s1
  
  s2= Sys.time()
  v2 = funHDDC( fdlist, K=ngrp, itermax=100, init="random", nb.rep=40 )
  f2 = Sys.time()
  t2 = f2-s2
  
  
  s3 = Sys.time()
  v3 = try.fkmeans(FD=gmfdat, metric="L2", n.cl=ngrp, p=10^8)
  f3 = Sys.time()
  t3 = f3-s3
  
  if(is.na(v3)[1]){ kclass = rep("nay", times=dim(mdat)[3]) }
  else{             kclass = v3$cluster                    }
  
  
  #Return results
  tclass = rep( seq_len(ngrp), each=dim(mdat)[3]/ngrp )
  
  list( "ARI"        = c(mixture::ARI(v1$class, tclass), 
                         mixture::ARI(v2$class, tclass),
                         mixture::ARI(kclass, tclass) ),
        "Times"      = c(t1, t2, t3),
        "Complexity" = list(v1$complexity, v2$complexity, v1$qdl, v2$d) 
  )
  
  
}






# Parameters --------------------------------------------------------------

m    = 100 #reps
G    = 2   #number of groups in data
qds1 = list( c(2,3), c(3,3), c(4,3) ) # hyperparameters for mfsf



# Test Simulation ---------------------------------------------------------

simsets = expand.grid(parnames, nobsnames, stringsAsFactors=FALSE)

stringr::str_detect(parnames, "n1$") |>
  ifelse(nobsnames[1], nobsnames[2]) |>
  (  \(u){ cbind(parnames, u) }  )() |>
  cbind(c(rep(TRUE, times=8), rep(FALSE, times=12))) ->
  simsets

hands = makeCluster( 120, type="FORK" )
registerDoParallel(hands)

res = foreach( u=seq_len(nrow(simsets)) ) %:%
  foreach( i=seq_len(m) ) %dopar%
  {
    sim_fun(dstr  = simsets[u,1],
            nstr  = simsets[u,2],
            pset  = pars,
            nset  = nobslist,
            ngrp  = G,
            qdset = qds1,
            sb    = as.logical(simsets[u,3])
    )
  }

stopCluster(hands)


save(res, simsets, file="results_d7.RData")

# //