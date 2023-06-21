#Parameter Recovery
#Generate Data for simulation and run sim
source("frankenclustAuto.R")
source("parrec_funs.R")
library(purrr)

load("trueparams4.RData") 

# II is with scale=2
# III is with scale=4 

# Helpers -----------------------------------------------------------------

scale.array <- function(x=NULL) {
  xx = apply(x, c(1,2), scale)
  xx = aperm(xx, c(2,3,1))
  return(xx)
}
try.AECM = possibly(AECMn, otherwise=NA)
try.init = possibly(rgzpar, otherwise=NA)

# Simulation --------------------------------------------------------------

#simulation parameters

m=5e2
seeds = 90053 + seq_len(m)
mdim = c(4,4) 
sdim = c(2,2)
rW = diag( mdim[2] )

start =  Sys.time()
cls = makeCluster( 110, type="FORK" )
registerDoParallel( cls )

#2 groups
sim_g25e1v1 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g25e1v1"]], nobslist[["n5e1_g2"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e1_g2"]])) )

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
}
save(sim_g25e1v1, file="sim_g25e1v1_III.RData")
rm(sim_g25e1v1)

sim_g21e2v1 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g21e2v1"]], nobslist[["n1e2_g2"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n1e2_g2"]])) )

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g21e2v1, file="sim_g21e2v1_III.RData")
rm(sim_g21e2v1)


sim_g25e2v1 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g25e2v1"]], nobslist[["n5e2_g2"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e2_g2"]])) )

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g25e2v1, file="sim_g25e2v1_III.RData")
rm(sim_g25e2v1)




sim_g25e1v2 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g25e1v2"]], nobslist[["n5e1_g2"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e1_g2"]])) )

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g25e1v2, file="sim_g25e1v2_III.RData")
rm(sim_g25e1v2)

sim_g21e2v2 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g21e2v2"]], nobslist[["n1e2_g2"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n1e2_g2"]])) )

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g21e2v2, file="sim_g21e2v2_III.RData")
rm(sim_g21e2v2)


sim_g25e2v2 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g25e2v2"]], nobslist[["n5e2_g2"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e2_g2"]])) )

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(40))
      {
        gpars = try.init(dat, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=200)

    }

    return( res )
  }
save(sim_g25e2v2, file="sim_g25e2v2_III.RData")
rm(sim_g25e2v2)





#4 Groups
sim_g45e1v1 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g45e1v1"]], nobslist[["n5e1_g4"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e1_g4"]])) )

    gee=length(pars[["g45e1v1"]]$pi)

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g45e1v1, file="sim_g45e1v1_III.RData")
rm(sim_g45e1v1)

sim_g41e2v1 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g41e2v1"]], nobslist[["n1e2_g4"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n1e2_g4"]])) )

    gee=length(pars[["g45e2v1"]]$pi)

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g41e2v1, file="sim_g41e2v1_III.RData")
rm(sim_g41e2v1)

sim_g45e2v1 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g45e2v1"]], nobslist[["n5e2_g4"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e2_g4"]])) )

    gee=length(pars[["g45e2v1"]]$pi)

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g45e2v1, file="sim_g45e2v1_III.RData")
rm(sim_g45e2v1)




sim_g45e1v2 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g45e1v2"]], nobslist[["n5e1_g4"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e1_g4"]])) )

    gee=length(pars[["g45e1v2"]]$pi)

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g45e1v2, file="sim_g45e1v2_III.RData")
rm(sim_g45e1v2)

sim_g41e2v2 = foreach( s=seeds ) %dopar%
{
    set.seed(s)
    mdat = gpardat(pars[["g41e2v2"]], nobslist[["n1e2_g4"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n1e2_g4"]])) )

    gee=length(pars[["g41e2v2"]]$pi)

    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;

    while(is.na(best) && its <20)
    {
      best = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {

      #mini em
      for(k in seq_len(30))
      {
        gpars = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)

        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }

      }

      #full EM
      res = try.AECM(dat, best, niter=100)

    }

    return( res )
  }
save(sim_g41e2v2, file="sim_g41e2v2_III.RData")
rm(sim_g41e2v2)

sim_g45e2v2 = foreach( s=seeds ) %dopar% 
{
    set.seed(s)
    mdat = gpardat(pars[["g45e2v2"]], nobslist[["n5e2_g4"]])
    dat = array( do.call(cbind, mdat), c(mdim, sum(nobslist[["n5e2_g4"]])) )
    
    gee=length(pars[["g45e2v2"]]$pi)
    
    #loop parameters
    bll=-Inf;  best=NA;  its=0;  fail=FALSE;
    
    while(is.na(best) && its <20)
    {
      best = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
      its = its+1
    }
    if(is.na(best)[1])
    {
      fail=TRUE
      res = -1
    }else
    {
      
      #mini em
      for(k in seq_len(40))
      {
        gpars = try.init(dat, g=gee, ddim=mdim, idim=sdim, hW=rW, auto=FALSE)
        res = try.AECM(dat, gpars, niter=20)
        
        if( !is.na(res)[1] )
        {
          if(max(res$lls) > bll) { bll=max(res$lls); best=res }
        }
        
      }
      
      #full EM
      res = try.AECM(dat, best, niter=200)
      
    }
    
    return( res )
  }
save(sim_g45e2v2, file="sim_g45e2v2_III.RData") 
rm(sim_g45e2v2)

stopImplicitCluster()

fin=Sys.time()

print(fin-start)



#