# Run simulation on their data
#Run targets and analyze the results
setwd("/u/ab2sharp/ResearchDocuments/Research/FunMatNorm/ComparisonStudy/github_version/")
library(fda)
library(funHDDC)
library(purrr)
library(pbapply)
source("R/comps_funs.R")
source("R/frankenclustAuto.R")


# Helpers -----------------------------------------------------------------

scale.array <- function(x=NULL) {
  xx = apply(x, c(1,2), scale)
  xx = aperm(xx, c(2,3,1))
  return(xx)
}
try.AECM = possibly(AECMn, otherwise=NA)
try.init = possibly(rgzpar, otherwise=NA)
try.funclust = possibly(funclust, otherwise=NA)

sim_tridat = function(s=NULL)
{

  #Generate data
  d1 = gentridat(seed=s)
  d2 = gentridat(scale=-1,seed=s)
  
  dtot = list(
    "fundat"=c(d1$fundat, d2$fundat), 
    "spars"=d1$spars
  )
  
  psim = length(dtot$fundat)
  bsim = 25
  nsim = 1e3
  
  
  pdat = fit_trifuns(
    fdat  = dtot$fundat,
    ts    = dtot$spars$time,
    nbase = bsim
  )
  
  
  
  # fit funhddc to data
  fhddcr = funHDDC(data=pdat,K=4, init="random", algo="EM", nb.rep=20)
  fhddck = funHDDC(data=pdat,K=4)
  convr = is.na(fhddcr$allCriteria$BIC)
  convk = is.na(fhddck$allCriteria$BIC)

  if( all(!convr, !convk) )
  {
    if( fhddcr$BIC > fhddck$BIC ) fhddcb = fhddcr
    else fhddcb = fhddck
    fBIC = fhddcb$BIC
    fclass = fhddcb$class
  }
  else if(!convr && convk )
  {
    fBIC = fhddcr$BIC
    fclass = fhddcr$class
  }
  else if(convr && !convk )
  {
    fBIC = fhddck$BIC
    fclass = fhddck$class
  } else
  {
    fBIC   = "nope"
    fclass = rep("nah", times=nsim)
  }
  if( fhddcr$BIC > fhddck$BIC ) fhddcb = fhddcr
  else fhddcb = fhddck
  
  #fit funclust to the data
  fclust = try.funclust(pdat, K=4, nbInit=20, increaseDimension=TRUE, hard=TRUE)
  
  if( is.na(fclust)[1] )
  {
    jBIC = "no"
    jclass=rep("nada", times=nsim)
  }else
  {
    jBIC=fclust$bic
    jclass=fclust$cls
  }
  
  
  # fit mfsf model to data
  fdat = fd2coef_v2(pdat)
  rW   = fd2inprod(pdat)$B


  bll=-Inf;  best=NA;  its=0;  fail=FALSE; #loop parameters

  while(is.na(best) && its <20)
  {
    best = try.init(
      fdat,g=4, ddim=c(psim, bsim), idim=c(2,2), hW=rW, auto=FALSE
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
    for(k in seq_len(20))
    {
      gpars = try.init(
        fdat,g=4, ddim=c(psim, bsim), idim=c(2,2), hW=rW, auto=FALSE
      )
      res   = try.AECM(fdat, gpars, niter=20)

      if( !is.na(res)[1] )
      {
        if(max(res$lls) > bll) { bll=max(res$lls); best=res }
      }

    }

    #full EM
    res = try.AECM( fdat, best, niter=100 )

  }
  if( is.na(res)[1] )
  {
    rBIC = "negatory"
    rclass = rep("nay", times=nsim)
  }
  else
  {
    rBIC = BICfc(fdat, res)
    rclass = res$Z
  }


  val = list(
    estclass = list( fclass, rclass, jclass ),
    BICs     = list( fBIC, rBIC, jBIC )
  )

  return(val)
}


# Test the functions ------------------------------------------------------
m=1000
seeds = 90053 +seq_len(m)

cls = makeCluster(120, type="FORK")
registerDoParallel(cls)

res = pblapply(cls, FUN=function(i){ sim_tridat(i)} )


stopCluster(cls)
rm(cls)

save(res, file="sim_tridat_results.RData")


















#test2 = sim_tridat(90053)

# #Test
# #Generate data
# d1 = genbouvdat(seed=90053)
# d2 = genbouvdat(scale=-1,seed=90053)
# 
# dtot = list(
#   "fundat"=c(d1$fundat, d2$fundat), 
#   "spars"=d1$spars
# )
# 
# psim = length(dtot$fundat)
# bsim = 25
# nsim = 1e3
# 
# plot_mat(
#   fdat    = dtot$fundat[[3]], 
#   ti      = dtot$spars$time, 
#   nobs    = dtot$spars$n, 
#   strname = "funs3"
# )
# plot_mat(
#   fdat    = dtot$fundat[[4]], 
#   ti      = dtot$spars$time, 
#   nobs    = dtot$spars$n, 
#   strname = "funs4"
# )
