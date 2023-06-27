#Cluster the pitch data
#
# SEM on FF, CH, and SL pitches. Choose q and d auto.
#
source("frankenclustAuto.R")
library(purrr)


# Helpers -----------------------------------------------------------------
scale.array <- function(x=NULL) {
  xx = apply(x, c(1,2), scale)
  xx = aperm(xx, c(2,3,1))
  return(xx)
}

try.AECM = possibly(AECMn, otherwise=NA)
try.init = possibly(rgzpar, otherwise=NA)


# Analysis ----------------------------------------------------------------

#load the data
load("p2c_nonoise.RData")


# Decide which pitches to have in data ------------------------------------

fst = which(keypees$pitch_type %in% c("FF", "CU", "SL", "CH", "SI"))
fstdf = keypees[fst,]
cdat = nfuns[,,fst]


set.seed(90053)
cdat = scale.array(cdat + rnorm(prod(dim(cdat)) ))
tlabs = unlist(fstdf[,2]) #true labels


# SEM runs ----------------------------------------------------------------

#Initialize the parameters
Gs = 5
Qs = c(3,4,5,6) 
qdims = length(Qs)
tqd = expand.grid(1:qdims, 1:qdims, 1:qdims, 1:qdims, 1:qdims)

m=20 #number of sem runs for given row of tqd

cls = makeCluster(120, type="FORK")
registerDoParallel(cls)

#Run some SEMs to get starting parameters
bicsearch = foreach( t=seq_len(nrow(tqd)) )%:%
  foreach( i=seq_len(m) ) %dopar% 
  {
    seed = 10000*t  + 10*i
    set.seed(seed)
    
    qset = Qs[as.numeric(tqd[t,])]
    init = try.init(cdat, g=Gs, ddim=dim(cdat)[1:2], qs=qset, hW=rW, iter=5)
    
    if(!is.na(init)[1])
    {
      bic.val = SEMz(data=cdat, gpar0=init, niter=100, burnin=3)$BIC
    }
    else
    { 
      bic.val = Inf
    }
    return(bic.val)
  }

# Find best q settings ----------------------------------------------------


bestbics = sapply(seq_along(bicsearch), FUN=function(i)
{
  min(unlist(bicsearch[[i]]))
})

best = which.min(bestbics)
bQs = Qs[as.numeric(tqd[best,])] #best setting according to bic

# Fit model on best q settings --------------------------------------------

#Number of initializations from best settings
N=50

#Fit the model N times
mods = foreach(i = seq_len(N)) %dopar% 
  {
    set.seed(90053 + i)  
    init = try.init(cdat, g=Gs, ddim=dim(cdat)[1:2], qs=bQs, hW=rW, iter=5)
    if( !is.na(init)[1] )
    {
      modfit = SEMz(data=cdat, gpar0=init, niter=300, burnin=3)
    }else
    {
      modfit = "nope"
    }
    modfit
  }

stopCluster(cls)
rm(cls)

#We save everything here, just in case. 
save.image(file="sem5pitches_aqad.RData")

#Additionally, keep only the fitted models
# save(mods, file="onlymods_sem5pitches_aqad.RData")