#Parameter Recovery
setwd("~/Documents/unidrive/ResearchDocuments/Research/FunMatNorm/ParameterRecovery")
#setwd("/u/ab2sharp/ResearchDocuments/Research/FunMatNorm/ParameterRecovery")
source("frankenclustAuto.R")
source("parrec_funs.R")


#Sample size and number of clusters and noise ???
#seed (set as desired)
seed1 = 90053
seed2 = 90053
scl = 4
n5e1_g2 = rep(50, 2)
n5e1_g4 = rep(50, 4)

n1e2_g2 = rep(1e2, 2)
n1e2_g4 = rep(1e2, 4)

n5e2_g2 = rep(5e2, 2)
n5e2_g4 = rep(5e2, 4)



#G=2

#Sample sizes and noise
set.seed(seed1)
g25e1v1 = genpars(g=2, ddim=c(4,4), idim=c(2,2), nobs=n5e1_g2)
set.seed(seed1)
g21e2v1 = genpars(g=2, ddim=c(4,4), idim=c(2,2), nobs=n1e2_g2)
set.seed(seed1)
g25e2v1 = genpars(g=2, ddim=c(4,4), idim=c(2,2), nobs=n5e2_g2)

set.seed(seed2)
g25e1v2 = genpars(g=2, ddim=c(4,4), idim=c(2,2), nobs=n5e1_g2, noise.scale=scl)
set.seed(seed2)
g21e2v2 = genpars(g=2, ddim=c(4,4), idim=c(2,2), nobs=n1e2_g2, noise.scale=scl)
set.seed(seed2)
g25e2v2 = genpars(g=2, ddim=c(4,4), idim=c(2,2), nobs=n5e2_g2, noise.scale=scl)



set.seed(seed1)
g45e1v1 = genpars(g=4, ddim=c(4,4), idim=c(2,2), nobs=n5e1_g4)
set.seed(seed1)
g41e2v1 = genpars(g=4, ddim=c(4,4), idim=c(2,2), nobs=n1e2_g4)
set.seed(seed1)
g45e2v1 = genpars(g=4, ddim=c(4,4), idim=c(2,2), nobs=n5e2_g4)

set.seed(seed2)
g45e1v2 = genpars(g=4, ddim=c(4,4), idim=c(2,2), nobs=n5e1_g4, noise.scale=scl)
set.seed(seed2)
g41e2v2 = genpars(g=4, ddim=c(4,4), idim=c(2,2), nobs=n1e2_g4, noise.scale=scl)
set.seed(seed2)
g45e2v2 = genpars(g=4, ddim=c(4,4), idim=c(2,2), nobs=n5e2_g4, noise.scale=scl)






# For each of these parameter sets, I need to generate x number of datasets
# fit them, get the parameter estimates, and check how close they are to the
# true values.

strings=ls()
param.names = stringr::str_detect(strings, "^(g2|g4)")
parnames = strings[param.names]

length(parnames) == 12

# make a list of the variables and set the names using this character vector.
# loop over the character vector and generate datasets for these parameters
# that we can then use for fitting the model. 

pars = vector( "list", length(parnames) )

for( i in seq_along(parnames) )
{
  pars[[i]] = get( parnames[i] )
}

names(pars) = parnames


#Save the sample sizes in a list as well

nobs.names = stringr::str_detect(strings, "(_g2|_g4)$")
nobsnames = strings[nobs.names]

nobslist = vector("list", length(nobsnames) )

for(i in seq_along(nobsnames))
{
  nobslist[[i]] = get(nobsnames[i])
}

names(nobslist) = nobsnames


save(pars, parnames, nobslist, nobsnames, file="params_fixxi.RData")
