#Generate comparison parameter sets

source("R/frankenclustAuto.R")
source("R/comps_funs.R")


#
# Purpose: Creates parameters using comps_funs.R for testing funHDDC on 
#          data generated according to our model. 
#
#
#
# Output: `params_d#.RData`,  which contains the parameter sets to be
#         used in the file `vs_sim.R`
#


#
# This script is run for different values of di to obtain different parameter
# sets. The values of di chosen are:
#
#     3.5, 5, 7
#
#
# This is controlled using the dm parameter. 
#




sd    = 90053
G     = 2
dimd1 = c(6,11)
dimd2 = c(18,33)
qds1  = c(2,3)
qds2  = c(2,3)
uni1  = c(5,5.5)
nset1 = rep(50, G)
#nset2 = rep(100, G)
nset3 = rep(250, G)

dm = 3.5

k1 =qds1[2]*dimd1[1] + qds1[1]*dimd1[2] - prod(qds1) 
k2 =qds2[2]*dimd2[1] + qds2[1]*dimd2[2] - prod(qds2) 

k1; k2

# Hard SB -----------------------------------------------------------------

set.seed(sd)
gparHSBl = grspar_co(G, ddim=dimd1, idim=qds1, accept=1, di=dm, unif1=uni1)

p1.d1n1 = genz4par(gpars=gparHSBl, nobs=nset1)
p1.d1n3 = genz4par(gpars=gparHSBl, nobs=nset3)

set.seed(sd)
gparHSBh = grspar_co(G, ddim=dimd2, idim=qds2, accept=1, di=dm, unif1=uni1)

p1.d2n1 = genz4par(gpars=gparHSBh, nobs=nset1)
p1.d2n3 = genz4par(gpars=gparHSBh, nobs=nset3)




# SB and Noisy Kronecker --------------------------------------------------

set.seed(sd)
gparSBNl = grspar_co(G, ddim=dimd1, idim=qds1, accept=0.5, di=dm, unif1=uni1)

p2.d1n1 = genz4par(gpars=gparSBNl, nobs=nset1)
p2.d1n3 = genz4par(gpars=gparSBNl, nobs=nset3)

set.seed(sd)
gparSBNh = grspar_co(G, ddim=dimd2, idim=qds2, accept=0.5, di=dm, unif1=uni1)

p2.d2n1 = genz4par(gpars=gparSBNh, nobs=nset1)
p2.d2n3 = genz4par(gpars=gparSBNh, nobs=nset3)





# SB and Kronecker --------------------------------------------------------


set.seed(sd)
#lowdim
gparSBKl = genpars_co(G, ddim=dimd1, idim=qds1, accept=0, di=dm, unif1=uni1)

p3.d1n1 = genz4par(gpars=gparSBKl, nobs=nset1)
p3.d1n3 = genz4par(gpars=gparSBKl, nobs=nset3)


set.seed(sd)
#high dim
gparSBKh = genpars_co(G, ddim=dimd2, idim=qds2, accept=0, di=dm, unif1=uni1)

p3.d2n1 = genz4par(gpars=gparSBKh, nobs=nset1)
p3.d2n3 = genz4par(gpars=gparSBKh, nobs=nset3)




# Noisy SB and Kronecker --------------------------------------------------

#middle

set.seed(sd)
#lowdim
gparNKRl = genpars_co(G, ddim=dimd1, idim=qds1, accept=0.5, di=dm, unif1=uni1)

p4.d1n1 = genz4par(gpars=gparNKRl, nobs=nset1)
p4.d1n3 = genz4par(gpars=gparNKRl, nobs=nset3)


set.seed(sd)
#high dim
gparNKRh = genpars_co(G, ddim=dimd2, idim=qds2, accept=0.5, di=dm,  unif1=uni1)

p4.d2n1 = genz4par(gpars=gparNKRh, nobs=nset1)
p4.d2n3 = genz4par(gpars=gparNKRh, nobs=nset3)




# Hard Kronecker ----------------------------------------------------------


#regular model

set.seed(sd)
#lowdim
gparHKRl = genpars_co(G, ddim=dimd1, idim=qds1, accept=1, di=dm, unif1=uni1)

p5.d1n1 = genz4par(gpars=gparHKRl, nobs=nset1)
p5.d1n3 = genz4par(gpars=gparHKRl, nobs=nset3)


set.seed(sd)
#high dim
gparHKRh = genpars_co(G, ddim=dimd2, idim=qds2, accept=1, di=dm, unif1=uni1)    

p5.d2n1 = genz4par(gpars=gparHKRh, nobs=nset1)
p5.d2n3 = genz4par(gpars=gparHKRh, nobs=nset3)




# Numerical Test ----------------------------------------------------------


#all should return same value as dm
(p1.d1n1$dpar[[1]]$mu - p1.d1n1$dpar[[2]]$mu)|> 
  (\(u){ sqrt(sum(u^2))  })()

(p1.d2n1$dpar[[1]]$mu - p1.d2n1$dpar[[2]]$mu)|> 
  (\(u){ sqrt(sum(u^2))  })()

(p5.d1n1$dpar[[1]]$mu - p5.d1n1$dpar[[2]]$mu)|> 
  (\(u){ sqrt(sum(u^2))  })()

(p5.d2n1$dpar[[1]]$mu - p5.d2n1$dpar[[2]]$mu)|> 
  (\(u){ sqrt(sum(u^2))  })()



# Checking that each set of parameters produced the same eigenvalues
# and means.
id=1
p1.d1n1$dpar[[id]]$eta2; p2.d1n1$dpar[[id]]$eta2; p3.d1n1$dpar[[id]]$eta2
p4.d1n1$dpar[[id]]$eta2; p5.d1n1$dpar[[id]]$eta2

p1.d1n1$dpar[[id]]$mu[2,2]; p2.d1n1$dpar[[id]]$mu[2,2]; 
p3.d1n1$dpar[[id]]$mu[2,2]; p4.d1n1$dpar[[id]]$mu[2,2]; 
p5.d1n1$dpar[[id]]$mu[2,2]



# Plotting the eigenvalues ------------------------------------------------


# Plot the eigenvalues to make sure subspace clustering assumptions 
# arent satisfied

#pard = par()
id=2
#par(mfrow=c(5,1), mar=c(1, 1, 1, 1))
plot(
  c( 
    p1.d1n1$dpar[[id]]$omg, 
    rep(
      p1.d1n1$dpar[[id]]$eta, 
      times=( prod(dimd1)-length(p1.d1n1$dpar[[id]]$omg) )
    )
), 
type="b")

plot(
  c(
    p2.d1n1$dpar[[id]]$omg, 
    rep(
      p2.d1n1$dpar[[id]]$eta, 
      times=( prod(dimd1)-length(p2.d1n1$dpar[[id]]$omg) )
    )
), 
type="b")

tlam1 = p3.d1n1$dpar[[id]]$lam1; txi1 = p3.d1n1$dpar[[id]]$xi1
tsig1 = tlam1 %*% t(tlam1) + diag(txi1)
tdel2 = diag( 
  c(
    p3.d1n1$dpar[[id]]$omg2, 
    rep(p3.d1n1$dpar[[id]]$eta2, times=(dimd1[2]-qds1[2])) 
  ) 
)
tsig2 = p3.d1n1$dpar[[id]]$gam2 %*% tdel2 %*% t(p3.d1n1$dpar[[id]]$gam2)
tsig = tsig2 %x% tsig1
plot(sort( eigen(tsig)$values, decreasing=TRUE ), type="b")
#plot( sort(eigen(tsig1)$values, decreasing=TRUE), type="b")



tlam1 = p4.d1n1$dpar[[id]]$lam1; txi1 = p4.d1n1$dpar[[id]]$xi1
tsig1 = tlam1 %*% t(tlam1) + diag(txi1)
tdel2 = diag( 
  c(
    p4.d1n1$dpar[[id]]$omg2, 
    rep( p4.d1n1$dpar[[id]]$eta2, times=(dimd1[2]-qds1[2]) ) 
  ) 
)
tsig2 = p4.d1n1$dpar[[id]]$gam2 %*% tdel2 %*% t(p4.d1n1$dpar[[id]]$gam2)
tsig = tsig2 %x% tsig1
plot(sort( eigen(tsig)$values, decreasing=TRUE ), type="b")
#plot( sort(eigen(tsig1)$values, decreasing=TRUE), type="b")


tlam1 = p5.d1n1$dpar[[id]]$lam1; txi1 = p5.d1n1$dpar[[id]]$xi1
tsig1 = tlam1 %*% t(tlam1) + diag(txi1)
tdel2 = diag( 
  c(
    p5.d1n1$dpar[[id]]$omg2, 
    rep(p5.d1n1$dpar[[id]]$eta2, times=(dimd1[2]-qds1[2])) 
  ) 
)
tsig2 = p5.d1n1$dpar[[id]]$gam2 %*% tdel2 %*% t(p5.d1n1$dpar[[id]]$gam2)
tsig = tsig2 %x% tsig1
plot(sort( eigen(tsig)$values, decreasing=TRUE ), type="b")
#plot( sort(eigen(tsig1)$values, decreasing=TRUE), type="b")


tlam1 = p3.d2n3$dpar[[id]]$lam1; txi1 = p3.d2n3$dpar[[id]]$xi1
tsig1 = tlam1 %*% t(tlam1) + diag(txi1)
v1 = sort(eigen(tsig1)$values, decreasing=TRUE)
#plot( sort(eigen(tsig1)$values, decreasing=TRUE), type="b", 
#      ylab="eigenvalue", xlab="index")

tlam1 = p4.d2n3$dpar[[id]]$lam1; txi1 = p4.d2n3$dpar[[id]]$xi1
tsig1 = tlam1 %*% t(tlam1) + diag(txi1)
v2 = sort(eigen(tsig1)$values, decreasing=TRUE)


tlam1 = p5.d2n3$dpar[[id]]$lam1; txi1 = p5.d2n3$dpar[[id]]$xi1
tsig1 = tlam1 %*% t(tlam1) + diag(txi1)
v3 = sort(eigen(tsig1)$values, decreasing=TRUE)

matrix(c(v1,v2,v3), ncol=3) |> 
  matplot(type="b", pch=1, ylab="eigenvalue", xlab="index")

matplot(v3, type="b", pch=1, col="green",
        ylab="eigenvalue", xlab="index", bty="n")
matlines(v2, type="b", pch=1, col="red")
matlines(v1, type="b", pch=1, col="black")
legend("topright", legend = c("   0", "0.5", "   1"), col=1:3, pch=1,
       title=expression(paste("Value of ", a)))

strings = ls()


param.names = stringr::str_detect(strings, "^(p1|p2|p3|p4|p5)\\.")
parnames = strings[param.names]


pars = vector( "list", length(parnames) )

for(i in seq_along(parnames))
{
  pars[[i]] = get(parnames[i])
}

names(pars) = parnames


nobs.names = stringr::str_detect(strings, "^nset")
nobsnames = strings[nobs.names]

nobslist = vector("list", length(nobsnames) )

for( i in seq_along(nobsnames) )
{
  nobslist[[i]] = get(nobsnames[i])
}

names(nobslist) = nobsnames


parnames
save(pars, nobslist, parnames, nobsnames, file="params_d3h.RData")