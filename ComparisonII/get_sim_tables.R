#Get tables of the simulation results
library(Rfast)
library(magrittr)


# d7 Results --------------------------------------------------------------

load("results_d7.RData")

#Get ARI values out of results and put into a matrix
val = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})

#complexities
cxdf = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})


#d5 Results ---------------------------------------------------------------

load("results_d5_PartI.RData")
val1 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf1 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)


load("results_d5_PartII.RData")
val2 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf2 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)


load("results_d5_PartIII.RData")
val3 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf3 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)


load("results_d5_PartIV.RData")
val4 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf4 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)


load("results_d5_PartV.RData")
val5 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf5 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)


load("results_d5_PartX.RData")
val10 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf10 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)

list(val1, val2, val3, val4, val5, val10) |> 
  unlist() |> 
  array(c(3,20,6)) -> val

#Check to make sure its accurate
val[,9,1]
val1[,9]
val[,7,2]
val2[,7]
val[,11,6]
val10[,11]

#Calculate the means
val |> apply(1:2, mean) |> round(4) -> means5

means5 |> t() -> means5

list(cxdf1, cxdf2, cxdf3, cxdf4, cxdf5, cxdf10) |> 
  unlist() |> 
  array(c(2,20,6)) -> cxdf

cxdf[,9,1]
cxdf1[,9]
cxdf[,7,2]
cxdf2[,7]
cxdf[,11,6]
cxdf10[,11]

cxdf |> apply(1:2, mean) |>   round(1) |> t() -> cxdf

cbind(means5, cxdf) |> View()


#d3h Results --------------------------------------------------------------

load("results_d3h_PartI.RData")
val1 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf1 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)

load("results_d3h_PartII.RData")
val2 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf2 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)

load("results_d3h_PartIII.RData")
val3 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf3 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)

load("results_d3h_PartIV.RData")
val4 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf4 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)

load("results_d3h_PartX.RData")
val10 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    res[[i]][[j]]$ARI 
  }) |> 
    rowmeans()
})
cxdf10 = sapply(seq_along(res), FUN=function(i)
{
  sapply( seq_along(res[[i]]), FUN=function(j)
  {
    c( res[[i]][[j]]$Complexity[[1]],
       res[[i]][[j]]$Complexity[[2]]
    )
  }) |> 
    rowmeans()
})
rm(res, simsets)

list(val1, val2, val3, val4, val10) |> 
  unlist() |> 
  array(c(3,20,5)) -> val

#Check to make sure its accurate
val[,9,1]
val1[,9]
val[,7,2]
val2[,7]


#Calculate the means
val |> apply(1:2, mean, na.rm=TRUE) |> round(3) |> t() -> means3h
#means3h |> t() |> View()

list(cxdf1, cxdf2, cxdf3, cxdf4, cxdf10) |> 
  unlist() |> 
  array(c(2,20,5)) -> cxdf

cxdf[,9,1]
cxdf1[,9]
cxdf[,7,2]
cxdf2[,7]
cxdf[,11,5]
cxdf10[,11]

cxdf |> apply(1:2, mean, na.rm=TRUE) |> round(1) |> t() |> 
  (\(u){ cbind(means3h, u) })() -> tab3h



# //
