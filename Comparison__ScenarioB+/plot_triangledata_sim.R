#Analyze tridata simulation
library(Rfast)
library(tibble)
library(magrittr)
library(ggthemes)

#Load data
load("sim_tridat_results.RData")

# Helpers -----------------------------------------------------------------


get_ari = function(alist=NULL, tclass=NULL)
{
  sq = seq_along(alist$estclass)
  val= vector("list", length=length(sq))
  for( i in sq )
  {
    val[[i]] = mixture::ARI(alist$estclass[[i]], tclass)
  }
  
  val
}


tclass = rep( c(1,2,3,4), each=250 )



# Get Results -------------------------------------------------------------

#Get ARI of each method
arilist = lapply( res, get_ari, tclass=tclass )

#get comp times
times = sapply( seq_along(res), FUN=function(i)
{
  res[[i]]$time  
}) 

#Add results of functional k-means
load("sim_tridat_kmeans.RData")

arikmeans = lapply( seq_along(res), FUN=function(u)
{
  mixture::ARI(res[[u]]$class, tclass)
})

timesk = sapply( seq_along(res), FUN=function(u)
{
  res[[u]]$time
}) |>
  (  \(u){ rbind(times,u) }  )() ->
  times

times |> as.numeric() |> matrix(nrow=4) |>  rowMeans(na.rm = TRUE)



for(i in seq_along(arilist))
{
  arilist[[i]][[length(arilist[[i]]) +1]] = arikmeans[[i]]
}




#Data
matrix( unlist(arilist), nrow=4 ) |> 
  t() |> 
  as_tibble() |>
  set_colnames( c("funHDDC", "MFSF", "Funclust", "L2FD kMeans") ) |> 
  pivot_longer(
    cols = everything(), names_to = "Algorithm", values_to = "ARI"
  ) ->
  plotdat

#Histograms
plotdat |> 
  ggplot( aes(x=ARI, group=Algorithm, fill=Algorithm) ) + 
  geom_histogram( color="black", bins=15 ) +
  facet_grid( cols=vars(Algorithm) ) +
  scale_fill_colorblind() +
  theme_hc()


#violins
plotdat |> 
  ggplot( aes(x=Algorithm, y=ARI, group=Algorithm, fill=Algorithm) ) + 
  geom_violin( adjust=1, scale="width" ) +
  #facet_wrap(vars(Algorithm)) +
  scale_fill_gdocs() +
  theme_hc()

#violins
plotdat |> 
  ggplot(aes(x=Algorithm, y=ARI, group=Algorithm, fill=Algorithm)) + 
  geom_violin(adjust=1, scale="width") +
  geom_jitter(height = 0, width = 0.03, alpha=0.3, cex=0.9) +
  #facet_wrap(vars(Algorithm)) +
  scale_fill_gdocs() +
  theme_hc()

#boxplots
plotdat |> 
  ggplot(aes(x=Algorithm, y=ARI, group=Algorithm, fill=Algorithm)) + 
  geom_boxplot(outlier.alpha = 0.3) +
  #facet_wrap(vars(Algorithm)) +
  scale_fill_gdocs() +
  theme_hc()



# Computation Times -------------------------------------------------------


times |> 
  t() |>
  as_tibble() |>
  set_colnames(c("funHDDC", "MFSF", "Funclust")) |> 
  pivot_longer(
    cols = everything(), 
    names_to = "Algorithm", 
    values_to = "Computation Time (mins)"
  ) |> 
  mutate("Computation Time (mins)" = as.numeric(`Computation Time (mins)`))->
  plotdat
  
#Histograms
plotdat |> 
  ggplot(aes(x=`Computation Time (mins)`, group=Algorithm, fill=Algorithm)) + 
  geom_histogram(
    color="black", binwidth = function(x) 10 * IQR(x) / (length(x)^(1/3))
  ) +
  facet_wrap( vars(Algorithm) ) +
  scale_fill_colorblind() +
  theme_hc()

#Histograms
plotdat |> 
  ggplot( aes(x=`Computation Time (mins)`, group=Algorithm, fill=Algorithm) ) + 
  geom_histogram(
    color="black", 
    alpha=0.7, 
    binwidth=function(x) 10 * IQR(x) / (length(x)^(1/3))
  ) +
  scale_fill_colorblind() +
  theme_hc()


#violins
plotdat |> 
  ggplot(
    aes(
      x=Algorithm, 
      y=`Computation Time (mins)`, 
      group=Algorithm, 
      fill=Algorithm
    )
  ) + 
  geom_violin( adjust=4, scale="width" ) +
  #facet_wrap(vars(Algorithm)) +
  scale_fill_gdocs() +
  theme_hc()

#violins
plotdat |> 
  ggplot(
    aes(
      x=Algorithm, 
      y=`Computation Time (mins)`, 
      group=Algorithm, 
      fill=Algorithm
    )
  ) + 
  geom_violin( adjust=1, scale="width" ) +
  geom_jitter(height = 0, width = 0.03, alpha=0.3, cex=0.9) +
  #facet_wrap(vars(Algorithm)) +
  scale_fill_gdocs() +
  theme_hc()

#boxplots
plotdat |> 
  ggplot(
    aes(
      x=Algorithm, 
      y=`Computation Time (mins)`, 
      group=Algorithm, 
      fill=Algorithm
    )
  ) + 
  geom_boxplot( outlier.alpha = 0.3 ) +
  #facet_wrap(vars(Algorithm)) +
  scale_fill_gdocs() +
  theme_hc()


#Mean times
sapply(
  seq_along(res), 
  FUN=function(i)
  {
    as.numeric(res[[i]]$time)  
  }
) |> 
  rowmeans()

#Mean/variance ARI
matrix( unlist(arilist), nrow=4 ) |> 
  rowmeans()

matrix( unlist(arilist), nrow=4 ) |> 
  rowVars(std=TRUE)

# //
