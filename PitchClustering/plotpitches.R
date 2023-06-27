#
# Plot the pitch functions
#
# This file produces the plot of Jacob DeGrom's fastball used in the paper.
#

library("fda")
library("rgl")
library(purrr)
library(Rfast)
library(mixture)
#Data
load("p2c_nonoise.RData")

#


#Creat plotable data
tt= c(0,seq_len(999))/999
fdat = eval.fd(pfd, tt)


#plotting function
plotpitchh = Vectorize(function(ind=NULL, leng=NULL)
{
  palette <- colorRampPalette(
    c("darkviolet", "dodgerblue", "gold", "indianred")
  )
  
  colors = palette(5)
  npnts = seq_len(leng)
  lines3d(
    fdat[npnts,ind,1:3], 
    color=colors[5], lwd=2.3*lwdth, alpha=apval, shininess=0
  )   #zone 1
  lines3d(
    fdat[npnts,ind,4:6], 
    color=colors[5], lwd=2.3*lwdth, alpha=apval, shininess=0
  )   #2
  lines3d(
    fdat[npnts,ind,7:9], 
    color=colors[5], lwd=2.3*lwdth, alpha=apval, shininess=0
  )   #3
  lines3d(fdat[npnts,ind,10:12], color=colors[4], lwd=1.5*lwdth, alpha=apval) #4
  lines3d(fdat[npnts,ind,13:15], color=colors[4], lwd=1.5*lwdth, alpha=apval) #5
  lines3d(fdat[npnts,ind,16:18], color=colors[4], lwd=1.5*lwdth, alpha=apval) #6
  lines3d(fdat[npnts,ind,19:21], color=colors[2], lwd=1.3*lwdth, alpha=apval) #7
  lines3d(fdat[npnts,ind,22:24], color=colors[2], lwd=1.3*lwdth, alpha=apval) #8
  lines3d(fdat[npnts,ind,25:27], color=colors[2], lwd=1.3*lwdth, alpha=apval) #9
  lines3d(fdat[npnts,ind,28:30], color=colors[3], lwd=lwdth, alpha=apval) #11
  lines3d(fdat[npnts,ind,31:33], color=colors[3], lwd=lwdth, alpha=apval) #12
  lines3d(
    fdat[npnts,ind,34:36], 
    color=colors[1], lwd=lwdth, alpha=apval, shininess=0
  ) #13
  lines3d(
    fdat[npnts,ind,37:39], 
    color=colors[1], lwd=lwdth, alpha=apval, shininess=0
  ) #14
})

plotpitchl = Vectorize(function(ind=NULL, leng=NULL)
{
  palette <- colorRampPalette(
    c("darkviolet", "dodgerblue", "gold", "indianred")
  )
  
  colors = palette(5)
  npnts = seq_len(leng)
  lines3d(
    fdat[npnts,ind,1:3], 
    color=colors[5], lwd=2.3*lwdth, alpha=apval, shininess=0
    )   #zone 1
  lines3d(
    fdat[npnts,ind,4:6], 
    color=colors[5], lwd=2.3*lwdth, alpha=apval, shininess=0
    )   #2
  lines3d(
    fdat[npnts,ind,7:9], 
    color=colors[5], lwd=2.3*lwdth, alpha=apval, shininess=0
    )   #3
  lines3d(fdat[npnts,ind,10:12], color=colors[4], lwd=2.3*lwdth, alpha=apval) #4
  lines3d(fdat[npnts,ind,13:15], color=colors[4], lwd=2.3*lwdth, alpha=apval) #5
  lines3d(fdat[npnts,ind,16:18], color=colors[4], lwd=2.3*lwdth, alpha=apval) #6
  lines3d(fdat[npnts,ind,19:21], color=colors[2], lwd=2.3*lwdth, alpha=apval) #7
  lines3d(fdat[npnts,ind,22:24], color=colors[2], lwd=2.3*lwdth, alpha=apval) #8
  lines3d(fdat[npnts,ind,25:27], color=colors[2], lwd=2.3*lwdth, alpha=apval) #9
  lines3d(fdat[npnts,ind,28:30], color=colors[1], lwd=2.3*lwdth, alpha=apval)#11
  lines3d(fdat[npnts,ind,31:33], color=colors[1], lwd=2.3*lwdth, alpha=apval)#12
  lines3d(
    fdat[npnts,ind,34:36], 
    color=colors[3], lwd=2.3*lwdth, alpha=apval, shininess=0
  ) #13
  lines3d(
    fdat[npnts,ind,37:39], 
    color=colors[3], lwd=2.3*lwdth, alpha=apval, shininess=0
  ) #14
})
#Chosen pitch
indpitch = 241 #deGrom FF

#Plot parameters
lwdth = 1.5
apval = 1
open3d()
plotpitchl(indpitch, leng=380)

#Add strikezone
aval = 1
lwdzone = 2.5
axes3d()
rgl.lines(c(-17/24, 17/24), c(17/12, 17/12), c(1.5, 1.5), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(-17/24, -17/24), c(17/12, 17/12), c(1.5, 3.5), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(-17/24, 17/24), c(17/12, 17/12), c(3.5, 3.5), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(17/24, 17/24), c(17/12, 17/12) ,c(1.5, 3.5), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(-17/24 + 34/72,-17/24 + 34/72),c(17/12, 17/12), c(1.5, 3.5), 
          color="black", alpha = aval, lwd=lwdzone)
rgl.lines(c(-17/24 + 2*(34/72),-17/24 + 2*(34/72)),c(17/12, 17/12), c(1.5, 3.5),
          color="black", alpha = aval, lwd=lwdzone )
rgl.lines(c(-17/24, 17/24), c(17/12, 17/12) ,c(1.5 + (2/3), 1.5 + (2/3)), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(-17/24, 17/24), c(17/12, 17/12) ,c(1.5 + (4/3), 1.5 + (4/3)), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(-1-17/24, -17/24), c(17/12, 17/12) ,c(10/4, 10/4), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(17/24, 1+17/24), c(17/12, 17/12) ,c(10/4, 10/4), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(0, 0), c(17/12, 17/12) ,c(3.5, 4.5), 
          color = "black", alpha = aval, lwd=lwdzone)
rgl.lines(c(0, 0), c(17/12, 17/12) ,c(0.5, 1.5), 
          color = "black", alpha = aval, lwd=lwdzone)
quads3d(c(-17/24, -17/24, 17/24, 17/24),
          c(17/12, 17/12, 17/12, 17/12),
          c(1.5, 3.5, 3.5, 1.5), lit=FALSE, shininess=0, color="gray")
aspect3d(1,1.5,1)

# //