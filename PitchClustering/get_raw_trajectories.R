setwd("G:/ResearchDocuments/Research/FunMatNorm/Analysis_Pitches")
library("tidyverse")
source("frankenclust.R")

#library(car)
#library(rgl)
#library(baseballr)
library(dbplyr)
library(knitr)
library(DBI)
library(RSQLite)
#library(pracma)
#library(tidyverse)
library(magrittr)
library(devtools)



# Helpers -----------------------------------------------------------------

#Time to reach given y position function
getTimeToReachYPosition <- function(vy0, ay, y0, y_p)
{
  return((-vy0 - sqrt(vy0^2 - 2*ay*(y0-y_p)))/ay)
}




# Get pitch data ----------------------------------------------------------


#Use this code to create a connection to the database
db <-  src_sqlite("scDB.sqlite3", create = FALSE)

#Get pitcher/pitch combos with more than 100 observations
#   Creates a tibble with two columns:
#   1: pitcher name
#   2: Name of pitch with 100+ uses
db%>%
  tbl("ccast2")%>%
  dplyr::select( pitcher_name, pitch_type)%>%
  collect()%>%
  group_by(pitcher_name, pitch_type)%>%
  summarise(useage = n()) %>%
  filter(useage >= 100) %>% 
  dplyr::select(pitcher_name, pitch_type)->
  keypees

nrow( keypees )

#Import the data into the current environment 
# CAUTION: this may be large
db %>% 
  tbl("ccast2") %>% 
  collect() ->
  pdat

#A list of dataframes, which we will bind together
lpxdat = foreach(i=seq_len(nrow(keypees))) %do% {
  
  pdat %>% 
    filter(pitcher_name==as.character(keypees[i,1]), 
           pitch_type==as.character(keypees[i,2])
    )->
    val
  return(val)
  
}
pxdat = bind_rows(lpxdat)


#Readability
pxdat %>% 
  rename(rpx = release_pos_x, 
         rpy = release_pos_y,
         rpz = release_pos_z)->
  pxdat


pxdat %>% 
  select(
    pitcher_name, pitch_type, zone, 
    rpx, rpy, rpz,
    vxr, vyr, vzr, 
    ax, ay, az
  ) %>% 
  group_by(pitcher_name, pitch_type, zone) %>% 
  summarise(across(everything(), mean)) %>% 
  ungroup() ->
  pmeans

# View(pmeans)

#Remove those without atleast one pitch in each zone
pmeans %>%
  group_by(pitcher_name, pitch_type) %>%
  summarise(zones = n()) %>%
  filter(zones==13) %>%
  select(pitcher_name, pitch_type) ->
  fullzones

lmdat = foreach(i=seq_len(nrow(fullzones))) %do% {
  pmeans %>%
    filter(pitcher_name==as.character(fullzones[i,1]),
           pitch_type==as.character(fullzones[i,2])
    )->
    val
  return(val)
}

pmdat = bind_rows(lmdat)

#Get new set of key-value pairs for pitchers
pmdat %>% 
  select(pitcher_name, pitch_type) %>% 
  unique() ->
  keypees

#Get value for time at home plate for each pitch
pmdat %>% 
  mutate(time_50i = getTimeToReachYPosition(vyr, ay, y0=rpy, y_p = 50 ),
         time_50f = getTimeToReachYPosition(vyr, ay, y0=rpy, y_p = -50 )
  )->
  pmdat


get_traj <- function(xr=NULL, yr=NULL, zr=NULL, 
                     vxr=NULL,vyr=NULL, vzr=NULL, 
                     ax=NULL, ay=NULL, az=NULL,
                     ti=NULL, tf=1)
{
  if(is.na(xr) || is.na(tf)){ val = rep(0,times = 3*1000) 
  }else{
    interval <- seq(ti, tf, length.out = 1000)
    val = matrix( c(xr + vxr*(interval) + (0.5)*ax*(interval)^2,
                    yr + vyr*(interval) + (0.5)*ay*(interval)^2,
                    zr + vzr*(interval) + (0.5)*az*(interval)^2),
                  nrow=3
    )
  }
  return(val)
}


zones = 1:14

lfdat = foreach(i = seq_len(nrow(keypees))) %do% {
  pmdat %>% 
    filter(pitcher_name == as.character(keypees[i,1]),
           pitch_type == as.character(keypees[i,2])
    ) %>% 
    select( -matches(c("pitcher_name", "pitch_type", 
                       "zone", "time_plate", "time_zeroy"))
    ) %>% 
    apply(1, FUN = function(r)
    {
      get_traj( xr=r[1],  yr=r[2],  zr=r[3],
                vxr=r[4], vyr=r[5], vzr=r[6],
                ax=r[7],  ay=r[8],  az=r[9],
                ti=r[10])
    })->
    val
  
  return(val)
  
} #correct

fdat = do.call(cbind,lfdat)

farray = array(fdat, dim=c(1000,39,nrow(keypees)))

fary = aperm(farray, c(1,3,2)) #ready for fda

save(fary, keypees, file="pitchtrajs_uniformtime.RData")

#
