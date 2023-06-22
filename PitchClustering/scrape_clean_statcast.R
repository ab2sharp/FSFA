
# Baseball Data -----------------------------------------------------------
source("frankenclust.R")

#library(car)
#library(rgl)
#library(baseballr)
library(dbplyr)
library(knitr)
library(DBI)
#library(pracma)
#library(tidyverse)
library(magrittr)
library(devtools)


# Helpers -----------------------------------------------------------------

get_statcast = function( start_date, end_date ) 
{
  if (!is.character(start_date) | !is.character(end_date)) {
    stop("Please wrap your dates in quotations in 'yyyy-mm-dd' format.")
  }
  if (as.Date(start_date) <= "2015-03-01") {
    warning(
      "Some metrics such as Exit Velocity and Batted Ball Events 
       have only been compiled since 2015."
    )
  }
  if (as.Date(start_date) <= "2008-03-25") {
    stop("The data are limited to the 2008 MLB season and after.")
  }
  if (as.Date(start_date) > as.Date(end_date)) {
    stop("The start date is later than the end date.")
  }
  year <- substr(start_date, 1, 4)
  days <- seq.Date(as.Date(start_date), as.Date(end_date), 
                   by = "day")
  start_days <- as.character(days[(1:length(days))%%7 == 1])
  end_days <- as.character(days[(1:length(days))%%7 == 0])
  res <- list()
  n <- max(length(start_days), length(end_days))
  res <- foreach(i = 1:n) %do% {
    if (i == n) 
      end_days[i] <- end_date
    url <- paste0(
      "https://baseballsavant.mlb.com/statcast_search/csv?all=true&hfPT=&",
      "hfAB=&hfBBT=&hfPR=&hfZ=&stadium=&hfBBL=&hfNewZones=&",
      "hfGT=R%7CPO%7CS%7C&hfC&hfSea=", 
      year, 
      "%7C&hfSit=&hfOuts=&opponent=&pitcher_throws=&batter_stands=&",
      "hfSA=&player_type=pitcher&hfInfield=&team=&position=&hfOutfield=&hfRO=&",
      "home_road=&game_date_gt=", 
      start_days[i], 
      "&game_date_lt=", 
      end_days[i], 
      "&hfFlag=&hfPull=&metric_1=&hfInn=&min_pitches=0&min_results=0&",
      "group_by=name&sort_col=pitches&player_event_sort=h_launch_speed&",
      "sort_order=desc&min_abs=0&type=details"
    )
    suppressMessages(
      suppressWarnings(
        readr::read_csv(url, na = "null")
      )
    ) %>% 
  dplyr::select(
    game_year, game_date, game_pk, pitcher_name=player_name, inning, 
    inning_topbot, strikes, balls, outs_when_up, p_throws, pitch_number, 
    pitch_type,pitch_name, release_speed, release_pos_x, release_pos_y, 
    release_pos_z, plate_x, plate_z, vx0, vy0, vz0, ax, ay, az, launch_speed, 
    launch_angle, effective_speed, release_spin_rate, release_extension,
    launch_speed_angle, zone, type, at_bat_number, stand, events, description, 
    bb_type, hit_location, hc_x, hc_y, hit_distance_sc
)
  }
  res_data <- do.call("rbind", res) %>% 
  arrange(
    game_year,game_date, game_pk, inning, desc(inning_topbot), at_bat_number, 
    pitch_number
  ) %>% 
    as.data.frame()

  return(res_data)
}

#Format decimal value with no leading 0
numformat <- function(val)
{ 
  sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) 
}

#scale a vector of data
scale_v <- function(x)
{
  ( x - mean(x, na.rm=TRUE) ) / sd(x, na.rm=TRUE)
}

#Time to reach given y position, function
getTimeToReachYPosition <- function( vy0, ay, y0, y_p )
{
  return(  ( -vy0 - sqrt(vy0^2 - 2*ay*(y0-y_p)) ) / ay  )
}

#Velocity at given time in trajectory
velocityAtTimet <- Vectorize( function( v0, a, t )
{
  vt <-  v0  + a*t
  return(vt)
})

#Position at time t function 
position_at_t <- function( p_initial, velo_initial, acceleration, t, t_initial )
{
  return(
    p_initial + velo_initial*(t-t_initial)+(0.5)*acceleration*(t-t_initial)^2
  )
}

#3D vector length function
vectorNorm <- Vectorize( function ( x, y, z )
{
  return( sqrt(x^2 + y^2 + z^2) )
} )

writeStatcastToDB <- function(date_init, date_fin, table = "statcast")
{
  csv_data = get_statcast(date_init, date_fin)
  dbWriteTable(db$con, 
               table, 
               csv_data, 
               append = TRUE, 
               overwrite = FALSE, 
               row.names = FALSE)
  return("fin")
  
}

#Moody_Mudskippers table writing function
# Writes a new table in the database "data", with the name "name".
create <- function( data, name )
{
  DBI::dbSendQuery(
    data$src$con,
    paste( "CREATE TABLE", name, "AS", dbplyr::sql_render(data) )
  )
  dplyr::tbl(data$src$con,name)
}


# Explore Data ------------------------------------------------------------

#Create a new sqlite database
db <- src_sqlite("scDB.sqlite3", create = TRUE)

initdat = get_statcast("2019-03-20", "2019-03-20") 

#Write the data to the database you created
dbWriteTable(db$con, "statcast", initdat, overwrite=TRUE)

#Generate the dates to be used in the baseballr pitcher_all function. 

#Set the first interval to be 9 days, including the inital date
dateInitial <- as.Date("2019-03-21")
endFirstDateInterval <- dateInitial + 8
#The last day of the season
dateFinal <- as.Date("2019-09-29")


# Generate the end points for each date interval. 
# Append the last day to finish_dates 
# because the final interval is less than 9 days
begin_dates <- seq(dateInitial, dateFinal, "9 days") %>%  as.character() 
finish_dates <- seq(endFirstDateInterval, dateFinal, "9 days") %>%
  append(dateFinal) %>% 
  as.character()

mapply(writeStatcastToDB, 
       date_init = begin_dates, 
       date_fin = finish_dates, 
       table = "statcast")

#Use this code to create a connection to the database
db <-  src_sqlite("scDB.sqlite3", create = FALSE)


db%>%
  tbl("statcast")%>%
  colnames()

keep = c(
  "pitcher_name", "pitch_type", "zone", "release_speed", "release_pos_x", 
  "release_pos_y", "release_pos_z", "vx0", "vy0", "vz0", "ax", "ay", "az", 
  "release_extension", "effective_speed", "release_spin_rate", "pitch_number", 
  "stand", "p_throws", "plate_x", "plate_z", "description"
)


db %>%
  tbl("statcast")%>%
  collect() %>% 
  select(all_of(keep))%>%
  dbWriteTable(db$con, "ccast2", ., overwrite=TRUE)


db %>% 
  tbl("cleancast") %>% 
  pull(pitch_type) %>% 
  table()

ptypes = c("SL", "FF", "KC", "CH", "CU", "FS", "SI", "FT", "FC")

db%>%
  tbl("cleancast")%>%
  tally() #732473

db %>% 
  tbl("cleancast") %>% 
  filter(pitch_type %in% ptypes) %>% 
  distinct(description) %>% 
  kable()

db %>% 
  tbl("cleancast") %>% 
  filter(pitch_type %in% ptypes) %>% 
  filter(description == "blocked_ball") %>% 
  tally() #17421, lol...

db %>% 
  tbl("ccast2") %>% 
  filter(pitch_type %in% ptypes) %>% 
  collect()->
  ccast

#Check reported speeds
ccast %>% 
  mutate(
    time_release = getTimeToReachYPosition(vy0, ay, y0=50, y_p = release_pos_y),
    time_home = getTimeToReachYPosition(vy0, ay, y0=50, y_p = 17/12 ), 
    flight_time = time_home - time_release 
  )%>%
  mutate(
    vxR = velocityAtTimet(vx0, ax, t = time_release),  #The release velocities
    vyR = velocityAtTimet(vy0, ay, t = time_release), 
    vzR = velocityAtTimet(vz0, az, t = time_release)
  )%>%
  #Fix incorrect release speeds
  mutate(
    calc_speed_release = vectorNorm(vxR, vyR, vzR)*(15/22),#15/22: ft/sec to mph
    diff_speed = release_speed - calc_speed_release
  ) %>% 
  filter( abs(diff_speed) > 1 ) %>% 
  select( 
    pitcher_name, 
    release_speed, 
    calc_speed_release, 
    release_speed 
  ) %>% 
  View() #works!!
  

#Check reported plate positions
ccast %>% 
  mutate(
    time_release = getTimeToReachYPosition(vy0, ay, y0=50, y_p = release_pos_y),
    time_home = getTimeToReachYPosition(vy0, ay, y0=50, y_p = 17/12 ), 
    flight_time = time_home - time_release 
  )%>%
  mutate(
    vxR = velocityAtTimet(vx0, ax, t = time_release),  #The release velocities
    vyR = velocityAtTimet(vy0, ay, t = time_release), 
    vzR = velocityAtTimet(vz0, az, t = time_release)
  )%>%
  mutate(
    p_z = position_at_t(release_pos_z, vzR, az, time_home, time_release), 
    p_x = position_at_t(release_pos_x, vxR, ax, time_home, time_release),
  ) %>% 
  filter( abs(plate_z - p_z) >1 ) %>% #both p_z and p_x are correct
  View()


# Cleaning and Calculations -----------------------------------------------

ccast%>%
  #Remove pitches with uniterpretable pitch type / no Statcast measurments
  filter( pitch_type %in% ptypes,
          p_throws=="R",
          # #Remove poorly thrown pitches, which may skew the data
          # !description %in% c("blocked_ball", "hit_by_pitch") 
  )%>%  
 #time _release: gives value of t that returns pos at release from traj eqn 
 #time_home: value of t that returns the pos at home plate from traj eqn
 #flight_time: total time from release to home plate, in seconds
  mutate(
    time_release = getTimeToReachYPosition(vy0, ay, y0=50, y_p = release_pos_y),
    time_home = getTimeToReachYPosition(vy0, ay, y0=50, y_p = 17/12 ), 
    flight_time = time_home - time_release 
  )%>%
  mutate(
    vxr = velocityAtTimet(vx0, ax, t = time_release),  #The release velocities
    vyr = velocityAtTimet(vy0, ay, t = time_release), 
    vzr = velocityAtTimet(vz0, az, t = time_release)
  )%>%
  #Fix incorrect release speeds
  mutate(
    calc_speed_release = vectorNorm(vxr, vyr, vzr)*(15/22),#15/22: ft/sec to mph
  )%>%
  #Get plate positions according to the given trajectories.
  mutate(
    p_z = position_at_t(release_pos_z, vzr, az, time_home, time_release), 
    p_x = position_at_t(release_pos_x, vxr, ax, time_home, time_release),
    zone = if_else(
      abs(plate_x - p_x) > 1,   #Mirror the zones to be correct
      case_when(
        zone == 1 ~3,
        zone == 3 ~1,
        zone == 4 ~6,
        zone == 6 ~4,
        zone == 7 ~9,
        zone == 9 ~7,
        zone == 11 ~12,
        zone == 12 ~11,
        zone == 13 ~14,
        zone == 14 ~13
      ), 
      as.numeric(zone)
    )
  )%>%
  #Drop old plate values to keep plate values consistent with trajectories.
  select(-c(plate_x, plate_z))%>%
  dbWriteTable(db$con, "ccast2", ., overwrite = TRUE)


db %>% 
  tbl("ccast2") %>% 
  distinct(pitch_type)

#Mirror pitches so that every pitcher is right handed.
db %>% 
  tbl("ccast2") %>% 
  collect() %>% 
  #Mirrors all lefty pitches to be righty pitches
  mutate(release_pos_x = case_when(p_throws =="L" ~ -release_pos_x,
                                   p_throws =="R" ~  release_pos_x),
         #vxr = case_when(p_throws =="L" ~ -vxr,
         #                p_throws =="R" ~  vxr),
         zone = case_when(p_throws =="L" & zone == 1 ~ 3,
                          p_throws =="L" & zone == 3 ~ 1,
                          p_throws =="L" & zone == 4 ~ 6,
                          p_throws =="L" & zone == 6 ~ 4,
                          p_throws =="L" & zone == 7 ~ 9,
                          p_throws =="L" & zone == 9 ~ 7,
                          p_throws =="L" & zone == 11 ~ 12,
                          p_throws =="L" & zone == 12 ~ 11,
                          p_throws =="L" & zone == 13 ~ 14,
                          p_throws =="L" & zone == 14 ~ 13,
                          TRUE                        ~ zone),
         stand = case_when(p_throws =="L" & stand =="L" ~ "R",
                           p_throws =="L" & stand =="R" ~ "L",
                           TRUE                         ~ stand),
         #p_throws = "R" #This is a dumb idea. Change asap.
  )%>%
  dbWriteTable(db$con, "ccast2", ., overwrite = TRUE)

#Remove NAs in zone field
db %>% 
  tbl("ccast2") %>% 
  collect() %>% 
  filter(!is.na(zone)) %>% 
  dbWriteTable(db$con, "ccast2", ., overwrite = TRUE)


# Plotting ----------------------------------------------------------------

#Want to plot, for each pitch type, its percent of usage by zone

#Use this code to create a connection to the database
db <-  src_sqlite("scDB.sqlite3", create = FALSE)


# Vectors used to specify text locations inside the zones
# They are used in the mutate/case_when below
zones_x <- c(1,4,7)
zones_y <- c(1,2,3)
corners <- c(11,12,13,14)

# Used to plot the actual zones in the Statcast coordinates
# using ggplot
test_zones <- tibble(
  x_lower = rep(seq(-17/24, 17/72, length.out =3), each = 3),
  x_upper = rep(seq(-17/72,17/24, length.out =3), each = 3),
  y_lower = rep( seq(3/2,17/6, length.out = 3), 3),
  y_upper = rep( seq(13/6,7/2, length.out = 3), 3),
)


db%>%
  tbl("calccast")%>%
  select( pitch_type, p_x, p_z, zone )%>%
  collect()%>%
  group_by(pitch_type, zone)%>%
  #Calculate whiff percent for each pitch type, in each zone
  summarise(useage = n()) %>% 
  mutate(useage = useage/sum(useage)) %>% 
  #Defines the coordinates for plotting the whiff percentage 
  #in the center of zone to which it corresponds
  mutate( zone_x = case_when(
    zone %in% zones_x      ~ -17/36,
    zone %in% (zones_x+1)  ~  0,
    zone %in% (zones_x+2)  ~  17/36,
    zone == 11           ~ -34/36,                 
    zone == 12           ~  34/36, 
    zone == 13           ~ -34/36,
    zone == 14           ~  34/36,
  ),
  zone_y = case_when( 
    zone %in% zones_y    ~  19/6,
    zone %in% (zones_y+3)  ~  15/6,
    zone %in% (zones_y+6)  ~  11/6,
    zone == 11           ~  23/6,                 
    zone == 12           ~  23/6, 
    zone == 13           ~  7/6,
    zone == 14           ~  7/6,
  )
  )->
  sumuse

db%>%
  tbl("calccast")%>%
  select( pitch_type, p_x, p_z, zone )%>%
  collect()%>%
  group_by(pitch_type, zone)%>%
  ggplot() +
  #A heatmap of the pitch locations
  stat_density_2d( aes( x = p_x,
                        y = p_z,
                        fill = ..density..),
                   geom = "raster",
                   contour = FALSE,
                   show.legend = FALSE
  ) +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  coord_cartesian(xlim = c(-2, 2),
                  ylim = c(0.75,4)
  ) +
  # #This plots the percentages
  # geom_text(data = sumuse, 
  #           aes(x=zone_x, 
  #               y=zone_y, 
  #               label = paste0(round(100*useage), "%"),
  #               #label=useage,
  #               fontface = "bold")
  # ) +
  geom_rect(data = test_zones,
            aes( xmin = x_lower,
                 xmax = x_upper, 
                 ymin = y_lower,
                 ymax = y_upper), 
            fill = NA, 
            color = "grey20",
            lwd = 1
  ) +
  geom_segment(aes(x = 0, y = 3.5 , xend =  0, yend = 3.5 + (2/3))) +
  geom_segment(aes(x = 0, xend = 0, y = 0.8333, yend = 1.5)) +
  geom_segment(aes(x = 17/24, xend = (17/24)+(17/36), y = 2.5, yend = 2.5)) +
  geom_segment(aes(x = (-17/24 - 17/36), xend = -17/24, y = 2.5, yend = 2.5)) +
  ggtitle("Pitch Useage By Zone") +
  theme(legend.position = "none") +
  theme_dark() + 
  facet_wrap(~pitch_type, nrow = 3)


#decide which zones to keep for each pitch, then keep them, then group them
#by pitcher_name, then get averages, then get their function values. 


db%>%
  tbl("calccast")%>%
  select( pitcher_name, pitch_type, zone)%>%
  collect()%>%
  group_by(pitcher_name,pitch_type, zone)%>%
  #Calculate whiff percent for each pitch type, in each zone
  summarise(useage = n()) %>% 
  View()


db%>%
  tbl("calccast")%>%
  select( pitcher_name)%>%
  collect()%>%
  group_by(pitcher_name)%>%
  summarise(useage = n()) %>% 
  filter(useage >= 750) %>% 
  pull(pitcher_name)->
  keypees


db%>%
  tbl("calccast")%>%
  select( pitcher_name, pitch_type, zone)%>%
  filter(pitcher_name %in% keypees) %>% 
  collect()%>%
  group_by(pitcher_name,pitch_type, zone)%>%
  summarise(useage = n()) %>% 
  filter(useage >10) %>% 
  select(pitcher_name, pitch_type, zone) ->
  keypitches

db %>% 
  tbl("calccast") %>% 
  collect() ->
  pdat

keypitches[1,]
pdat %>% 
  filter(pitcher_name==as.character(keypitches[1,1]), 
         pitch_type==as.character(keypitches[1,2]), 
         zone==as.numeric(keypitches[1,3]) ) %>% 
  View()

#Foreach over this. 