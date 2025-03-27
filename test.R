library(celltrackR)
library(trajr)
library(ggplot2)
#library(xlsx)
library(car)
library(rgl)
library(dplyr)
library(ggpubr)
library(pracma)
library(fractaldim)
library(here)
library(GGally)
library(tibble)
library(plyr)
library(dunn.test)
library(DescTools)
library(FactoMineR)
library(factoextra)
library(Hmisc)
library(gplots)
library(ggsignif)
library(stringi)
library(PMCMRplus)
#library(pca3d)
# DF2_Control_pas.8
traj_analysis_DF2_Control_pas.8 <- function(input) {
  
  data <- read.csv(input,  skip = 1, header = FALSE, dec = ".")
  
  traj_params <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                          c("track",
                            "length",
                            "distance",
                            "straight",
                            "square_displacement",
                            "mean_speed",
                            "sd_speed",
                            "max_speed",
                            "min_speed",
                            "sinuosity",
                            "emax",
                            "DC",
                            "SDDC",
                            "mean_angle",
                            "probe"))
  
  traj_params$probe <- as.character(traj_params$probe)
  
  for (i in unique(data$V1)) {
    
    # Define x, y, and time coordinates
    coords <- data.frame(x = data$V3[data$V1 == i],
                         y = data$V4[data$V1 == i],
                         # times = c(1:96))
                         timeCol = data$V2[data$V1 == i],
                         spatialUnits = "pixels", timeUnits = "hours")
    
    trj <- TrajFromCoords(coords, spatialUnits = "pixels", timeUnits = "seconds", fps = 96/3600/24) #fps frames/3600/hours = 1/120
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    # A 0.6579 object had length 1 pixels in the video, scale to micrometres
    trj <- TrajScale(trj, 0.6525 / 1, "micrometer") # 7.5 px = 10 microns
    TrajGetUnits(trj) # Returns the spatial units of a trajectory
    TrajGetTimeUnits(trj)	#Returns the temporal units of a trajectory
    TrajStepLengths(trj)	#Returns the lengths of each step within the trajectory
    # Rediscretization
    # The function TrajResampleTime linearly interpolates points along a trajectory
    # to create a new trajectory with fixed step time intervals.
    # trj <- TrajResampleTime(trj, 901)
    TrajDuration(trj) # Returns the temporal duration of the trajectory (or a portion)
    TrajGetFPS(trj)
    par(mar=c(5,5,5,5))
    # Plot it
    plot(trj, lwd = 2)
    points(trj, draw.start.pt = FALSE, pch = 21, col = "black", cex = 1.2)
    # Trajectory analysis
    # The TrajDerivatives function calculates linear speed and acceleration along a Trajectory
    derivs <- TrajDerivatives(trj)
    
    traj_params <- add_row(traj_params,
                           track = i,
                           # total length of the trajectory
                           length = TrajLength(trj),
                           # straight-line distance from the start to the end of the trajectory
                           distance = TrajDistance(trj),
                           # Straightness index
                           straight = TrajStraightness(trj), # D/L ratio
                           # expected square displacement of a correlated random walk
                           square_displacement = TrajExpectedSquareDisplacement(trj),
                           # Measures of speed
                           mean_speed = mean(derivs$speed),
                           sd_speed = sd(derivs$speed),
                           max_speed = max(derivs$speed),
                           min_speed = min(derivs$speed),
                           # Measures of straightness
                           sinuosity = TrajSinuosity2(trj),
                           emax = TrajEmax(trj),
                           SDDC  =  sd(TrajDirectionalChange(trj)),
                           DC = mean(TrajDirectionalChange(trj)),
                           #mean_angle = TrajMeanVectorOfTurningAngles(trj),
                           mean_angle = 0,
                           probe = 'DF2_Control_pas.8'
                           #probe = 1
    )
    
    head(traj_params)
    
  }
  
  # print(traj_params)
  write.csv(traj_params, file = 'DF2_Control_pas.8.csv')
  #tracks <<- traj_params
  return(traj_params) } # 97 frames ok

###

# Set names
mynames <- c("track",
             "length",
             "distance",
             "straight",
             "square_displacement",
             "mean_speed",
             "sd_speed",
             "max_speed",
             "min_speed",
             "sinuosity",
             "emax",
             "DC",
             "SDDC",
             "mean_angle",
             "probe")
# set constants
hundred = 100

# Set names
alltracks <- setNames(data.frame(matrix(ncol = 15, nrow = 0)), mynames)

alltracks$probe <- as.character(alltracks$probe)

#tracks <- alltracks
tracks_DF2_Control_pas.8 <- alltracks # 

# Choose dir ---------
file_list_DF2_Control_pas.8 <- list.files(path = , choose.dir(),
                                          pattern = "csv",
                                          all.files = FALSE,
                                          full.names = TRUE, recursive = TRUE,
                                          ignore.case = FALSE, include.dirs = FALSE,
                                          no.. = FALSE)

############# Call to function ##################
#################################################

# Start scan for Control
for (file_name in file_list_DF2_Control_pas.8) {
  tracks_DF2_Control_pas.8 <- rbind(tracks_DF2_Control_pas.8, 
                                    traj_analysis_DF2_Control_pas.8(file_name))}

# Merge all tracks-----------------------
alltracks <- rbind(tracks_DF2_Control_pas.8) # collect all tracks

# Set time to hours, remove tracks and mean_angle columns
all.h <- alltracks
head(all.h)
data <- cbind(all.h[,c(6,7,8,9)]*3600, all.h[,c(2,3,4,5,10,11,12,13,15)])
data$time <- "24h"
head(data)
tail(data)
write.csv(data, file = 'data_DF2_Control_pas.8.csv')




