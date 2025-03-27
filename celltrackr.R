# FetMSC movement trajectory analysis created by Dan Bobkov, 2021 # dan.bobkov@gmail.com
# Images of cell movement were obtained over a 24 h period using CQ1 confocal cytometer. 
# Coordinate values were extracted from images using the Manual Tracking plugin of the ImageJ.
# Images were collected from wells 1,4, (control), 2, 5 (LPA), 3, 6 (Y27632)
# If you use this script in your work, please cite: https://doi.org/10.1007/s11033-020-05476-6
# This script based on celltrackR package: Wortel, Inge MN, et al. "CelltrackR: an R package for 
# fast and flexible analysis of immune cell migration data." BioRxiv (2019): 670505.
##################################################
source(here::here('1_libraries.R'))
source(here::here('3_functions.R'))
# Read manuals
browseVignettes('celltrackR')
?TrackMeasures
#
# mean(sapply( t, speed ))
# sd(sapply( t, speed ))
# hist( sapply( t, speed ))

# Read cell tracks obtained from ImageJ and
# assign a unique number to each track



tracks <- list()
hundreds <- 1
rm(tmp)

for (file_name in file_list_control) {
  
  tmp <- read.csv(file_name, skip = 1, header = FALSE)[1:4]
  tmp$V1 <- tmp$V1 + hundreds * 100 # rename tracks
  tracks <- rbind(tracks, tmp)
  hundreds <- hundreds + 1
  
}

colnames(tracks)<- c("Track","Slice","X","Y")
tracks$Slice <- (tracks$Slice - 1)*.25 # time to hours: 15 min = 0.25 h
tracks$X <- tracks$X*1.3333 # xy calibration:         
tracks$Y <- tracks$Y*1.3333 # 
length(unique(tracks$Track))
write.csv(tracks, 'control.csv', row.names=FALSE) #
##
#
#


tracks <- list()
hundreds <- 1
rm(tmp)

for (file_name in file_list_inhib) {
  
  tmp <- read.csv(file_name, skip = 1, header = FALSE)[1:4]
  tmp$V1 <- tmp$V1 + hundreds * 100 # rename tracks
  tracks <- rbind(tracks, tmp)
  hundreds <- hundreds + 1
  
}

colnames(tracks)<- c("Track","Slice","X","Y")
tracks$Slice <- (tracks$Slice - 1)*.25 # time to hours: 15 min = 0.25 h
tracks$X <- tracks$X*1.3333 # xy calibration:         
tracks$Y <- tracks$Y*1.3333 # 
length(unique(tracks$Track))
write.csv(tracks, 'inhib.csv', row.names=FALSE) #
##
#



tracks <- list()
hundreds <- 1
rm(tmp)




#
# MSD
#
#
# Now, let's compare the MSD plot for the three different types of cells:
# Combine into a single dataframe with one column indicating the celltype
# To truly compare them, report subtrack length not in number of steps but
# in their duration (which may differ between different datasets)
#
t_control.msd <- aggregate( t_control, squareDisplacement, FUN = "mean.se" )
t_control.msd$cells <- "control"
t_control.msd$dt <- t_control.msd$i * timeStep( t_control )

t_inhib.msd <- aggregate( t_HAM_KC, squareDisplacement, FUN = "mean.se" )
t_HAM_KC.msd$cells <- "FC_KC"
t_HAM_KC.msd$dt <- t_HAM_KC.msd$i * timeStep( t_HAM_KC )

t_TCPS_KC.msd <- aggregate( t_TCPS_KC, squareDisplacement, FUN = "mean.se" )
t_TCPS_KC.msd$cells <- "FC_KC"
t_TCPS_KC.msd$dt <- t_TCPS_KC.msd$i * timeStep( t_TCPS_KC )





msddata8p <- rbind( t_FC_KC.msd, t_HAM_KC.msd, t_TCPS_KC.msd )
#msddata8p <- rbind( t_FC_KC.msd )
head(msddata8p)

# Plot
p8 <- ggplot( msddata8p, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (hours)") ),
        y = "MSD") +
  theme_classic() + ggtitle("DF2 cells at passage 8")



p8







































# Also make a zoomed version to look only at first part of the plot
#pzoom <- p1 + coord_cartesian( xlim = c(0,50), ylim = c(0,5000) ) 
#gridExtra::grid.arrange( p1, pzoom, ncol = 2 ) 
# Now, let's compare the MSD plot for the three different types of cells:
# Combine into a single dataframe with one column indicating the celltype
# To truly compare them, report subtrack length not in number of steps but
# in their duration (which may differ between different datasets)
# t_15p
t_15p.msd <- aggregate( t_15p, squareDisplacement, FUN = "mean.se" )
t_15p.msd$cells <- "Control"
t_15p.msd$dt <- t_15p.msd$i * timeStep( t_15p )


msddata15p <- rbind( t_15p.msd, t_8p.msd )
head(msddata15p)

# Plot
p15 <- ggplot( msddata15p, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (hours)") ),
        y = "MSD") +
  theme_classic() + ggtitle("FetMSC cells at passage 15")
p15
# Also make a zoomed version to look only at first part of the plot
#pzoom <- p1 + coord_cartesian( xlim = c(0,50), ylim = c(0,5000) ) 
#gridExtra::grid.arrange( p1, pzoom, ncol = 2 ) 

gridExtra::grid.arrange( p9, p15, nrow = 2 )

# Autocorrelation plots
# autocovariance

t_8p.acov <- aggregate( t_8p, overallDot, 
                         FUN = "mean.se" )
t_8p.acov$dt <- t_8p.acov$i * timeStep(t_8p)
t_8p.acov$cells <- "Pas.9 Control"
# Normalized autocovariance
t_8p.acov[,2:4] <- t_8p.acov[,2:4] / t_8p.acov$mean[1]



t_L_8p.acov <- aggregate( t_L_8p, overallDot, 
                          FUN = "mean.se" )
t_L_8p.acov$dt <- t_L_8p.acov$i * timeStep(t_L_8p)
t_L_8p.acov$cells <- "Pas.9 LPA"
# Normalized autocovariance
t_L_8p.acov[,2:4] <- t_L_8p.acov[,2:4] / t_L_8p.acov$mean[1]


t_Y_8p.acov <- aggregate( t_Y_8p, overallDot, 
                          FUN = "mean.se" )
t_Y_8p.acov$dt <- t_Y_8p.acov$i * timeStep(t_Y_8p)
t_Y_8p.acov$cells <- "Pas.9 Y27632"
# Normalized autocovariance
t_Y_8p.acov[,2:4] <- t_Y_8p.acov[,2:4] / t_Y_8p.acov$mean[1]

# Compare
acovdata_8p <- rbind( t_8p.acov, t_L_8p.acov, t_Y_8p.acov )

# compare
p9ac <- ggplot( acovdata_8p, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (hours)") ),
        y = "autocovariance" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

p9ac

p9zoom <- p9ac + coord_cartesian( xlim = c(0,5), ylim = c(0,.5) )


# Autocorrelation plots
# autocovariance

t_15p.acov <- aggregate( t_15p, overallDot, 
                           FUN = "mean.se" )
t_15p.acov$dt <- t_15p.acov$i * timeStep(t_15p)
t_15p.acov$cells <- "Pas.15 Control"
# Normalized autocovariance
t_15p.acov[,2:4] <- t_15p.acov[,2:4] / t_15p.acov$mean[1]



t_L_15p.acov <- aggregate( t_L_15p, overallDot, 
                           FUN = "mean.se" )
t_L_15p.acov$dt <- t_L_15p.acov$i * timeStep(t_L_15p)
t_L_15p.acov$cells <- "Pas.15 LPA"
# Normalized autocovariance
t_L_15p.acov[,2:4] <- t_L_15p.acov[,2:4] / t_L_15p.acov$mean[1]


t_Y_15p.acov <- aggregate( t_Y_15p, overallDot, 
                           FUN = "mean.se" )
t_Y_15p.acov$dt <- t_Y_15p.acov$i * timeStep(t_Y_15p)
t_Y_15p.acov$cells <- "Pas.15 Y27632"
# Normalized autocovariance
t_Y_15p.acov[,2:4] <- t_Y_15p.acov[,2:4] / t_Y_15p.acov$mean[1]

# Compare
acovdata_15p <- rbind( t_15p.acov, t_L_15p.acov, t_Y_15p.acov )

# compare
p15ac <- ggplot( acovdata_15p, aes( x = dt , y = mean, color = cells, fill = cells ) ) +
  geom_hline( yintercept = 0 ) +
  geom_ribbon( aes( ymin = lower, ymax = upper) , alpha = 0.2 ,color=NA ) +
  geom_line( ) +
  labs( x = expression( paste(Delta,"t (hours)") ),
        y = "autocovariance" ) +
  theme_classic() + 
  theme( axis.line.x = element_blank() )

p15ac

p15zoom <- p15ac + coord_cartesian( xlim = c(0,5), ylim = c(0,.5) )


gridExtra::grid.arrange(p9zoom, p15zoom, ncol = 2 )


#######
#######
#track.lengths <- sapply( t_8p, nrow ) - 1
#hist( track.lengths, xlab = "Track length (#steps)" )
# Detecting global directionality: hotellingsTest (Texor et al., 2011)
par( mfrow=c(1,2) )
plot( t_8p )
hotellingsTest( t_8p, plot = TRUE, col = "gray" )
#there seems to be a significant directionality in the  y
# direction (as shown by the blue ellipse).

# However, when we only consider steps that are some 
# distance apart using the step.spacing argument:
hotellingsTest( t_8p, plot = TRUE, 
                col = "gray", step.spacing = 3 ) # step.spacing = 5 
# there is no longer evidence for directionality. 
# (Note, however, that the power of the test is also 
# reduced because the larger step spacing reduces the total number of steps).


# Beltman et al (2009) proposed an analysis of angles 
# versus distance between cell pairs to detect global 
# directionality in a dataset. Use the function 
# analyzeCellPairs to get a dataframe with for each pair of tracks 
# in the dataset the angle (between their displacement vectors) 
# and the distance (min distance between the tracks at any timepoint).

library(pracma)
# Detecting double tracking: angle versus distance between cell pairs

df <- analyzeCellPairs( t_8p )

# label cellpairs that have both angle and distance below threshold
angle.thresh <- 90 # in degrees
dist.thresh <- 10 # this should be the expected cell radius
#df$id <- paste0( df$cell1,"-",df$cell2 )
#df$id[ !(df$angle < angle.thresh & df$dist < dist.thresh) ] <- ""
df$id[ !(df$angle < angle.thresh & df$dist < dist.thresh) ] <- "" 



# Plot
ggplot( df, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  geom_text( aes( label = id ), color = "red" ) +
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  geom_hline( yintercept = angle.thresh, col = "blue",lty=2 ) +
  geom_vline( xintercept = dist.thresh, col = "blue", lty=2) +
  theme_classic()


# Function to remove all data after given timepoint
# x must be a single track matrix, which is what this function will
# receive from lapply
removeAfterT <- function( x, time.cutoff ){
  
  # Filter out later timepoints
  x2 <- x[ x[,"t"] <= time.cutoff, ]
  
  # Return the new matrix, or NULL if there are no timepoints before the cutoff
  if( nrow(x2) == 0 ){
    return(NULL)
  } else {
    return(x2)
  }
}

# Call function on each track using lapply
filtered.t <- lapply( t_8p, function(x) removeAfterT( x, 6 ) )
# Remove any tracks where NULL was returned
filtered.t <- filtered.t[ !sapply( filtered.t, is.null )]
filtered.t <- as.tracks( filtered.t )
is.tracks( filtered.t )
#par(mfrow=c(1,2))
plot( t, main = "Unfiltered data")
plot( filtered.t, main = "Filtered on timepoints < 20" )


# The filtering function must return TRUE or FALSE for each track given to it
my.filter <- function(x){
  return( nrow(x) > 15 )
}

# Filter with this function using filterTracks
long.tracks <- filterTracks( my.filter, t )

# Plot the result
#par(mfrow=c(1,2))
plot( t, main = "All tracks")
plot( long.tracks, main = "Long tracks only" )

# The function selectTracks() selects tracks 
# based on upper and lower bounds of a certain measure. 
# For example, we can get the fastest half of the T cells:

# Filter with this function using filterTracks
median.speed <- median( sapply( filtered.t, speed ) )
fast.tracks <- selectTracks( filtered.t, speed, median.speed, Inf )

# Plot the result
par(mfrow=c(1,2))
plot( filtered.t, main = "All tracks")
plot( fast.tracks, main = "Fastest half" )


# Converted to dataframe
t.df <- as.data.frame(t)
str( t.df )

# Simple track visualization

#par( mfrow = c(2,2) )
plot( t_8p, main = "Control FetMSC cells" )
plot( t_L_8p, main = "FetMSC after 1 h LPA" )
plot( t_Y_8p, main = "FetMSC after 1 h Y27632" )
plot( normalizeTracks(t_Y_8p), main = "Normalizes Y tracks" )


#par( mfrow = c(2,2) )
plot( normalizeTracks(t_8p), main = "C" )
plot( normalizeTracks(t_L_8p), main = "L" )
plot( normalizeTracks(t_Y_8p), main = "Ys" )

# Now let's repeat for the other two celltypes and plot them next to each other for comparison

t_8p.speeds <- sapply( t_8p, speed )
t_L_8p.speeds <- sapply( t_L_8p, speed )
t_Y_8p.speeds <- sapply( t_Y_8p, speed )

# Create a dataframe of all data
dt_8p <- data.frame( cells = "t_8p", speed = t_8p.speeds )
dt_L_8p <- data.frame( cells = "t_L_8p", speed = t_L_8p.speeds )
dt_Y_8p <- data.frame( cells = "t_Y_8p", speed = t_Y_8p.speeds )
d <- rbind( dt_8p, dt_L_8p, dt_Y_8p )

# Compare:
ggplot( d, aes( x = cells, y = speed ) ) +
  ggbeeswarm::geom_quasirandom() +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


t_8p.step.speeds <- sapply( subtracks(t_8p,1), speed )
t_L_8p.step.speeds <- sapply( subtracks(t_L_8p,1), speed )
t_Y_8p.step.speeds <- sapply( subtracks(t_Y_8p,1), speed )

# Create a dataframe of all data
d_C_8p <- data.frame( cells = "t_8p", speed = t_8p.step.speeds )
d_L_8p <- data.frame( cells = "t_L_8p", speed = t_L_8p.step.speeds )
d_Y_8p <- data.frame( cells = "t_Y_8p", speed = t_Y_8p.step.speeds )
dstep <- rbind( d_C_8p, d_L_8p, d_Y_8p )

# Compare:
ggplot( dstep, aes( x = cells, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()


#Step-based analyses
#It has been suggested that cell-based analyses 
#introduce biases, and that analyses should be performed 
#on individual ?steps? from all tracks combined instead (Beltman et al (2009)).
t_8p.steps <- subtracks( t_8p, i = 1 )
# Obtain instantaneous speeds for each step using sapply
t_8p.step.speeds <- sapply( t_8pl.steps, speed )

summary( t_8p.step.speeds )
hist( t_8p.step.speeds )
# Compare cell-based versus step-based:
d.step <- data.frame( method = "step-based", speed = t_8p.step.speeds )
d.cell <- data.frame( method = "cell-based", speed = t_8p.speeds )
d.method <- rbind( d.step, d.cell )

ggplot( d.method, aes( x = method, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()

# And again compare between the different celltypes:

t_8p.step.speeds <- sapply( subtracks(t_8p,1), speed )
t_L_8p.step.speeds <- sapply( subtracks(t_L_8p,1), speed )

# Create a dataframe of all data
dT <- data.frame( cells = "t_8p", speed = t_8p.step.speeds )
dB <- data.frame( cells = "t_L_8p", speed = t_L_8p.step.speeds )
dN <- data.frame( cells = "t_Y_8p", speed = t_Y_8p.step.speeds )
dstep <- rbind( dT, dB, dN )

# Compare:
ggplot( dstep, aes( x = cells, y = speed ) ) +
  ggbeeswarm::geom_quasirandom( size = 0.3 ) +
  scale_y_continuous( limits = c(0,NA) ) +
  theme_classic()

# Comparing data with models
# get displacement vectors
step.displacements <- t( sapply( subtracks(t_8p,1), displacementVector ) )

# get mean and sd of displacement in each dimension
step.means <- apply( step.displacements, 2, mean )
step.sd <- apply( step.displacements, 2, sd )


#
# The Beauchemin model
beauchemin.tracks <- simulateTracks( 10, beaucheminTrack(sim.time=20) )
plot( beauchemin.tracks )

# Bootstrapping method for simulating migration
bootstrap.tracks <- simulateTracks( 10, bootstrapTrack( nsteps = 20, t_8p ) )
plot( bootstrap.tracks )



#
# Simulate more tracks to reduce noice
bootstrap.tracks <- simulateTracks( 100, bootstrapTrack( nsteps = 96, t_8p ) )

# Compare step speeds in real data to those in bootstrap data
real.speeds <- sapply( subtracks( t_8p,1 ), speed )
bootstrap.speeds <- sapply( subtracks( bootstrap.tracks,1), speed )
dspeed <- data.frame( tracks = c( rep( "data", length( real.speeds ) ),
                                  rep( "bootstrap", length( bootstrap.speeds ) ) ),
                      speed = c( real.speeds, bootstrap.speeds ) )

# Same for turning angles
real.angles <- sapply( subtracks( t_8p,2 ), overallAngle, degrees = TRUE )
bootstrap.angles <- sapply( subtracks( bootstrap.tracks,2), overallAngle, degrees = TRUE )
dangle <- data.frame( tracks = c( rep( "data", length( real.angles ) ),
                                  rep( "bootstrap", length( bootstrap.angles ) ) ),
                      angle = c( real.angles, bootstrap.angles ) )

# plot
pspeed <- ggplot( dspeed, aes( x = tracks, y = speed ) ) +
  geom_violin( color = NA, fill = "gray" ) +
  geom_boxplot( width = 0.3 ) +
  theme_classic()

pangle <- ggplot( dangle, aes( x = tracks, y = angle ) ) +
  geom_violin( color = NA, fill = "gray" ) +
  geom_boxplot( width = 0.3 ) +
  theme_classic()

gridExtra::grid.arrange( pspeed, pangle, ncol = 2 )



#Example: Comparing data with models
#Mean square displacement plot
# Simulate more tracks
brownian.tracks <- simulateTracks( 100, brownianTrack( nsteps = 96, dim = 3,
                                                       mean = step.means,
                                                       sd = step.sd ) )

bootstrap.tracks <- simulateTracks( 100, bootstrapTrack( nsteps = 96, t_8p ) )

msd.data <- aggregate( t_8p, squareDisplacement, FUN = "mean.se" )
msd.data$data <- "data"
msd.brownian <- aggregate( brownian.tracks, squareDisplacement, FUN = "mean.se" )
msd.brownian$data <- "brownian"
msd.bootstrap <- aggregate( bootstrap.tracks, squareDisplacement, FUN = "mean.se" )
msd.bootstrap$data <-"bootstrap"

msd <- rbind( msd.data, msd.brownian, msd.bootstrap )
ggplot( msd, aes( x = i, y = mean, 
                  ymin = lower, 
                  ymax = upper, 
                  color = data, 
                  fill = data ) ) +
  geom_ribbon( color= NA, alpha  = 0.2 ) +
  geom_line() +
  labs( x = "t_8p (steps)",
        y = "square displacement" ) +
  scale_x_continuous(limits= c(0, 95) ) +
  theme_bw()


#random walk model overestimate the MSD as observed in real data 
#This suggests there is less directional persistence in the real data 
#than those models assume.

# To check for directional persistence, we generate an autocovariance plot:
# compute autocorrelation
acor.data <- aggregate( t_8p, overallDot, FUN = "mean.se" )
acor.data$data <- "data"
acor.brownian <- aggregate( brownian.tracks, overallDot, FUN = "mean.se" )
acor.brownian$data <- "brownian"
acor.bootstrap <- aggregate( bootstrap.tracks, overallDot, FUN = "mean.se" )
acor.bootstrap$data <-"bootstrap"

acor <- rbind( acor.data, acor.brownian, acor.bootstrap )
ggplot( acor, aes( x = i, y = mean, 
                   ymin = lower, ymax = upper, 
                   color = data, fill = data ) ) +
  geom_ribbon( color= NA, alpha  = 0.2 ) +
  geom_line() +
  labs( x = "dt (steps)",
        y = "autocovariance" ) +
  scale_x_continuous(limits= c(0,20) ) +
  theme_bw()

#Indeed, the autocovariance drops same steeply for the real 
#data, which indicates that there is no directional persistence 
#in the cell data that is not captured by the brownian and bootstrapping models. 
#Thus, even when a model captures some aspects of the walk statistics in the data, 
#it may still behave differently in other respects.

# Combine them in a single dataset, where labels also indicate celltype:

T2 <- t_8p
names(T2) <- paste0( "C", names(T2) )
tlab <- rep( "C", length(T2) )

B2 <- t_L_8p
names(B2) <- paste0( "L", names(B2) )
blab <- rep( "L", length(B2) )

N2 <- t_Y_8p
names(N2) <- paste0( "Y", names(Neutrophils) )
nlab <- rep( "Y", length( N2) )

all.tracks <- c( T2, B2, N2 )
real.celltype <- c( tlab, blab, nlab )

# Extracting a feature matrix

m <- getFeatureMatrix( all.tracks, 
                       c(speed, meanTurningAngle, 
                         outreachRatio, squareDisplacement) )

# We get a matrix with a row per track and one column for each metric:
head(m)

plot( m, xlab = "speed", ylab = "mean turning angle" )

# PCA
pca <- trackFeatureMap( all.tracks, 
                        c(speed,meanTurningAngle,squareDisplacement,
                          maxDisplacement,outreachRatio ), method = "PCA", 
                        labels = real.celltype, return.mapping = TRUE )


pc1 <- pca[,1]
pc2 <- pca[,2]
track.speed <- sapply( all.tracks, speed )
cor.test( pc1, track.speed )

# multidimensional scaling (MDS)
trackFeatureMap( all.tracks,
                 c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                   outreachRatio ), method = "MDS",
                 labels = real.celltype )

# Uniform manifold approximate and projection (UMAP) 
trackFeatureMap( all.tracks,
                 c(speed,meanTurningAngle,squareDisplacement,
                   maxDisplacement,outreachRatio ), method = "UMAP",
                 labels = real.celltype )


# Hierarchical clustering
clusterTracks( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "hclust", labels = real.celltype )


# K-means clustering
clusterTracks( all.tracks,
               c(speed,meanTurningAngle,squareDisplacement,maxDisplacement,
                 outreachRatio ), method = "kmeans", 
               labels = real.celltype, centers = 3 )







#
#
# STAGGERED analysis
#
# compute metrics on all possible subtracks in a single track
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( applyStaggered( t_8p[[15]], speed, matrix = TRUE ),
                plot.title = title(main = "FetMSC, pas. 9\nControl"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 25)

filled.contour( applyStaggered( t_L_8p[[25]], speed, matrix = TRUE ),
                plot.title = title(main = "FetMSC, pas. 9\nLPA"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 25)

filled.contour( applyStaggered( t_Y_8p[[23]], speed, matrix = TRUE ),
                plot.title = title(main = "FetMSC, pas. 9\nY27632"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#nlevels = 25)

filled.contour( applyStaggered( t_15p[[44]], speed, matrix = TRUE ),
                plot.title = title(main = "FetMSC, pas. 15\nControl"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 24)

filled.contour( applyStaggered( t_L_15p[[74]], speed, matrix = TRUE ),
                plot.title = title(main = "FetMSC, pas. 15\nLPA"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 24)

filled.contour( applyStaggered( t_Y_15p[[45]], speed, matrix = TRUE ),
                plot.title = title(main = "FetMSC, pas. 15\nY27632"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 24)
#
# STAGGERED
#
# compute metrics on all possible subtracks in a track
# average step matrices (time consuming)
stag_C_8p <- matrix( rep( 0, len=9409), nrow = 97) 
for (i in 1:length(t_8p)) {
  stag_C_8p <- stag_C_8p + applyStaggered( t_8p[[i]], speed, matrix = TRUE )
}
mean_stag_C_8p <- stag_C_8p / length(t_8p)

#
stag_L_8p <- matrix( rep( 0, len=9409), nrow = 97) 
for (i in 1:length(t_L_8p)) {
  stag_L_8p <- stag_L_8p + applyStaggered( t_L_8p[[i]], speed, matrix = TRUE )
}
mean_stag_L_8p <- stag_L_8p / length(t_L_8p)

#
stag_Y_8p <- matrix( rep( 0, len=9409), nrow = 97)
for (i in 1:length(t_Y_8p)) {
  stag_Y_8p <- stag_Y_8p + applyStaggered( t_Y_8p[[i]], speed, matrix = TRUE )
}
mean_stag_Y_8p <- stag_Y_8p / length(t_Y_8p)

#
stag_C_15p <- matrix( rep( 0, len=9409), nrow = 97) 
for (i in 1:length(t_15p)) {
  stag_C_15p <- stag_C_15p + applyStaggered( t_15p[[i]], speed, matrix = TRUE )
}
mean_stag_C_15p <- stag_C_15p / length(t_15p)
#

stag_L_15p <- matrix( rep( 0, len=9409), nrow = 97) 
for (i in 1:length(t_L_15p)) {
  stag_L_15p <- stag_L_15p + applyStaggered( t_L_15p[[i]], speed, matrix = TRUE )
}
mean_stag_L_15p <- stag_L_15p / length(t_L_15p)
#
stag_Y_15p <- matrix( rep( 0, len=9409), nrow = 97) 
for (i in 1:length(t_Y_15p)) {
  stag_Y_15p <- stag_Y_15p + applyStaggered( t_Y_15p[[i]], speed, matrix = TRUE )
}
mean_stag_Y_15p <- stag_Y_15p / length(t_Y_15p)
#
# 
# Set Color Palettes
#library("viridis") # viridis(n), magma(n), inferno(n), plasma(n)
#library("ggsci") # pal_npg(), pal_aaas(), pal_lancet(), pal_jco()
#palette <- colorRampPalette(c("darkblue", "blue", "lightblue1",
#                              "green", "yellow", "red", "darkred"))
#palette <- colorRampPalette(c('darkred','red','blue','lightblue'))
#palette <- colorRampPalette(c("blue","yellow","red"))
#palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
#palette <- viridis(35)
#cols = palette
#cols = plasma(24)
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_C_8p,
                plot.title = title(main = "FetMSC, pas. 9\nControl"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 25)
                

cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_L_8p,
                plot.title = title(main = "FetMSC, pas. 9\nLPA"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 20)
                

cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_Y_8p,
                plot.title = title(main = "FetMSC, pas. 9\nY27632"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 20)
                


# Passage 15
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_C_15p,
                plot.title = title(main = "FetMSC, pas. 15\nControl"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 20)
                

cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_L_15p,
                plot.title = title(main = "FetMSC, pas. 15\nLPA"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 20)

cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_Y_15p,
                plot.title = title(main = "FetMSC, pas. 15\nY27632"),
                #xlab = "X",ylab = "Y"),
                key.title = title(main="Speed\n(μm/h)"),
                col = cols)#,nlevels = 20)

#
# ROSE
# PLOT
# 
plot(t_8p[1:10], main = "DF2 at passage 9\n(N = 10)")
# Overlay starting points to view directionality
plot( normalizeTracks(t_8p[1:10]), 
      main = "DF2 at passage 9\n Normalized data (N = 10)" )


# Function to remove all data after given timepoint
# x must be a single track matrix, which is what this function will
# receive from lapply
removeAfterT <- function( x, time.cutoff ){
  
  # Filter out later timepoints
  x2 <- x[ x[,"t"] <= time.cutoff, ]
  
  # Return the new matrix, or NULL if there are no timepoints before the cutoff
  if( nrow(x2) == 0 ){
    return(NULL)
  } else {
    return(x2)
  }
}

# Call function on each track using lapply
filtered.t <- lapply( t_8p, function(x) removeAfterT( x, 6 ) )

# Remove any tracks where NULL was returned
filtered.t <- filtered.t[ !sapply( filtered.t, is.null )]
filtered.t <- as.tracks( filtered.t )
is.tracks( filtered.t )


stag_8p_f6h <- matrix(rep(0, len=25^2), nrow = 25)

system.time(
  for (i in 1:length(filtered.t)) {
    stag_8p_f6h <- stag_8p_f6h + applyStaggered(filtered.t[[i]], speed, matrix = TRUE)
  }
)


### CLEAN TRACK 9

# Dropping all tracks with gaps
split.gap <- repairGaps( t_8p, how = "split" )
# Detecting double tracking: angle versus distance between cell pairs
df <- analyzeCellPairs( t_8p )
# label cellpairs that have both angle and distance below threshold
angle.thresh <- 9 # in degrees
dist.thresh <- 10 # this should be the expected cell radius
df$id <- paste0( df$cell1,"-",df$cell2 )
df$id[ !(df$angle < angle.thresh & df$dist < dist.thresh) ] <- "" 
# Plot
ggplot( df, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  geom_text( aes( label = id ), color = "red" ) +
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  geom_hline( yintercept = angle.thresh, col = "blue",lty=2 ) +
  geom_vline( xintercept = dist.thresh, col = "blue", lty=2) +
  theme_classic()
# Plot these tracks specifically to check:
plot( t_8p[c("2201","2209")], main = "Duplicates")
# Drop duplicates
t_8p$`2201` <- NULL

#

### CLEAN TRACK 36

# Dropping all tracks with gaps
split.gap <- repairGaps( t_36p, how = "split" )
# Detecting double tracking: angle versus distance between cell pairs
df <- analyzeCellPairs( t_36p )
# label cellpairs that have both angle and distance below threshold
angle.thresh <- 9 # in degrees
dist.thresh <- 9 # this should be the expected cell radius
df$id <- paste0( df$cell1,"-",df$cell2 )
df$id[ !(df$angle < angle.thresh & df$dist < dist.thresh) ] <- "" 
# Plot
ggplot( df, aes( x = dist, y = angle ) ) +
  geom_point( color = "gray40" ) +
  geom_text( aes( label = id ), color = "red" ) +
  labs( x = "distance between cell pairs",
        y = "angle between cell pairs" ) +
  geom_hline( yintercept = angle.thresh, col = "blue",lty=2 ) +
  geom_vline( xintercept = dist.thresh, col = "blue", lty=2) +
  theme_classic()
# Plot these tracks specifically to check:
plot( t_36p[c("1608","1616")], main = "Duplicates")
plot( t_36p[c("1011","1019")], main = "Duplicates")
# Drop duplicates
t_36p$`1608` <- NULL
t_36p$`1505` <- NULL
t_36p$`1011` <- NULL
#


# STAGGERED analysis
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))

# compute metrics on all possible subtracks in a single track
filled.contour( applyStaggered( t_8p[[19]], meanTurningAngle, matrix = TRUE ),
                plot.title = title(main = "DF2, pas. 15\nTrack #15"),
                key.title = title(main="Mean\nturning\nangle"),
                col = cols)

# compute metrics on all possible subtracks in a single track
filled.contour( applyStaggered( t_38p[[100]], speed, matrix = TRUE ),
                plot.title = title(main = "DF2, pas. 38\nTrack #100"),
                key.title = title(main="Speed\n(??m/h)"),
                col = cols)

# trackLength(x)   duration(x)   speed(x)   displacementRatio(x)
# displacement(x, from = 1, to = nrow(x)) outreachRatio(x)
# squareDisplacement(x, from = 1, to = nrow(x)) 
# displacementVector(x) asphericity(x) straightness(x)
# maxDisplacement(x) hurstExponent(x) fractalDimension(x)




























# STAGGERED analysis

# trackLength(x)   duration(x)   speed(x)   displacementRatio(x)
# displacement(x, from = 1, to = nrow(x)) outreachRatio(x)
# squareDisplacement(x, from = 1, to = nrow(x)) 
# displacementVector(x) asphericity(x) straightness(x)
# maxDisplacement(x) hurstExponent(x) fractalDimension(x)

# compute metrics on all possible subtracks in all tracks
# (time consuming)
stag_8p_C <- matrix(rep( 0, len=97^2), nrow = 97)

for (i in 1:length(t_8p_C)) {
  stag_8p_C <- stag_8p_C + applyStaggered( t_8p_C[[i]], speed, matrix = TRUE )
}

mean_stag_8p_C <- stag_8p_C / length(t_8p_C)
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
filled.contour( mean_stag_8p_C,
                plot.title = title(main = "DF2 at passage 9\n24 h, Control"),
                key.title = title(main="Speed"),
                col = cols)







# DF2 Passage 8 # 96 frames # 1 px = 0.6666667 micrometers
file_list_8p_L <- list.files(path = choose.dir(), pattern = "csv",
                             all.files = FALSE, full.names = TRUE,
                             recursive = TRUE, ignore.case = FALSE,
                             include.dirs = FALSE, no.. = FALSE)


tracks <- list()
hundreds <- 1
rm(tmp)

for (file_name in file_list_8p_L) {
  
  tmp <- read.csv(file_name, skip = 1, header = FALSE)[1:4]
  tmp$V1 <- tmp$V1 + hundreds * 100 # rename tracks
  tracks <- rbind(tracks, tmp)
  hundreds <- hundreds + 1
  
}

colnames(tracks)<- c("Track","Slice","X","Y")
tracks$Slice <- (tracks$Slice - 1)*.25 # time to hours: 15 min = 0.25 h
tracks$X <- tracks$X*0.6666667 # xy calibration:         
tracks$Y <- tracks$Y*0.6666667 # 0.6666667 mkm = 1 px
length(unique(tracks$Track))
write.csv(tracks, 'DF2_p8_L.csv', row.names=FALSE) # N=100
# READ
# TRACKS
t_8p_L <- read.tracks.csv('DF2_p8_L.csv',
                          id.column = 1,
                          time.column = 2,
                          pos.columns = c(3, 4),
                          scale.t = 1, # a value by which to multiply each time point
                          scale.pos = 1, # a value, by which to multiply each spatial position
                          header = TRUE,
                          sep = ",",
                          track.sep.blankline = FALSE)

# STAGGERED analysis

# trackLength(x)   duration(x)   speed(x)   displacementRatio(x)
# displacement(x, from = 1, to = nrow(x)) outreachRatio(x)
# squareDisplacement(x, from = 1, to = nrow(x)) 
# displacementVector(x) asphericity(x) straightness(x)
# maxDisplacement(x) hurstExponent(x) fractalDimension(x)

# compute metrics on all possible subtracks in all tracks
# (time consuming)
stag_8p_L <- matrix(rep( 0, len=97^2), nrow = 97)
for (i in 1:length(t_8p_L)) {
  stag_8p_L <- stag_8p_L + applyStaggered( t_8p_L[[i]], speed, matrix = TRUE )
}
mean_stag_8p_L <- stag_8p_L / length(t_8p_L)
cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
stag8l <-  filled.contour( mean_stag_8p_L,
                           plot.title = title(main = "DF2 at passage 8\n24 h, LPA"),
                           key.title = title(main="Speed"),
                           col = cols)
stag8l


# DF2 Passage 8 # 96 frames # 1 px = 0.6666667 micrometers
file_list_8p_Y <- list.files(path = choose.dir(), pattern = "csv",
                             all.files = FALSE, full.names = TRUE,
                             recursive = TRUE, ignore.case = FALSE,
                             include.dirs = FALSE, no.. = FALSE)


tracks <- list()
hundreds <- 1
rm(tmp)

for (file_name in file_list_8p_Y) {
  
  tmp <- read.csv(file_name, skip = 1, header = FALSE)[1:4]
  tmp$V1 <- tmp$V1 + hundreds * 100 # rename tracks
  tracks <- rbind(tracks, tmp)
  hundreds <- hundreds + 1
  
}

colnames(tracks)<- c("Track","Slice","X","Y")
tracks$Slice <- (tracks$Slice - 1)*.25 # time to hours: 15 min = 0.25 h
tracks$X <- tracks$X*0.6666667 # xy calibration:         
tracks$Y <- tracks$Y*0.6666667 # 0.6666667 mkm = 1 px
length(unique(tracks$Track))
write.csv(tracks, 'DF2_p8_Y.csv', row.names=FALSE) # N=100
# READ
# TRACKS
t_8p_Y <- read.tracks.csv('DF2_p8_Y.csv',
                          id.column = 1,
                          time.column = 2,
                          pos.columns = c(3, 4),
                          scale.t = 1, # a value by which to multiply each time point
                          scale.pos = 1, # a value, by which to multiply each spatial position
                          header = TRUE,
                          sep = ",",
                          track.sep.blankline = FALSE)

# STAGGERED analysis
# compute metrics on all possible subtracks in all tracks
# (time consuming)
stag_8p_Y <- matrix(rep( 0, len=97^2), nrow = 97)

for (i in 1:length(t_8p_Y)) {
  stag_8p_Y <- stag_8p_Y + applyStaggered( t_8p_Y[[i]], speed, matrix = TRUE )
}

mean_stag_8p_Y <- stag_8p_Y / length(t_8p_Y)

cols = rev(colorRampPalette(c('darkred','red','blue','lightblue'))(30))
stag8y <- filled.contour( mean_stag_8p_Y,
                          plot.title = title(main = "DF2 pas. 8\n24 h, Y27632"),
                          key.title = title(main="Speed"),
                          col = cols)






mean_stag_8p <- (mean_stag_8p_C + mean_stag_8p_L + mean_stag_8p_Y)/3

