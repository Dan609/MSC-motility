# MSC cells movement trajectory analysis created by Dan Bobkov, 2024 # dan.bobkov@gmail.com
# Images of cell movement were obtained over a 24 h period using CQ1 confocal cytometer. 
# Coordinate values were extracted from images using the Manual Tracking plugin of the ImageJ.
# Images were collected from wells 
# If you use this script in your work, please cite https://doi.org/10.1007/s11033-020-05476-6
# This script based on trajr package
# If you use trajr in your work, please cite https://doi.org/10.1111/eth.12739.
##################################################
# par( mfrow = c(1,1) )
# Load libraries
# library(trajr)
#source(here::here('1_libraries.R'))
# Load functions 
#source(here::here('3_functions.R'))
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


#alltracks_p12C <- setNames(data.frame(matrix(ncol = 15, nrow = 0)), mynames)
#alltracks_p12L <- setNames(data.frame(matrix(ncol = 15, nrow = 0)), mynames)
#alltracks_p12Y <- setNames(data.frame(matrix(ncol = 15, nrow = 0)), mynames)


#alltracks_p12C$probe <- as.character(alltracks_p12C$probe)
#alltracks_p12L$probe <- as.character(alltracks_p12L$probe)
#alltracks_p12Y$probe <- as.character(alltracks_p12Y$probe)

alltracks$probe <- as.character(alltracks$probe)

#tracks <- alltracks
tracks_DF2_Control_pas.8 <- alltracks # 
#tracks_p12L <- alltracks # 
#tracks_p12Y <- alltracks # 



#################################################
# Choose dir function for mac <<<<<<<<<<<<<<<<<<
choose.dir <- function() {
  system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
         intern = FALSE, ignore.stderr = TRUE)
  p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
  return(ifelse(length(p), p, NA))
}









##### Tracks importing from local dirs ##########
#################################################
# Choose dir ---------
file_list_DF2_Control_pas.8 <- list.files(path = , choose.dir(),
                           pattern = "csv",
                           all.files = FALSE,
                           full.names = TRUE, recursive = TRUE,
                           ignore.case = FALSE, include.dirs = FALSE,
                           no.. = FALSE)

file_list_p12L <- list.files(path = , choose.dir(),
                             pattern = "csv",
                             all.files = FALSE,
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = FALSE, include.dirs = FALSE,
                             no.. = FALSE)

file_list_p12Y <- list.files(path = , choose.dir(),
                             pattern = "csv",
                             all.files = FALSE,
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = FALSE, include.dirs = FALSE,
                             no.. = FALSE)


#################################################
############# Call to function ##################
#################################################

# Start scan for Control
for (file_name in file_list_DF2_Control_pas.8) {
  tracks_DF2_Control_pas.8 <- rbind(tracks_DF2_Control_pas.8, 
                       traj_analysis_DF2_Control_pas.8(file_name))}


# Start scan for LPA
for (file_name in file_list_p12L) {
  tracks_p12L <- rbind(tracks_p12L, 
                          traj_analysis_p12L(file_name))}


# Start scan for Y27632
for (file_name in file_list_p12Y) {
  tracks_p12Y <- rbind(tracks_p12Y, 
                          traj_analysis_p12Y(file_name))}



#################################################
############# Prepare data ######################
#################################################
# Merge all tracks-----------------------
alltracks <- rbind(tracks_DF2_Control_pas.8) # collect all tracks





# Set time to hours, remove tracks and mean_angle columns
all.h <- alltracks
head(all.h)
data <- cbind(all.h[,c(6,7,8,9)]*3600, all.h[,c(2,3,4,5,10,11,12,13,15)])
data$time <- "24h"
head(data)
tail(data)
write.csv(data, file = 'data_0.csv')










#summary(alltracks)

#write.csv(alltracks, file = 'alltracks.csv') #save results

# Order probe levels
alltracks$probe <- as.factor(alltracks$probe)

alltracks$probe <- ordered(alltracks$probe,
                           levels = c("control", "inhib"))

#################################################
######## Plot all correlations ##################
#################################################
#ggpairs(data)
chart.Correlation(data[,1:12], histogram=TRUE, pch=19, legend = TRUE)
ggcorr(data, palette = "RdBu", label = TRUE)

# mean_speed
# Shapiro test for normality
shapiro.test(data[data$probe=='FC_KC',]$mean_speed) # data is not normally distributed
shapiro.test(data[data$probe=='HAM_KC',]$mean_speed) # data is not normally distributed
shapiro.test(data[data$probe=='TCPS_KC',]$mean_speed)

# Shapiro test for normality
shapiro.test(data[data$probe=='FC_KC',]$sinuosity) # data is not normally distributed
shapiro.test(data[data$probe=='HAM_KC',]$sinuosity) # data is not normally distributed
shapiro.test(data[data$probe=='TCPS_KC',]$sinuosity)

# Shapiro test for normality
shapiro.test(data[data$probe=='FC_KC',]$straight) # data is not normally distributed
shapiro.test(data[data$probe=='HAM_KC',]$straight) # data is not normally distributed
shapiro.test(data[data$probe=='TCPS_KC',]$straight)# data is not normally distributed



# QQ plot 
ggqqplot(data[data$probe=='FC_KC',]$mean_speed, main = 'mean_speed FC_KC')
ggqqplot(data[data$probe=='HAM_KC',]$mean_speed, main = 'mean_speed HAM_KC')
ggqqplot(data[data$probe=='TCPS_KC',]$mean_speed, main = 'mean_speed TCPS_KC')

# QQ plot 
ggqqplot(data[data$probe=='FC_KC',]$sinuosity, main = 'sinuosity FC_KC')
ggqqplot(data[data$probe=='HAM_KC',]$sinuosity, main = 'sinuosity HAM_KC')
ggqqplot(data[data$probe=='TCPS_KC',]$sinuosity, main = 'sinuosity TCPS_KC')

# QQ plot 
ggqqplot(data[data$probe=='FC_KC',]$straight, main = 'straight FC_KC')
ggqqplot(data[data$probe=='HAM_KC',]$straight, main = 'straight HAM_KC')
ggqqplot(data[data$probe=='TCPS_KC',]$straight, main = 'straight TCPS_KC')

# Mean and sd 

mean(data[data$probe=='FC_KC',]$mean_speed)
sd(data[data$probe=='FC_KC',]$mean_speed)

mean(data[data$probe=='HAM_KC',]$mean_speed)
sd(data[data$probe=='HAM_KC',]$mean_speed)

mean(data[data$probe=='TCPS_KC',]$mean_speed)
sd(data[data$probe=='TCPS_KC',]$mean_speed)

# Mean and sd 

mean(data[data$probe=='FC_KC',]$sinuosity)
sd(data[data$probe=='FC_KC',]$sinuosity)

mean(data[data$probe=='HAM_KC',]$sinuosity)
sd(data[data$probe=='HAM_KC',]$sinuosity)

mean(data[data$probe=='TCPS_KC',]$sinuosity)
sd(data[data$probe=='TCPS_KC',]$sinuosity)


# Mean and sd 

mean(data[data$probe=='FC_KC',]$straight)
sd(data[data$probe=='FC_KC',]$straight)

mean(data[data$probe=='HAM_KC',]$straight)
sd(data[data$probe=='HAM_KC',]$straight)

mean(data[data$probe=='TCPS_KC',]$straight)
sd(data[data$probe=='TCPS_KC',]$straight)


# Kruskal-Wallis rank sum test
kruskal.test(mean_speed ~ probe, data = data_KC)
kruskal.test(sinuosity ~ probe, data = data_KC)
kruskal.test(straight ~ probe, data = data_KC)

#fit the one-way ANOVA model
#res.aov <- aov(mean_speed ~ probe, data = data)
#summary(res.aov)
# Pairwise comparisons

compare_means(mean_speed ~ probe,
              data = data_KC,
              method = "wilcox.test")

compare_means(sinuosity ~ probe,
              data = data_KC,
              method = "wilcox.test")

compare_means(straight ~ probe,
              data = data_KC,
              method = "wilcox.test")



# Set X - axis names
CellSciGuylabs <- c("FC_KC", "HAM_KC", "TCPS_KC") #


# Summarize the data :
df_summary <- data_summary(data_KC, varname = "mean_speed", 
                    groupnames = c("probe"))

# Convert dose to a factor variable
df_summary$probe=as.factor(df_summary$probe)
head(df_summary)


########### MEAN SPEED bar plot 
mean_speed <- ggplot(df_summary, aes(x=probe, y=mean_speed, fill="black")) + 
  
  geom_bar(stat="identity",
           position=position_dodge(),
           fill = 'gray') +
  
  geom_errorbar(aes(ymin=mean_speed-sd, ymax=mean_speed+sd), width=.2,
                position=position_dodge(.9)) +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  scale_y_continuous(name = "Mean speed, micron/h",
                     labels = function(x) format(x, scientific = FALSE),
                     breaks = c(seq(0, 25, 10)), 
                     limits = c(-3, 70)) +
  
  labs(y = 'Mean speed, micron/h',
       x = "24 h motility")  +
  
  theme_classic(base_size=14) +
  
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1,     
                                   size = 12, face="bold",
                                   colour="black"))
  
mean_speed
















########### LSC @@@ #############

# Mean and sd 

# mean_speed
# Shapiro test for normality
shapiro.test(data[data$probe=='FC_LSC',]$mean_speed) # data is not normally distributed
shapiro.test(data[data$probe=='HAM_LSC',]$mean_speed) # data is not normally distributed
shapiro.test(data[data$probe=='TCPS_LSC',]$mean_speed)

# Shapiro test for normality
shapiro.test(data[data$probe=='FC_LSC',]$sinuosity) # data is not normally distributed
shapiro.test(data[data$probe=='HAM_LSC',]$sinuosity) # data is not normally distributed
shapiro.test(data[data$probe=='TCPS_LSC',]$sinuosity)

# Shapiro test for normality
shapiro.test(data[data$probe=='FC_LSC',]$straight) # data is not normally distributed
shapiro.test(data[data$probe=='HAM_LSC',]$straight) # data is not normally distributed
shapiro.test(data[data$probe=='TCPS_LSC',]$straight)# data is not normally distributed



# QQ plot 
ggqqplot(data[data$probe=='FC_LSC',]$mean_speed, main = 'mean_speed FC_LSC')
ggqqplot(data[data$probe=='HAM_LSC',]$mean_speed, main = 'mean_speed HAM_LSC')
ggqqplot(data[data$probe=='TCPS_LSC',]$mean_speed, main = 'mean_speed TCPS_LSC')

# QQ plot 
ggqqplot(data[data$probe=='FC_LSC',]$sinuosity, main = 'sinuosity FC_LSC')
ggqqplot(data[data$probe=='HAM_LSC',]$sinuosity, main = 'sinuosity HAM_LSC')
ggqqplot(data[data$probe=='TCPS_LSC',]$sinuosity, main = 'sinuosity TCPS_LSC')

# QQ plot 
ggqqplot(data[data$probe=='FC_LSC',]$straight, main = 'straight FC_LSC')
ggqqplot(data[data$probe=='HAM_LSC',]$straight, main = 'straight HAM_LSC')
ggqqplot(data[data$probe=='TCPS_LSC',]$straight, main = 'straight TCPS_LSC')

# Mean and sd 

mean(data[data$probe=='FC_LSC',]$mean_speed)
sd(data[data$probe=='FC_LSC',]$mean_speed)

mean(data[data$probe=='HAM_LSC',]$mean_speed)
sd(data[data$probe=='HAM_LSC',]$mean_speed)

mean(data[data$probe=='TCPS_LSC',]$mean_speed)
sd(data[data$probe=='TCPS_LSC',]$mean_speed)

# Mean and sd 

mean(data[data$probe=='FC_LSC',]$sinuosity)
sd(data[data$probe=='FC_LSC',]$sinuosity)

mean(data[data$probe=='HAM_LSC',]$sinuosity)
sd(data[data$probe=='HAM_LSC',]$sinuosity)

mean(data[data$probe=='TCPS_LSC',]$sinuosity)
sd(data[data$probe=='TCPS_LSC',]$sinuosity)


# Mean and sd 

mean(data[data$probe=='FC_LSC',]$straight)
sd(data[data$probe=='FC_LSC',]$straight)

mean(data[data$probe=='HAM_LSC',]$straight)
sd(data[data$probe=='HAM_LSC',]$straight)

mean(data[data$probe=='TCPS_LSC',]$straight)
sd(data[data$probe=='TCPS_LSC',]$straight)


# Kruskal-Wallis rank sum test
kruskal.test(mean_speed ~ probe, data = data_LSC)
kruskal.test(sinuosity ~ probe, data = data_LSC)
kruskal.test(straight ~ probe, data = data_LSC)

#fit the one-way ANOVA model
#res.aov <- aov(mean_speed ~ probe, data = data)
#summary(res.aov)
# Pairwise comparisons

compare_means(mean_speed ~ probe,
              data = data_LSC,
              method = "wilcox.test")

compare_means(sinuosity ~ probe,
              data = data_LSC,
              method = "wilcox.test")

compare_means(straight ~ probe,
              data = data_LSC,
              method = "wilcox.test")



# Set X - axis names
CellSciGuylabs <- c("FC_LSC", "HAM_LSC", "TCPS_LSC") #


# Summarize the data :
df_summary <- data_summary(data_LSC, varname = "mean_speed", 
                           groupnames = c("probe"))

# Convert dose to a factor variable
df_summary$probe=as.factor(df_summary$probe)
head(df_summary)


########### MEAN SPEED bar plot 
mean_speed <- ggplot(df_summary, aes(x=probe, y=mean_speed, fill="black")) + 
  
  geom_bar(stat="identity",
           position=position_dodge(),
           fill = 'gray') +
  
  geom_errorbar(aes(ymin=mean_speed-sd, ymax=mean_speed+sd), width=.2,
                position=position_dodge(.9)) +
  
  scale_x_discrete(labels= CellSciGuylabs) +
  
  scale_y_continuous(name = "Mean speed, micron/h",
                     labels = function(x) format(x, scientific = FALSE),
                     breaks = c(seq(0, 25, 10)), 
                     limits = c(-3, 70)) +
  
  labs(y = 'Mean speed, micron/h',
       x = "24 h motility")  +
  
  theme_classic(base_size=14) +
  
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1,     
                                   size = 12, face="bold",
                                   colour="black"))

mean_speed













































#geom_text(size=3.5,aes(y=-2, 
#                       label=c("n = 102", "n = 101", "n = 100"))) # "n = 102", "n = 107", "n = 102"


scatter3d(x = data_KC$mean_speed, y = data_KC$sinuosity, z = data_KC$straight, 
          groups = data_KC$probe,
          xlab=deparse(substitute(mean_speed)), 
          ylab=deparse(substitute(sinuosity)),
          zlab=deparse(substitute(straight)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE,
          ellipsoid = TRUE,
          level=.9, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)



qplot(probe, mean_speed, data = data_KC,
      geom = c("jitter", "boxplot"), alpha = I(0.3), fill = probe,
      main = "Mean speed") +
  labs(y = 'Micrometers per hour',
       x = "") +
  theme_classic(base_size=14)















  # xmin / xmax positions should match the x-axis labels' positions
  geom_signif(y_position = c(101),
              xmin = c(1),
              xmax = c(2),
              annotation = "****",
              tip_length = 0.04)
  
  
  
  
  
  
  
  
  
  # ggplots ----------------------------
  ggplot(data = data_KC, aes(x = probe)) +
    geom_line(aes(y = mean_speed, colour = "mean_speed"), lwd=2, alpha = .99, lty=2) +
    scale_colour_manual("", 
                        breaks = c("mean_speed"),
                        values = c( 
                          "mean_speed"="red"
                        )) +
    xlab("Passages") +
    scale_y_continuous("Micrometers per hour", limits = c(0,50)) + 
   # scale_x_continuous(breaks = c(9,15,28), labels = c("9", "15", "28")) +
    labs(title="") + 
    geom_vline(xintercept = 15, lwd=1, lty=3) +
    geom_vline(xintercept = 28, lwd=1, lty=3) +
    theme_classic2()
  
  
  
  
  
  
  ggplot(data = all, aes(x = probe)) +
    geom_line(aes(y = sinuosity, colour = "Sinuosity"), lwd=2, alpha = .99, lty=2) +
    scale_colour_manual("", 
                        breaks = c("Sinuosity"),
                        values = c( 
                          "Sinuosity"="blue"
                        )) +
    xlab("Passages") +
    scale_y_continuous("Sinuosity", limits = c(0.2,.35)) + 
    scale_x_continuous(breaks = c(9,15,28,36), labels = c("9", "15", "28", "36")) +
    labs(title="") + 
    geom_vline(xintercept = 15, lwd=1, lty=3) +
    geom_vline(xintercept = 28, lwd=1, lty=3) +
    theme_classic2()
  






########### SINUOSITY bar plot 
  
  
  # Summarize the data :
  df_summary <- data_summary(data, varname = "sinuosity", 
                             groupnames = c("probe"))
  
  # Convert dose to a factor variable
  df_summary$probe=as.factor(df_summary$probe)
  head(df_summary)
  
    
  sinuosity <-ggplot(df_summary, aes(x=probe, y=sinuosity, fill="black")) + 
    
    geom_bar(stat="identity",
             position=position_dodge(),
             fill = 'gray') +
    
    geom_errorbar(aes(ymin=sinuosity-sd, ymax=sinuosity+sd), width=.2,
                  position=position_dodge(.9)) +
    
    scale_x_discrete(labels= CellSciGuylabs) +
    
    scale_y_continuous(name = "sinuosity",
                       labels = function(x) format(x, scientific = FALSE),
                       breaks = c(seq(0, 1, .10)), 
                       limits = c(-.03, 1 )) +
    
    labs(y = 'sinuosity',
         x = "DF2 24 h")  +
    
    theme_classic(base_size=14) +
    
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1,     
                                     size = 12, face="bold",
                                     colour="black")) +
    
    geom_text(size=3.5,aes(y=-2, 
                           label=c("n = 102", "n = 101", "n = 100"))) # "n = 102", "n = 107", "n = 102"
  
  sinuosity 
  
  geom_signif(y_position = c(25),
              xmin = c(1),
              xmax = c(2),
              annotation = "*", 
              tip_length = 0.04)
  
  
  ########### emax bar plot 
  
  
  # Summarize the data :
  df_summary <- data_summary(data, varname = "emax", 
                             groupnames = c("probe"))
  
  # Convert dose to a factor variable
  df_summary$probe=as.factor(df_summary$probe)
  head(df_summary)
  
  
  
emax <-  ggplot(df_summary, aes(x=probe, y=emax, fill="black")) + 
    
    geom_bar(stat="identity",
             position=position_dodge(),
             fill = 'gray') +
    
    geom_errorbar(aes(ymin=emax-sd, ymax=emax+sd), width=.2,
                  position=position_dodge(.9)) +
    
    scale_x_discrete(labels= CellSciGuylabs) +
    
    scale_y_continuous(name = "emax",
                       labels = function(x) format(x, scientific = FALSE),
                       breaks = c(seq(0, 2, .20)), 
                       limits = c(-.03, 2.2)) +
    
    labs(y = 'emax',
         x = "FetMSC 24 h")  +
    
    theme_classic(base_size=14) +
    
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1,     
                                     size = 12, face="bold",
                                     colour="black")) +
    
    geom_text(size=3.5,aes(y=-2, 
                           label=c("n = 102", "n = 101", "n = 100"))) # "n = 102", "n = 107", "n = 102"

emax

gridExtra::grid.arrange(mean_speed, sinuosity, emax, ncol = 3 )


#options(rgl.printRglwidget = TRUE)
## Plot data in 3D factor space
scatter3d(x = data$mean_speed, y = data$sinuosity, z = data$emax, 
          groups = data$probe,
          xlab=deparse(substitute(mean_speed)), 
          ylab=deparse(substitute(sinuosity)),
          zlab=deparse(substitute(emax)), 
          axis.scales=FALSE, axis.ticks=TRUE,
          surface=FALSE,
          ellipsoid = TRUE,
          level=0.7, ellipsoid.alpha=0.5, id=FALSE,
          model.summary=FALSE)
