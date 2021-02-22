# Authors: Nayela Zeba, Timothy D. Berry, Thea L. Whitman
# Modified from: An open-source, automated, gas sampling peripheral 
# for laboratory incubation experiments; Authors: Timothy D. Berry, Thea L. Whitman


#Load necessary packages
library (tidyverse)
library (zoo)
library (stats)
library (RColorBrewer)
library (broom)

#The working directory will need to be set by the user to point to where their
#data is stored
#setwd("C:/Users/nayel/Desktop/Pyronema_incubation")

#Changes overflow options to suppress automatic conversion to scientific notation
options(scipen=99)

################################################################################
#Parameter setup
################################################################################

#When smoothing data, how many points wide should the sliding window be? 
#Larger values = smoother
smooth_window = 50

#When checking for outliers, how many points should be analyzed (nslice) and what
#threshold constitutes an outlier (outliercutoff) smaller values will find more 
#outliers, but increase chances of false positives. It is recommended these values
#not be changed unless you have very choppy data even after smoothing
nslice = 10
outliercutoff = 4   ####changed to 4

#Designates the valves used for flushing manifolds (flush_positions)
#and for moving gas through system without sampling jars (idle_positions).
#Anything not listed here will be treated like a sample!
flush_positions <-  c(1,17)
idle_positions <-  c(16, 32)

#Which valves are connected to samples that are your soil treatment 
#(or your first C source) and which have had PyOM (or your second C source) added?
treatments_fungi_biochar <- c(6,9,12,15,20)
treatments_control <- c(4,7,11,13,21)
treatments_fungi_water <- c(5,8,14,18,19)
blanks <- c(2,10,22)

#Fraction of system volume that comes from the sample jar. In a system with no-deadspace 
#this would be 1. In reality, this depends on jar volume and
#must be determined for each specific sample container that is used
fj <- 0.85

#Parameters specific to incubation experimental design:
#Jar volume (L), dry mass of soil (grams), pyom added (mg), and
#PyOM d13C (d13C per mille vPDB)
jarvol_L <- 0.12
pyom_d13C <- 833.11

#Sets the times during which standards are run. Measurements within these windows
#will not be analyzed as if they were samples
#standards_begin <- 1565634800
#standards_end <- 1565638000
  
################################################################################
#Data preprocessing (skip this section if using previously created merged file)
################################################################################

#Assembles a list of log files and Picarro data files from the chosen directory and combines them
logfile_list <- list.files(path="./raw_data", pattern = 'log_*', recursive = T, full.names = TRUE)
logfile_list_read = lapply(logfile_list, read.csv, header = TRUE)
relay_log = do.call(rbind, logfile_list_read)

picarro_files <- list.files(path = "./raw_data", pattern='.*[.]dat$', recursive = T, full.names = TRUE)
picarro_files

all.pic.data <- as.data.frame(c())
total.pic.files <- length(picarro_files)

for(i in 1:total.pic.files){
  temp <- read_table(picarro_files[i])
  all.pic.data <- rbind.data.frame(temp, all.pic.data)
}

rm(i, picarro_files)

#Keep only data columns that are needed for this analysis
columns_to_keep <- c("TIME", "EPOCH_TIME", "12CO2_dry", "Delta_Raw_iCO2","Delta_5min_iCO2")

cleaned.pic.data <- all.pic.data[,columns_to_keep]
rm(all.pic.data, temp)

#Removes unneeded columns and times from data that are from 
#before or after the relay was active (since it's not experiment data)
#Also renames the CO2 concentration column to not start with
#a number so that GG plot doesn't encounter issues later
relay_log= relay_log[order(relay_log$Epoch_time,decreasing = FALSE),]
pic.trimmed <-  cleaned.pic.data[cleaned.pic.data$EPOCH_TIME > relay_log$Epoch_time[1],]
pic.trimmed= pic.trimmed[order(pic.trimmed$EPOCH_TIME),]
pic.trimmed <- pic.trimmed[pic.trimmed$EPOCH_TIME < relay_log$Epoch_time[length(relay_log$Epoch_time)],]
names(pic.trimmed)[3] <- "dry_12CO2"


#Coverts epoch times into integers and removes duplicate timestamps caused by
#the Picarro and relay board not measuring at exact 1 second intervals. 
#This results in a loss of some data points, but there should still be plenty 
#left to redraw the curves and do our calculations with
pic.trimmed <- mutate(pic.trimmed, EPOCH_TIME = as.integer(EPOCH_TIME))
pic.trimmed.deduped <- distinct(pic.trimmed, EPOCH_TIME, .keep_all = TRUE)

#Renames the Epoch_time in the relay file to EPOCH_TIME,
#removes duplicate time stamps, and creates a data frame for the
#relay logs wherein the epoch time matches with a picarro time
names(relay_log)[1] <- "EPOCH_TIME"
relay.deduped <- distinct(relay_log, EPOCH_TIME, .keep_all = TRUE)
relay.matched <- filter(relay.deduped, EPOCH_TIME %in% pic.trimmed.deduped$EPOCH_TIME)

merged <-  merge.data.frame(relay.matched, pic.trimmed.deduped, by.x = "EPOCH_TIME")

#Cleanup to salvage some memory
rm(cleaned.pic.data, pic.trimmed, relay.deduped, relay_log, logfile_list_read, logfile_list,
   pic.trimmed.deduped, relay.matched)

#saveRDS(merged, "merged_r2")

#########################################
### Start here if using previously merged file
#########################################
merged = readRDS("merged_r2")

#Smoothing using the ZOO package. Takes a moving average (determined by the smooth_window parameter, set above). 
# For this dataset, 300 points roughly corresponds to a 5 minute moving average;
# 50 points corresponds to a 50 second moving average.
CO2_smooth <-rollapply(data = merged$dry_12CO2, FUN = "mean", width = smooth_window, partial = TRUE, align = "right")
merged <- cbind(merged,CO2_smooth)

iCO2_smooth <-rollapply(data = merged$Delta_Raw_iCO2, FUN = "mean", width = smooth_window, partial = TRUE, align = "right")
merged <- cbind(merged,iCO2_smooth)

rm(CO2_smooth, iCO2_smooth)

#Removes standard for calibration from merged data frame and creates a second 
#"standards dataframe". 

#standards <- filter(merged, EPOCH_TIME > standards_begin & EPOCH_TIME < standards_end)

#merged <- filter(merged, EPOCH_TIME > standards_end)


################################################################################
#Data Annotation (adds status and cycle info)
################################################################################

#Remove the redundant Time column in the merged dataframe and creates a status, 
#sample number,and cycle number, and cycle_bound columns which are then assigned 
#dummy values for now and are populated in the next steps.

merged$Time <- NULL
merged$status <- 0
merged$sample_num <- 0
merged$cycle <- NA
merged$cycle_bound <- FALSE

#Determines the status of the analysis (stored in the merged dataframe) at each
#time point and assigns sample numbers to jar-specific steps. The status column is 
#mutated to replace the dummy value with a status value assigned based on
#the combination of valves open. For the sake of this process, it is assumed that the default 
#manifold configuration is used and the first position on each manifold is used to 
#flush while the last position on each is an "idle" position in which no jar is sampled.

merged <- merged %>% 
  mutate(status = replace(status, is.element(Active_relay1, flush_positions) & is.element(Active_relay2, idle_positions),"flushing lines")) %>%
  mutate(status = replace(status, is.element(Active_relay1, flush_positions) & is.na(Active_relay2),"flushing lines")) %>%
  mutate(status = replace(status, is.element(Active_relay1,flush_positions) & !is.element(Active_relay2,idle_positions) & !is.na(Active_relay2),"flushing jar")) %>% 
  mutate(status = replace(status, !is.element(Active_relay1,flush_positions) & !is.element(Active_relay1,idle_positions),"measuring jar")) %>%
  mutate(status = replace(status, is.element(Active_relay1,idle_positions),"idle")) %>%
  mutate(status = replace(status, status != lag(status,1),"boundary")) %>%
  mutate(status = replace(status, is.na(lag(Date,1)),"boundary")) %>%
  mutate(sample_num = case_when(.$status == "measuring jar" ~ .$Active_relay1,
                                .$status == "flushing jar" ~ .$Active_relay2,
                                .$status == "boundary" ~ as.integer(0),
                                .$status == "flushing lines" ~ as.integer(0),
                                .$status == "idle" ~ as.integer(0))) %>%
  mutate(cycle_bound = 
           replace(cycle_bound,status == "boundary" &
                     lag(status) == "flushing jar" &
                     lead(status) == "idle", TRUE))

#!!!NOTE: This function assigns cycle numbers based on boundaries between each 
#measurement cycle. When using this function it is assumed that cycles were run 
#to completetion and following the suggested analysis timing template. If this
#is not the case, this function must be modified to accurately label cycle numbers

cycle_counter <- function(x){
                  boundary_vec <- which(x$cycle_bound == TRUE)
                  boundary_vec <- c(1, boundary_vec)
                  for(i in 2:length(boundary_vec)){
                   x[c(boundary_vec[i-1]:boundary_vec[i]),]$cycle <- i -1 
                  }
                  return(x)
}

merged_complete<-cycle_counter(merged)

rm(merged)

################################################################################
#Summarize data and correct for dilution
################################################################################

#Remove boundary and idle statuses
jars <- merged_complete %>% group_by(sample_num, cycle) %>% 
  filter(cycle != 0) %>% filter(status != "boundary") %>% filter(status != "idle")

rm(merged_complete)

#This version of the script has be modified to use the minimum time (aka first data point in measurement) as the "dilute value" so that it can account
#for "negative peaks". To do that the minimum_time column in the min_time df has been redefined to be min(EPOCH_TIME) for a given sample_num/cycle 

most_dilute = jars %>% filter(status == "measuring jar") %>% 
  group_by(sample_num, cycle) %>% summarize(cycle_minimum = min(CO2_smooth))

min_time = jars %>% filter(status == "measuring jar") %>%
  group_by(sample_num, cycle) %>% summarize(minimum_time = min(EPOCH_TIME))

iCO2_dilute = jars%>% 
  # Grab the datapoints right when the jar is being measured
  filter(status == "measuring jar") %>%
  # Create a series of columns identifying the difference between neighboring values,
  # This is not really necessary for the smoothed data but is not computationally
  #intensive so there is little reason to remove it yet.
  select(sample_num,cycle,iCO2_smooth,CO2_smooth, EPOCH_TIME,status) %>%
  mutate(NeighborDist1 = c(NA,diff(iCO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(iCO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(iCO2_smooth,lag=3)),
         NeighborDist4 = c(diff(iCO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(iCO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(iCO2_smooth,lead=3),NA)) %>%
  # Ungroup everything
  group_by(EPOCH_TIME)%>%
  # Create an outlier flag column, that pings when the distance between nearby points is above the outlier cutoff
  mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6),na.rm=TRUE)>outliercutoff,1,0))%>%
  # Group them up again
  group_by(sample_num,cycle) %>%  merge(.,min_time) %>%  
  mutate(difference = minimum_time - EPOCH_TIME) %>% 
  filter(abs(difference) <= nslice/2) %>%
  filter(Flag != 1) %>% group_by(sample_num,cycle) %>%
  # Take the mean value of that raw iCO2 input
  # End result is a table with sample number, cycle number, and mean iCO2
  summarize(dilute_iCO2=mean(iCO2_smooth),dilute_CO2 = mean(CO2_smooth))

#Do the same thing at the end of the flush step to get the purge values (concentration)
#of CO2 in system after flushing. This is also where the summarized sample's time 
#entry comes from - the last point in the purge step is when accumulation begins.
iCO2_purged = jars%>%
  filter(status == "flushing jar") %>%
  select(sample_num,cycle,iCO2_smooth,CO2_smooth, EPOCH_TIME,status) %>%
  mutate(NeighborDist1 = c(NA,diff(iCO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(iCO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(iCO2_smooth,lag=3)),
         NeighborDist4 = c(diff(iCO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(iCO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(iCO2_smooth,lead=3),NA)) %>%
  group_by(EPOCH_TIME)%>%
  mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),
                           abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6)
                           ,na.rm=TRUE)>outliercutoff,1,0)) %>%
  group_by(sample_num,cycle)%>%
  slice((n()-nslice+1):n())   %>%
  filter(Flag != 1) %>%
  summarize(purged_iCO2=mean(iCO2_smooth),
            purged_CO2 = mean(CO2_smooth), time = last(EPOCH_TIME))

#Peak selection takes the last point of the "measuring jar" status and points preceding
#this point by 10 seconds to be representative of the peak value: this tends 
#to represent the most stable part of the "peak" and makes the measurement
#much less likely to be erroneous when the previous sample was very concentrated

meas_end <- jars %>% group_by(sample_num, cycle) %>% 
  filter(status == "measuring jar") %>% summarize(last_point = last(EPOCH_TIME))

iCO2_peak = jars%>%
  filter(status == "measuring jar") %>%
  select(sample_num,cycle,iCO2_smooth,CO2_smooth, EPOCH_TIME,status) %>%
  mutate(NeighborDist1 = c(NA,diff(iCO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(iCO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(iCO2_smooth,lag=3)),
         NeighborDist4 = c(diff(iCO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(iCO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(iCO2_smooth,lead=3),NA)) %>%
  group_by(EPOCH_TIME) %>%
  mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),
                           abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6),
                           na.rm=TRUE)>outliercutoff,1,0)) %>% 
  group_by(sample_num, cycle) %>% merge(.,meas_end) %>%
  mutate(difference = last_point - EPOCH_TIME) %>% filter(difference <= nslice) %>%
  filter(Flag != 1) %>% group_by(sample_num,cycle) %>%
  summarize(peak_iCO2 = mean(iCO2_smooth), peak_CO2 = mean(CO2_smooth))

####Check your cycle assignments####
test <-jars %>% filter(cycle == 15)

ggplot(data = test, aes(x = EPOCH_TIME, y = CO2_smooth, color = status)) +
  geom_point() + facet_wrap(~sample_num, scale = "free")
ggplot(data = test, aes(x = EPOCH_TIME, y = iCO2_smooth, color = status)) +
  geom_point() + facet_wrap(~sample_num, scale = "free")

#Merge these 3 sets of values together and remove unnecessary intermediate files  
summary <- merge(iCO2_dilute, iCO2_peak) %>%
  merge(.,iCO2_purged) %>%
  arrange(sample_num,cycle)

summary_long1 <- summary %>% rename("dilute" = dilute_CO2, "peak" = peak_CO2, "purged" = purged_CO2) %>%
  gather(key = "landmark", value = "CO2", dilute, peak, purged) %>% 
  select(-dilute_iCO2, -peak_iCO2, -purged_iCO2)

summary_long2 <- summary %>% rename("dilute" = dilute_iCO2, "peak" = peak_iCO2, "purged" = purged_iCO2) %>%
  gather(key = "landmark", value = "iCO2", dilute, peak, purged) %>% group_by(sample_num, landmark) %>%
  select(-dilute_CO2, -peak_CO2, -purged_CO2)

summary_long <- merge(summary_long1, summary_long2)

#Plot figures to see the CO2 concentration and isotope profiles acrossed cycles
#This is an easy way to visualize trends in respiration of samples and to assess
#whether your relay command file is adequate. Ideally the dilute and purge values
#should remain consistent throughout the course of the experiment. Changes in
#these values over the course of the experiment likely mean insufficient time is
#being given for the "flushing jars" step. The script can account for this deviations
#but to maximize data quality the flush time should be increased.
respiration_profiles <- ggplot(data = summary_long) + 
  geom_col(aes(x = cycle, y = CO2, fill = landmark), position = "dodge") + 
  facet_wrap(~sample_num, scale = "free")

isotope_profiles <- ggplot(data = summary_long) + 
  geom_col(aes(x = cycle, y = iCO2, fill = landmark), position = "dodge") + 
  facet_wrap(~sample_num, scale = "free")

respiration_profiles
isotope_profiles

#Remove unnecessary intermediate files to free up memory
rm(summary_long1, summary_long2)
rm(meas_end, min_time, most_dilute)
rm(iCO2_dilute, iCO2_peak, iCO2_purged, jars)

##########################################
#Calculations
##########################################

#Dilution parameters are determined by injecting known volumnes of pure CO2 into
#purged vessels to determine the fraction of the contribution of the sample jar's gas (fj) to observed gas


#The total CO2 measured at the peak represents a combination of what's in the system
#before you start measuring and what's in the jar. The total volume of the system is jar
#+ gas lines, instrument, manifolds, etc. Thus, the CO2 in an equilibrated system will 
#be whatever  CO2 came from the system (dilute_CO2) multiplied by its fraction (1-fj),
#plus the CO2 that was in the jar multiplied by its fraction (fj).
#That results in the following equation:
corrected <- summary %>% 
  mutate(peak_CO2_cor = (peak_CO2 - dilute_CO2 + (dilute_CO2 * fj))/fj)

#Next, we need to correct the iCO2 at the peak measurement to account for the fact
#that  our measurement includes mostly  CO2 from our sample, but also some
#from whatever gas remained in the system before initiating the measurement.
#Fortunately, we have those values. First, we calculate what fraction of all the CO2 
#molecules come from our actual sample (pj), vs. what was initially in the system (1-pj). 
#We can get to this by taking the amount of each source of CO2, multiplied by its 
#relative volume. Once we have that value, pj, we can use the same approach as we did
#for the total CO2, above.
i_corrected <- corrected %>% 
  mutate(pj = (peak_CO2_cor * fj)/((peak_CO2_cor * fj + (1-fj) * dilute_CO2))) %>%
  mutate(final_peak_iCO2 = (peak_iCO2 - dilute_iCO2 + (dilute_iCO2 * pj))/pj)

# Then, we just remove the columns we do not want to keep
final <- i_corrected %>% 
  select(-pj, -dilute_iCO2, -peak_iCO2, -peak_CO2, -peak_iCO2, -peak_iCO2)


#Now, we have the estimated values for the volume of CO2 in each jar at measurement (peak_CO2_cor)
#and the isotopic signature of that CO2 (final_peak_iCO2).
#Plus, we have all our baseline values of what the whole system was at when
#it was last purged: purged_CO2 and purged_iCO2_cor. Note, purged_CO2 does not need to
#be corrected, because the whole system was open and mixing at the time of its 
#measurement, so it should be a good estimate for what is in the jars when their
#valves close. The only correction applied to the purged_iCO2_cor was adjusting for
#any machine drift in isotopic measurement.

#Add and populate a treatment column in the finalized data summary
final$treatment <- NA
final <- final %>% mutate(treatment = case_when(
  is.element(sample_num, blanks) ~ "blank",
  is.element(sample_num, treatments_control) ~ "control",
  is.element(sample_num, treatments_fungi_biochar)~ "fungi_biochar",
  is.element(sample_num, treatments_fungi_water)~ "fungi_agar"))

### We will need an estimate of the agar_d13C.
# We can do this from the fungi_agar jars.

agar_13C_calcs <- final %>%
  filter(treatment == "fungi_agar")%>%
  group_by(sample_num)%>%
  mutate(prev_purged_CO2 = lag(purged_CO2))%>%
  mutate(prev_purged_iCO2 = lag(purged_iCO2))

# So, we expect that:
# final_peak_iCO2 = f_prev_purged_CO2*prev_purged_iCO2 + f_new_CO2*agar_d13C
# We have each of these, except agar_d13C, which we will solve for.
# First, adding in the parameters not previously calculated:
agar_13C_calcs = agar_13C_calcs %>%
  mutate(f_prev_purged_CO2 = prev_purged_CO2 / peak_CO2_cor)%>%
  mutate(f_new_CO2 = (peak_CO2_cor-prev_purged_CO2) / peak_CO2_cor)

# Then, let's solve for agar_d13C:
# (final_peak_iCO2 - f_prev_purged_CO2*prev_purged_iCO2) / f_new_CO2 = agar_d13C
# We can do this in the dataset now, too, to get agar_d13C estimates for each timepoint
agar_13C_calcs = agar_13C_calcs %>%
  mutate(agar_d13C = (final_peak_iCO2 - f_prev_purged_CO2*prev_purged_iCO2) / f_new_CO2)
p = ggplot(data=agar_13C_calcs, aes(y=agar_d13C,x=cycle,color=as.factor(sample_num)))+geom_point()
p

# Our estimates vary across the dataset, but are generally in the range of what we'd expect.
# There are a few outliers, especially for the earliest cycles where there just wasn't much CO2
# to speak of. Since there's not a super clear trend, we'll just cut the first cycle.
agar_13C_calcs = agar_13C_calcs %>%
    filter(cycle > 2)
mean(agar_13C_calcs$agar_d13C, na.rm=TRUE)

# Setting the isotopic end-member for non-PyOM-derived emissions
agar_d13C = mean(agar_13C_calcs$agar_d13C, na.rm=TRUE)

####################################################
# Code to estimate abiotic CO2 sorption and release
# E.g., by photolysis or pH equilibration
####################################################

# First, we want to calculate how much CO2 is sorbed by the PyOM, and how much is released each cycle
# We are assuming that no agar C is released abiotically, but that makes more sense, as 
# we wouldn't expect there to be inorganic C in the agar to start with.
sorprelease <- final %>%
  filter(treatment == "control")%>%
  group_by(sample_num)%>%
  mutate(prev_purged_CO2 = lag(purged_CO2))%>%
  mutate(prev_purged_iCO2 = lag(purged_iCO2))%>%
  mutate(frac_new_pyom = ((final_peak_iCO2-prev_purged_iCO2) / (pyom_d13C-prev_purged_iCO2)))%>%
  mutate(net_new_pyom_CO2 = peak_CO2_cor*frac_new_pyom)%>%
  mutate(non_pyom_CO2 = peak_CO2_cor-net_new_pyom_CO2)%>%
  mutate(net_non_pyom_CO2 = non_pyom_CO2-prev_purged_CO2)%>%
  mutate(net_total_sorprelease_CO2 = peak_CO2_cor-prev_purged_CO2)%>%
  mutate(fractional_change = net_total_sorprelease_CO2/peak_CO2_cor)%>%
  mutate(net_effective_sorprelease_iCO2 = final_peak_iCO2/fractional_change - prev_purged_iCO2/fractional_change + prev_purged_iCO2)

# Note: net_non_pyom_CO2 will be negative if there is net sorption, positive if net release.

# Is there a clear relationship to CO2 status at the time the jar closes and how much (nonPyOM) CO2 was sorbed?
ggplot(sorprelease,aes(y=net_non_pyom_CO2,x=prev_purged_CO2,color=as.factor(cycle)))+geom_point()
# I.e., we might be interested if high initial CO2 led to more sorbed CO2.
# If anything, we see the opposite - sorption was greater when there was less CO2 in the jar.
# The purge effects only really vary between ~385-400 - 15ppm - but the differences in net sorption
# are much greater - 200-300 ppm span. Thus, I'd expect we're more seeing the effect of time here
# The purge efficacy varies slightly over time (maybe due to small changes in room conditions, flushing efficiency as net CO2 release gets higher),
# and the sorption varies a lot over time (presumably due to saturation of whatever processes are active)
# The fact that it saturates may mean we don't need to worry too much about whether
# higher levels of CO2 would lead to higher sorption (presumably to some extent, but not over the long term) - 
# it might just saturate a bit sooner.
ggplot(sorprelease,aes(y=net_non_pyom_CO2,x=cycle,color=as.factor(sample_num)))+geom_point()
# We see some sorption trends over time, 
# with more sorption at the beginning, levelling off over time
# Checking also the jar IDs
ggplot(sorprelease,aes(y=net_non_pyom_CO2,x=cycle,color=as.factor(sample_num)))+geom_point()
ggplot(sorprelease,aes(y=net_non_pyom_CO2,x=prev_purged_CO2,color=as.factor(sample_num)))+geom_point()

# As we observed in the lab, jar 13 didn't purge well; possible structural issue with jar seal; will exclude from blanks

# How much CO2 from PyOM do we see being released?
ggplot(sorprelease,aes(y=net_new_pyom_CO2,x=cycle,color=as.factor(sample_num)))+geom_point()
# Abiotic PyOM-CO2 release follows the general trends for the control jars
ggplot(sorprelease,aes(y=net_new_pyom_CO2,x=prev_purged_CO2,color=as.factor(cycle)))+geom_point()
# We see a range of net PyOM release, and also see grouping by cycle (as we saw in the total flux data)
ggplot(sorprelease,aes(y=net_new_pyom_CO2,x=prev_purged_CO2,color=as.factor(sample_num)))+geom_point()
# Again, we see issues with jar 13.

# We will match timepoints between control jars and sample jars

# Now, we are ready to make adjustments to the sample jars to account for sorption and release.
# Since the individual fungal and control jars aren't paired, 
# we will work with mean values from the controls
# for a given cycle.
# Also note here we exclude jar 13, as noted above.
sorprelease.means <- sorprelease %>%
  filter(sample_num != 13)%>%
  group_by(cycle)%>%
  summarize(mean_net_new_abiotic_pyom_CO2 = mean(net_new_pyom_CO2), 
            mean_net_non_pyom_CO2 =mean(net_non_pyom_CO2),
            mean_net_total_sorprelease_CO2=mean(net_total_sorprelease_CO2),
            mean_net_effective_sorprelease_iCO2=mean(net_effective_sorprelease_iCO2))

# We will join these data to our fungal data
final.fungal <- final %>%
  filter(treatment == "fungi_biochar")
final.fungal = merge(final.fungal, sorprelease.means, by="cycle")

# Now we are ready to adjust our total CO2 values in the sample jars with fungi and biochar
# to subtract the net effect of any abiotic sorption and release of CO2.
# First, we take the reference purged CO2 (from the cycle before - what was in the jar at the start of the cycle)
# Then, we adjust its CO2 and iCO2 values to account for other processes we want to subtract: 
# abiotic PyOM release and CO2 sorption.
# After that, we can use these new biotic purged values to calculate the net biotic CO2

final.fungal <- final.fungal %>%
  group_by(sample_num) %>%
  mutate(prev_purged_CO2 = lag(purged_CO2)) %>%
  mutate(prev_purged_iCO2 = lag(purged_iCO2)) %>%
  mutate(sample_biotic_CO2 = peak_CO2_cor - prev_purged_CO2 - mean_net_total_sorprelease_CO2) %>%
  mutate(sample_biotic_iCO2 = (final_peak_iCO2 - prev_purged_iCO2*(prev_purged_CO2/peak_CO2_cor)
                               -pyom_d13C*(mean_net_new_abiotic_pyom_CO2/peak_CO2_cor)
                               -prev_purged_iCO2*(mean_net_non_pyom_CO2/peak_CO2_cor))
         /(sample_biotic_CO2/peak_CO2_cor))

#plot(x=final.fungal2$sample_biotic_iCO2,y=final.fungal$sample_biotic_iCO2)
#plot(x=final.fungal2$sample_biotic_CO2,y=final.fungal$sample_biotic_CO2)
#plot(x=final.fungal2$sample_biotic_CO2,y=final.fungal2$sample_biotic_iCO2)

# Let's make sure - are these values reasonable?
ggplot(final.fungal,aes(x=cycle,y=sample_biotic_CO2,color=as.factor(cycle)))+geom_point()
ggplot(final.fungal,aes(x=cycle,y=sample_biotic_iCO2,color=as.factor(cycle)))+geom_point()
# We see the total CO2 follow the previous pattern
# We see initial high iCO2 contributions, decreasing somewhat
# That generally makes sense.

# We now have the total CO2 that was biotically derived and its 13CO2 signature
# That CO2 is a combination of the PyOM-C and any non-PyOM-C.

# Thus, we will now partition the total CO2 between PyOM-derived and non-PyOM-derived CO2 ("agar_d13C")
final.fungal <- final.fungal %>%
  mutate(frac_pyom =((sample_biotic_iCO2-agar_d13C) / (pyom_d13C-agar_d13C)))%>%
  mutate(sample_biotic_pyom_CO2 = sample_biotic_CO2*frac_pyom)%>%
  mutate(sample_biotic_agar_CO2 = sample_biotic_CO2*(1-frac_pyom))

# In addition to the partitioned data, 
# We are also interested in the totals for the inoculated jars with biochar
final.inoculated <- final %>%
  filter(treatment == "fungi_biochar")%>%
  group_by(sample_num) %>% 
  mutate(sample_CO2 = peak_CO2_cor - lag(purged_CO2)) %>%
  mutate(sample_iCO2 = (final_peak_iCO2 - lag(purged_iCO2)*(lag(purged_CO2)/(lag(purged_CO2) + sample_CO2)))
         /(sample_CO2/(sample_CO2 + lag(purged_CO2))))
colnames(final.inoculated)[c(10,11)] = c("inoculated_CO2","inoculated_iCO2")

# Covert CO2 values from ppm to a mass basis for control CO2 data.
# We make some simple assumptions that the 
# sample is at standard temperature and pressure
# We also adjust time into time since incubations started and convenient units
mineralization.inoculated <- final.inoculated %>%
  dplyr::group_by(sample_num) %>%
  dplyr::mutate(time_s = time-min(time)) %>%
  dplyr::mutate(time_min = time_s/60) %>%
  dplyr::mutate(time_hr = ((time_min/60)+72)) %>% ### 1st measurement was done 3 days after inoculation 
  dplyr::arrange(time_hr) %>%
  dplyr::mutate(inoculated_CO2C_mg = inoculated_CO2/1000000/22.4*jarvol_L*12.01*1000) %>% 
  dplyr::filter(!is.na(inoculated_CO2))%>%
  select(sample_num, cycle, time_s, time_min, time_hr,inoculated_CO2C_mg, inoculated_iCO2)

#Calculating cumulative CO2-C mineralized over each measurement time for control CO2 data
mineralization.inoculated = mineralization.inoculated%>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(Total_inoculated_CO2C_mg.cumulative = cumsum(inoculated_CO2C_mg))

# We are also interested in the uninoculated jars with biochar
# Select the uninoculated jars only
final.uninoculated <- final %>%
  filter(treatment == "control")%>%
  group_by(sample_num) %>% 
  mutate(sample_CO2 = peak_CO2_cor - lag(purged_CO2)) %>%
  mutate(sample_iCO2 = (final_peak_iCO2 - lag(purged_iCO2)*(lag(purged_CO2)/(lag(purged_CO2) + sample_CO2)))
         /(sample_CO2/(sample_CO2 + lag(purged_CO2))))
colnames(final.uninoculated)[c(10,11)] = c("uninoculated_CO2","uninoculated_iCO2")

# Covert CO2 values from ppm to a mass basis for control CO2 data.
# We make some simple assumptions that the 
# sample is at standard temperature and pressure
# We also adjust time into time since incubations started and convenient units
mineralization.uninoculated <- final.uninoculated %>%
  dplyr::group_by(sample_num) %>%
  dplyr::mutate(time_s = time-min(time)) %>%
  dplyr::mutate(time_min = time_s/60) %>%
  dplyr::mutate(time_hr = ((time_min/60)+72)) %>% ### 1st measurement was done 3 days after inoculation 
  dplyr::arrange(time_hr) %>%
  dplyr::mutate(uninoculated_CO2C_mg = uninoculated_CO2/1000000/22.4*jarvol_L*12.01*1000) %>% 
  dplyr::filter(!is.na(uninoculated_CO2))%>%
  select(sample_num, cycle, time_s, time_min, time_hr,uninoculated_CO2C_mg, uninoculated_iCO2)

#Calculating cumulative CO2-C mineralized over each measurement time for control CO2 data
mineralization.uninoculated = mineralization.uninoculated%>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(Total_uninoculated_CO2C_mg.cumulative = cumsum(uninoculated_CO2C_mg))

# We will do the same for the partitioned inoculated data
final_summary <- final.fungal %>% 
  select(sample_num, cycle, sample_biotic_iCO2, sample_biotic_pyom_CO2, sample_biotic_agar_CO2,  time)

#As above, we covert CO2 values from ppm to a mass basis,
# We make some simple assumptions that the 
#sample is at standard temperature and pressure
mineralization <- final_summary %>%
  dplyr::group_by(sample_num) %>%
  dplyr::mutate(time_s = time-min(time)) %>%
  dplyr::mutate(time_min = time_s/60) %>%
  dplyr::mutate(time_hr = ((time_min/60)+72)) %>% ### 1st measurement was done 3 days after inoculation 
  dplyr::arrange(time_hr) %>%
  dplyr::mutate(sample_biotic_pyom_CO2C_mg = sample_biotic_pyom_CO2/1000000/22.4*jarvol_L*12.01*1000) %>% 
  dplyr::mutate(sample_biotic_agar_CO2C_mg = sample_biotic_agar_CO2/1000000/22.4*jarvol_L*12.01*1000) %>% 
  dplyr::filter(!is.na(sample_biotic_pyom_CO2))

#As above, we will calculate cumulative CO2-C mineralized for each measurement time
df.c = mineralization%>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(Biotic_pyom_CO2C_mg.cumulative = cumsum(sample_biotic_pyom_CO2C_mg))%>%
  dplyr::mutate(Biotic_agar_CO2C_mg.cumulative = cumsum(sample_biotic_agar_CO2C_mg))

#################################################################################
#Plotting
#################################################################################

#Creating standard error functions
lower_se <- function(mean, se){
  res <- mean - se
  return(res)
}
upper_se <- function(mean, se){
  res <- mean + se
  return(res)
}

# Calculating the means and standard errors for corrected CO2-C values at each time point 
# for the biotic partitioned CO2 from PyOM+fungi jars
df.means = df.c %>%
  dplyr::group_by(cycle) %>%
  dplyr::mutate(
    mean_biotic_pyom_CO2C_mg.cumulative = mean(Biotic_pyom_CO2C_mg.cumulative),
    sd_biotic_pyom_CO2C_mg.cumulative=sd(Biotic_pyom_CO2C_mg.cumulative),
    mean_biotic_agar_CO2C_mg.cumulative = mean(Biotic_agar_CO2C_mg.cumulative),
    sd_biotic_agar_CO2C_mg.cumulative=sd(Biotic_agar_CO2C_mg.cumulative),
    Mean_time_hr = mean(time_hr),
    n_CO2C_mg.cumulative= n()) %>%
  dplyr::mutate(SE_biotic_CO2C_mg.cumulative= sd_biotic_pyom_CO2C_mg.cumulative/sqrt(n_CO2C_mg.cumulative),
                lower_se_biotic_pyom_CO2C_mg.cumulative=lower_se(mean_biotic_pyom_CO2C_mg.cumulative, SE_biotic_CO2C_mg.cumulative),
                upper_se_biotic_pyom_CO2C_mg.cumulative=upper_se(mean_biotic_pyom_CO2C_mg.cumulative, SE_biotic_CO2C_mg.cumulative),
                SE_biotic_agar_CO2C_mg.cumulative= sd_biotic_agar_CO2C_mg.cumulative/sqrt(n_CO2C_mg.cumulative),
                lower_se_biotic_agar_CO2C_mg.cumulative=lower_se(mean_biotic_agar_CO2C_mg.cumulative, SE_biotic_agar_CO2C_mg.cumulative),
                upper_se_biotic_agar_CO2C_mg.cumulative=upper_se(mean_biotic_agar_CO2C_mg.cumulative, SE_biotic_agar_CO2C_mg.cumulative))

# We also want this for the non-partitioned PyOM+fungi jars.
# Total inoculated
mineralization.inoculated.means = mineralization.inoculated %>%
  dplyr::group_by(cycle) %>%
  dplyr::mutate(
    Mean_inoculated_CO2C_mg.cumulative = mean(Total_inoculated_CO2C_mg.cumulative),
    sd_inoculated_CO2C_mg.cumulative=sd(Total_inoculated_CO2C_mg.cumulative),
    Mean_time_hr = mean(time_hr),
    n_CO2C_mg.cumulative= n()) %>%
  dplyr::mutate(SE_CO2C_mg.cumulative= sd_inoculated_CO2C_mg.cumulative/sqrt(n_CO2C_mg.cumulative),
                lower_se_inoculated_CO2C_mg.cumulative=lower_se(Mean_inoculated_CO2C_mg.cumulative, SE_CO2C_mg.cumulative),
                upper_se_inoculated_CO2C_mg.cumulative=upper_se(Mean_inoculated_CO2C_mg.cumulative, SE_CO2C_mg.cumulative))

# Perform the same summary statistics with the uninoculated jars
# Total uninoculated
mineralization.uninoculated.means = mineralization.uninoculated %>%
  dplyr::group_by(cycle) %>%
  dplyr::mutate(
    Mean_uninoculated_CO2C_mg.cumulative = mean(Total_uninoculated_CO2C_mg.cumulative),
    sd_uninoculated_CO2C_mg.cumulative=sd(Total_uninoculated_CO2C_mg.cumulative),
    Mean_time_hr = mean(time_hr),
    n_CO2C_mg.cumulative= n()) %>%
  dplyr::mutate(SE_CO2C_mg.cumulative= sd_uninoculated_CO2C_mg.cumulative/sqrt(n_CO2C_mg.cumulative),
                lower_se_uninoculated_CO2C_mg.cumulative=lower_se(Mean_uninoculated_CO2C_mg.cumulative, SE_CO2C_mg.cumulative),
                upper_se_uninoculated_CO2C_mg.cumulative=upper_se(Mean_uninoculated_CO2C_mg.cumulative, SE_CO2C_mg.cumulative))

############## Figures for paper ##################
# Plot the final data, with points and error bars of SE, joined by lines
# Panel A - Cumulative total emissions for both PyOM jars - one with fungi, one without

# Create plot
p = ggplot()

# Add total inoculated PyOM plates data
p = p + geom_errorbar(data=mineralization.inoculated.means,aes(x=Mean_time_hr,
                                                               ymin=lower_se_inoculated_CO2C_mg.cumulative,
                                                               ymax=upper_se_inoculated_CO2C_mg.cumulative, 
                                                               width=5))
p = p + geom_point(data=mineralization.inoculated.means,aes(x=Mean_time_hr,y=Mean_inoculated_CO2C_mg.cumulative),color="grey42",shape=15, size=2)
p = p + geom_line (data=mineralization.inoculated.means,aes(x=Mean_time_hr,y=Mean_inoculated_CO2C_mg.cumulative))



# Add abiotic (uninoculated) data
p = p + geom_errorbar(data=mineralization.uninoculated.means,aes(x=Mean_time_hr,
                                                                 ymin=lower_se_uninoculated_CO2C_mg.cumulative,
                                                                 ymax=upper_se_uninoculated_CO2C_mg.cumulative, 
                                                                 width=5))
p = p + geom_point(data=mineralization.uninoculated.means,aes(x=Mean_time_hr,y=Mean_uninoculated_CO2C_mg.cumulative),color="red",shape=18, size=3)
p = p + geom_line (data=mineralization.uninoculated.means,aes(x=Mean_time_hr,y=Mean_uninoculated_CO2C_mg.cumulative))

# Formatting
p = p + geom_abline(slope=0,intercept=0)
p = p + xlab("Time since inoculation (hours)") + ylab(expression(Cumulative~total~CO[2]-C~(mg)))
p = p + theme_bw()
p = p + theme(legend.title = element_blank())
p = p + theme(axis.text=element_text(size=12),axis.title = element_text(size=14))
p = p + ylim(min=-0.06,max=0.31)
p

# Panel B - Cumulative total emissions for biotic PyOM and biotic non-PyOM C in inoculated PyOM jars

# Create plot
p = ggplot()

# Add Biotic non-PyOM data
p = p + geom_errorbar(data=df.means,aes(x=Mean_time_hr,
                                        ymin=lower_se_biotic_agar_CO2C_mg.cumulative,
                                        ymax=upper_se_biotic_agar_CO2C_mg.cumulative, 
                                        width=5))
p = p + geom_point(data=df.means,aes(x=Mean_time_hr,y=mean_biotic_agar_CO2C_mg.cumulative),color="grey84", shape=15,size=2)
p = p + geom_line (data=df.means,aes(x=Mean_time_hr,y=mean_biotic_agar_CO2C_mg.cumulative))

# Add biotic PyOM data
p = p + geom_errorbar(data=df.means,aes(x=Mean_time_hr,
                                        ymin=lower_se_biotic_pyom_CO2C_mg.cumulative,
                                        ymax=upper_se_biotic_pyom_CO2C_mg.cumulative, 
                                        width=5))
p = p + geom_point(data=df.means,aes(x=Mean_time_hr,y=mean_biotic_pyom_CO2C_mg.cumulative),color="Black", shape=15,size=2)
p = p + geom_line(data=df.means,aes(x=Mean_time_hr,y=mean_biotic_pyom_CO2C_mg.cumulative))

# Formatting
p = p + geom_abline(slope=0,intercept=0)
p = p + xlab("Time since inoculation (hours)") + ylab(expression(Cumulative~fungal~CO[2]-C~(mg)))
p = p + theme_bw()
p = p + theme(legend.title = element_text())
p = p + theme(axis.text=element_text(size=12),axis.title = element_text(size=14))
p = p + ylim(min=-0.06,max=0.31)
p
