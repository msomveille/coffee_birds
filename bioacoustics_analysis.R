library(tidyverse)
library(elevatr)
library(lubridate)
library(pheatmap)
library(ggpubr)
library(sp)
library(sf)
library(terra)
library(rnaturalearth)
library(httr)
library(tidyterra)
library(ggspatial)


# Load detection data

bosque <- read_csv("resources/Buggs/Detections/Bugg_detections_Bosque_fullyear.csv") %>%
  subset(Confidence > 0.8) %>%
  mutate(File = as_datetime(Date) - hours(6)) %>%
  mutate(Date = as_datetime(Date) - hours(6) + seconds(Begin.Time..s.)) %>%
  separate(Date, into = c("Date", "Time"), sep = " ") %>%
  separate(Time, into = c("Hour", "Min", "Sec"), sep = ":") %>%
  dplyr::mutate(Hour = as.numeric(Hour)) %>%
  dplyr::mutate(Min = as.numeric(Min)) %>%
  dplyr::mutate(Sec = as.numeric(Sec)) %>%
  mutate(site = "bosque") %>%
  filter(Date > "2022-02-20" & Date <= "2023-02-20") %>%
  mutate(jday = yday(as.Date(Date))) %>%
  select(site, File, Date, jday, Hour, Min, Sec, Species.Code, Common.Name, Confidence)

interior <- read_csv("resources/Buggs/Detections/Bugg_detections_IrlandaInterior_fullyear.csv") %>%
  subset(Confidence > 0.8) %>%
  mutate(File = as_datetime(Date) - hours(6)) %>%
  mutate(Date = as_datetime(Date) - hours(6) + seconds(Begin.Time..s.)) %>%
  separate(Date, into = c("Date", "Time"), sep = " ") %>%
  separate(Time, into = c("Hour", "Min", "Sec"), sep = ":") %>%
  dplyr::mutate(Hour = as.numeric(Hour)) %>%
  dplyr::mutate(Min = as.numeric(Min)) %>%
  dplyr::mutate(Sec = as.numeric(Sec)) %>%
  mutate(site = "interior") %>%
  filter(Date > "2022-02-20" & Date <= "2023-02-20") %>%
  mutate(jday = yday(as.Date(Date))) %>%
  select(site, File, Date, jday, Hour, Min, Sec, Species.Code, Common.Name, Confidence)

hamburgo <- read_csv("resources/Buggs/Detections/Bugg_detections_Hamburgo_fullyear.csv") %>%
  subset(Confidence > 0.8) %>%
  mutate(File = as_datetime(Date) - hours(6)) %>%
  mutate(Date = as_datetime(Date) - hours(6) + seconds(Begin.Time..s.)) %>%
  separate(Date, into = c("Date", "Time"), sep = " ") %>%
  separate(Time, into = c("Hour", "Min", "Sec"), sep = ":") %>%
  dplyr::mutate(Hour = as.numeric(Hour)) %>%
  dplyr::mutate(Min = as.numeric(Min)) %>%
  dplyr::mutate(Sec = as.numeric(Sec)) %>%
  mutate(site = "hamburgo") %>%
  filter(Date > "2022-02-20" & Date <= "2023-02-20") %>%
  mutate(jday = yday(as.Date(Date))) %>%
  select(site, File, Date, jday, Hour, Min, Sec, Species.Code, Common.Name, Confidence)


# BirdNET calibration

validation_data <- read_csv("resources/Buggs/Validation/Checking_Chiapas_detections.csv")
validation_species <- unique(validation_data$Species)

validation_data_spp <- vector()
for(i in 1:length(validation_species)){
  validation_data_spp[i] <- validation_data %>% filter(Species == validation_species[i]) %>%
    filter(BirdNET_correct == "Yes") %>%
    filter(Confidence %in% c("High", "high", "Medium")) %>%
    count()/50
}
validation_results <- data.frame(species = validation_species, 
                                 precision = as.numeric(unlist(validation_data_spp)))
species_validated <- validation_results$species[which(validation_results$precision > 0.8)] # 50 species

bosque <- bosque %>%
  filter(Common.Name %in% species_validated)
interior <- interior %>%
  filter(Common.Name %in% species_validated)
hamburgo <- hamburgo %>%
  filter(Common.Name %in% species_validated)

# Total number of detections
nrow(bosque) + nrow(interior) + nrow(hamburgo) # 123355


# Estimate sampling effort

## Load audio data
bosque_audio <- read_csv("resources/Buggs/Detections/file_dates_Bosque_fullyear.csv") %>%
  rename(Date = x) %>% mutate(Date = as_datetime(Date) - hours(6)) %>% separate(Date, into = c("Day", "Time"), sep = " ", remove=F) %>%
  separate(Time, into = c("Hour", "Min", "Sec"), sep = ":", remove=F) %>%
  dplyr::mutate(Hour = as.numeric(Hour)) %>% 
  mutate(site = "bosque") %>%
  filter(Day > "2022-02-20" & Day <= "2023-02-20") %>%
  mutate(jday = yday(as.Date(Day)))
interior_audio <- read_csv("resources/Buggs/Detections/file_dates_IrlandaInterior_fullyear.csv") %>%
  rename(Date = x) %>% mutate(Date = as_datetime(Date) - hours(6)) %>% separate(Date, into = c("Day", "Time"), sep = " ", remove=F) %>%
  separate(Time, into = c("Hour", "Min", "Sec"), sep = ":", remove=F) %>%
  dplyr::mutate(Hour = as.numeric(Hour)) %>% 
  mutate(site = "interior") %>%
  filter(Day > "2022-02-20" & Day <= "2023-02-20") %>%
  mutate(jday = yday(as.Date(Day)))
hamburgo_audio <- read_csv("resources/Buggs/Detections/file_dates_Hamburgo_fullyear.csv") %>%
  rename(Date = x) %>% mutate(Date = as_datetime(Date) - hours(6)) %>% separate(Date, into = c("Day", "Time"), sep = " ", remove=F) %>%
  separate(Time, into = c("Hour", "Min", "Sec"), sep = ":", remove=F) %>%
  dplyr::mutate(Hour = as.numeric(Hour)) %>% 
  mutate(site = "hamburgo") %>%
  filter(Day > "2022-02-20" & Day <= "2023-02-20") %>%
  mutate(jday = yday(as.Date(Day)))

# Total sampling effort (number of 5-min files) in hours
(nrow(bosque_audio) + nrow(interior_audio) + nrow(hamburgo_audio)) * 5/60 # 17,654 hours


# Sampling effort per 5-day period
first_day <- seq(yday(as.Date("2022-02-21")), yday(as.Date("2022-02-21")) + 365, 5)
first_day <- if_else(first_day > 365, first_day - 365, first_day)
bosque_audio_period <- rep(0, nrow(bosque_audio))
interior_audio_period <- rep(0, nrow(interior_audio))
hamburgo_audio_period <- rep(0, nrow(hamburgo_audio))
bosque_detection_period <- rep(0, nrow(bosque))
interior_detection_period <- rep(0, nrow(interior))
hamburgo_detection_period <- rep(0, nrow(hamburgo))
for(i in 1:(length(first_day)-1)){
  if(first_day[i] + 5 > 365){
    bosque_audio_period[which(bosque_audio$jday >= first_day[i] | bosque_audio$jday < first_day[i+1])] <- i
    interior_audio_period[which(interior_audio$jday >= first_day[i] | interior_audio$jday < first_day[i+1])] <- i
    hamburgo_audio_period[which(hamburgo_audio$jday >= first_day[i] | hamburgo_audio$jday < first_day[i+1])] <- i
    bosque_detection_period[which(bosque$jday >= first_day[i] | bosque$jday < first_day[i+1])] <- i
    interior_detection_period[which(interior$jday >= first_day[i] | interior$jday < first_day[i+1])] <- i
    hamburgo_detection_period[which(hamburgo$jday >= first_day[i] | hamburgo$jday < first_day[i+1])] <- i
  }else{
    bosque_audio_period[which(bosque_audio$jday >= first_day[i] & bosque_audio$jday < first_day[i+1])] <- i
    interior_audio_period[which(interior_audio$jday >= first_day[i] & interior_audio$jday < first_day[i+1])] <- i
    hamburgo_audio_period[which(hamburgo_audio$jday >= first_day[i] & hamburgo_audio$jday < first_day[i+1])] <- i
    bosque_detection_period[which(bosque$jday >= first_day[i] & bosque$jday < first_day[i+1])] <- i
    interior_detection_period[which(interior$jday >= first_day[i] & interior$jday < first_day[i+1])] <- i
    hamburgo_detection_period[which(hamburgo$jday >= first_day[i] & hamburgo$jday < first_day[i+1])] <- i
  }
}
bosque_audio$period <- bosque_audio_period
interior_audio$period <- interior_audio_period
hamburgo_audio$period <- hamburgo_audio_period
bosque$period <- bosque_detection_period
interior$period <- interior_detection_period
hamburgo$period <- hamburgo_detection_period


# Figure 3: plots of recording effort and detections across 5-day periods. 
audio_all <- rbind(hamburgo_audio, interior_audio, bosque_audio)
audio_heat <- table(audio_all$site, audio_all$period)
audio_heat <- audio_heat[c(2,3,1),]
row.names(audio_heat) <- c("sun coffee", "shade coffee", "forest")
detections_all <- rbind(bosque, interior, hamburgo)
detections_heat <- table(detections_all$site, detections_all$period)
detections_heat <- detections_heat[c(2,3,1),]
row.names(detections_heat) <- c("sun coffee", "shade coffee", "forest")

data_for_plot <- data.frame(site = rep(row.names(audio_heat), each=73),
                            period = rep(seq.Date(from=as.Date("2022-02-23"), to=as.Date("2023-02-18"), by=5), 6),
                            recording_effort = c(audio_heat[1,], audio_heat[2,], audio_heat[3,]),
                            detections = c(detections_heat[1,], detections_heat[2,], detections_heat[3,])) %>%
  mutate(detection_rate = detections/recording_effort)

cols <- c("forest" = "purple2", "shade coffee" = "darkgreen", "sun coffee" = "chocolate1")
g_audio <- ggplot(data=data_for_plot) +
  geom_point(aes(x=period, y=recording_effort*5, col=site), size=1) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("recording per 5-day period (in minutes)") + theme_light() + ggtitle("(a)") +
  theme(legend.title = element_text(size=13), legend.text = element_text(size=12))
g_det <- ggplot(data=data_for_plot) +
  geom_point(aes(x=period, y=detections, col=site), size=1) +
  geom_smooth(aes(x=period, y=detections, col=site), se=F, span=0.25) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Number of detections per 5-day period") + theme_light() + ggtitle("(b)") +
  theme(legend.title = element_text(size=13), legend.text = element_text(size=12))
g_det_rate <- ggplot(data=data_for_plot) +
  geom_point(aes(x=period, y=detections/(recording_effort*5), col=site), size=1) +
  geom_smooth(aes(x=period, y=detections/(recording_effort*5), col=site), se=F, span=0.25) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Detection rate (detections per minute)") + theme_light() + ggtitle("(c)") +
  theme(legend.title = element_text(size=13), legend.text = element_text(size=12))

pdf(file = "results/Figures/Figure_3.pdf", width = 6, height = 10)
ggarrange(g_audio, g_det, g_det_rate, nrow=3, ncol=1, common.legend = TRUE, legend="bottom")
dev.off()



# Quantify species temporal occurence

##  Hamburgo (sun coffee)  ##

spp_list_hamburgo <- names(table(hamburgo$Common.Name))[which(table(hamburgo$Common.Name) > 20)] # species with at least 20 detections
spp_detections_hamburgo <- matrix(0, ncol=73, nrow=length(spp_list_hamburgo))
for(k in 1:length(spp_list_hamburgo)){
  hamburgo_spp <- hamburgo %>%
    filter(Common.Name == spp_list_hamburgo[k])
  files_periods <- table(hamburgo_spp$File, hamburgo_spp$period)
  spp_det <- apply(files_periods, 2, function(x) length(which(x > 0)))
  spp_detections_hamburgo[k,][match(as.numeric(names(spp_det)), 1:73)] <- spp_det
}
spp_detection_rate_hamburgo <- t(apply(spp_detections_hamburgo, 1, function(x) x/audio_heat[1,]))
spp_detection_rate_hamburgo <- t(apply(spp_detection_rate_hamburgo, 1, function(x) x/max(x, na.rm=T)))
spp_detection_rate_hamburgo2 <- cbind(spp_detection_rate_hamburgo, spp_detection_rate_hamburgo, spp_detection_rate_hamburgo)
# Estimate species' temporal range
spp_ranges_hamburgo <- matrix(0, ncol=73, nrow=length(spp_list_hamburgo))
for(k in 1:length(spp_list_hamburgo)){
  pres1 <- which(spp_detection_rate_hamburgo2[k,] > 0.05)
  pres_spp <- vector()
  for(i in 1:ncol(spp_detection_rate_hamburgo2)){
    aa <- i - pres1
    if(min(aa[aa>=0]) + min(abs(aa[aa<=0])) > 8){
      pres_spp[i] <- 0
    }else{
      pres_spp[i] <- 1
    }
  }
  pres2 <- which(pres_spp > 0)
  for(i in 1:length(pres_spp)){
    aa <- i - pres2
    if(min(min(aa[aa>0]), min(abs(aa[aa<0]))) > 4){
      pres_spp[i] <- 0
    }
  }
  spp_ranges_hamburgo[k,] <- pres_spp[74:(73*2)]
}
spp_ranges_hamburgo <- spp_ranges_hamburgo[which(apply(spp_ranges_hamburgo, 1, sum) > 0),]


##  Interior (shade coffee)  ##

spp_list_interior <- names(table(interior$Common.Name))[which(table(interior$Common.Name) > 20)] # species with at least 20 detections
spp_detections_interior <- matrix(0, ncol=73, nrow=length(spp_list_interior))
for(k in 1:length(spp_list_interior)){
  interior_spp <- interior %>%
    filter(Common.Name == spp_list_interior[k])
  files_periods <- table(interior_spp$File, interior_spp$period)
  spp_det <- apply(files_periods, 2, function(x) length(which(x > 0)))
  spp_detections_interior[k,][match(as.numeric(names(spp_det)), 1:73)] <- spp_det
}
spp_detection_rate_interior <- t(apply(spp_detections_interior, 1, function(x) x/audio_heat[2,]))
spp_detection_rate_interior <- t(apply(spp_detection_rate_interior, 1, function(x) x/max(x, na.rm=T)))
spp_detection_rate_interior2 <- cbind(spp_detection_rate_interior, spp_detection_rate_interior, spp_detection_rate_interior)
# Estimate species' temporal range
spp_ranges_interior <- matrix(0, ncol=73, nrow=length(spp_list_interior))
for(k in 1:length(spp_list_interior)){
  pres1 <- which(spp_detection_rate_interior2[k,] > 0.05)
  pres_spp <- vector()
  for(i in 1:ncol(spp_detection_rate_interior2)){
    aa <- i - pres1
    if(min(aa[aa>=0]) + min(abs(aa[aa<=0])) > 8){
      pres_spp[i] <- 0
    }else{
      pres_spp[i] <- 1
    }
  }
  pres2 <- which(pres_spp > 0)
  for(i in 1:length(pres_spp)){
    aa <- i - pres2
    if(min(min(aa[aa>0]), min(abs(aa[aa<0]))) > 4){
      pres_spp[i] <- 0
    }
  }
  spp_ranges_interior[k,] <- pres_spp[74:(73*2)]
}
spp_ranges_interior <- spp_ranges_interior[which(apply(spp_ranges_interior, 1, sum) > 0),]


##  Bosque (forest remnant)  ##

spp_list_bosque <- names(table(bosque$Common.Name))[which(table(bosque$Common.Name) > 20)] # species with at least 20 detections
spp_detections_bosque <- matrix(0, ncol=73, nrow=length(spp_list_bosque))
for(k in 1:length(spp_list_bosque)){
  bosque_spp <- bosque %>%
    filter(Common.Name == spp_list_bosque[k])
  files_periods <- table(bosque_spp$File, bosque_spp$period)
  spp_det <- apply(files_periods, 2, function(x) length(which(x > 0)))
  spp_detections_bosque[k,][match(as.numeric(names(spp_det)), 1:73)] <- spp_det
}
spp_detection_rate_bosque <- t(apply(spp_detections_bosque, 1, function(x) x/audio_heat[3,]))
spp_detection_rate_bosque <- t(apply(spp_detection_rate_bosque, 1, function(x) x/max(x, na.rm=T)))
spp_detection_rate_bosque2 <- cbind(spp_detection_rate_bosque, spp_detection_rate_bosque, spp_detection_rate_bosque)
# Estimate species' temporal range
spp_ranges_bosque <- matrix(0, ncol=73, nrow=length(spp_list_bosque))
for(k in 1:length(spp_list_bosque)){
  pres1 <- which(spp_detection_rate_bosque2[k,] > 0.05)
  pres_spp <- vector()
  for(i in 1:ncol(spp_detection_rate_bosque2)){
    aa <- i - pres1
    if(min(aa[aa>=0]) + min(abs(aa[aa<=0])) > 8){
      pres_spp[i] <- 0
    }else{
      pres_spp[i] <- 1
    }
  }
  pres2 <- which(pres_spp > 0)
  for(i in 1:length(pres_spp)){
    aa <- i - pres2
    if(min(min(aa[aa>0]), min(abs(aa[aa<0]))) > 4){
      pres_spp[i] <- 0
    }
  }
  spp_ranges_bosque[k,] <- pres_spp[74:(73*2)]
}
spp_ranges_bosque <- spp_ranges_bosque[which(apply(spp_ranges_bosque, 1, sum) > 0),]

# List of unique species across zones
spp_list_all <- unique(c(spp_list_bosque, spp_list_interior, spp_list_hamburgo))
#write.csv(spp_list_all, "resources/species_list2.csv")


## Put species functional trait data together ##

spp_data <- data.frame(Species = spp_list_all)
#spp_data <- read_csv("resources/species_list2.csv")[,-1] %>%
#  rename(Species = x)
# Get foraging strata data
foraging_strata <- read_csv("resources/foraging_strata.csv")
# Get trait data from AVONET
ebird.taxonomy <- read_csv("resources/eBird_Taxonomy_v2019.csv")
spp_data$sci_name <- ebird.taxonomy$SCI_NAME[match(spp_data$Species, ebird.taxonomy$PRIMARY_COM_NAME)]
avonet.data <- read_csv("resources/AVONET2_ebird.csv") %>%
  filter(Species2 %in% spp_data$sci_name) %>%
  rename(sci_name = Species2)
spp_data <- spp_data %>% left_join(avonet.data) %>% left_join(foraging_strata)


## Figure 4: species temporal ranges

spp_ranges_plot <- list()
for(k in 1:length(spp_list_bosque)){
  if(spp_ranges_bosque[k,][1] == 1 & spp_ranges_bosque[k,][73] == 1){
    aa <- spp_ranges_bosque[k,][1:72] - spp_ranges_bosque[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="forest", Species=spp_list_bosque[k], start=c(1, which(aa == -1)+1), end=c(which(aa == 1), 73))
  }else if(spp_ranges_bosque[k,][1] == 1 & spp_ranges_bosque[k,][73] == 0){
    aa <- spp_ranges_bosque[k,][1:72] - spp_ranges_bosque[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="forest", Species=spp_list_bosque[k], start=c(1, which(aa == -1)+1), end=which(aa == 1))
  }else if(spp_ranges_bosque[k,][1] == 0 & spp_ranges_bosque[k,][73] == 1){
    aa <- spp_ranges_bosque[k,][1:72] - spp_ranges_bosque[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="forest", Species=spp_list_bosque[k], start=which(aa == -1)+1, end=c(which(aa == 1), 73))
  }else if(spp_ranges_bosque[k,][1] == 0 & spp_ranges_bosque[k,][73] == 0){
    aa <- spp_ranges_bosque[k,][1:72] - spp_ranges_bosque[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="forest", Species=spp_list_bosque[k], start=which(aa == -1)+1, end=which(aa == 1))
  }
}
spp_ranges_bosque_plot <- do.call(rbind, spp_ranges_plot)

spp_ranges_plot <- list()
for(k in 1:length(spp_list_interior)){
  if(spp_ranges_interior[k,][1] == 1 & spp_ranges_interior[k,][73] == 1){
    aa <- spp_ranges_interior[k,][1:72] - spp_ranges_interior[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="shade coffee", Species=spp_list_interior[k], start=c(1, which(aa == -1)+1), end=c(which(aa == 1), 73))
  }else if(spp_ranges_interior[k,][1] == 1 & spp_ranges_interior[k,][73] == 0){
    aa <- spp_ranges_interior[k,][1:72] - spp_ranges_interior[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="shade coffee", Species=spp_list_interior[k], start=c(1, which(aa == -1)+1), end=which(aa == 1))
  }else if(spp_ranges_interior[k,][1] == 0 & spp_ranges_interior[k,][73] == 1){
    aa <- spp_ranges_interior[k,][1:72] - spp_ranges_interior[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="shade coffee", Species=spp_list_interior[k], start=which(aa == -1)+1, end=c(which(aa == 1), 73))
  }else if(spp_ranges_interior[k,][1] == 0 & spp_ranges_interior[k,][73] == 0){
    aa <- spp_ranges_interior[k,][1:72] - spp_ranges_interior[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="shade coffee", Species=spp_list_interior[k], start=which(aa == -1)+1, end=which(aa == 1))
  }
}
spp_ranges_interior_plot <- do.call(rbind, spp_ranges_plot)

spp_ranges_plot <- list()
for(k in 1:length(spp_list_hamburgo)){
  if(spp_ranges_hamburgo[k,][1] == 1 & spp_ranges_hamburgo[k,][73] == 1){
    aa <- spp_ranges_hamburgo[k,][1:72] - spp_ranges_hamburgo[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="sun coffee", Species=spp_list_hamburgo[k], start=c(1, which(aa == -1)+1), end=c(which(aa == 1), 73))
  }else if(spp_ranges_hamburgo[k,][1] == 1 & spp_ranges_hamburgo[k,][73] == 0){
    aa <- spp_ranges_hamburgo[k,][1:72] - spp_ranges_hamburgo[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="sun coffee", Species=spp_list_hamburgo[k], start=c(1, which(aa == -1)+1), end=which(aa == 1))
  }else if(spp_ranges_hamburgo[k,][1] == 0 & spp_ranges_hamburgo[k,][73] == 1){
    aa <- spp_ranges_hamburgo[k,][1:72] - spp_ranges_hamburgo[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="sun coffee", Species=spp_list_hamburgo[k], start=which(aa == -1)+1, end=c(which(aa == 1), 73))
  }else if(spp_ranges_hamburgo[k,][1] == 0 & spp_ranges_hamburgo[k,][73] == 0){
    aa <- spp_ranges_hamburgo[k,][1:72] - spp_ranges_hamburgo[k,][2:73]
    spp_ranges_plot[[k]] <- data.frame(site="sun coffee", Species=spp_list_hamburgo[k], start=which(aa == -1)+1, end=which(aa == 1))
  }
}
spp_ranges_hamburgo_plot <- do.call(rbind, spp_ranges_plot)

spp_ranges_plot <- rbind(spp_ranges_bosque_plot, spp_ranges_interior_plot, spp_ranges_hamburgo_plot)

# Convert start and end times to date
period_dates <- seq.Date(as.Date("2022-02-21"), as.Date("2023-02-20"), 5)
spp_ranges_plot$start <- period_dates[spp_ranges_plot$start]
spp_ranges_plot$end <- period_dates[spp_ranges_plot$end]

spp_ranges_plot <- spp_ranges_plot %>%
  left_join(spp_data)

# forest specialists
forest_spp <- as.data.frame(spp_data)[which(spp_data$Habitat == "Forest" & spp_data$Habitat.Density == 1),]
forest_spp <- forest_spp$Species[-which(forest_spp$Species %in% c("Wilson's Warbler", "Magnolia Warbler", "Black-and-white Warbler", "Western Tanager", "Gray Hawk", "Roadside Hawk", "Rufous-tailed Hummingbird"))]
spp_ranges_plot$forest_spe <- if_else(spp_ranges_plot$Species %in% forest_spp, 1, 0)

# Plot figure
cols <- c("forest" = "purple2", "shade coffee" = "darkgreen", "sun coffee" = "chocolate1")
g_ranges <- ggplot(spp_ranges_plot, aes(x = start, y = reorder(reorder(Species, Mass), forest_spe))) + ylab("") + xlab("") +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_linerange(aes(xmin = start, xmax = end, col=site, group=site), size=1, position = position_dodge(width=0.7)) + theme_bw() +
  scale_color_manual(values = cols, name="Site") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=20), legend.position="bottom", legend.title = element_text(size=22))



##  Functional signatures  ##

spp_data$forest_specialist <- if_else(spp_data$Species %in% forest_spp, 1, 0)
spp_data$understory_specialist <- if_else(spp_data$Understory == 1 & spp_data$Mid_vegetation == 0 & spp_data$Canopy == 0, 1, 0)
spp_data$strata_generalist <- if_else(spp_data$Understory == 1 & spp_data$Mid_vegetation == 1 & spp_data$Canopy == 1, 1, 0)

# Bosque
spp_ranges_bosque_period <- apply(spp_ranges_bosque, 2, function(x) spp_list_bosque[which(x==1)])
functional_signatures_bosque <- vector()
for(i in 1:length(spp_ranges_bosque_period)){
  spp_data_subset <- spp_data[match(spp_ranges_bosque_period[[i]], spp_data$Species),] %>%
    select(Species, sci_name, Habitat, Habitat.Density, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle, forest_specialist, understory_specialist, strata_generalist)
  functional_signatures_bosque <- rbind(functional_signatures_bosque, 
                                        c(sum(spp_data_subset$forest_specialist) / nrow(spp_data_subset),
                                          sum(spp_data_subset$understory_specialist, na.rm=T) / nrow(spp_data_subset),
                                          sum(spp_data_subset$strata_generalist, na.rm=T) / nrow(spp_data_subset),
                                          length(which(spp_data_subset$Migration == 3)) / nrow(spp_data_subset),
                                          length(which(spp_data_subset$Trophic.Niche == "Invertivore")) / nrow(spp_data_subset),
                                          length(which(spp_data_subset$Trophic.Niche == "Omnivore")) / nrow(spp_data_subset),
                                          length(which(spp_data_subset$Migration == 3 & spp_data_subset$Trophic.Niche == "Invertivore")) / nrow(spp_data_subset)))
}
functional_signatures_bosque <- as.data.frame(functional_signatures_bosque)
colnames(functional_signatures_bosque) <- c("prop_forest_specialist", "prop_understory_specialist", "prop_strata_generalist", "prop_migrants", "prop_invertivores", "prop_omnivores", "prop_migrants_invertivores")
functional_signatures_bosque$site <- "forest"
functional_signatures_bosque$date <- seq.Date(from=as.Date("2022-02-23"), to=as.Date("2023-02-18"), by=5)

# Interior
spp_ranges_interior_period <- apply(spp_ranges_interior, 2, function(x) spp_list_interior[which(x==1)])
functional_signatures_interior <- vector()
for(i in 1:length(spp_ranges_interior_period)){
  spp_data_subset <- spp_data[match(spp_ranges_interior_period[[i]], spp_data$Species),] %>%
    select(Species, sci_name, Habitat, Habitat.Density, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle, forest_specialist, understory_specialist, strata_generalist)
  functional_signatures_interior <- rbind(functional_signatures_interior, 
                                          c(sum(spp_data_subset$forest_specialist) / nrow(spp_data_subset),
                                            sum(spp_data_subset$understory_specialist, na.rm=T) / nrow(spp_data_subset),
                                            sum(spp_data_subset$strata_generalist, na.rm=T) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Migration == 3)) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Trophic.Niche == "Invertivore")) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Trophic.Niche == "Omnivore")) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Migration == 3 & spp_data_subset$Trophic.Niche == "Invertivore")) / nrow(spp_data_subset)))
}
functional_signatures_interior <- as.data.frame(functional_signatures_interior)
colnames(functional_signatures_interior) <- c("prop_forest_specialist", "prop_understory_specialist", "prop_strata_generalist", "prop_migrants", "prop_invertivores", "prop_omnivores", "prop_migrants_invertivores")
functional_signatures_interior$site <- "shade coffee"
functional_signatures_interior$date <- seq.Date(from=as.Date("2022-02-23"), to=as.Date("2023-02-18"), by=5)

# Hamburgo
spp_ranges_hamburgo_period <- apply(spp_ranges_hamburgo, 2, function(x) spp_list_hamburgo[which(x==1)])
functional_signatures_hamburgo <- vector()
for(i in 1:length(spp_ranges_hamburgo_period)){
  spp_data_subset <- spp_data[match(spp_ranges_hamburgo_period[[i]], spp_data$Species),] %>%
    select(Species, sci_name, Habitat, Habitat.Density, Migration, Trophic.Level, Trophic.Niche, Primary.Lifestyle, forest_specialist, understory_specialist, strata_generalist)
  functional_signatures_hamburgo <- rbind(functional_signatures_hamburgo, 
                                          c(sum(spp_data_subset$forest_specialist) / nrow(spp_data_subset),
                                            sum(spp_data_subset$understory_specialist, na.rm=T) / nrow(spp_data_subset),
                                            sum(spp_data_subset$strata_generalist, na.rm=T) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Migration == 3)) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Trophic.Niche == "Invertivore")) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Trophic.Niche == "Omnivore")) / nrow(spp_data_subset),
                                            length(which(spp_data_subset$Migration == 3 & spp_data_subset$Trophic.Niche == "Invertivore")) / nrow(spp_data_subset)))
}
functional_signatures_hamburgo <- as.data.frame(functional_signatures_hamburgo)
colnames(functional_signatures_hamburgo) <- c("prop_forest_specialist", "prop_understory_specialist", "prop_strata_generalist", "prop_migrants", "prop_invertivores", "prop_omnivores", "prop_migrants_invertivores")
functional_signatures_hamburgo$site <- "sun coffee"
functional_signatures_hamburgo$date <- seq.Date(from=as.Date("2022-02-23"), to=as.Date("2023-02-18"), by=5)

# Figure 5: temporal trends in functional signatures 
functional_signatures <- rbind(functional_signatures_bosque, functional_signatures_interior, functional_signatures_hamburgo)
functional_signatures$season <- if_else(functional_signatures$date >= as.Date("2022-05-10") & functional_signatures$date < as.Date("2022-09-15"), "rain", if_else(functional_signatures$date >= as.Date("2022-04-01") & functional_signatures$date < as.Date("2022-05-10"), "spring", if_else(functional_signatures$date >= as.Date("2022-09-15") & functional_signatures$date < as.Date("2022-11-01"), "fall", "dry")))

functional_signatures_boxplot <- functional_signatures %>%
  filter(season %in% c("dry", "rain"))

cols <- c("forest" = "purple2", "shade coffee" = "darkgreen", "sun coffee" = "chocolate1")
g_forest <- ggplot(data=functional_signatures) +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_point(aes(x=date, y=prop_forest_specialist, col=site), size=1) +
  geom_smooth(aes(x=date, y=prop_forest_specialist, col=site), se=F, span=0.3) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Proportion of forest specialists") + theme_light() +
  theme(legend.position="none")

g_forest_boxplot <- ggplot(data=functional_signatures_boxplot) +
  geom_boxplot(aes(x=site, y=prop_forest_specialist, fill = season)) +
  xlab("") + ylab("Proportion of forest specialists") + theme_light() +
  theme(legend.position="none")

g_under <- ggplot(data=functional_signatures) +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_point(aes(x=date, y=prop_understory_specialist, col=site), size=1) +
  geom_smooth(aes(x=date, y=prop_understory_specialist, col=site), se=F, span=0.3) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Proportion of understory specialists") + theme_light() +
  theme(legend.position="none")

g_under_boxplot <- ggplot(data=functional_signatures_boxplot) +
  geom_boxplot(aes(x=site, y=prop_understory_specialist, fill = season)) +
  xlab("") + ylab("Proportion of understory specialists") + theme_light() +
  theme(legend.position="none")

g_strata <- ggplot(data=functional_signatures) +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_point(aes(x=date, y=prop_strata_generalist, col=site), size=1) +
  geom_smooth(aes(x=date, y=prop_strata_generalist, col=site), se=F, span=0.3) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Proportion of strata generalists") + theme_light() +
  theme(legend.position="none")

g_strata_boxplot <- ggplot(data=functional_signatures_boxplot) +
  geom_boxplot(aes(x=site, y=prop_strata_generalist, fill = season)) +
  xlab("") + ylab("Proportion of strata generalists") + theme_light() +
  theme(legend.position="none")

g_migr <- ggplot(data=functional_signatures) +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_point(aes(x=date, y=prop_migrants, col=site), size=1) +
  geom_smooth(aes(x=date, y=prop_migrants, col=site), se=F, span=0.25) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Proportion of migrants") + theme_light() +
  theme(legend.position="none")

g_migr_boxplot <- ggplot(data=functional_signatures_boxplot) +
  geom_boxplot(aes(x=site, y=prop_migrants, fill = season)) +
  xlab("") + ylab("Proportion of migrants") + theme_light() +
  theme(legend.position="none")

g_inv <- ggplot(data=functional_signatures) +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_point(aes(x=date, y=prop_invertivores, col=site), size=1) +
  geom_smooth(aes(x=date, y=prop_invertivores, col=site), se=F, span=0.25) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Proportion of invertivores") + theme_light() +
  theme(legend.position="none")

g_inv_boxplot <- ggplot(data=functional_signatures_boxplot) +
  geom_boxplot(aes(x=site, y=prop_invertivores, fill = season)) +
  xlab("") + ylab("Proportion of invertivores") + theme_light() +
  theme(legend.position="none")

g_omn <- ggplot(data=functional_signatures) +
  geom_rect(xmin=as.Date("2022-05-10"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="azure2", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-04-01"), xmax=as.Date("2022-05-10"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_rect(xmin=as.Date("2022-09-15"), xmax=as.Date("2022-10-31"), ymin=0, ymax=Inf, fill="antiquewhite", alpha = 0.5) +
  geom_point(aes(x=date, y=prop_omnivores, col=site), size=1) +
  geom_smooth(aes(x=date, y=prop_omnivores, col=site), se=F, span=0.25) +
  scale_color_manual(values = cols) +
  xlab("") + ylab("Proportion of omnivores") + theme_light() +
  theme(legend.position="bottom")

g_omn_boxplot <- ggplot(data=functional_signatures_boxplot) +
  geom_boxplot(aes(x=site, y=prop_omnivores, fill = season)) +
  xlab("") + ylab("Proportion of omnivores") + theme_light() +
  theme(legend.position="bottom")


png(filename = "results/Figures/Figure_5a.png", width = 400, height = 1200, units = "px", pointsize = 12, bg = "white")
ggarrange(g_forest, g_under, g_strata, g_migr, g_inv, g_omn, nrow=6, ncol=1, legend="bottom", common.legend=T)
dev.off()
png(filename = "results/Figures/Figure_5b.png", width = 250, height = 1200, units = "px", pointsize = 12, bg = "white")
ggarrange(g_forest_boxplot, g_under_boxplot, g_strata_boxplot, g_migr_boxplot, g_inv_boxplot, g_omn_boxplot, nrow=6, ncol=1, legend="bottom", common.legend=T)
dev.off()


# T-tests
functional_signatures_boxplot_forest <- functional_signatures_boxplot %>% filter(site=="forest")
functional_signatures_boxplot_shade <- functional_signatures_boxplot %>% filter(site=="shade coffee")
functional_signatures_boxplot_sun <- functional_signatures_boxplot %>% filter(site=="sun coffee")

t.test(prop_forest_specialist ~ season, data=functional_signatures_boxplot_forest)
t.test(prop_forest_specialist ~ season, data=functional_signatures_boxplot_shade)

t.test(prop_strata_generalist ~ season, data=functional_signatures_boxplot_forest)
t.test(prop_strata_generalist ~ season, data=functional_signatures_boxplot_shade)
t.test(prop_strata_generalist ~ season, data=functional_signatures_boxplot_sun)

t.test(prop_invertivores ~ season, data=functional_signatures_boxplot_forest)
t.test(prop_invertivores ~ season, data=functional_signatures_boxplot_shade)
t.test(prop_invertivores ~ season, data=functional_signatures_boxplot_sun)


