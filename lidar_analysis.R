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
library(lidR)
library(leafR)
library(rGEDI)


# Figure 1 — map of study area

#study_area <- st_bbox(c(xmin = -92.35, xmax = -92.321, ymax = 15.182, ymin = 15.161), crs = "+proj=longlat +datum=WGS84 +no_defs")
study_area1 <- st_bbox(c(xmin = -92.37, xmax = -92.31, ymax = 15.19, ymin = 15.15), crs = "+proj=longlat +datum=WGS84 +no_defs")
study_area2 <- st_bbox(c(xmin = -93.0, xmax = -92.0, ymax = 15.65, ymin = 15.0), crs = "+proj=longlat +datum=WGS84 +no_defs")

# Get elevation data
srtm <- terra::rast("resources/elevation/N15W093.hgt")
#srtm_bb <- st_bbox(ext(srtm), crs = "+proj=longlat +datum=WGS84 +no_defs")

# Get canopy height
canopy_height <- terra::rast("resources/canopy_height/ETH_GlobalCanopyHeight_10m_2020_N15W093_Map.tif")
canopy_height <- terra::crop(canopy_height, terra::ext(study_area1))

Buggs <- data.frame(
  site = c("shade coffee", "edge", "restoration", "forest", "sun coffee"),
  lon = c(-92.34371, -92.33772, -92.33588, -92.33578, -92.33044),
  lat = c(15.16741, 15.16983, 15.17079, 15.16844, 15.17186)
) %>%
  filter(site %in% c("shade coffee", "forest", "sun coffee"))
Buggs_sf <- st_as_sf(Buggs, coords = c("lon","lat"), crs="+proj=longlat +datum=WGS84 +no_defs") 

mexico_guat <- terra::vect(ne_countries(country=c("Mexico", "Guatemala", "Belize", "El Salvador", "Honduras", "Nicaragua"), returnclass = "sf") %>% st_combine)
crs(mexico_guat) <- crs(srtm)
g_mexico <- ggplot() + geom_spatvector(data = mexico_guat) +
  theme_void() + xlim(c(-105, -87.5)) + ylim(c(13, 21.5)) +
  geom_sf(data=st_as_sfc(study_area2), col="red", fill=NA, lwd=1.5)

pdf(file = "results/Figures/Fig1_b.pdf", width = 10, height = 5)
g_mexico
dev.off()

mexico_guat_elev <- terra::mask(srtm, mexico_guat)
mexico_guat_elev <- terra::crop(mexico_guat_elev, ext(study_area2))
g_elev <- ggplot() + geom_spatraster(data = mexico_guat_elev) +
  scale_fill_viridis_c(option = "cividis", name = "Elevation") +
  theme_void() +
  geom_sf(data=st_as_sfc(study_area1), col="red", fill=NA, lwd=1.5) +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=13.5))

pdf(file = "results/Figures/Fig1_c.pdf", width = 8, height = 4.7)
g_elev
dev.off()

cols <- c("forest" = "purple2", "shade coffee" = "darkgreen", "sun coffee" = "chocolate1")
g_map <- ggplot() +
  geom_spatraster(data = canopy_height) +
  scale_fill_viridis_c(option = "viridis", name = "Canopy height") +
  geom_sf(data = Buggs_sf, aes(col=site), shape=1, size=4, stroke = 2) +
  scale_color_manual(values = cols, name = "Acoustic monitoring site") +
  theme_void() + annotation_scale(location = "tl", pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in")) +
  theme(legend.title = element_text(size=13), legend.text = element_text(size=12))

pdf(file = "results/Figures/Fig1_d.pdf", width = 8, height = 5)
g_map
dev.off()





# Analysis of handheld lidar

sites <- c("bosque", "irlanda_interior", "hamburgo")

## Handheld lidar data obtained from the location where the acoustic devices were deployed

# Height normalisation
mycsf <- csf(sloop_smooth = TRUE)
for(i in 1:length(sites)){
  # 2022
  las_file_path <- list.files(paste0("resources/lidar/Bugg/", sites[i]), pattern=".las")
  PC <- readLAS(paste0("resources/lidar/Bugg/", sites[i], "/", las_file_path))
  PC <- classify_ground(PC, mycsf)     
  PC <- normalize_height(PC, algorithm = knnidw(k = 10, p = 2))
  PC <- PC[which(PC@data$Classification != 2)] # remove ground points
  PC <- clip_circle(PC, radius = 25, xcenter = 0, ycenter = 0) # Clip to 25m
  writeLAS(PC, paste0("resources/lidar/Bugg_heightNorm/", sites[i], "_heightNorm_noGround.las"))
  # 2023
  las_file_path <- list.files(paste0("resources/lidar/Bugg_2023/", sites[i]), pattern=".las")
  PC <- readLAS(paste0("resources/lidar/Bugg_2023/", sites[i], "/", las_file_path))
  PC <- classify_ground(PC, mycsf)     
  PC <- normalize_height(PC, algorithm = knnidw(k = 10, p = 2))
  PC <- PC[which(PC@data$Classification != 2)] # remove ground points
  PC <- clip_circle(PC, radius = 25, xcenter = 0, ycenter = 0) # Clip to 25m
  writeLAS(PC, paste0("resources/lidar/Bugg_2023_heightNorm/", sites[i], "_2023_heightNorm_noGround.las"))
  print(i)
}

# Calculate vegetation metrics
lai_buggs_2022 <- fhd_buggs_2022 <- gini_buggs_2022 <- vector()
lai_buggs_2023 <- fhd_buggs_2023 <- gini_buggs_2023 <- vector()
LAI_rasters_buggs_2022 <- LAI_rasters_buggs_2023 <- list()
for(i in 1:length(sites)){
  # 2022
  VOXELS_LAD <- lad.voxels(paste0("Resources/lidar/Bugg_heightNorm/",sites[i],"_heightNorm_noGround.las"), grain.size=1)
  LAI_rasters_buggs_2022[[i]] <- lai.raster(VOXELS_LAD)
  lad_profile <- lad.profile(VOXELS_LAD)
  fhd_buggs_2022[i] <- FHD(lad_profile)
  gini_buggs_2022[i] <- GC(paste0("Resources/lidar/Bugg_heightNorm/",sites[i],"_heightNorm_noGround.las"))
  lai_0_5 <- lai(lad_profile, min = 0, max = 5)
  lai_5_10 <- lai(lad_profile, min = 5, max = 10)
  lai_10_20 <- lai(lad_profile, min = 10, max = 20)
  lai_20_50 <- lai(lad_profile, min = 20, max = 50)
  lai_buggs_2022 <- rbind(lai_buggs_2022, c(lai_0_5,lai_5_10,lai_10_20,lai_20_50)) 
  # 2023
  VOXELS_LAD <- lad.voxels(paste0("Resources/lidar/Bugg_2023_heightNorm/",sites[i],"_2023_heightNorm_noGround.las"), grain.size=1)
  LAI_rasters_buggs_2023[[i]] <- lai.raster(VOXELS_LAD)
  lad_profile <- lad.profile(VOXELS_LAD)
  fhd_buggs_2023[i] <- FHD(lad_profile)
  gini_buggs_2023[i] <- GC(paste0("Resources/lidar/Bugg_2023_heightNorm/",sites[i],"_2023_heightNorm_noGround.las"))
  lai_0_5 <- lai(lad_profile, min = 0, max = 5)
  lai_5_10 <- lai(lad_profile, min = 5, max = 10)
  lai_10_20 <- lai(lad_profile, min = 10, max = 20)
  lai_20_50 <- lai(lad_profile, min = 20, max = 50)
  lai_buggs_2023 <- rbind(lai_buggs_2023, c(lai_0_5,lai_5_10,lai_10_20,lai_20_50)) 
  print(i)
}

vegetation_metrics_buggs <- data.frame(site = rep(rep(sites,4), 2),
                                       year = rep(c("2022", "2023"), each=nrow(lai_buggs_2022)*4),
                                       height = rep(c("understory\n(0-5)", "mid height\n(5-10)", "canopy\n(10-20)", "high canopy\n(20+)"), each=nrow(lai_buggs_2022)),
                                       lai = c(lai_buggs_2022[,1], lai_buggs_2022[,2], lai_buggs_2022[,3], lai_buggs_2022[,4], lai_buggs_2023[,1], lai_buggs_2023[,2], lai_buggs_2023[,3], lai_buggs_2023[,4])
)
vegetation_metrics_buggs$site[vegetation_metrics_buggs$site == "bosque"] <- "forest"
vegetation_metrics_buggs$site[vegetation_metrics_buggs$site == "irlanda_interior"] <- "shade coffee"
vegetation_metrics_buggs$site[vegetation_metrics_buggs$site == "hamburgo"] <- "sun coffee"


## Handheld lidar data obtained from the location where the acoustic devices were deployed

# Height normalisation
mycsf <- csf(sloop_smooth = TRUE)
for(k in 1:45){
  lf <- list.files(paste0("resources/lidar/point_counts/PC", k))
  las_file_path <- lf[unlist(lapply(strsplit(lf, "[.]"), function(x) x[2])) == "las"]
  PC <- readLAS(paste0("resources/lidar/point_counts/PC", k, "/", las_file_path))
  PC <- classify_ground(PC, mycsf)     
  PC <- normalize_height(PC, algorithm = knnidw(k = 10, p = 2))
  PC <- PC[which(PC@data$Classification != 2)] # remove ground points
  PC <- clip_circle(PC, radius = 25, xcenter = 0, ycenter = 0) # Clip to 25m
  writeLAS(PC, paste0("resources/lidar/point_counts_heightNorm/PC", k, "_heightNorm_noGround.las"))
  print(k)
}

# Calculate vegetation metrics
PCs <- paste0("PC", 1:45)
lai_PCs <- fhd_PCs <- gini_PCs <- ch_PCs <- vector()
LAI_rasters_PCs <- list()
for(i in 1:length(PCs)){
  VOXELS_LAD <- lad.voxels(paste0("resources/lidar/point_counts_heightNorm/",PCs[i],"_heightNorm_noGround.las"), grain.size=1)
  LAI_rasters_PCs[[i]] <- lai.raster(VOXELS_LAD)
  lad_profile <- lad.profile(VOXELS_LAD)
  ch_PCs[i] <- max(lad_profile$height[which(lad_profile$lad > 0.01)])
  fhd_PCs[i] <- FHD(lad_profile)
  gini_PCs[i] <- GC(paste0("resources/lidar/point_counts_heightNorm/PC", i, "_heightNorm_noGround.las"))
  lai_0_5 <- lai(lad_profile, min = 0, max = 5)
  lai_5_10 <- lai(lad_profile, min = 5, max = 10)
  lai_10_20 <- lai(lad_profile, min = 10, max = 20)
  lai_20_50 <- lai(lad_profile, min = 20, max = 50)
  lai_PCs <- rbind(lai_PCs, c(lai_0_5,lai_5_10,lai_10_20,lai_20_50)) 
  print(i)
}

PC_data <- read.csv("resources/Point_counts.csv")
PC_data <- PC_data[duplicated(PC_data$Point_count) == F,] %>%
  select(Point_count, Zone, Lat, Lon, Canopy_cover) %>%
  filter(Point_count <= 45) %>%
  rename(site = Point_count)

vegetation_metrics_PC <- data.frame(site = rep(1:nrow(lai_PCs),4),
                                    height = rep(c("understory\n(0-5)", "mid height\n(5-10)", "canopy\n(10-20)", "high canopy\n(20+)"), each=nrow(lai_PCs)),
                                    lai = c(lai_PCs[,1], lai_PCs[,2], lai_PCs[,3], lai_PCs[,4]),
                                    fhd = rep(fhd_PCs, 4),
                                    ch = rep(ch_PCs, 4)
) %>%
  left_join(PC_data)

vegetation_metrics_PC$Zone[vegetation_metrics_PC$Zone == "bosque" | vegetation_metrics_PC$Zone == "restoration"] <- "forest"
vegetation_metrics_PC$Zone[vegetation_metrics_PC$Zone == "hamburgo"] <- "sun coffee"
vegetation_metrics_PC$Zone[vegetation_metrics_PC$Zone == "irlanda" ] <- "shade coffee"

# Plot the figures (Figure 2)

#vegetation_metrics_PC <- read.csv("results/vegetation_metrics_PC.csv")
#vegetation_metrics_buggs <- read.csv("results/vegetation_metrics_buggs.csv")

cols <- c("forest" = "purple2", "shade coffee" = "darkgreen", "sun coffee" = "chocolate1")
g_lai_lidar <- ggplot() +
  geom_boxplot(data = vegetation_metrics_PC, aes(x=lai, y=reorder(height,-lai), fill=Zone)) +
  scale_fill_manual(values = cols) +
  geom_boxplot(data=vegetation_metrics_buggs, aes(x=lai, y=reorder(height,-lai), group=interaction(site,height)), col="grey45", alpha=0) +
  labs(y="Canopy height (m)", x=bquote('LAI'~(m^2/m^2))) +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(a)")

g_fhd_PC <- ggplot() +
  geom_boxplot(data = vegetation_metrics_PC, aes(y=fhd, x=Zone, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="Ground LiDAR foliage height diversity", x="") +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(c)")

g_cover_obs <- ggplot() +
  geom_boxplot(data = vegetation_metrics_PC, aes(y=Canopy_cover/100, x=Zone, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="Canopy cover visual estimate", x="") +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(d)")

g_ch_PC <- ggplot() +
  geom_boxplot(data = vegetation_metrics_PC, aes(y=ch*100, x=Zone, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="Ground LiDAR canopy height (cm)", x="") +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(e)")




## GEDI analysis

# Study area boundary box coordinates
ul_lat <- 15.182
lr_lat <- 15.161
ul_lon <- -92.35
lr_lon <- -92.321

# Specifying the date range
daterange=c("2019-02-20","2024-02-20")

# Get path to GEDI data
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="002",daterange=daterange)

# Set output dir for downloading the files
outdir=paste0(getwd(),"/resources/GEDI")

# Downloading GEDI data
gediDownload(filepath=gLevel2B,outdir=outdir)

study_area <- st_bbox(c(xmin = -92.35, xmax = -92.321, ymax = 15.182, ymin = 15.161), crs = "+proj=longlat +datum=WGS84 +no_defs")

gedi_mexico_pai <- gedi_mexico_fhd <- list()
for(k in 1:length(gLevel2B)){
  
  # Reading GEDI data
  gedilevel2b <- readLevel2B(level2Bpath = paste0(outdir, "/", list.files(outdir))[k])
  
  # Process GEDI data
  level2B_pai <- getLevel2BPAIProfile(gedilevel2b) %>% left_join(getLevel2BPAVDProfile(gedilevel2b))
  level2B_fhd <- getLevel2BVPM(gedilevel2b)
  
  # Converting shot_number as "integer64" to "character"
  level2B_pai$shot_number<-as.character(level2B_pai$shot_number)
  level2B_fhd$shot_number<-as.character(level2B_fhd$shot_number)
  
  # Converting GEDI Vegetation Profile Biophysical Variables as data.table to SpatialPointsDataFrame
  level2B_pai_sf <- st_as_sf(SpatialPointsDataFrame(cbind(level2B_pai$lon_lowestmode, level2B_pai$lat_lowestmode), data=level2B_pai))
  level2B_fhd_sf <- st_as_sf(SpatialPointsDataFrame(cbind(level2B_fhd$longitude_bin0, level2B_fhd$latitude_bin0), data=level2B_fhd))
  st_crs(level2B_pai_sf) <- st_crs(level2B_fhd_sf) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # Extracting values just within the study area
  gedi_mex_pai <- st_within(level2B_pai_sf, st_as_sfc(study_area), sparse=F)
  gedi_mex_fhd <- st_within(level2B_fhd_sf, st_as_sfc(study_area), sparse=F)
  gedi_mexico_pai[[k]] <- level2B_pai_sf[which(gedi_mex_pai[,1] == T),]
  gedi_mexico_fhd[[k]] <- level2B_fhd_sf[which(gedi_mex_fhd[,1] == T),]
}
gedi_mexico_pai <- bind_rows(gedi_mexico_pai)
gedi_mexico_pai <- gedi_mexico_pai %>% filter(height_bin0 > 0)
gedi_mexico_fhd <- bind_rows(gedi_mexico_fhd)
gedi_mexico_fhd <- gedi_mexico_fhd %>% filter(cover >= 0)



# Figure 1 — map of lidar and gedi

PC_data$Zone[PC_data$Zone == "bosque" | PC_data$Zone == "restoration"] <- "forest"
PC_data$Zone[PC_data$Zone == "hamburgo"] <- "sun coffee"
PC_data$Zone[PC_data$Zone == "irlanda" ] <- "shade coffee"
PC_data_sf = st_as_sf(PC_data, coords = c("Lon","Lat"), remove = FALSE, crs="+proj=longlat +datum=WGS84 +no_defs")

hulls <- PC_data_sf %>%
  group_by(Zone) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_convex_hull()

gedi_mexico_pai_hulls <- st_within(gedi_mexico_pai, hulls, sparse=F)
gedi_mexico_pai_hulls_forest <- gedi_mexico_pai[which(gedi_mexico_pai_hulls[,1] == T),] %>% mutate(Zone = "forest")
gedi_mexico_pai_hulls_shade <- gedi_mexico_pai[which(gedi_mexico_pai_hulls[,2] == T),] %>% mutate(Zone = "shade coffee")
gedi_mexico_pai_hulls_sun <- gedi_mexico_pai[which(gedi_mexico_pai_hulls[,3] == T),] %>% mutate(Zone = "sun coffee")
gedi_mexico_pai_hulls <- rbind(gedi_mexico_pai_hulls_forest, gedi_mexico_pai_hulls_shade, gedi_mexico_pai_hulls_sun)

gedi_mexico_fhd_hulls <- st_within(gedi_mexico_fhd, hulls, sparse=F)
gedi_mexico_fhd_hulls_forest <- gedi_mexico_fhd[which(gedi_mexico_fhd_hulls[,1] == T),] %>% mutate(Zone = "forest")
gedi_mexico_fhd_hulls_shade <- gedi_mexico_fhd[which(gedi_mexico_fhd_hulls[,2] == T),] %>% mutate(Zone = "shade coffee")
gedi_mexico_fhd_hulls_sun <- gedi_mexico_fhd[which(gedi_mexico_fhd_hulls[,3] == T),] %>% mutate(Zone = "sun coffee")
gedi_mexico_fhd_hulls <- rbind(gedi_mexico_fhd_hulls_forest, gedi_mexico_fhd_hulls_shade, gedi_mexico_fhd_hulls_sun)

#cols <- c("forest" = "darkgreen", "shade coffee" = "olivedrab3", "sun coffee" = "chocolate1")
cols <- c("forest" = "purple2", "shade coffee" = "darkgreen", "sun coffee" = "chocolate1")
g_map_lidar <- ggplot() +
  geom_spatraster(data = canopy_height) +
  geom_sf(data=hulls, col=cols, fill=NA, lwd=1) +
  geom_sf(data = gedi_mexico_pai_hulls, col="gray2") +
  geom_sf(data = PC_data_sf, aes(col=Zone), size = 1.2, fill = NA, shape = 21, stroke = 1.5) +
  geom_sf(data = Buggs_sf, shape=18, size=6, stroke = 1.5) +
  labs( x = "Longitude", y = "Latitude") +
  scale_color_manual(values = cols) +
  scale_fill_viridis_c(option = "viridis", name = "Canopy height") +
  theme_void() + annotation_scale(location = "tl", pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in")) +
  theme(text=element_text(size=14), plot.title=element_text(size=20), legend.position = "right") 

pdf(file = "results/Figures/Figure_1_d_new.pdf", width = 7.5, height = 5)
g_map_lidar
dev.off()


# Plots for figure 2
vegetation_metrics_gedi <- data.frame(Zone = rep(gedi_mexico_pai_hulls$Zone,4),
                                      height = rep(c("understory\n(0-5)", "mid height\n(5-10)", "canopy\n(10-20)", "high canopy\n(20+)"), each=nrow(gedi_mexico_pai_hulls)),
                                      pai = c(gedi_mexico_pai_hulls$pai_z0_5m, gedi_mexico_pai_hulls$pai_z5_10m, (gedi_mexico_pai_hulls$pai_z10_15m + gedi_mexico_pai_hulls$pai_z15_20m)/2, (gedi_mexico_pai_hulls$pai_z20_25m + gedi_mexico_pai_hulls$pai_z25_30m + gedi_mexico_pai_hulls$pai_z30_35m + gedi_mexico_pai_hulls$pai_z35_40m + gedi_mexico_pai_hulls$pai_z40_45m + gedi_mexico_pai_hulls$pai_z45_50m)/6)
)

g_pai_gedi <- ggplot() +
  geom_boxplot(data = vegetation_metrics_gedi, aes(x=pai, y=reorder(height,-pai), fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="", x=bquote('PAI'~(m^2/m^2))) +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "right") +
  ggtitle("(b)")

pdf(file = "results/Figures/Figure_2_a.pdf", width = 15, height = 7)
ggarrange(g_lai_lidar, g_pai_gedi, nrow=1, ncol=2, common.legend = TRUE)
dev.off()

vegetation_metrics_gedi_fhd <- gedi_mexico_fhd_hulls %>%
  mutate(Zone = gedi_mexico_fhd_hulls$Zone)

g_fhd_gedi <- ggplot() +
  geom_boxplot(data = vegetation_metrics_gedi_fhd, aes(y=fhd_normal * pai, x=Zone, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="GEDI foliage height diversity", x="") +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(f)")

g_cover_gedi <- ggplot() +
  geom_boxplot(data = vegetation_metrics_gedi_fhd, aes(y=cover, x=Zone, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="GEDI canopy cover", x="") +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(g)")

g_rh100_gedi <- ggplot() +
  geom_boxplot(data = vegetation_metrics_gedi_fhd, aes(y=rh100, x=Zone, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="GEDI canopy height (cm)", x="") +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  ggtitle("(h)")

pdf(file = "results/Figures/Figure_2_b.pdf", width = 15, height = 10)
ggarrange(g_fhd_PC, g_cover_obs, g_ch_PC,
          g_fhd_gedi, g_cover_gedi, g_rh100_gedi, nrow=2, ncol=3)
dev.off()


# Pairwise t tests for handheld
vegetation_metrics_PC %>% filter(height == "understory\n(0-5)") %>%
  pairwise_t_test(lai ~ Zone, pool.sd = FALSE, p.adjust.method = "bonferroni")
vegetation_metrics_PC %>% filter(height == "mid height\n(5-10)") %>%
  pairwise_t_test(lai ~ Zone, pool.sd = FALSE, p.adjust.method = "bonferroni")
vegetation_metrics_PC %>% filter(height == "canopy\n(10-20)") %>%
  pairwise_t_test(lai ~ Zone, pool.sd = FALSE, p.adjust.method = "bonferroni")

# Pairwise t tests for GEDI
vegetation_metrics_gedi %>% filter(height == "understory\n(0-5)") %>%
  pairwise_t_test(pai ~ Zone, pool.sd = FALSE, p.adjust.method = "bonferroni")
vegetation_metrics_gedi %>% filter(height == "mid height\n(5-10)") %>%
  pairwise_t_test(pai ~ Zone, pool.sd = FALSE, p.adjust.method = "bonferroni")
vegetation_metrics_gedi %>% filter(height == "canopy\n(10-20)") %>%
  pairwise_t_test(pai ~ Zone, pool.sd = FALSE, p.adjust.method = "bonferroni")

# Correlation with canopy cover
vegetation_metrics_PC %>% filter(height == "understory\n(0-5)") %>%
  ggplot() + geom_point(aes(x=Canopy_cover, y=lai, col=Zone))
vegetation_metrics_PC %>% filter(height != "understory\n(0-5)") %>%
  ggplot() + geom_point(aes(x=Canopy_cover, y=lai, col=Zone))

ggplot() +
  geom_boxplot(data = vegetation_metrics_PC, aes(x=Canopy_cover, fill=Zone)) +
  scale_fill_manual(values = cols) +
  labs(y="", x=bquote('PAI'~(m^2/m^2))) +
  theme_classic() + theme(text=element_text(size=20), plot.title=element_text(size=20), legend.position = "none") +
  
  
  plot(vegetation_metrics_PC_mid$Canopy_cover, vegetation_metrics_PC_mid$lai)































