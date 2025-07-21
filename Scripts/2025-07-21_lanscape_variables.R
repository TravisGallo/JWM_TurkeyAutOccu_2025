source("./Scripts/TurkeyMultiOccu_package_load.R")

## Calculating landscape variables


# Create spatial sites ----------------------------------------------------


# load sites
sites <- read.csv("./Data/TurkeyOccupancyReport.csv") |>
  select(Site, Latitude, Longitude) |>
  distinct(Latitude, Longitude, .keep_all = TRUE)

# turn sites spatial and convert to UTM
points <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(crs = 26918)
mapview(points)


# Read in Land Cover ------------------------------------------------------
# https://www.chesapeakeconservancy.org/projects/cbp-land-use-land-cover-data-project

# # make a vector of file names
# files <- list.files("./Data/GIS/LandCover/", recursive = TRUE, pattern = ".tif$")
# 
# # create a list object to store the rasters
# rast_list <- vector("list", length(files))
# 
# for(i in 1:length(rast_list)){
#   rast_list[[i]] <- rast(paste0("./Data/GIS/",files[i]))
# }
# 
# # make a spatrast collection
# rsrc <- sprc(rast_list)
# 
# # merge the raster
# lc_rast <- merge(rsrc)
# 
# # save the raster
# writeRaster(lc_rast, filename = "./Data/GIS/study_area_landcover.tif")

# load study area raster
lc <- rast("./Data/GIS/study_area_landcover.tif")


# Distance to Variables ---------------------------------------------------


## Distance to Road, Water, and Wetland -----------------------------------

# reproject points to match raster
points_rp <- st_transform(points, st_crs(lc))

# create empty list to save items in
dist2road <- rep(NA, nrow(points_rp))
dist2water <- rep(NA, nrow(points_rp))
dist2wetland <- rep(NA, nrow(points_rp))


# loop through each site. Remember sites (i.e. features) are in rows. So this is
# doing the processing one site at a time, and the results for that site are 
# stored as one element 
for(i in 1:nrow(points_rp)){
  # create a buffer around the site
  site_buff <- st_buffer(points_rp[i,], 5000)
  # crop the big raster down to the same size as the buffer
  r_reduced <- crop(lc, site_buff)
  # calculate the distance from the point to all cells in the reduced raster
  d <- distance(r_reduced, points_rp[i,])
  # calculate the minimum distance to each category
  # these results are in meters
  min_dist <- zonal(d, r_reduced, fun='min')
  # find the minimum to roads
  # Appendix C: Land Use/Land Cover
  # https://cicwebresources.blob.core.windows.net/docs/LU_Classification_Methods_2017_2018.pdf
  dist2road[i] <- min_dist[which(min_dist[,1] == "Impervious Roads"), 2]
  # find the minimum to wetlands
  dist2wetland[i] <- min(min_dist[which(min_dist[,1] == "Wetlands, Riverine Non-forested" |
                                          min_dist[,1] == "Wetlands, Terrene Non-forested" |
                                          min_dist[,1] == "Wetlands, Tidal Non-forested"), 2],
                         na.rm = TRUE)
  # find minimum distance to water
  dist2water[i] <- min(min_dist[which(min_dist[,1] == "Water"), 2],
                       na.rm = TRUE)
  
}


## Distance to Linear Water (stream) ------------------------------------------------

# load in rivers and stream layer
#https://hub.arcgis.com/datasets/esri::usa-rivers-and-streams/explore
rivers_streams <- st_read("./Data/GIS/USA_Rivers_and_Streams-shp/",
                layer = "9ae73184-d43c-4ab8-940a-c8687f61952f2020328-1-r9gw71.0odx9")

points_rp <- st_transform(points, crs = st_crs(rivers_streams))

near_tmp = st_nearest_feature(points_rp, rivers_streams)

dist2stream = st_distance(points_rp, rivers_streams[near_tmp,], by_element=TRUE) |>
  drop_units()


## Distance to Trails ------------------------------------------------------

# https://download.geofabrik.de/north-america/us.html

dc_osm <- st_read("./Data/GIS/district-of-columbia-latest-free.shp/",
                  layer = "gis_osm_roads_free_1") |>
  filter(fclass == "cycleway" | fclass == "footway")

md_osm <- st_read("./Data/GIS/maryland-latest-free.shp/",
                     layer = "gis_osm_roads_free_1") |>
  filter(fclass == "cycleway" | fclass == "footway")

va_osm <- st_read("./Data/GIS/virginia-latest-free.shp/",
                  layer = "gis_osm_roads_free_1") |>
  filter(fclass == "cycleway" | fclass == "footway")

comb_osm <- rbind(dc_osm, md_osm, va_osm)

points_rp <- st_transform(points, crs = st_crs(comb_osm))

nearest = st_nearest_feature(points_rp, comb_osm)

dist2trail = st_distance(points_rp, comb_osm[nearest,], by_element=TRUE) |>
  drop_units()


# Multi-scale variables ---------------------------------------------------

# create a list of buffers at three different scales.
buff_list <- list(buff_4000 = st_buffer(points, dist = 4000),
                  buff_2000 = st_buffer(points, dist = 2000),
                  buff_1000 = st_buffer(points, dist = 1000))

## Entropy -----------------------------------------------------------------

buffs_rp <- lapply(buff_list, function(x){
  st_transform(x, crs = st_crs(lc))
})

extractEntropy <- function(x){     
  
  entropy <- rep(NA, nrow(x))
  
  for(i in 1:nrow(x)){
    # crop the big raster down to the same size as the buffer
    r_reduced <- crop(lc, x[i,])
    
    # reclassify natural areas into their own category and everything else as NA
    # FOREST
    r_reduced[r_reduced == 41 | r_reduced == 65 | r_reduced == 75 | r_reduced == 95] <- 99
    # Natural Succession
    r_reduced[r_reduced == 16 | r_reduced == 54 | r_reduced == 55 | r_reduced == 56] <- 98
    # Wetlands, Riverine Non-Forested
    r_reduced[r_reduced == 61 | r_reduced == 62 | r_reduced == 63] <- 97
    # Wetlands, Terrene Non-forested
    r_reduced[r_reduced == 71 | r_reduced == 22 | r_reduced == 73] <- 96
    # Tidal wetlands
    r_reduced[r_reduced == 91 | r_reduced == 92 | r_reduced == 93] <- 95
    
    # change everything else to NA
    r_reduced[r_reduced != 99 & 
                r_reduced != 98 &
                r_reduced != 97 &
                r_reduced != 96 &
                r_reduced != 95] <- NA
    
    entropy[i] <- lsm_l_ent(r_reduced)$value
    
  }
  
  return(entropy)
  
}

entropy_4000 <- extractEntropy(buffs_rp[[1]])
entropy_2000 <- extractEntropy(buffs_rp[[2]])
entropy_1000 <- extractEntropy(buffs_rp[[3]])


## Available Habitat and Impervious -----------------------------------------

# this could be parallelized to go faster, but the landcover raster might blow up mem
extractLandCover <- function(x, buff_size){
  
  lc_mat <- matrix(NA, nrow = nrow(x), ncol = 2)
  
  for(i in 1:nrow(x)){
    
  lc_extract <- terra::extract(lc, x[i,])

  lc_prop <- prop.table(table(lc_extract))
  
  # available habitat
  lc_mat[i,1] <- sum(lc_prop[,c("Forest", "Natural Succession", 
                                "Tree Canopy, Other",
                                "Wetlands, Tidal Non-forested",
                                "Wetlands, Riverine Non-forested",
                                "Wetlands, Terrene Non-forested")], na.rm = TRUE)
  # impervious
  lc_mat[i, 2] <- sum(lc_prop[,c("Impervious Roads", "Impervious Structures", 
                                "Impervious, Other", "Tree Canopy over Impervious")],
                      na.rm = TRUE)
  
  }
  
  colnames(lc_mat) <- c(paste0("habitat_", buff_size),
                        paste0("impervious_", buff_size))
  
  return(lc_mat)
  
}

lc_4000 <- extractLandCover(buffs_rp[[1]], buff_size = 4000)
lc_2000 <- extractLandCover(buffs_rp[[2]], buff_size = 2000)
lc_1000 <- extractLandCover(buffs_rp[[3]], buff_size = 1000)


## Elevation & Ruggedness ----------------------------------------------------

elev_files <- list.files("./Data/GIS/Elevation/", recursive = TRUE, pattern = ".img$")

# create a list object to store the rasters
elev_list <- vector("list", length(elev_files))

for(i in 1:length(elev_list)){
  elev_list[[i]] <- rast(paste0("./Data/GIS/Elevation/",elev_files[i]))
}

# make a spatrast collection
rsrc_tmp <- sprc(elev_list)

# merge the raster
elev_rast <- merge(rsrc_tmp)

# reproject buffers
buffs_rp <- lapply(buff_list, function(x){
  st_transform(x, crs = st_crs(elev_rast))
})

# extract the mean elevation for each scale
elev_4000 <- terra::extract(elev_rast, buffs_rp[[1]], mean, na.rm = TRUE)
elev_2000 <- terra::extract(elev_rast, buffs_rp[[2]], mean, na.rm = TRUE)
elev_1000 <- terra::extract(elev_rast, buffs_rp[[3]], mean, na.rm = TRUE)

# calculate ruggedness
ruggedness <- terrain(elev_rast, v = "TRI")

# extract the mean ruggedness for each scale

rugged_4000 <- terra::extract(ruggedness, buffs_rp[[1]], mean, na.rm = TRUE)
rugged_2000 <- terra::extract(ruggedness, buffs_rp[[2]], mean, na.rm = TRUE)
rugged_1000 <- terra::extract(ruggedness, buffs_rp[[3]], mean, na.rm = TRUE)


## Canopy Height -----------------------------------------------------------

# https://glad.umd.edu/dataset/gedi
# load raster
cp_tmp <- rast("./Data/GIS/ForestHeight_GLAD/Forest_height_2019_NAM.tiff")

# crop down to elevation extent
cp <- crop(cp_tmp, ruggedness)

#reclassify canopy raster - no data value
cp[cp == 103] <- NA

# convert everything else greater than 60 to a 0 - not real values
cp[cp > 60] <- 0 

# extract the mean from a 1km radius
# projection happens to be in the same as ruggedness
canopy_4000 <- terra::extract(cp, buffs_rp[[1]], mean, na.rm = TRUE)
canopy_2000 <- terra::extract(cp, buffs_rp[[2]], mean, na.rm = TRUE)
canopy_1000 <- terra::extract(cp, buffs_rp[[3]], mean, na.rm = TRUE)


## Population Density --------------------------------------------------------

pop_2022 <- get_acs(
  geography = "cbg",
  variables = "B02001_001", # total pop
  year = 2022,
  state = c("DC","VA","MD"),
  geometry = TRUE,
  output = "wide"
) |>
  select(-ends_with("M"))

pop_2022 <- pop_2022 |>
  mutate(area = st_area(pop_2022))

buffs_rp <- lapply(buff_list, function(x){
  st_transform(x, crs = st_crs(pop_2022))
})

getPop <- function(x){
  
  pop_dens <- rep(NA, nrow(x))
  
  for(i in 1:length(pop_dens)){
    
    # create a sf that is just the census data that intersects with a buffer
    pop_adjusted <- st_intersection(pop_2022, x[i,]) %>%
      # calculate overlap area and proportion of coverage
      mutate(intersect_area = st_area(.),
             intersect_prop = intersect_area/area,
             total_adjusted = B02001_001E * intersect_prop) |>
      pull(total_adjusted) |>
      sum(na.rm = TRUE)
    
    # calculte density - person/km^2
    pop_dens[i] <- pop_adjusted / (st_area(x)[1]/1e6) |>
      drop_units()
  }
  
  return(pop_dens)
  
}

# use function to get population density then convert to person/km^2
pop_dens_4000 <- getPop(buffs_rp[[1]])
pop_dens_2000 <- getPop(buffs_rp[[2]])
pop_dens_1000 <- getPop(buffs_rp[[3]])


# Combine all the variables -----------------------------------------------


covs <- data.frame(dist2road = dist2road,
                   dist2wetland = dist2wetland,
                   dist2stream = dist2stream,
                   dist2water = dist2water,
                   dist2trail = dist2trail,
                   entropy_4000 = entropy_4000,
                   entropy_2000 = entropy_2000,
                   entropy_1000 = entropy_1000,
                   habitat_4000 = lc_4000[,1],
                   habitat_2000 = lc_2000[,1],
                   habitat_1000 = lc_1000[,1],
                   imp_4000 = lc_4000[,2],
                   imp_2000 = lc_2000[,2],
                   imp_1000 = lc_1000[,2],
                   elev_4000 = elev_4000[,2],
                   elev_2000 = elev_2000[,2],
                   elev_1000 = elev_1000[,2],
                   rugged_4000 = rugged_4000[,2],
                   rugged_2000 = rugged_2000[,2],
                   rugged_1000 = rugged_1000[,2],
                   canopy_4000 = canopy_4000[,2],
                   canopy_2000 = canopy_2000[,2],
                   canopy_1000 = canopy_1000[,2],
                   pop_dens_4000 = pop_dens_4000,
                   pop_dens_2000 = pop_dens_2000,
                   pop_dens_1000 = pop_dens_1000
                   )

# to take a look at correlations
cor_table <- cor(covs)
#write.csv(cor_table, "./Results/cor_table_covs.csv")

for(i in 1:ncol(covs)){
  hist(covs[,i], main = colnames(covs)[i])
}

# distance to trail needed to be logged
covs$log_dist_trail <- log(covs$dist2trail)

# remove the non-logged distance to trail variable
covs <- covs[,!(names(covs) == "dist2trail")]

write.csv(covs, "./Data/2024-12-19_sitecovariates.csv")

# Scale covariates for modeling --------------------------------------------

# scale the covariates
scaled_covs <- apply(covs, 2, function(x) {
  (x - mean(x))/sd(x)
})

# give them unique column names
colnames(scaled_covs) <- paste0(colnames(scaled_covs), "_scaled")
write.csv(scaled_covs, "./Data/2024-12-19_scaled_covariates.csv", 
          row.names = FALSE)

