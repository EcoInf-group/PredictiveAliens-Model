# rationale: Prima et al. (2024) https://doi.org/10.1111/2041-210X.14444
# also consider Eichenberg et al (2021) https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.15447 for environmental variables

# use Copernicus data (Corine Land Cover, Small Woody Features, water bodies)
# merraclim for environmental variables
# merraclim legend:
# BIO1: Annual Mean Temperature
# BIO2: Mean Diurnal Range
# BIO3: Isothermality
# BIO4: Temperature Seasonality
# BIO5: Max Temperature of Warmest Month
# BIO6: Min Temperature of Coldest Month
# BIO7: Temperature Annual Range
# BIO8: Mean Temperature of Wettest Quarter
# BIO9: Mean Temperature of Driest Quarter
# BIO10: Mean Temperature of Warmest Quarter
# BIO11: Mean Temperature of Coldest Quarter
# BIO12: Annual Precipitation
# BIO13: Precipitation of Wettest Month
# BIO14: Precipitation of Driest Month
# BIO15: Precipitation Seasonality
# BIO16: Precipitation of Wettest Quarter
# BIO17: Precipitation of Driest Quarter
# BIO18: Precipitation of Warmest Quarter
# BIO19: Precipitation of Coldest Quarter

library(terra)
library(biomod2)
library(PresenceAbsence)
library(dplyr)
library(sf)

setwd("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/")

# Senecio inaequidens ====
ger <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_0") %>%
  dplyr::filter(COUNTRY == "Germany")

sen.ina <- read.csv("data/species occurrence data/senecio inaequidens/gbif.22.01.2026.csv", sep = "\t")
add <- read.csv("data/species occurrence data/senecio inaequidens/SenecioSpread_HegerBoehmerCaspianSuppl.csv",
                sep = ";")

sen.ina <- tibble(
  senecio = 1,
  x = c(sen.ina$decimalLongitude, add$Long),
  y = c(sen.ina$decimalLatitude, add$Lat)
) %>%
  na.omit()

sen.ina.XY <- sen.ina[, c("x", "y")]
traffic <- rast(
  "data/traffic data/eu.traffic.volume.raster.german.glm.predict.complete.tiff"
) %>%
  crop(vect(ger)) %>%
  mask(vect(ger), inverse = FALSE)

mcl <- list.files("data/environmental data/merraclim 2_5m_mean_00s",
                  full.names = TRUE)
mc.stack <- rast(mcl[c(1, 4, 17, 7)]) %>%
  crop(vect(ger)) %>%
  mask(vect(ger), inverse = FALSE) %>%
  resample(traffic)

nitr <- rast("data/environmental data/soil nitrogen/out.tif") %>%
  crop(vect(ger)) %>%
  mask(vect(ger), inverse = FALSE) %>%
  resample(traffic)

# important for LC: first reclassify THEN resample, otherwise it makes up new values during resampling (bilinear) and these cannot be reclassified. or make categorical first, then resample with other mode than bilinear, see ?resample
LC <- rast(
  "data/environmental data/copernicus landcover classes germany/U2018_CLC2018_V2020_20u1.tif"
)

LC <- subst(LC, 1:11, 1, raw = TRUE) # built-up
LC <- subst(LC, 12:22, 12, raw = TRUE) # agricultural land
LC <- subst(LC, 23:25, 23, raw = TRUE) # forests
LC <- subst(LC, 26:29, 26, raw = TRUE) # scrubs and herbs
LC <- subst(LC, 30:34, 30, raw = TRUE) # open or little or no vegetation
LC <- subst(LC, 35:39, 35, raw = TRUE) # wetlands
LC <- subst(LC, 40:44, 40, raw = TRUE) # water bodies
LC <- subst(LC, 48:255, 48, raw = TRUE) # NA

LC.legend <-
  readxl::read_xls(
    "data/environmental data/copernicus landcover classes germany/Documentation U2018_CLC2018_V2020_20u1_raster100m_tiled_doc/Info/Legend/Vector/clc_legend.xls"
  )

LC.legend <- data.frame(id = LC.legend$GRID_CODE,
                        lc.cover = c(LC.legend$LABEL2)) #IMPORTANT: choose LABEL1 or LABEL2 defines the level of detail!

LC.legend <- LC.legend[c(1, 12, 23, 26, 30, 35, 40, 48), ]
LC.r <- LC
levels(LC.r) <- LC.legend
plot(LC.r)

LC.r <- LC.r %>%
  project(traffic) %>%
  crop(vect(ger)) %>%
  mask(vect(ger), inverse = FALSE) %>%
  resample(traffic)
plot(LC.r)

expl.var <- c(mc.stack$"2_5m_mean_00s_bio1",
              mc.stack$"2_5m_mean_00s_bio12",
              nitr,
              traffic,
              LC.r)

names(expl.var) <- c("bio1", "bio12", "nitrogen", "traffic.glm", "LC")

terra::layerCor(
  expl.var,
  fun = "cor",
  w,
  asSample = TRUE,
  use = "everything",
  maxcell = Inf
)$correlation

one.rec.per.grid.cell <- function(ref.raster, occurences) {
  cells <- cellFromXY(ref.raster, occurences) %>% unique()
  coords <- xyFromCell(ref.raster, cells)
}

sen.ina.red <- one.rec.per.grid.cell(
  ref.raster = expl.var[[1]],
  occurences = cbind(
    lon = sen.ina$x,
    lat = sen.ina$y
  )
) %>%
  na.omit()

sen.ina <- tibble(senecio = 1,
                  x = sen.ina.red[,1],
                  y = sen.ina.red[,2])

sen.ina.data <- BIOMOD_FormatingData(
  resp.var = sen.ina$senecio,
  expl.var = expl.var,
  resp.xy = sen.ina[, c("x", "y")],
  resp.name = "senecio",
  PA.nb.rep = 2,
  PA.nb.absences = nrow(sen.ina),
  PA.strategy = "random",
  filter.raster = TRUE,
  dir.name = "data/biomod working directory/senecio inaequidens/"
)

sen.ina.models <- BIOMOD_Modeling(
  bm.format = sen.ina.data,
  models = c("GBM"),
  #OPT.user = mod.opt, # returns an error if no formula provided in mod.opt.
  CV.strategy = "random",
  CV.nb.rep = 2,
  CV.perc = .7,
  nb.cpu = 8,
  modeling.id = "sen.ina",
  var.import = 2
)

sen.ina.models

# Get evaluation scores & variables importance
eval <- get_evaluations(sen.ina.models)
eval %>% dplyr::filter(run %in% c("RUN1", "RUN2") &
                         metric.eval %in% c("ROC")) %>% summary()
eval %>% dplyr::filter(run %in% c("RUN1", "RUN2") &
                         metric.eval %in% c("TSS")) %>% summary()



# Represent evaluation scores
bm_PlotEvalMean(bm.out = sen.ina.models, dataset = 'calibration')
bm_PlotEvalMean(bm.out = sen.ina.models, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = sen.ina.models, group.by = c('algo', 'run'))

# Represent variables importance
bm_PlotVarImpBoxplot(bm.out = sen.ina.models,
                     group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = sen.ina.models,
                     group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = sen.ina.models,
                     group.by = c('algo', 'expl.var', 'run'))

# Response curves:
mods <- get_built_models(sen.ina.models, run = 'RUN1')
bm_PlotResponseCurves(bm.out = sen.ina.models,
                      models.chosen = mods,
                      fixed.var = 'median')

bm_PlotResponseCurves(bm.out = sen.ina.models,
                      models.chosen = mods,
                      fixed.var = 'min')

mods <- get_built_models(sen.ina.models, full.name = "senecio_PA2_allRun_GBM")
bm_PlotResponseCurves(
  bm.out = sen.ina.models,
  models.chosen = mods,
  fixed.var = 'median',
  do.bivariate = TRUE
)

sen.ina.proj <- BIOMOD_Projection(
  bm.mod = sen.ina.models,
  proj.name = 'Current.traffic',
  new.env = expl.var,
  models.chosen = 'all'
)


plot(sen.ina.proj)
sen.ina.proj@proj.out # get the link
sen.ina.proj.rast <- rast("./senecio/proj_Current.traffic/proj_Current.traffic_senecio.tif") # raster file from the link
sen.ina.proj.rast
names(sen.ina.proj.rast)

sen.ina.GBM <- sen.ina.proj.rast[[c(grep("GBM", names(sen.ina.proj.rast)))]]
sen.ina.GBM <- sen.ina.GBM[[c(grep("PA", names(sen.ina.GBM)))]]
sen.ina.GBM <- sen.ina.GBM[[c(grep("RUN", names(sen.ina.GBM)))]]
sen.ina.GBM <- mean(sen.ina.GBM, na.rm = TRUE)
plot(sen.ina.proj.rast["senecio_allData_allRun_GBM"])
plot(sen.ina.GBM)

# write for use in resistant kernel:
writeRaster(
  sen.ina.GBM,
  "data/simulation input data/senecio inaequidens/biomod2_GBM.2PA.average.bio1.bio12.nitrogen.traffic.LC.tif"
)

# Tapinoma magnum: ====
library(terra)
library(biomod2)
library(PresenceAbsence)
library(dplyr)
library(sf)
setwd("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/")

eur <- c(
  "Albania",
  "Austria",
  "Belarus",
  "Belgium",
  "Bulgaria",
  "Croatia",
  "Czechia",
  "Denmark",
  "Finland",
  "France",
  "Germany",
  "Greece",
  "Hungary",
  "Ireland",
  "Italy",
  "Latvia",
  "Luxembourg",
  "Netherlands",
  "Norway",
  "Portugal",
  "Poland",
  "Romania",
  "San Marino",
  "Sardinia",
  "Serbia",
  "Sicily",
  "Slovakia",
  "Slovenia",
  "Spain",
  "Sweden",
  "Switzerland",
  "Slovenia",
  "Turkey",
  "Ukraine",
  "United Kingdom",
  "Estonia",
  "Lithuania",
  "Bosnia and Herzegovina",
  "Montenegro",
  "Kosovo",
  "North Macedonia",
  "Liechtenstein",
  "Andorra"
)

gadm.0 <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_0") %>%
  st_make_valid() %>%
  dplyr::filter(COUNTRY %in% eur)

# distributional data:
dist.rec.am <- readxl::read_xlsx("data/species occurrence data/tapinoma magnum/records tapinoma magnum [antmaps].xlsx")
dist.rec.am <- tibble(
  tap.mag = 1,
  x = dist.rec.am$dec_lon,
  y = dist.rec.am$dec_lat,
  year = dist.rec.am$`Year Collection`
) %>%
  na.omit()

dist.rec.am$x <- as.numeric(dist.rec.am$x)
dist.rec.am$y <- as.numeric(dist.rec.am$y)
dist.rec.am[370, 2] <- -3.36354

# seiferts data:
dist.rec.s1 <- readxl::read_xlsx("data/species occurrence data/tapinoma magnum/TAPINO_C [Seifert].xlsx") %>%
  dplyr::filter(HYP == "magn")
dist.rec.s2 <- readxl::read_xlsx("data/species occurrence data/tapinoma magnum/TAPINTRO [Seifert].xlsx") %>%
  dplyr::filter(HYP == "magn")
dist.rec.s <- tibble(
  tap.mag = 1,
  x = c(dist.rec.s1$LON, dist.rec.s2$LON),
  y = c(dist.rec.s1$LAT, dist.rec.s2$LAT),
  year = c(dist.rec.s1$YEAR, dist.rec.s2$YEAR)
)

dist.rec <- rbind(dist.rec.am, dist.rec.s) %>%
  na.omit() %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

dist.rec.bbox <- dist.rec %>%
  st_bbox() %>%
  st_as_sfc()

dist.rec.XY <- as_tibble(st_coordinates(dist.rec))

# environmental data:
traffic <- rast(
  "data/traffic data/eu.traffic.volume.raster.german.glm.predict.complete.tiff"
) %>%
  crop(vect(gadm.0)) %>%
  mask(vect(gadm.0), inverse = FALSE)
plot(traffic)

mcl <- list.files("data/environmental data/merraclim 2_5m_mean_00s",
                  full.names = TRUE)

mc.stack <- rast(mcl[c(2, 3, 5, 6)])
names(mc.stack) <- c(#"bio1",
  "bio10", "bio11", "bio13", "bio14")
mc.stack <- mc.stack %>%
  crop(vect(st_transform(gadm.0, st_crs(mc.stack)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(mc.stack))), inverse = FALSE) %>%
  project(traffic) %>%
  resample(traffic)
plot(mc.stack)
# important for LC: first reclassify THEN resample, otherwise it makes up new values during resampling (bilinear) and these cannot be reclassified. or make categorical first, then resample with other mode than bilinear, see ?resample

LC <- rast(
  "data/environmental data/copernicus landcover classes europe/U2018_CLC2018_V2020_20u1_aggregated_LVL1.tif"
)
LC <- LC %>%
  crop(vect(st_transform(gadm.0, st_crs(LC)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(LC))), inverse = FALSE) %>%
  clamp(lower = 1,
        upper = 40,
        values = FALSE)
# clamp to restrict values to this range, everything else (i.e. everything that has no legend entry) is then NA

# this is how values in LC were subsituted for level 1 aggregation:
#LC <- subst(LC, 1:11, 1)
#LC <- subst(LC, 12:22, 12, raw = TRUE)
#LC <- subst(LC, 23:34, 23, raw = TRUE)
#LC <- subst(LC, 35:39, 35, raw = TRUE)
#LC <- subst(LC, c(40:44, 50), 40, raw = TRUE)

# this is how values in LC were subsituted for level 2 aggregation:
#LC <- subst(LC, 2, 1)
#LC <- subst(LC, 3:6, 3, raw = TRUE)
#LC <- subst(LC, 7:9, 7, raw = TRUE)
#LC <- subst(LC, 10:11, 10, raw = TRUE)
#LC <- subst(LC, 12:14, 12, raw = TRUE)
#LC <- subst(LC, 15:17, 15, raw = TRUE)
#LC <- subst(LC, 19:22, 19, raw = TRUE)
#LC <- subst(LC, 23:25, 23, raw = TRUE)
#LC <- subst(LC, 26:29, 26, raw = TRUE)
#LC <- subst(LC, 30:34, 30, raw = TRUE)
#LC <- subst(LC, 35:36, 35, raw = TRUE)
#LC <- subst(LC, 37:39, 37, raw = TRUE)
#LC <- subst(LC, 40:41, 40, raw = TRUE)
#LC <- subst(LC, 42:44, 42, raw = TRUE)
#LC <- subst(LC, 49:50, 49, raw = TRUE)
#LC.legend <- 
#  readxl::read_xls("SDMs/environmental data/copernicus landcover classes europe/Documentation U2018_CLC2018_V2020_20u1_raster100m_tiled_doc/Info/Legend/Vector/clc_legend.xls")
#LC.legend <- data.frame(id = LC.legend$GRID_CODE, lc.cover = c(LC.legend$LABEL1)) #IMPORTANT: choose LABEL1 or LABEL2 defines the level of detail!
#
#LC.legend <- LC.legend[c(1,3,7,10,12,15,18,19,23,26,30,35,37,40,42,45,46,48),] # level 2 legend
LC.legend <- tibble(
  id = c(1, 12, 23, 35, 40),
  lc.cover = c(
    "artificial.surfaces",
    "agricultural.areas",
    "forest.and.semi.natural.areas",
    "wetlands",
    "water.bodies"
  )
)
# 48 = NA

levels(LC) <- LC.legend
traffic.training <- traffic %>%
  crop(vect(st_transform(dist.rec.bbox, st_crs(traffic)))) %>%
  mask(vect(st_transform(dist.rec.bbox, st_crs(traffic))), inverse = FALSE)
plot(traffic.training)

LC.training <- LC %>%
  crop(vect(st_transform(dist.rec.bbox, st_crs(LC)))) %>%
  mask(vect(st_transform(dist.rec.bbox, st_crs(LC))), inverse = FALSE) %>%
  project(traffic.training) %>%
  resample(traffic.training)
plot(LC.training)

mc.stack.training <- mc.stack %>%
  crop(vect(st_transform(dist.rec.bbox, st_crs(mc.stack)))) %>%
  mask(vect(st_transform(dist.rec.bbox, st_crs(mc.stack))), inverse = FALSE) %>%
  project(traffic.training) %>%
  resample(traffic.training)
plot(mc.stack.training)

# note: population probably heavily biased because reflects sampling effort/ probability of detection
pop <- rast(
  "data/environmental data/GHS_POP_E2030_GLOBE_R2023A_54009_100_V1_0/GHS_POP_E2030_GLOBE_R2023A_54009_100_V1_0.tif"
)
pop <- pop %>%
  crop(vect(st_transform(gadm.0, st_crs(pop)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(pop))), inverse = FALSE)
plot(pop)

pop.training <- pop %>%
  crop(vect(st_transform(dist.rec.bbox, st_crs(pop)))) %>%
  mask(vect(st_transform(dist.rec.bbox, st_crs(pop))), inverse = FALSE) %>%
  project(traffic.training) %>%
  resample(traffic.training)
plot(pop.training)

expl.var.training <- c(
  mc.stack.training$bio10,
  mc.stack.training$bio13,
  mc.stack.training$bio14,
  LC.training,
  pop.training
)

names(expl.var.training) <- c("bio10", "bio13", "bio14", "LC", "pop")

### check for autocorrelation --------------------------------------------------
terra::layerCor(
  expl.var.training,
  fun = "cor",
  w,
  asSample = TRUE,
  use = "everything",
  maxcell = Inf
)$correlation

## model building --------------------------------------------------------------
gc()

tap.mag.data <- BIOMOD_FormatingData(
  resp.var = dist.rec$tap.mag,
  expl.var = expl.var.training,
  resp.xy = dist.rec.XY[, c("X", "Y")],
  resp.name = "tap.mag",
  PA.nb.rep = 2,
  PA.nb.absences = nrow(dist.rec),
  PA.strategy = "disk",
  # usually "random"
  filter.raster = TRUE,
  PA.dist.min = 50000, # used for PA.strategy = "disk"
  #PA.dist.max = 100000 # used for PA.strategy = "disk"
  dir.name = "data/biomod working directory/tapinoma magnum/"
)

used.vars <- sort(unique(tap.mag.data@data.env.var$LC))
used.vars <- as.character(droplevels(used.vars))

tap.mag.models <- BIOMOD_Modeling(
  bm.format = tap.mag.data,
  models = c("GBM"),
  CV.strategy = "random",
  CV.nb.rep = 2,
  CV.perc = .7,
  nb.cpu = 8,
  modeling.id = "tap.mag",
  var.import = 2
)

tap.mag.models

# Get evaluation scores & variables importance
evalw14 <- get_evaluations(tap.mag.models)
evalw14 %>% dplyr::filter(run %in% c("RUN1", "RUN2") &
                            metric.eval %in% c("ROC")) %>% summary()
evalw14 %>% dplyr::filter(run %in% c("RUN1", "RUN2") &
                            metric.eval %in% c("TSS")) %>% summary()

get_variables_importance(tap.mag.models)

# Represent evaluation scores
bm_PlotEvalMean(bm.out = tap.mag.models, dataset = 'calibration')
bm_PlotEvalMean(bm.out = tap.mag.models, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = tap.mag.models, group.by = c('algo', 'run'))

# Represent variables importance
bm_PlotVarImpBoxplot(bm.out = tap.mag.models,
                     group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = tap.mag.models,
                     group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = tap.mag.models,
                     group.by = c('algo', 'expl.var', 'run'))

# Response curves:
mods <- get_built_models(tap.mag.models, algo = "GBM")
bm_PlotResponseCurves(bm.out = tap.mag.models,
                      models.chosen = mods,
                      fixed.var = 'median')

bm_PlotResponseCurves(bm.out = tap.mag.models,
                      models.chosen = mods,
                      fixed.var = 'min')

mods <- get_built_models(tap.mag.models, full.name = "tap.mag_PA2_allRun_GBM")
bm_PlotResponseCurves(
  bm.out = tap.mag.models,
  models.chosen = mods,
  fixed.var = 'median',
  do.bivariate = TRUE
)

mc.stack.proj <- mc.stack
LC.proj <- LC %>%
  project(mc.stack.proj$bio10) %>%
  resample(mc.stack.proj$bio10)
pop.proj <- pop %>%
  project(mc.stack.proj$bio10) %>%
  resample(mc.stack.proj$bio10)

expl.var.proj <- c(mc.stack.proj$bio10,
                   mc.stack.proj$bio13,
                   mc.stack.proj$bio14,
                   LC.proj,
                   pop.proj)

names(expl.var.proj) <- c("bio10", "bio13", "bio14", "LC", "pop")

tap.mag.proj <- BIOMOD_Projection(
  bm.mod = tap.mag.models,
  proj.name = 'Current.no.traffic',
  new.env = expl.var.proj,
  models.chosen = 'all'
)

plot(tap.mag.proj)
tap.mag.proj@proj.out # get the link
tap.mag.proj.rast <- rast("./tap.mag/proj_Current.no.traffic/proj_Current.no.traffic_tap.mag.tif") # raster file from the link
names(tap.mag.proj.rast)

tap.mag.GBM <- tap.mag.proj.rast[[c(grep("GBM", names(tap.mag.proj.rast)))]]
tap.mag.GBM <- tap.mag.GBM[[c(grep("PA", names(tap.mag.GBM)))]]
tap.mag.GBM <- tap.mag.GBM[[c(grep("RUN", names(tap.mag.GBM)))]]
tap.mag.GBM <- mean(tap.mag.GBM, na.rm = TRUE)
plot(tap.mag.proj.rast["tap.mag_allData_allRun_GBM"])
plot(tap.mag.GBM)
points(
  x = dist.rec.XY[, 1],
  y = dist.rec.XY[, 2],
  col = "red",
  pch = 19
)

# write for use in resistant kernel:
writeRaster(
  tap.mag.GBM,
  "data/simulation input data/tapinoma magnum/biomod2_GBM.bio10.bio13.bio14.LC.pop.density.tif",
  overwrite = TRUE
)

# Myocastor coypus -------------------------------------------------------------
setwd("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/")
library(terra)
library(biomod2)
library(PresenceAbsence)
library(dplyr)
library(sf)
library(spThin)

eur <- c(
  "Albania",
  "Austria",
  "Belarus",
  "Belgium",
  "Bulgaria",
  "Croatia",
  "Czechia",
  "Denmark",
  "Finland",
  "France",
  "Germany",
  "Greece",
  "Hungary",
  "Ireland",
  "Italy",
  "Latvia",
  "Luxembourg",
  "Netherlands",
  "Norway",
  "Portugal",
  "Poland",
  "Romania",
  "San Marino",
  "Sardinia",
  "Serbia",
  "Sicily",
  "Slovakia",
  "Slovenia",
  "Spain",
  "Sweden",
  "Switzerland",
  "Slovenia",
  "Turkey",
  "Ukraine",
  "United Kingdom",
  "Estonia",
  "Lithuania",
  "Bosnia and Herzegovina",
  "Montenegro",
  "Kosovo",
  "North Macedonia",
  "Liechtenstein",
  "Andorra"
)

gadm.0 <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_0") %>%
  st_make_valid() %>%
  dplyr::filter(COUNTRY %in% eur)

myo.coy <- read.csv("data/species occurrence data/myocastor coypus/M coypus gbif 13082025.csv", sep = "\t")
myo.coy <- tibble(
  myocastor = 1,
  longitude = myo.coy$decimalLongitude,
  latitude = myo.coy$decimalLatitude,
  year = myo.coy$year
) %>%
  na.omit()

#
#
#
#
# with Anna's data:
a.s.rec <- st_read("data/species occurrence data/myocastor coypus/nutria_occurrence_records_1980_2018 anna schertler.gpkg")###

points <- st_centroid(a.s.rec) %>%
  st_transform("epsg:4326")

plot(points["year"])

points$dec <- points$year

points$dec[points$dec <= 1990] <- 1990
points$dec[points$dec <= 2000 & points$dec > 1990] <- 2000
points$dec[points$dec <= 2010 & points$dec > 2000] <- 2010
points$dec[points$dec <= 2020 & points$dec > 2010] <- 2020

plot(points["dec"])

points.c <- st_coordinates(points)
points <- tibble(
  myocastor = 1,
  longitude = points.c[, 1],
  latitude = points.c[, 2],
  year = points$year
)
#
#
#
#
#

myo.coy <- bind_rows(points, myo.coy)

one.rec.per.grid.cell <- function(ref.raster, occurences) {
  cells <- cellFromXY(ref.raster, occurences) %>% unique()
  coords <- xyFromCell(ref.raster, cells)
}

coords.matrix <- one.rec.per.grid.cell(
  ref.raster = rast(
    "data/traffic data/eu.traffic.volume.raster.german.glm.predict.complete.tiff"
  ),
  occurences = cbind(lon = myo.coy$longitude, lat = myo.coy$latitude)
) %>%
  na.omit()

myo.coy.xy <- st_as_sf(as_tibble(coords.matrix),
                       coords = c("x", "y"),
                       crs = 4326)

myo.coy <- myo.coy %>% dplyr::filter(year > 1961)
myo.coy.70 <- myo.coy %>% dplyr::filter(year <= 1970)
myo.coy.80 <- myo.coy %>% dplyr::filter(year > 1970 & year <= 1980)
myo.coy.90 <- myo.coy %>% dplyr::filter(year > 1980 & year <= 1990)
myo.coy.00 <- myo.coy %>% dplyr::filter(year > 1990 & year <= 2000)
myo.coy.10 <- myo.coy %>% dplyr::filter(year > 2000 & year <= 2010)
myo.coy.20 <- myo.coy %>% dplyr::filter(year > 2010 & year <= 2020)
myo.coy.25 <- myo.coy %>% dplyr::filter(year > 2020 & year <= 2025)

dist.rec.bbox <- myo.coy.xy %>%
  st_bbox() %>%
  st_as_sfc()

# environmental data:
traffic <- rast(
  "data/traffic data/eu.traffic.volume.raster.german.glm.predict.complete.tiff"
) #%>%
mcl <- list.files("data/environmental data/merraclim 2_5m_mean_00s",
                  full.names = TRUE)
# according to Schertler et al 2020:
mc.stack <- rast(mcl[c(12, 16, 2, 7, 9)]) %>%
  crop(vect(st_transform(gadm.0, st_crs(mc.stack)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(mc.stack))), inverse = FALSE)
names(mc.stack) <- c("bio2", "bio6", "bio10", "bio15", "bio17") #%>%

# important for LC: first reclassify THEN resample, otherwise it makes up new values during resampling (bilinear) and these cannot be reclassified. or make categorical first, then resample with other mode than bilinear, see ?resample
LC <- rast(
  "data/environmental data/copernicus landcover classes europe/U2018_CLC2018_V2020_20u1_aggregated_LVL1.tif"
)
LC <- LC %>%
  crop(vect(st_transform(gadm.0, st_crs(LC)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(LC))), inverse = FALSE) %>%
  clamp(lower = 1,
        upper = 40,
        values = FALSE)
# clamp to restrict values to this range, everything else (i.e. everything that has no legend entry) is then NA

LC.legend <- tibble(
  id = c(1, 12, 23, 35, 40),
  lc.cover = c(
    "artificial.surfaces",
    "agricultural.areas",
    "forest.and.semi.natural.areas",
    "wetlands",
    "water.bodies"
  )
)
levels(LC) <- LC.legend

plot(st_transform(dist.rec.bbox, st_crs(LC)), add = TRUE)
plot(
  st_transform(myo.coy.xy, st_crs(LC)),
  add = TRUE,
  pch = 19,
  cex = .5,
  col = "red"
)

traffic <- traffic %>%
  crop(vect(st_transform(gadm.0, st_crs(traffic)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(traffic))), inverse = FALSE)
plot(traffic)

LC.r <- LC
LC.r <- LC.r %>%
  project(traffic) %>%
  resample(traffic)
plot(LC.r)

mc.stack <- mc.stack %>%
  project(traffic) %>%
  resample(traffic) %>%
  mask(LC.r, inverse = FALSE)
plot(mc.stack)

pop <- rast(
  "data/environmental data/GHS_POP_E2030_GLOBE_R2023A_54009_100_V1_0/GHS_POP_E2030_GLOBE_R2023A_54009_100_V1_0.tif"
)
pop.r <- pop %>%
  crop(vect(st_transform(gadm.0, st_crs(pop)))) %>%
  mask(vect(st_transform(gadm.0, st_crs(pop))), inverse = FALSE) %>%
  project(traffic) %>%
  resample(traffic) %>%
  mask(LC.r, inverse = FALSE)
plot(pop.r)

expl.var <- c(mc.stack, LC.r, pop.r)

names(expl.var) <- c("bio2", "bio6", "bio10", "bio15", "bio17", "LC", "pop")

### check for autocorrelation ====
terra::layerCor(
  expl.var,
  fun = "cor",
  w,
  asSample = TRUE,
  use = "everything",
  maxcell = Inf
)$correlation
# bio6 and bio17 highly autocorrelated
# 6 = min temperature of coldest quarter, 17 = precipitation of driest quarter

expl.var$bio17 = NULL

plot(expl.var$bio2)
plot(dist.rec.bbox, add = TRUE)
expl.var.training <- expl.var %>%
  crop(vect(dist.rec.bbox)) %>%
  mask(vect(dist.rec.bbox), inverse = FALSE)

myo.coy <- tibble(myocastor = 1, x = coords.matrix[, 1], y = coords.matrix[, 2])

myo.coy.data <- BIOMOD_FormatingData(
  resp.var = myo.coy$myocastor,
  expl.var = expl.var.training,
  resp.xy = myo.coy[, c("x", "y")],
  resp.name = "myocastor",
  PA.nb.rep = 2,
  PA.nb.absences = nrow(myo.coy),
  PA.strategy = "random",
  # usually "random"
  filter.raster = TRUE,
  dir.name = "data/biomod working directory/myocastor coypus/"
)

used.vars <- sort(unique(myo.coy.data@data.env.var$LC))
used.vars <- as.character(droplevels(used.vars))

myo.coy.models <- BIOMOD_Modeling(
  bm.format = myo.coy.data,
  models = c(#"GLM", "GAM",
    #"RF",
    "GBM"),
  # GBM usually best, so the others muted here to save time
  #OPT.user = mod.opt, # returns an error if no formula provided in mod.opt.
  CV.strategy = "random",
  CV.nb.rep = 2,
  CV.perc = .7,
  nb.cpu = 4,
  modeling.id = "myo.coy",
  var.import = 2
)

myo.coy.models

eval <- get_evaluations(myo.coy.models)
eval %>% dplyr::filter(run %in% c("RUN1", "RUN2") &
                         metric.eval %in% c("ROC")) %>% summary()
eval %>% dplyr::filter(run %in% c("RUN1", "RUN2") &
                         metric.eval %in% c("TSS")) %>% summary()


# Represent evaluation scores
bm_PlotEvalMean(bm.out = myo.coy.models, dataset = 'calibration')
bm_PlotEvalMean(bm.out = myo.coy.models, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myo.coy.models, group.by = c('algo', 'run'))

# Represent variables importance
bm_PlotVarImpBoxplot(bm.out = myo.coy.models,
                     group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myo.coy.models,
                     group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myo.coy.models,
                     group.by = c('algo', 'expl.var', 'run'))

# Response curves:
mods <- get_built_models(myo.coy.models, run = 'RUN1')
bm_PlotResponseCurves(bm.out = myo.coy.models,
                      models.chosen = mods,
                      fixed.var = 'median')

bm_PlotResponseCurves(bm.out = myo.coy.models,
                      models.chosen = mods,
                      fixed.var = 'min')

mods <- get_built_models(myo.coy.models, full.name = "myocastor_PA2_allRun_GBM")
bm_PlotResponseCurves(
  bm.out = myo.coy.models,
  models.chosen = mods,
  fixed.var = 'median',
  do.bivariate = TRUE
)

myo.coy.proj <- BIOMOD_Projection(
  bm.mod = myo.coy.models,
  proj.name = 'Current.no.traffic',
  new.env = expl.var,
  models.chosen = 'all'
)

plot(myo.coy.proj)
myo.coy.proj@proj.out # get the link
myo.coy.proj.rast <- rast("./myocastor/proj_Current.no.traffic/proj_Current.no.traffic_myocastor.tif") # raster file from the link
names(myo.coy.proj.rast)

myo.coy.GBM <- myo.coy.proj.rast[[c(grep("GBM", names(myo.coy.proj.rast)))]]
myo.coy.GBM <- myo.coy.GBM[[c(grep("PA", names(myo.coy.GBM)))]]
myo.coy.GBM <- myo.coy.GBM[[c(grep("RUN", names(myo.coy.GBM)))]]
myo.coy.GBM <- mean(myo.coy.GBM, na.rm = TRUE)
plot(myo.coy.proj.rast["myocastor_allData_allRun_GBM"])
plot(myo.coy.GBM)
points(
  x = myo.coy$x,
  y = myo.coy$y,
  col = "red",
  pch = 19
)

# write for use in resistant kernel:
writeRaster(
  myo.coy.GBM,
  "data/simulation input data/myocastor coypus/biomod2_GBM.update.bio02.bio6.bio10.bio15.LC.pop.gbif and anna.tif",
  overwrite = TRUE
)


# Plot SDM maps ----------------------------------------------------------------
sen.ina.GBM <- rast("data/simulation input data/senecio inaequidens/biomod2_GBM.2PA.average.bio1.bio12.nitrogen.traffic.LC.tif")
sen.ina.gbm.r <- 1 / max(values(sen.ina.GBM), na.rm = TRUE) * sen.ina.GBM
plot(sen.ina.gbm.r)

tap.mag.GBM <- rast("data/simulation input data/tapinoma magnum/biomod2_GBM.bio10.bio13.bio14.LC.pop.density.tif")
tap.mag.gbm.r <- 1 / max(values(tap.mag.GBM), na.rm = TRUE) * tap.mag.GBM
plot(tap.mag.gbm.r)

myo.coy.GBM <- rast("data/simulation input data/myocastor coypus/biomod2_GBM.update.bio02.bio6.bio10.bio15.LC.pop.gbif and anna.tif")
myo.coy.gbm.r <- 1 / max(values(myo.coy.GBM), na.rm = TRUE) * myo.coy.GBM # with this step it is set to a scale of 0 to 1
plot(myo.coy.gbm.r)
