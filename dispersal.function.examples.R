# DISPERSAL SIMULATION FUNCTION ------------------------------------------------

setwd("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/")

# load packages:
library(terra)
library(dplyr)
library(sf)
library(gridprocess) # devtools::install_github("ethanplunkett/gridprocess")
library(mltools)
#library(foreach)
#library(doParallel)
library(purrr)
library(mirai)

# dispersal() functions are stored here:
source("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/code/PredictiveAliens-Model_git/dispersal function.R")
source("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/code/PredictiveAliens-Model_git/parallel dispersal function.R") 

# dispersal() allows to provide several named landscape resistance matrices for rawspread, e.g. "power.2" to compare different transformations in the same run.
# par.dispersal() cannot do that but parallelizes the land.spread by using a defined number of cores.


#
#
#
#
#
# SENECIO INAEQUIDENS ----------------------------------------------------------
# define inputs needed for function:
## data preparation ------------------------------------------------------------
### geographic reference -------------------------------------------------------
ref.raster <- rast("data/simulation input data/senecio inaequidens/biomod2_GBM.2PA.average.bio1.bio12.nitrogen.traffic.LC.tif") %>%
  terra::aggregate(5)
names(ref.raster) <- "layer"
empty.r <- ref.raster
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything

ger <- st_read("data/environmental data/gadm41_DEU.gpkg", # to crop to germany
               layer = "ADM_ADM_0") %>%
  st_transform(st_crs(ref.raster))


### occurrence data (CASPIAN & gbif) -------------------------------------------
sen.spread <- data.table::fread("data/species occurrence data/senecio inaequidens/SenecioSpread_HegerBoehmerCaspianSuppl.csv") # data from Hanno for the Caspian suppl map
sen.spread <- sen.spread[, -c(2, 3)]
sen.gbif <- read.csv("data/species occurrence data/senecio inaequidens/gbif.22.01.2026.csv", sep = "\t")
sen.gbif <- tibble(
  Long = sen.gbif$decimalLongitude,
  Lat = sen.gbif$decimalLatitude,
  year = sen.gbif$year
) %>%
  na.omit()
sen.spread <- bind_rows(sen.spread, sen.gbif)

sen.spread <- st_as_sf(sen.spread, coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(st_crs(ref.raster))
sen.spread <- st_intersection(sen.spread, ger$geom) 


sen.spread$decade <- sen.spread$year
# convert to decades instead of years:
sen.spread$decade[sen.spread$decade < 1979] <- 1979
sen.spread$decade[between(sen.spread$decade, 1980, 1989)] <- 1989
sen.spread$decade[between(sen.spread$decade, 1990, 1999)] <- 1999
sen.spread$decade[between(sen.spread$decade, 2000, 2009)] <- 2009
sen.spread$decade[between(sen.spread$decade, 2010, 2019)] <- 2019
sen.spread$decade[between(sen.spread$decade, 2020, 2025)] <- 2025

sen.ini <- sen.spread %>% # filter based on record day which occurrence records are used as initial distribution (i.e. starting points) in the simulation
  dplyr::filter(year <= 1979)

ini.dist.r <- empty.r %>%
  terra::mask(vect(st_buffer(sen.ini, 10000)), updatevalue = 1, inverse = TRUE) %>%
  mask(ref.raster) # has to be masked again because otherwise the buffered occurrence points extent beyond the country borders
plot(ini.dist.r)

# final reference distribution:
ref.p <- sen.spread %>%
  dplyr::filter(year <= 2009)

ref.dist.r <- c(terra::mask(ini.dist.r, vect(st_buffer(dplyr::filter(sen.spread, year <= 1989), 10000)), updatevalue = 1, inverse = TRUE),
                terra::mask(ini.dist.r, vect(st_buffer(dplyr::filter(sen.spread, year <= 1999), 10000)), updatevalue = 1, inverse = TRUE),
                terra::mask(ini.dist.r, vect(st_buffer(dplyr::filter(sen.spread, year <= 2009), 10000)), updatevalue = 1, inverse = TRUE)
)
names(ref.dist.r) <- c("ref.1989", "ref.1999", "ref.2009")
plot(ref.dist.r)
#
#
#
#
#


### read SDM map & make resistance raster for gridprocess::rawpsread() ---------
gbm.r <- rast("data/simulation input data/senecio inaequidens/biomod2_GBM.2PA.average.bio1.bio12.nitrogen.traffic.LC.tif") %>%  # biomod SDM output from boosted regression trees (gbm)
  project(ref.raster)
names(gbm.r) <- "layer"
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r <- terra::mask(gbm.r, ref.raster)
gbm.r <- subst(gbm.r, NA, 0)

mask.thresh <- 0.35
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(ger)) %>%                     # crops it to the area of interest -> CHECK IF NECESSARY # deleted because if I do this, then the grid allows to cross borders over time
  cells()                                         # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r ^ 1                                # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r[mask] <- 0
#x11()
#plot(gbm.r)

gbm.r.inv <- gbm.r * -1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

gbm.r.inv <- subst(gbm.r.inv, # must not have NAs for the function below, so replace with 9999 to make these areas not crossable
                   c(1,NA), 9999) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network because of color-scale. reduce spread.limit to make it visible again.

gbm.r.inv <- asgrid( # convert to grid for spread function
  gbm.r.inv, 
  xll = xmin(gbm.r.inv),
  yll = ymin(gbm.r.inv),
  cellsize = 1000
) # update cellsize with the aggregate factor * 1000m

# rename for parametrization loop:
power.1 <- gbm.r.inv


gbm.r <- rast("data/simulation input data/senecio inaequidens/biomod2_GBM.2PA.average.bio1.bio12.nitrogen.traffic.LC.tif") %>%  # biomod SDM output from boosted regression trees (gbm)
  project(ref.raster)
names(gbm.r) <- "layer"
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r <- terra::mask(gbm.r, ref.raster)
gbm.r <- subst(gbm.r, NA, 0)

mask.thresh <- 0.35
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(ger)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY # deleted because if I do this, then the grid allows to cross borders over time
  cells()                               # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r ^ 1.5 # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r[mask] <- 0
#x11()
#plot(gbm.r)

gbm.r.inv <- gbm.r * -1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

gbm.r.inv <- subst(gbm.r.inv, # must not have NAs for the function below, so replace with 9999 to make these areas not crossable
                   c(1,NA), 9999) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network becaus of color-scale. reduce spread.limit to make it visible again.
gbm.r.inv <- asgrid( # convert to grid for spread function
  gbm.r.inv,
  xll = xmin(gbm.r.inv),
  yll = ymin(gbm.r.inv),
  cellsize = 1000
) # update cellsize with the aggregate factor * 1000m

# rename for parametrization loop:
power.1.5 <- gbm.r.inv


gbm.r <- rast("data/simulation input data/senecio inaequidens/biomod2_GBM.2PA.average.bio1.bio12.nitrogen.traffic.LC.tif") %>%  # biomod SDM output from boosted regression trees (gbm)
  project(ref.raster)
names(gbm.r) <- "layer"
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r <- terra::mask(gbm.r, ref.raster)
gbm.r <- subst(gbm.r, NA, 0)

mask.thresh <- 0.35
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(ger)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY # deleted because if I do this, then the grid allows to cross borders over time
  cells()                               # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r ^ 2 # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r[mask] <- 0
#x11()
#plot(gbm.r)

gbm.r.inv <- gbm.r * -1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

gbm.r.inv <- subst(gbm.r.inv, # must not have NAs for the function below, so replace with 9999 to make these areas not crossable
                   c(1,NA), 9999) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network becaus of color-scale. reduce spread.limit to make it visible again.
gbm.r.inv <- asgrid(
  gbm.r.inv, # convert to grid for spread function
  xll = xmin(gbm.r.inv),
  yll = ymin(gbm.r.inv),
  cellsize = 1000
) # update cellsize with the aggregate factor * 1000m

# rename for parametrization loop:
power.2 <- gbm.r.inv

#
#
#
#
#


#
#
#
### plot to check if all data align and see what they look like. ---------------
plot(gbm.r)

ref.p %>% dplyr::filter(decade <= 2025) %>% 
  plot(add = TRUE, col = "yellow",
       pch = 19)
ref.p %>% dplyr::filter(decade <= 1999) %>% 
  plot(add = TRUE, col = "orange",
       pch = 19)
plot(sen.ini[, 1],
     add = TRUE,
     col = "red",
     pch = 19)
#
plot(ref.p[, 1],
     add = TRUE ,
     col = "yellow",
     pch = 19)

#
#


### read in traffic network ----------------------------------------------------
eu.links <- st_read("data/traffic data/1031.1165ua.GHS.800km.gpkg", # this is a shapefile with all least-cost paths (i.e. open street map routes) between all pairs of urban areas
                    layer = "1031.1165ua.GHS.800km.exp.predict.w.original.data")
distances <- eu.links$length
eu.links$link.id <- 1:nrow(eu.links) # paths need an ID for easier use later
eu.links <- eu.links %>%
  rename(o.ID = start.ID)
eu.links <- eu.links %>%
  rename(d.ID = dest.ID)

# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters:
eu.links <- eu.links[lengths(st_intersects(eu.links, ger)) > 0, ]
eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links) # geometries not needed for use in the dispersal() function.
#
#
#


#
#
#
### read in urban areas --------------------------------------------------------
nodes <- st_read("data/environmental data/eur.urban.areas.gpkg", layer = "GHS.29.10.2025.points")
nodes <- nodes[lengths(st_intersects(nodes, ger)) > 0, ]# crop the nodes to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters

nodes$occ <- terra::extract(ini.dist.r, vect(nodes))[, 2]
ini.nodes <- nodes %>% dplyr::filter(occ == 1)
ini.nodes <- ini.nodes$ID
#
#
#

spread.val <- 1
thresh.disp.factor <- 0.5
time.steps <- 24
agg.acc.fact <- 1
acc.vect <- st_union(st_buffer(ref.p, 30000))
min.tr <- 0.6011832
min.tr = quantile(eu.links$predicted, probs = 0.05, na.rm = TRUE)[[1]]
max.dist <- quantile(distances, probs = 0.5, na.rm = TRUE)


ext(ref.dist.r) == ext(empty.r)

## run function ----------------------------------------------------------------
out <- par.dispersal(
  land.spread = TRUE,
  net.spread = TRUE,
  spread.val = spread.val,
  thresh.disp.factor = thresh.disp.factor, 
  time.steps = time.steps,
  dist.ini = sen.ini,
  ini.nodes = ini.nodes,
  gbm.r.inv = get("power.1.5"),
  name.gbm.r.inv = "power.1.5",
  
  ref.raster = empty.r,
  result.r = empty.r,
  ref.dist.r = ref.dist.r,
  
  plot.result = TRUE, 
  sample.nodes.from.raster = TRUE,
  unsuitability.mask = mask,
  acc.vect = acc.vect,
  min.tr = min.tr,
  max.dist = max.dist,
  cores = 6
)

x11()
par(mfrow = c(2,2))
plot(ini.dist.r, 
     main = "initial distribution", 
     background = "darkgrey")
plot(out[[1]], 
     main = paste(
       "s.v =", spread.val, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "transformation ^5"), 
     background = "darkgrey")
plot(out[[1]], 
     main = paste(
       "s.v =", spread.val, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "ini.points = green ;",
       "\nreached ua (sim.) = red"), 
     background = "darkgrey")
plot(sen.ini, add = TRUE, col = "darkgreen", pch = 19) # starting points
plot(dplyr::filter(nodes, ID %in% out[[2]])[,1], add = TRUE, col = "red", pch = 19) # reached ua
plot(ref.dist.r[[nlyr(ref.dist.r)]], 
     main = "final reference distribution", 
     background = "darkgrey")
out[[3]]


#
#
#
#
#
## parameter estimation --------------------------------------------------------
### without network to copy to other core ====
parameters <- tidyr::crossing(
  spread.val = c(1, 1.5, 2),
  thresh.disp.factor = c(0.5, 0.75, 1),
  #min.tr = c(0.1, 0.5, 0.95), # irrelevant with net.spread = FALSE, greatly reduces the number of parameter combinations if off
  transformation = c("power.1", "power.1.5", "power.2"),
  net.spread = c(FALSE) # TRUE on other core in parallel
)

for(i.p in 1:nrow(parameters)){
  print(i.p)
  spread.val <- parameters[i.p,]$spread.val
  thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
  time.steps <- 40
  #min.tr <- quantile(eu.links$predicted, probs = c(parameters[i.p,]$min.tr), na.rm = TRUE)[[1]] 
  min.tr <- quantile(eu.links$predicted, probs = 0.99, na.rm = TRUE)[[1]] # traffic not relevant with net.spread = FALSE
  max.dist <- 1000000
  gbm.r.inv <- get(parameters[i.p,]$transformation)
  
  out <- dispersal(
    land.spread = TRUE,
    net.spread = parameters[i.p,]$net.spread,
    spread.val = spread.val,
    thresh.disp.factor = thresh.disp.factor, 
    time.steps = time.steps,
    dist.ini = sen.ini,
    ini.nodes = ini.nodes,
    ref.raster = ref.raster,
    result.r = empty.r,
    ref.dist.r = ref.dist.r,
    plot.result = TRUE, 
    sample.nodes.from.raster = TRUE,
    unsuitability.mask = mask,
    acc.vect = acc.vect,
    min.tr = min.tr,
    max.dist = max.dist
  )
  
  if(i.p == 1){
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    optim.output <- out$accuracies
      } else {
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    optim.output <- bind_rows(optim.output, out$accuracies)
    }
}

optim.out.no.network <- optim.output

### with network (to copy to other core) ---------------------------------------
parameters <- tidyr::crossing(
  spread.val = c(1, 1.5, 2),
  thresh.disp.factor = c(0.5, 0.75, 1),
  min.tr = c(0.1, 0.5, 0.95), # irrelevant with net.spread = FALSE, greatly reduces the number of parameter combinations if off
  transformation = c("power.1", "power.1.5", "power.2"),
  net.spread = c(TRUE)
)

for(i.p in 1:nrow(parameters)){
  print(i.p)
  spread.val <- parameters[i.p,]$spread.val
  thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
  time.steps <- 40
  min.tr <- quantile(eu.links$predicted, probs = c(parameters[i.p,]$min.tr), na.rm = TRUE)[[1]] 
  #min.tr <- quantile(eu.links$predicted, probs = 0.99, na.rm = TRUE)[[1]] # traffic not relevant with net.spread = FALSE
  max.dist <- 1000000
  gbm.r.inv <- get(parameters[i.p,]$transformation)
  
  out <- dispersal(
    land.spread = TRUE,
    net.spread = parameters[i.p,]$net.spread,
    spread.val = spread.val,
    thresh.disp.factor = thresh.disp.factor, 
    time.steps = time.steps,
    dist.ini = sen.ini,
    ini.nodes = ini.nodes,
    ref.raster = ref.raster,
    result.r = empty.r,
    ref.dist.r = ref.dist.r,
    plot.result = TRUE, 
    sample.nodes.from.raster = TRUE,
    unsuitability.mask = mask,
    acc.vect = acc.vect,
    min.tr = min.tr,
    max.dist = max.dist
  )
  
  if(i.p == 1){
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    optim.output <- out$accuracies
    } else {
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    optim.output <- bind_rows(optim.output, out$accuracies)
    }
}

optim.out.with.network <- optim.output

head(optim.out.with.network)
head(optim.out.no.network)

out.full <- bind_rows(optim.out.no.network,
                      optim.out.with.network)
#write.csv(out.full,
#          "data/simulation output/senecio inaequidens/out.full.csv",
#          row.names = FALSE)

# show results with and without network spread ---------------------------------
library(ggplot2)

out.full<- data.table::fread("data/simulation output/senecio inaequidens/out.full.csv")
out.full$net.spread <- as.character(out.full$net.spread)

crit.F <- .975
#
#
## results with network --------------------------------------------------------
out.with.network <- out.full %>% 
  dplyr::filter(net.spread == "TRUE")
out.no.network <- out.full %>% 
  dplyr::filter(net.spread == "FALSE")

par(mfrow = c(1,2))
plot(F.ref.1989 ~ time.step, 
     data = out.with.network,
     col = alpha("black", .1),
     ylim = c(0,1),
     main = "with traffic network",
     ylab = "F-score")
points(F.ref.1999 ~ time.step, 
       data = out.with.network,
       col = alpha("orange",.1))
points(F.ref.2009 ~ time.step, 
       data = out.with.network,
       col = alpha("red",.1))

legend("topright", c("1989", "1999", "2009"), border="black", fill = c("black", "orange", "red"))

for(n.r in 1:81) {
  if(n.r == 1){
    r.id <- paste(rep(paste0("run.", n.r), 40))
  } else {
    x <- paste(rep(paste0("run.", n.r), 40))
    r.id <- c(r.id, x)
  }
}
out.with.network$run <- r.id

summary.output <- data.frame(F.ref.1989.ts = 1:81,
                             F.ref.1999.ts = 1:81,
                             F.ref.2009.ts = 1:81,
                             F.ref.1989.F = 1:81,
                             F.ref.1999.F = 1:81,
                             F.ref.2009.F = 1:81,
                             run = paste0("run.", 1:81))
for(i.d in 1:length(unique(out.with.network$run))){
  slice <- out.with.network %>% 
    dplyr::filter(run == unique(out.with.network$run)[i.d])
  summary.output$F.ref.1989.ts[i.d] <- slice[min(which(slice$F.ref.1989 >= crit.F * max(slice$F.ref.1989))),]$time.step
  summary.output$F.ref.1999.ts[i.d] <- slice[min(which(slice$F.ref.1999 >= crit.F * max(slice$F.ref.1999))),]$time.step
  summary.output$F.ref.2009.ts[i.d] <- slice[min(which(slice$F.ref.2009 >= crit.F * max(slice$F.ref.2009))),]$time.step
  
  summary.output$F.ref.1989.F[i.d] <- crit.F * max(slice$F.ref.1989)
  summary.output$F.ref.1999.F[i.d] <- crit.F * max(slice$F.ref.1999)
  summary.output$F.ref.2009.F[i.d] <- crit.F * max(slice$F.ref.2009)
  
  summary.output$run[i.d] <- unique(out.with.network$run)[i.d]
}
points(x = mean(summary.output$F.ref.1989.ts),
       y = mean(summary.output$F.ref.1989.F),
       col = "black",
       pch = 19,
       lwd = 2)
lines(x = c(mean(summary.output$F.ref.1989.ts) - sd(summary.output$F.ref.1989.ts),
            mean(summary.output$F.ref.1989.ts) + sd(summary.output$F.ref.1989.ts)),
      y = c(mean(summary.output$F.ref.1989.F), mean(summary.output$F.ref.1989.F)),
      col = "black",
      lwd = 2)
lines(y = c(mean(summary.output$F.ref.1989.F) - sd(summary.output$F.ref.1989.F),
            mean(summary.output$F.ref.1989.F) + sd(summary.output$F.ref.1989.F)),
      x = c(mean(summary.output$F.ref.1989.ts), mean(summary.output$F.ref.1989.ts)),
      col = "black",
      lwd = 2)

points(x = mean(summary.output$F.ref.1999.ts),
       y = mean(summary.output$F.ref.1999.F),
       col = "orange",
       pch = 19,
       lwd = 2)
lines(x = c(mean(summary.output$F.ref.1999.ts) - sd(summary.output$F.ref.1999.ts),
            mean(summary.output$F.ref.1999.ts) + sd(summary.output$F.ref.1999.ts)),
      y = c(mean(summary.output$F.ref.1999.F), mean(summary.output$F.ref.1999.F)),
      col = "orange",
      lwd = 2)
lines(y = c(mean(summary.output$F.ref.1999.F) - sd(summary.output$F.ref.1999.F),
            mean(summary.output$F.ref.1999.F) + sd(summary.output$F.ref.1999.F)),
      x = c(mean(summary.output$F.ref.1999.ts), mean(summary.output$F.ref.1999.ts)),
      col = "orange",
      lwd = 2)

points(x = mean(summary.output$F.ref.2009.ts),
       y = mean(summary.output$F.ref.2009.F),
       col = "red",
       pch = 19,
       lwd = 2)
lines(x = c(mean(summary.output$F.ref.2009.ts) - sd(summary.output$F.ref.2009.ts),
            mean(summary.output$F.ref.2009.ts) + sd(summary.output$F.ref.2009.ts)),
      y = c(mean(summary.output$F.ref.2009.F), mean(summary.output$F.ref.2009.F)),
      col = "red",
      lwd = 2)
lines(y = c(mean(summary.output$F.ref.2009.F) - sd(summary.output$F.ref.2009.F),
            mean(summary.output$F.ref.2009.F) + sd(summary.output$F.ref.2009.F)),
      x = c(mean(summary.output$F.ref.2009.ts), mean(summary.output$F.ref.2009.ts)),
      col = "red",
      lwd = 2)

t.test(summary.output$F.ref.1999.ts, 
       summary.output$F.ref.2009.ts, 
       paired = TRUE)

#
#
#
## results without network -----------------------------------------------------
plot(F.ref.1989 ~ time.step, 
     data = out.no.network,
     col = alpha("black", .1),
     ylim = c(0,1),
     main = "without traffic network",
     ylab = "F-score")
points(F.ref.1999 ~ time.step, 
       data = out.no.network,
       col = alpha("orange", .1))
points(F.ref.2009 ~ time.step, 
       data = out.no.network,
       col = alpha("red",.1))

legend("topright", c("1989", "1999", "2009"), border="black", fill = c("black", "orange", "red"))

for(n.r in 1:27) {
  if(n.r == 1){
    r.id <- paste(rep(paste0("run.", n.r), 40))
  } else {
    x <- paste(rep(paste0("run.", n.r), 40))
    r.id <- c(r.id, x)
  }
}
out.no.network$run <- r.id

summary.output <- data.frame(F.ref.1989.ts = 1:27,
                             F.ref.1999.ts = 1:27,
                             F.ref.2009.ts = 1:27,
                             F.ref.1989.F = 1:27,
                             F.ref.1999.F = 1:27,
                             F.ref.2009.F = 1:27,
                             run = paste0("run.", 1:27))
for(i.d in 1:length(unique(out.no.network$run))){
  slice <- out.no.network %>% 
    dplyr::filter(run == unique(out.no.network$run)[i.d])
  summary.output$F.ref.1989.ts[i.d] <- slice[min(which(slice$F.ref.1989 >= crit.F * max(slice$F.ref.1989))),]$time.step
  summary.output$F.ref.1999.ts[i.d] <- slice[min(which(slice$F.ref.1999 >= crit.F * max(slice$F.ref.1999))),]$time.step
  summary.output$F.ref.2009.ts[i.d] <- slice[min(which(slice$F.ref.2009 >= crit.F * max(slice$F.ref.2009))),]$time.step
  
  summary.output$F.ref.1989.F[i.d] <- crit.F * max(slice$F.ref.1989)
  summary.output$F.ref.1999.F[i.d] <- crit.F * max(slice$F.ref.1999)
  summary.output$F.ref.2009.F[i.d] <- crit.F * max(slice$F.ref.2009)
  
  summary.output$run[i.d] <- unique(out.no.network$run)[i.d]
}
points(x = mean(summary.output$F.ref.1989.ts),
       y = mean(summary.output$F.ref.1989.F),
       col = "black",
       pch = 19,
       lwd = 2)
lines(x = c(mean(summary.output$F.ref.1989.ts) - sd(summary.output$F.ref.1989.ts),
            mean(summary.output$F.ref.1989.ts) + sd(summary.output$F.ref.1989.ts)),
      y = c(mean(summary.output$F.ref.1989.F), mean(summary.output$F.ref.1989.F)),
      col = "black",
      lwd = 2)
lines(y = c(mean(summary.output$F.ref.1989.F) - sd(summary.output$F.ref.1989.F),
            mean(summary.output$F.ref.1989.F) + sd(summary.output$F.ref.1989.F)),
      x = c(mean(summary.output$F.ref.1989.ts), mean(summary.output$F.ref.1989.ts)),
      col = "black",
      lwd = 2)

points(x = mean(summary.output$F.ref.1999.ts),
       y = mean(summary.output$F.ref.1999.F),
       col = "orange",
       pch = 19,
       lwd = 2)
lines(x = c(mean(summary.output$F.ref.1999.ts) - sd(summary.output$F.ref.1999.ts),
            mean(summary.output$F.ref.1999.ts) + sd(summary.output$F.ref.1999.ts)),
      y = c(mean(summary.output$F.ref.1999.F), mean(summary.output$F.ref.1999.F)),
      col = "orange",
      lwd = 2)
lines(y = c(mean(summary.output$F.ref.1999.F) - sd(summary.output$F.ref.1999.F),
            mean(summary.output$F.ref.1999.F) + sd(summary.output$F.ref.1999.F)),
      x = c(mean(summary.output$F.ref.1999.ts), mean(summary.output$F.ref.1999.ts)),
      col = "orange",
      lwd = 2)

points(x = mean(summary.output$F.ref.2009.ts),
       y = mean(summary.output$F.ref.2009.F),
       col = "red",
       pch = 19,
       lwd = 2)
lines(x = c(mean(summary.output$F.ref.2009.ts) - sd(summary.output$F.ref.2009.ts),
            mean(summary.output$F.ref.2009.ts) + sd(summary.output$F.ref.2009.ts)),
      y = c(mean(summary.output$F.ref.2009.F), mean(summary.output$F.ref.2009.F)),
      col = "red",
      lwd = 2)
lines(y = c(mean(summary.output$F.ref.2009.F) - sd(summary.output$F.ref.2009.F),
            mean(summary.output$F.ref.2009.F) + sd(summary.output$F.ref.2009.F)),
      x = c(mean(summary.output$F.ref.2009.ts), mean(summary.output$F.ref.2009.ts)),
      col = "red",
      lwd = 2)


t.test(summary.output$F.ref.1999.ts, 
       summary.output$F.ref.2009.ts, 
       paired = TRUE)

# the F-score stays constant over time after it reached its maximum because the F-score calculation does not include True Negatives

out.full %>% dplyr::filter()

###
# end senecio
###
###
###
###
###


# TAPINOMA MAGNUM --------------------------------------------------------------
## data preparation ------------------------------------------------------------
# area of interest (the smaller the faster):
e <- ext(c(
  -11.8167128633999,
  25.9369450325551,
  35.5922736627615,
  54.2843855437946
))

### native and non.native areas ------------------------------------------------
# Info from Manuela's Database:
native <- c("Corse", "Italy", "Sardegna", "Sicily", "Spain") # these are known non.native countries; Sardegna = Sardinia, Corse = Corsica; "Morocco", "Tunisia","Algeria" left out because here only europe
non.native <- c("Belgium",
                "France",
                "Germany",
                "Netherlands",
                "Slovenia",
                "Switzerland") 
pot.countries <- c( # these are potential countries
  "Portugal",
  "United Kingdom",
  "Ireland",
  "Norway",
  "Sweden",
  "Finland",
  "Estonia",
  "Latvia",
  "Lithuania",
  "Poland",
  "Czechia",
  "Austria",
  "Croatia",
  "Greece",
  "Slovenia",
  "San Marino",
  "Bosnia and Herzegovina",
  "Serbia",
  "Montenegro",
  "Albania",
  "Bulgaria",
  "Romania",
  "Ukraine",
  "Moldova",
  "Belarus",
  "Luxembourg",
  "Denmark",
  "Turkey",
  "Hungary",
  "Iceland",
  "Sicily",
  "Slovakia",
  "Kosovo",
  "North Macedonia",
  "Montenegro"
)
# "Azerbaijan" left out because no data available, majority in Europe

# filter with GADM instead of years => most records are from after 2010, even in the native region so I try here with regional separation, not temporal
gadm.0 <- sf::st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_0")
gadm.0.native <- gadm.0 %>%
  dplyr::filter(COUNTRY %in% native)
gadm.0.n.native <- gadm.0 %>%
  dplyr::filter(COUNTRY %in% non.native)
gadm.0.pot <- gadm.0 %>%
  dplyr::filter(COUNTRY %in% pot.countries)


gadm.1 <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_1")
gadm.1.native <- gadm.1 %>%
  dplyr::filter(NAME_1 %in% native)
gadm.1.n.native <- gadm.1 %>%
  dplyr::filter(NAME_1 %in% non.native)
gadm.1.pot <- gadm.1 %>%
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1.add <- gadm.1 %>% dplyr::filter(
  COUNTRY %in% c("North Macedonia" , "Moldova", "Montenegro", "Cyprus") | # are not in gadm.2 level.....
    NAME_1 == "Kaliningrad"
) 

gadm.2 <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_2") %>%
  dplyr::filter(COUNTRY %in% c(native, non.native, pot.countries))
gadm.2 <- dplyr::bind_rows(gadm.2, gadm.1.add)
gadm.2.native <- gadm.2 %>%
  dplyr::filter(COUNTRY %in% native)
gadm.2.n.native <- gadm.2 %>%
  dplyr::filter(COUNTRY %in% non.native)
gadm.2.ref <- gadm.2
gadm.2.pot <- gadm.2 %>%
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.3 <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_3") %>%
  dplyr::filter(COUNTRY %in% c(native, non.native, pot.countries))
gadm.3 <- dplyr::bind_rows(gadm.3, gadm.1.add)


### occurrence data from B. Seifert --------------------------------------------
# transform excel data into a shapefile.
tapintro <- readxl::read_excel("data/species occurrence data/tapinoma magnum/TAPINTRO [Seifert].xlsx")
tapinoc <- readxl::read_excel("data/species occurrence data/tapinoma magnum/TAPINO_C [Seifert].xlsx") %>%
  dplyr::filter(HYP == "magn")

tap <- bind_rows(tapinoc, tapintro)

comp <- tibble(
  tap.mag = "tapinoma magnum",
  x = tap$LON,
  y = tap$LAT,
  year = tap$YEAR
) %>%
  na.omit()
comp$x <- as.numeric(comp$x)
comp$y <- as.numeric(comp$y)

### occurrence data from Destour et al 2024 ------------------------------------
des <- data.table::fread("data/species occurrence data/tapinoma magnum/tapinoma [Destour et al 2024].csv") |>
  dplyr::filter(Species == "T. magnum")
des <- tibble(
  tap.mag = "tapinoma magnum",
  x = des$longitude,
  y = des$latitude,
  year = des$`earliest sample date`)
comp <- rbind(comp, des)
comp <- st_as_sf(comp, coords = c("x", "y"), crs = 4326)


tap.mag.ini <- comp %>% # filter based on record day which occurrence records are used as initial distribution (i.e. starting points) in the simulation
  dplyr::filter(year <= 2007) #%>%
  #st_intersection(st_union(gadm.0.native, gadm.1.native)) # restricting it to native means that french atlantic coast will not be reached

st_erase = function(x, y){st_difference(x, st_union(st_combine(y)))}

tap.mag.n.native <- comp #%>% # create a subset of occurrence non.native records.
#st_as_sf() %>%
#st_intersection(gadm.0.n.native) %>%
#st_erase(gadm.0.native) %>%
#st_erase(gadm.1.native)
#
#
#
### read SDM map & make resistance raster for gridprocess::rawpsread() ---------
# biomod2, gbm = boosted regression trees, variables are in the file name (merraclim database):
# BIO1: Annual Mean Temperature
# BIO10: Mean Temperature of Warmest Quarter
# BIO11: Mean Temperature of Coldest Quarter
# BIO13: Precipitation of Wettest Month
# BIO14: Precipitation of Driest Month
# LC = copernicus land-cover classes
# pop.density = population density raster (copernicus GHS)
gbm.r <- rast(
  "data/simulation input data/tapinoma magnum/biomod2_GBM.bio10.bio13.bio14.LC.pop.density.tif"
)
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1

gbm.r <- aggregate(gbm.r, 5, # making the raster more coarse to avoid overly high resolution and reduce computing time; INCLUDE IN ASGRID() BELOW!!
                   fun = mean, na.rm = TRUE)

mask.thresh <- 0.35
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>% # gets a spatraster that has only cells where the suitability is equal or below the mask.thresh, this will be used as unsuitability mask in the dispersal function
  cells() # gets the cell numbers of these cells; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r ^ 1 # exponential conversion from habitat suitability to resistance instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 (again) irrespective of the transformation

gbm.r <- terra::mask(gbm.r, vect(gadm.0))
gbm.r.inv <- gbm.r * -1 + max(values(gbm.r[[1]]), na.rm = TRUE) # transform (i.e. invert) suitability raster to resistance raster.

gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below or it will crash, so replace with 9999 to make these areas not crossable
                   from = c(1,NA), 
                   to = 9999)

gbm.r.inv <- asgrid(
  gbm.r.inv,
  # convert raster to grid format for spread function
  xll = xmin(gbm.r.inf),
  yll = ymin(gbm.r.inv),
  cellsize = 5000
) # update cellsize with the aggregate factor * 1000m
power.1 <- gbm.r.inv
# creation of resistance layer for gridprocess::rawspread() done.
#
#
#

#
#
#
gbm.r <- rast(
  "data/simulation input data/tapinoma magnum/biomod2_GBM.bio10.bio13.bio14.LC.pop.density.tif"
)
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1

gbm.r <- aggregate(gbm.r, 5, # making the raster more coarse to avoid overly high resolution and reduce computing time; INCLUDE IN ASGRID() BELOW!!
                   fun = mean, na.rm = TRUE)

mask.thresh <- 0.35
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>% # gets a spatraster that has only cells where the suitability is equal or below the mask.thresh, this will be used as unsuitability mask in the dispersal function
  cells() # gets the cell numbers of these cells; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r ^ 2 # exponential conversion from habitat suitability to resistance instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 (again) irrespective of the transformation

gbm.r <- terra::mask(gbm.r, vect(gadm.0))
gbm.r.inv <- gbm.r * -1 + max(values(gbm.r[[1]]), na.rm = TRUE) # transform (i.e. invert) suitability raster to resistance raster.

gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below or it will crash, so replace with 9999 to make these areas not crossable
                   from = c(1,NA), 
                   to = 9999)

gbm.r.inv <- asgrid(
  gbm.r.inv,
  # convert raster to grid format for spread function
  xll = xmin(gbm.r.inf),
  yll = ymin(gbm.r.inv),
  cellsize = 5000
) # update cellsize with the aggregate factor * 1000m
power.2 <- gbm.r.inv
# creation of resistance layer for gridprocess::rawspread() done.
#
#
#

#
#
#
gbm.r <- rast(
  "data/simulation input data/tapinoma magnum/biomod2_GBM.bio10.bio13.bio14.LC.pop.density.tif"
)
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1

gbm.r <- aggregate(gbm.r, 5, # making the raster more coarse to avoid overly high resolution and reduce computing time; INCLUDE IN ASGRID() BELOW!!
                   fun = mean, na.rm = TRUE)

mask.thresh <- 0.35
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>% # gets a spatraster that has only cells where the suitability is equal or below the mask.thresh, this will be used as unsuitability mask in the dispersal function
  cells() # gets the cell numbers of these cells; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r ^ 3 # exponential conversion from habitat suitability to resistance instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 (again) irrespective of the transformation

gbm.r <- terra::mask(gbm.r, vect(gadm.0))
gbm.r.inv <- gbm.r * -1 + max(values(gbm.r[[1]]), na.rm = TRUE) # transform (i.e. invert) suitability raster to resistance raster.

gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below or it will crash, so replace with 9999 to make these areas not crossable
                   from = c(1,NA), 
                   to = 9999)

gbm.r.inv <- asgrid(
  gbm.r.inv,
  # convert raster to grid format for spread function
  xll = xmin(gbm.r.inf),
  yll = ymin(gbm.r.inv),
  cellsize = 5000
) # update cellsize with the aggregate factor * 1000m
power.3 <- gbm.r.inv
# creation of resistance layer for gridprocess::rawspread() done.
#
#
#


#
#
#
### read in traffic network ----------------------------------------------------
eu.links <- st_read("data/traffic data/1031.1165ua.GHS.800km.gpkg", # this is a shapefile with all least-cost paths (i.e. open street map routes) between all pairs of urban areas
                    layer = "1031.1165ua.GHS.800km.exp.predict.w.original.data")

distances <- eu.links$length

eu.links$link.id <- 1:nrow(eu.links) # paths need an ID for easier use later
eu.links <- eu.links %>%
  rename(o.ID = start.ID)
eu.links <- eu.links %>%
  rename(d.ID = dest.ID)

# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters:
eu.links <- eu.links[lengths(st_intersects(eu.links, 
                                           st_transform(st_as_sf(vect(e, crs = "EPSG:4326")),st_crs(eu.links)))) > 0, ]
eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links) # geometries not needed for use in the dispersal() function.
#
#
#

#
#
#
### read in urban areas --------------------------------------------------------
nodes <- st_read("data/environmental data/eur.urban.areas.gpkg", 
                 layer = "GHS.29.10.2025.points")
nodes <- nodes[lengths(st_intersects(nodes, st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), st_crs(nodes)))) > 0, ]# crop the nodes to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters
#
#
#

#
#
#
### define input ---------------------------------------------------------------
empty.r <- gbm.r
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything, is used in dispersal function to make the result raster

native.dist.r <- empty.r %>%
  terra::mask(vect(st_buffer(tap.mag.ini, 10000)), updatevalue = 1, inverse = TRUE)

ref.dist.r <- c(terra::mask(native.dist.r, vect(st_buffer(comp, 40000)), updatevalue = 1, inverse = TRUE),
                terra::mask(native.dist.r, vect(gadm.3[lengths(st_intersects(gadm.3, comp)) > 0, ]), updatevalue = 1, inverse = TRUE),
                terra::mask(native.dist.r, vect(gadm.2[lengths(st_intersects(gadm.2, comp)) > 0, ]), updatevalue = 1, inverse = TRUE)
)
names(ref.dist.r) <- c("points.40km", "gadm.3", "gadm.2")

ini.nodes <- nodes %>%
  dplyr::filter(gc_ucn_mai_2025 %in% c(nodes[lengths(st_intersects(nodes, tap.mag.ini)) > 0, ]$gc_ucn_mai_2025)) # nodes which overlap with the known initial distribution are used as ini.nodes (i.e. initial nodes)

add.nodes <- tibble(
  ID = nodes$ID,
  is.occupied = terra::extract(native.dist.r, vect(nodes), ID = FALSE)[[1]]
) # nodes which lie within the known distributional range are added (number of added nodes depends on the buffer used above to make this raster)
add.nodes <- add.nodes %>%
  dplyr::filter(is.occupied > 0)
ini.nodes <- unique(c(ini.nodes$ID, add.nodes$ID))

agg.acc.fact <-  1 # how much should the result and reference raster be aggregated to calculate the accuracy?
acc.vect <- st_union(st_buffer(comp[, 1], 300000)) #acc.vect = accuracy vector, in meter, inside this area, accuracy is calculated.
ext(ref.dist.r) == ext(empty.r) # check if data match.

spread.val <- 1 # budget for spread used in gridprocess::rawspread()
thresh.disp.factor <- 0.95 # threshold which has to be reached in a cell to be treated as occupied/presence
time.steps <- 60 # how many consecutive iterations?
min.tr <- quantile(eu.links$predicted, probs = c(0.85), na.rm = TRUE)[[1]] # probs defines which quantile of the traffic volumes is used as minimum traffic value that filter or paths which are used in the traffic network.
max.dist <- quantile(distances, probs = c(0.9), na.rm = TRUE)[[1]] # paths in the network longer than this will not be used, in meter.
#
#
#

#
#
#
## run function ----------------------------------------------------------------
out <- par.dispersal(
  land.spread = TRUE,
  net.spread = TRUE,
  dist.ini = tap.mag.ini,
  spread.val = spread.val,
  thresh.disp.factor = thresh.disp.factor,
  ini.nodes = ini.nodes,
  ref.raster = empty.r,
  result.r = empty.r,
  
  gbm.r.inv = get("power.1"),
  name.gbm.r.inv = "power.1",
  
  time.steps = time.steps,
  ref.dist.r = ref.dist.r,
  
  sample.nodes.from.raster = TRUE,
  unsuitability.mask = mask,
  #acc.vect = acc.vect,
  min.tr = min.tr,
  max.dist = max.dist,
  
  plot.result = TRUE,
  cores = 15
)


out[[4]] <- out[[4]] %>% 
  mask(empty.r) 

out[[1]] <- mask(out[[1]], vect(gadm.0))

par(mfrow = c(2,2))
plot(native.dist.r, 
     main = "initial distribution", 
     background = "darkgrey")
plot(out[[1]], 
     main = paste0(
       "s.v:", spread.val, "; ",
       "init: 1;",
       "\nmin.tr:", round(min.tr,2), "; ",
       "max.dist:", round(max.dist/1000,2), "km; ",
       "t.s:", time.steps, ";",
       "\nt.d.f:", thresh.disp.factor, "; ",
       "transformation ²"), 
     background = "darkgrey")
plot(out[[1]], 
     main = paste(
       "s.v =", spread.val, ";",
       "init = 1;",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "ini.points = green ;",
       "\nreached ua (sim.) = red"), 
     background = "darkgrey")
plot(tap.mag.ini, add = TRUE, col = "green") # starting points
plot(dplyr::filter(nodes, ID %in% c(out[[2]])), add = TRUE, col = "red")
plot(ref.dist.r, 
     main = "final reference distribution", 
     background = "darkgrey")
out[[3]]
x11()
plot(out[[4]][[52]], background = "darkgrey")
writeRaster(out[[4]], 
            "C:/Users/JLU-SU/Nextcloud/Predictive Aliens/data/simulation output/tapinoma magnum/02062026.tiff")
#
#
#


#
#
#
#
#
## parameter estimation --------------------------------------------------------

parameters <- tidyr::crossing(
  spread.val = c(1), # , 1.5, 2
  thresh.disp.factor = c(.9), #0.5, 0.75, 1
  min.tr = c(0.05, 0.5, 0.95),
  max.dist = c(0.05, 0.5, 0.95),
  transformation = c("power.1"), # , "power.2", "power.3"
  net.spread = c(TRUE) # with TRUE on other core
)


for(i.p in 1:nrow(parameters)){
  print(i.p)
  spread.val <- parameters[i.p,]$spread.val
  thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
  time.steps <- 21
  min.tr <- quantile(eu.links$predicted, probs = c(parameters[i.p,]$min.tr), na.rm = TRUE)[[1]]
  #min.tr <- 0.5 # derived from the reference distribution (see below)
  max.dist <- quantile(eu.links$length, probs = c(parameters[i.p,]$length), na.rm = TRUE)[[1]]
  gbm.r.inv <- get(parameters[i.p,]$transformation)
  
  
  out <- dispersal(
    land.spread = TRUE,
    net.spread = parameters[i.p,]$net.spread,
    spread.val = spread.val,
    thresh.disp.factor = thresh.disp.factor, 
    time.steps = time.steps,
    dist.ini = tap.mag.ini,
    ini.nodes = ini.nodes,
    ref.raster = empty.r,
    result.r = empty.r,
    ref.dist.r = ref.dist.r,
    plot.result = TRUE, 
    sample.nodes.from.raster = TRUE,
    unsuitability.mask = mask,
    acc.vect = acc.vect,
    min.tr = min.tr,
    max.dist = max.dist
  )
  
  if(i.p == 1){
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    optim.output <- out$accuracies
    # add power here
  } else {
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    optim.output <- bind_rows(optim.output, out$accuracies)
    # add power here
  }
}

#write.csv(optim.output,
#          "data/simulation output/tapinoma magnum/optim.output.net.spread.FALSE.csv", 
#          row.names = FALSE)

no.net <- data.table::fread("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/data/simulation output/tapinoma magnum/optim.output.net.spread.FALSE.csv")
net <- data.table::fread("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/data/simulation output/tapinoma magnum/optim.output.net.spread.TRUE.csv")

all <- rbind(no.net, net)
View(all)
#
#
#
#
###
# end tapinoma
###
###
###
###


# MYOCASTOR COYPUS -------------------------------------------------------------
## data preparation ------------------------------------------------------------
### native and non.native areas ------------------------------------------------
# Info from Manuela's Database:
native <- c("Argentina",
            "Bolivia",
            "Brazil",
            "Chile",
            "Paraguay",
            "Uruguay") # Sardegna = Sardinia, Corse = Corsica; "Morocco", "Tunisia","Algeria" left out because here only europe
non.native <- c(
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
  st_make_valid()

gadm.0 <- gadm.0 %>%
  dplyr::filter(COUNTRY %in% c(non.native, native))
gadm.0.native <- gadm.0 %>%
  dplyr::filter(COUNTRY %in% native)
gadm.0.n.native <- gadm.0 %>%
  dplyr::filter(COUNTRY %in% non.native)

sf_use_s2(FALSE)
gadm.2 <- st_read("data/environmental data/world_gadm_410-levels.gpkg", layer = "ADM_2") %>%
  dplyr::filter(COUNTRY %in% c(native, non.native)) %>%
  st_make_valid() %>%
  st_buffer(0)
sf_use_s2(TRUE)

gadm.2.native <- gadm.2 %>%
  dplyr::filter(COUNTRY %in% native)
gadm.2.n.native <- gadm.2 %>%
  dplyr::filter(COUNTRY %in% non.native)

waterways.add <- rast("data/environmental data/waterways.raster.tiff")
gbm.r <- rast(
  "data/simulation input data/myocastor coypus/biomod2_GBM.update.bio02.bio6.bio10.bio15.LC.pop.gbif and anna.tif"
)
### occurrence data (anna & gbif) ----------------------------------------------
myo.coy <- read.csv("data/species occurrence data/myocastor coypus/M coypus gbif 13082025.csv",
                    sep = "\t")
myo.coy <- tibble(myocastor = 1,
                  longitude = myo.coy$decimalLongitude,
                  latitude = myo.coy$decimalLatitude,
                  year = myo.coy$year
) %>%
  na.omit()

a.s.rec <- st_read("data/species occurrence data/myocastor coypus/nutria_occurrence_records_1980_2018 anna schertler.gpkg")###
points <- st_centroid(a.s.rec) %>% 
  st_transform("epsg:4326")

points.add <- tibble(myocastor = 1,
                     longitude = st_coordinates(points)[,1],
                     latitude = st_coordinates(points)[,2],
                     year = points$year) %>% 
  na.omit()

myo.coy.all <- dplyr::bind_rows(myo.coy, points.add)
myo.coy.xy <- st_as_sf(myo.coy.all,
                       coords = c("longitude", "latitude"),
                       crs = 4326)
myo.coy.xy.ini <- myo.coy.xy %>%
  dplyr::filter(year <= 1980)

one.rec.per.grid.cell <- function(ref.raster, occurences) {
  cells <- cellFromXY(ref.raster, occurences) %>% unique()
  coords <- xyFromCell(ref.raster, cells)
}

coords.matrix.ini <- one.rec.per.grid.cell(
  ref.raster = gbm.r,
  occurences = cbind(
    lon = st_coordinates(myo.coy.xy.ini)[, 1],
    lat = st_coordinates(myo.coy.xy.ini)[, 2]
  )
) %>%
  na.omit()

myo.coy.xy.ini <- st_as_sf(as_tibble(coords.matrix.ini),
                           coords = c("x", "y"),
                           crs = 4326)

coords.matrix <- one.rec.per.grid.cell(
  ref.raster = gbm.r,
  occurences = cbind(lon = st_coordinates(myo.coy.xy)[, 1], lat = st_coordinates(myo.coy.xy)[, 2])
) %>%
  na.omit()
myo.coy.xy <- st_as_sf(as_tibble(coords.matrix),
                       coords = c("x", "y"),
                       crs = 4326)

### read SDM map & make resistance raster for gridprocess::rawpsread() ---------
names(gbm.r) <- "layer"
gbm.r <- terra::mask(gbm.r, vect(gadm.0))

i.mask <- gbm.r
i.mask[i.mask != 0] <- 1

waterways.add <- mask(waterways.add, i.mask, inverse = FALSE)

gbm.r <- subst(gbm.r, NA, 0)
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
waterways.add[waterways.add != 0] <- 0.4 # to adjust if needed, also higher than the cut off threshold below
gbm.r <- gbm.r + waterways.add
gbm.r <- subst(gbm.r, NA, 0)
gbm.r <- aggregate(gbm.r, 5, fun = mean, na.rm = TRUE) # making the raster more coarse; INCLUDE IN ASGRID() BELOW!!

mask.thresh <- 0.5 
# gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)


gbm.r[gbm.r <= mask.thresh] <- 0
# with the extreme cut.off value of 0.5 it is probably not necessary to make a nonlinear conversion
# to reconstruct the final reference distribution
gbm.r <- gbm.r^1 # exponential conversion instead of linear.

mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  # terra::mask(vect(gadm.0)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY
  cells()      

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation

gbm.r.inv <- gbm.r*-1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix
#spread.limit <- limit*100 # set value for areas which are not crossable
gbm.r.inv <- subst(gbm.r.inv, # must not have NAs for the function below, so replace with 9999 to make these areas not crossable
                   c(1,NA), 9999) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network because of color-scale. reduce spread.limit to make it visible again.
gbm.r.inv <- asgrid( # convert to grid for spread function
  gbm.r.inv,
  xll = xmin(gbm.r.inv),
  yll = ymin(gbm.r.inv),
  cellsize = 1000
) # update cellsize with the aggregate factor * 1000m
power.1 <- gbm.r.inv
#
#
#
### plot to check if all data align and see what they look like. ---------------
plot(gbm.r)
plot(myo.coy.xy.ini[, 1],
     add = TRUE ,
     col = "yellow",
     pch = 19)
plot(myo.coy.xy[, 1],
     add = TRUE,
     col = "red",
     pch = 19)
#
#
#

#
#
#
### read in traffic network ----------------------------------------------------
eu.links <- st_read("data/traffic data/1031.1165ua.GHS.800km.gpkg", # this is a shapefile with all least-cost paths (i.e. open street map routes) between all pairs of urban areas
                    layer = "1031.1165ua.GHS.800km.exp.predict.w.original.data")
distances <- eu.links$length
eu.links$link.id <- 1:nrow(eu.links) # paths need an ID for easier use later
eu.links <- eu.links %>%
  rename(o.ID = start.ID)
eu.links <- eu.links %>%
  rename(d.ID = dest.ID)

# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters:
eu.link.IDs <- eu.links[lengths(st_intersects(eu.links, gadm.0)) > 0, ]$link.id # no need to use gadm.2 as well, all are here
eu.links <- eu.links %>% dplyr::filter(link.id %in% eu.link.IDs)

eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links) # geometries not needed for use in the dispersal() function.
#
#
#

#
#
#
### read in urban areas --------------------------------------------------------
nodes <- st_read("data/environmental data/eur.urban.areas.gpkg", 
                 layer = "GHS.29.10.2025.points")
nodes <- nodes[lengths(st_intersects(nodes, gadm.0)) > 0, ] # crop the nodes to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters
# no need to use gadm.2, all are here
#
#
#


#
#
#
### define input ---------------------------------------------------------------
empty.r <- gbm.r 
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything

ini.dist.r <- empty.r %>% 
  terra::mask(vect(st_buffer(myo.coy.xy.ini, 10000)),
              updatevalue = 1,
              inverse = TRUE)

ref.dist.r <- ini.dist.r %>% # used as reference to calculate the accuracy of the dispersal()-output.
  terra::mask(vect(st_buffer(myo.coy.xy, 10000)),
              updatevalue = 1,
              inverse = TRUE) 

nodes$occ <- terra::extract(ini.dist.r, vect(nodes))[,2]
ini.nodes <- nodes %>% 
  dplyr::filter(occ == 1)
ini.nodes <- ini.nodes$ID
#
#
#

### run function ---------------------------------------------------------------
ext(ref.dist.r) == ext(empty.r)

ini.nodes <-  ini.nodes
spread.val <- 1
nodes.cut.off <- 0.1
time.steps <- 40
thresh.disp.factor <- 0.25
agg.acc.fact <- 40
acc.vect <- gadm.0
min.tr <- quantile(eu.links$predicted, probs = c(0.50), na.rm = TRUE)[[1]] # probs defines which quantile of the traffic volumes is used as minimum traffic value that filter or paths which are used in the traffic network.
max.dist <- quantile(distances, probs = c(0.99), na.rm = TRUE)[[1]] # irrelevant anyways if net.spread = FALSE, paths in the network longer than this will not be used, in meter.
#

## run function ----------------------------------------------------------------
#out <- dispersal(
#  land.spread = TRUE,
#  net.spread = FALSE,
#  dist.ini = myo.coy.xy.ini,
#  spread.val = spread.val,
#  thresh.disp.factor = thresh.disp.factor,
#  ini.nodes = ini.nodes,
#  ref.raster = empty.r,
#  result.r = empty.r, 
#  gbm.r.inv = "power.1",
#  
#  time.steps = time.steps,
#  ref.dist.r = ref.dist.r,
#  
#  sample.nodes.from.raster = TRUE,
#  unsuitability.mask = mask,
#  
#  acc.vect = acc.vect,
#  min.tr = min.tr,
#  max.dist = max.dist,
#  plot.result = TRUE
#)

out <- par.dispersal(
  land.spread = TRUE,
  net.spread = FALSE,
  dist.ini = myo.coy.xy.ini,
  spread.val = spread.val,
  thresh.disp.factor = thresh.disp.factor,
  ini.nodes = ini.nodes,
  ref.raster = empty.r,
  result.r = empty.r, 
  gbm.r.inv = get("power.1"),
  name.gbm.r.inv = "power.1",
  
  time.steps = time.steps,
  ref.dist.r = ref.dist.r,
  
  sample.nodes.from.raster = TRUE,
  unsuitability.mask = mask,
  
  acc.vect = acc.vect,
  min.tr = min.tr,
  max.dist = max.dist,
  plot.result = TRUE,
  cores = 10
)

par(mfrow = c(2,2))
plot(mask(ini.dist.r, vect(gadm.0)), 
     main = "initial distribution", 
     background = "darkgrey")
plot(mask(out[[1]],vect(gadm.0)), 
     main = paste(
       "s.v =", spread.val, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "transformation ²"), 
     background = "darkgrey")
plot(mask(out[[1]],vect(gadm.0)), 
     main = paste(
       "s.v =", spread.val, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "ini.points = green ;",
       "\nreached ua (sim.) = red"), 
     background = "darkgrey")
plot(myo.coy.xy.ini, add = TRUE, col = "green") # starting points
plot(out[[2]], add = TRUE, col = "red") # reached ua
plot(mask(ref.dist.r, vect(gadm.0)), 
     main = "final reference distribution", 
     background = "darkgrey")
# plot(out[[2]], add = TRUE, col = "red", pch = 19) # reached ua
#summary(out[[2]]$tot.prop.traff)
out[[3]]
plot(mask(out[[4]][[31]],vect(gadm.0)),
     background = "darkgrey")

#
#
#
#
#
###
# end myocastor
###
###
###
###


#
#
#
#
# ADDITIONAL STEPS -------------------------------------------------------------
## plotting with tmap ----------------------------------------------------------
library(tmap)
out.pol <- as.polygons(out[[4]][[52]]) %>%
  st_as_sf()
plot(st_as_sf(out.pol))

reference.pol <- as.polygons(ref.dist.r$points.10km) %>%
  st_as_sf()

out.pol$layer <- as.character(out.pol$layer)
reference.pol$layer <- as.character(reference.pol$layer)

tm_shape(reference.pol) +
  tm_polygons(fill = "points.10km", 
              col = "black",
              fill.scale = tm_scale(
                values = c("0" = "lightgrey", "1" = "#005AB5")),
              fill_alpha = 1,
              fill.legend = tm_legend_hide()) +
  tm_shape(out.pol) +
  tm_polygons(fill = "layer", 
              col = "black",
              fill.scale = tm_scale(
                values = c("0" = rgb(0, 0, 0, alpha = 0), "1" = "#DC3220")),
              fill_alpha = 0.4,
              fill.legend = tm_legend_hide()) +
  tm_layout(bg.color = "lightblue")

tm_shape(st_as_sf(as.polygons(native.dist.r))) +
  tm_polygons(fill = "layer", 
              col = "black",
              fill.scale = tm_scale(
                values = c("0" = "lightgrey", "1" = "#005AB5")),
              fill_alpha = 1,
              fill.legend = tm_legend_hide()) +
  tm_layout(bg.color = "lightblue")
#
#
#
