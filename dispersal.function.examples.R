#
#
#
#
#
# set WD:
setwd("C:/Users/JLU-SU/Nextcloud/Predictive Aliens/")

# dispersal simulation function ====

# load packages: 
library(terra)
library(dplyr)
library(sf)
library(gridprocess) # devtools::install_github("ethanplunkett/gridprocess")

dispersal <- function(land.spread = TRUE, # is spread through the landscape allowed?
                      net.spread = FALSE, # is spread through the network allowed?
                      
                      dist.ini, # initial distribution coordinates from where the species spreads through the landscape
                      spread.val, # how far can the species spread in gridprocess::rawspread?
                      thresh.disp.factor, # how much of the spread.val must a pixel receive to be treated as occupied?
                      ini.nodes, # initial urban areas from which the species can spread through the traffic network AND the landscape
                      #thresh.nodes, # not used anymore: threshold for which cells from spread are treated as origin in following step? 
                      ref.raster, # reference raster (resolution etc.) for the output
                      result.r, # empty raster that will be updated after each time.step and thus turned into the result raster.
                      initiation = 1, # the initial traffic budget that each used node gets per time.step. this will then be distributed among all outgoing paths relative to the traffic flow on each path.
                      # I tried to scale this with gdp of the respecive node, but that did not improve the output. OPEN FOR DISCUSSION.
                      
                      nodes.cut.off, # not used anymore, perhaps OPEN FOR DISCUSSION:
                      # states how much traffic urban area "y" has to receive from urban area "x" to be treated as invaded. 
                      # each urban area gets an outgoing traffic budget of 1 and this budget is then distributed among all links going from this urban 
                      # area according to the traffic volume on this link. if the total traffic volume on all links going from urban area x is 1000 and the 
                      # link to urban area y has a volume of 10, then urban area y receives 1/1000*10 = 0.01 and it is checked whether or not this is above 
                      # the threshold for being treated as invaded. 
                      
                      time.steps, # how many steps should the simulation run?
                      dist.red = FALSE, 
                      # in each time.step an output raster is generated with areas which are occupied and from these areas the border cells are usedas starting points for further landscape spread in the following step.
                      # if dist.red = TRUE then only every xth cell is used as starting point with x = dist.red.param (next argument). greatly reduces computing time BUT MIGHT CHANGE OUTPUT PATTERNS!
                      dist.red.param, # define how strongly the initial cells should be reduced, numeric; 
                      # e.g. if 2, then every second cell/coordinate pair is chosen from the ini.update matrix
                      
                      ref.dist.r, # reference raster with the final known distribution with cells being present (1) or absent (0), used for accuracy calculation
                      # (accuracy <- (True Positives + True Negatives) / (True Positives + False Positives + False Negatives + True Negatives))
                      # the result.r output raster will be compared with this and the values needed for accuracy calculation derived from the comparison.
                      
                      plot.result = TRUE, # does not automatically plot anything, just makes an output raster if TRUE
                      
                      sample.nodes.from.raster = TRUE, # if TRUE then urban areas which are within the occupied areas on the result raster will be treated as 
                      # occupied and can serve as starting point in the traffic network. this allows a species to enter the traffic network from the landscape.
                      
                      regional.approach = FALSE, # not recommended. here calculations are based on gadm-0 region levels, but gadm will have to be provided and it takes ages...
                      regional.approach.acc = FALSE, # calculates accuracy not on raster cells but based on gadm regions. however, i recommend to provided an aggregated
                      # raster to use for calculation since it is much faster and ecologically more meaningful than usein administrative boundaries
                      # left here in case there is room for improvement to make this approach usable.
                      
                      unsuitability.mask = NULL, # if provided then every cell which is provided will be set to unoccupied at end of each timestep. can be used
                      # to make sure that unsuitable cells will not become occupied. they can still be crossed though.
                      
                      acc.vect, # accuracy vector; area inside which the accuracy will be calulated.
                      agg.acc.fact = 1, # how much should the output be aggregated (too avoid overly fine-scaled output and accuracy calculations)
                      min.tr, # traffic network: minimum traffic volume that paths must have to be used at all.
                      max.dist # traffic network: maximum length that a path can have to be used (longer distances will not be traveled)
) {
  a.int <- Sys.time() # just to keep track of how much time the simulation needs.
  
  eu.links <- eu.links %>% # filter for all paths which are connected with the chosen start.node
    dplyr::filter(predicted >= min.tr) %>% # filter all paths that have at least the given traffic volume
    dplyr::filter(length <= max.dist) # filter all paths which are up to the max.dist length
  
  for(t.s in 1:time.steps){ # iterate over time-steps
    
    # start spread through landscape:
    if(land.spread == TRUE){
      thresh.disp = spread.val * thresh.disp.factor # the threshold which has to be reached with the provided movement budget. depends on budget, resistance and distance to starting point.
      if(t.s == 1){
        id.ini <- rowColFromCell(ref.raster, cellFromXY(ref.raster, st_coordinates(dist.ini))) # gets the row and column id of the start sites, needed for gridprocess::spread function
      } else {
        id.ini <- id.ini.update # if it is a later time.step than 1, not the initial but the reached points from the previous time.step are used as starting points.
        id.ini <- na.omit(id.ini)} # somehow there is always a very small number of cells which have NA as rows and column numbers. maybe it is because of the reduction or because there is a NA value in the cell? I don't know, but this na.omit()-call prevents the function from crashing.
      
      empty.r <- ref.raster
      empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything
      
      for(i in 1:nrow(id.ini)){ # loop through all start-sites and determine where the species spreads to from here through the resistance layer.
        if(i == 1){spread.m <- rawspread(x = gbm.r.inv$m,
                                         spread.value = spread.val,
                                         row = id.ini[i,1],
                                         col = id.ini[i,2]#,
                                         #sd = sd # sd = bandwidth, not used here. "In the standard Gaussian kernel, the “bandwidth” which controls the spread of the kernel is equal to one standard deviation and accounts for 39% of the kernel volume." from: doi:10.1007/s10980-018-0653-9 
        )}else{
          spread.m.c <- rawspread(x = gbm.r.inv$m, # if it is not the first starting point then the output of this one will be added to the previous one, so that all reached cells are collected in one output raster. thus, the threshold can also be reached if a cell is reached just so from many starting locations.
                                  spread.value = spread.val,
                                  row = id.ini[i,1],
                                  col = id.ini[i,2]#,
                                  #sd = sd # sd = bandwidth, not used here. "In the standard Gaussian kernel, the “bandwidth” which controls the spread of the kernel is equal to one standard deviation and accounts for 39% of the kernel volume." from: doi:10.1007/s10980-018-0653-9 
          )
          
          spread.m <- spread.m + spread.m.c # summarizes all consecutive steps (a cell can be reached from different source cells)
        }}
      
      spread.m.r <- rast(spread.m, # the output matrix is converted to raster for plotting; might be possible to delete this step to speed up process?
                         extent = ext(ref.raster))
      
      spread.m.r[spread.m.r < thresh.disp] <- NA # below chosen threshold are considered as UNOCCUPIED cells
      spread.m.r[spread.m.r >= thresh.disp] <- 1 # below chosen threshold are considered as OCCUPIED cells
      crs(spread.m.r) <- crs(ref.raster) # has no crs after conversion from matrix to raster, so needs the reference crs.

      if(t.s == 1){ # all reached raster-cells are marked as such in the final result.r. this is updated at the end of each landscape spread iteration.
        dist.r <- mask(empty.r, spread.m.r, updatevalue = 1, inverse = TRUE)
        result.r <- mask(result.r, spread.m.r, updatevalue = 1, inverse = TRUE)
      } else {
        dist.r <- mask(dist.r, spread.m.r, updatevalue = 1, inverse = TRUE)
        result.r <- mask(result.r, spread.m.r, updatevalue = 1, inverse = TRUE)
      }
      # end spread through landscape
      #
      #
      #
      
      #
      #
      #
      # here the boundaries of the invaded areas are determined so that in a consecutive step, spread is only from here onwards to avoid that the function unnecessarily calculates spread from center areas again.
      dist.r <- result.r
      dist.r[is.na(dist.r)] <- 0 # sets the NA cells to 0 so that terra::boundaries can find the edges of the occupied patches below, otherwise it also identifies german border as edge (identifes all class differences between c(NA, 0, 1))
      b <- boundaries(dist.r, classes = TRUE, inner = FALSE)
      b.c <- cells(b, 1)[[1]] # identify boundary cells
      id.ini.update <- rowColFromCell(ref.raster, b.c) # in matrix: x,y ; lon, lat without resetting
      #
      #
      #
      
      #
      #
      #
      # optional: reduce number of starting-points for consecutive time-step and update initial coordinates for raster spread:
      # IMPORTANT: I HAD CASES WHERE THE REDUCTION DID NOT CHANGE THE RESULT AND OTHERS WHERE THE OUTPUT WAS VERY DIFFERENT.
      if(dist.red == TRUE){id.ini.update <- id.ini.update[c(seq(from = 0, 
                                                                to = nrow(id.ini.update), 
                                                                by = dist.red.param)), ]} else {} # takes only each nth cell (e.g. dist.red.param = 2 then each second cell)
      # end number reduction
    } else {} # better do nothing here.
    #
    #
    #
    
    #
    #
    #
    # start dispersal through the traffic network:
    if(net.spread == TRUE) { # net.spread = spread through traffic network (i.e. IDs of urban areas and the traffic flows between each pair)
      # start network spread:
      if(t.s == 1){} else {
        ini.nodes <- updated.ua.i} # updated.ua.i are the ones which were reached in a previous time.step.
      
      for(s.n in 1:length(ini.nodes)){ # loop over all ini.nodes
        start.node <- ini.nodes[s.n] # select the starting node of ini.nodes (s.n obviously has to be in starting.points table)
        dest.eu.links <- eu.links %>% 
          dplyr::filter(o.ID == start.node) # filter for all paths which start at the start.node
        
        dest.nodes <- nodes %>% 
          dplyr::filter(ID %in%
                          dest.eu.links$d.ID) # identify all nodes which are connected to start.node via the paths (i.e. via dest.eu.links).
        
        updated.ua.i <- dest.nodes$ID
        
        #
        #
        #
        # result generation for consecutive steps:
        # put into if statement above to use only if t.s > 1 to check if optim can better deal with it then
        id.ini.update <- rbind(id.ini.update, # these are used in consecutive landspread.
                               rowColFromCell(ref.raster, 
                                              cellFromXY(ref.raster, st_coordinates(dplyr::filter(nodes, ID %in% updated.ua.i)))
                               )
        )
        #
        #
        #
        
      }
      
      
      #
      #
      #
      #
      #
      # optional: sample nodes from the distributional landscape raster which lie in invaded territory - optional to do this. allows the species to switch from the landscape to the traffic network. the opposite way is always active but can be limited with spread.value and thresh.disp.fact.
      if(sample.nodes.from.raster == TRUE) {
        additional.nodes <- tibble(ID = nodes$ID, 
                                   is.occupied = terra::extract(result.r, vect(nodes), 
                                                                ID = FALSE)[[1]])
        additional.nodes <- additional.nodes[which(additional.nodes[,2] > 0), ]
        updated.ua.i <- sort(unique(c(updated.ua.i, additional.nodes$ID)))
        id.ini.update <- rbind(id.ini.update, 
                               rowColFromCell(ref.raster, 
                                              cellFromXY(ref.raster, st_coordinates(dplyr::filter(nodes, ID %in% updated.ua.i)))))
        
      } else {}
      
    } else {updated.ua.i <- nodes[ini.nodes,]} # end of if net.spread == TRUE
    #
    #
    #
    #
    #
    
    # accuracy measurement with random sub-samples of cells instead of whole raster to decrease the influence of the huge amount of 0s, analogue to pseudo-absence sampling:
    # currently not used but may be another option for accuracy calculation.
    #AUC.t <- tibble(cell = p.a.cells,
    #                ref = ref.dist.r[[1]][p.a.cells],
    #                pred = result.r[[1]][p.a.cells])
    #TP <- AUC.t %>% dplyr::filter(AUC.t[,2] == 1 &
    #                                AUC.t[,3] == 1) %>%
    #  nrow()
    #TN <- AUC.t %>% dplyr::filter(AUC.t[,2] == 0 &
    #                                AUC.t[,3] == 0) %>%
    #  nrow()
    # 
    #FN <- AUC.t %>% dplyr::filter(AUC.t[,2] == 1 &
    #                                AUC.t[,3] == 0) %>%
    #  nrow()
    #FP <- AUC.t %>% dplyr::filter(AUC.t[,2] == 0 &
    #                                AUC.t[,3] == 1) %>%
    #  nrow()
    #accuracy <- 1 - (TP + TN) / (TP + FP + FN + TN)
    

    #
    #
    #
    # optional: mask output with unsuitability mask (i.e. set cells with given IDs to 0)
    if(!is.null(unsuitability.mask)){ # if there is an unsuitability mask provided all the cells which have a habitat suitability below a certain threshold are set to unoccupied at the end of each iteration. this way, a species might cross an unsuitable habitat but it can not establich there an each population will be deleted at the end of each time.step.
      result.r[mask] <- 0 # the unsuitability mask has IDs of each cells which were determined as being unsuitable (making an unsuitability mask is a preparatory step and not part of this function).
    } else {}
    #
    #
    #
    
    #
    #
    #
    # option how to calculate accuracy:
    if(regional.approach.acc == TRUE){ # i tested this regional approach to determine the accuracy on a regional level instead of on a pixel level. However, the result is not convincing and the process is slow. is kept here just in case....
        ref.pol <- st_as_sf(as.polygons(result.r)) %>%
          dplyr::filter(layer == 1)
        
        gadm.2.ref.occ <- gadm.2.ref[lengths(st_intersects(gadm.2.ref, ref.pol)) > 0, ] # including this step first greatly reduces computing time by subsetting to one which are overlapping with predicted occupied area.
        gadm.2.ref.occ <- as_tibble(st_intersection(ref.pol, gadm.2.ref.occ)) 
        
        gadm.2.ref.occ$occ.area <- units::drop_units(st_area(gadm.2.ref.occ$geometry))
        gadm.2.ref.occ <- gadm.2.ref.occ %>% 
          dplyr::select(id, occ.area, area)
        
        gadm.2.ref.occ <- gadm.2.ref.occ[gadm.2.ref.occ$pred.occ <- 100/gadm.2.ref.occ$area * gadm.2.ref.occ$occ.area > 10,]
        if(length(gadm.2.ref.occ) > 0){
          gadm.2.ref[c(unique(gadm.2.ref.occ$id)), ]$pred.occ <- 1
        } else {}
        
        gadm.2.ref.acc <- gadm.2.ref %>% 
          dplyr::filter(accu.calc == 1) %>%
          dplyr::filter(COUNTRY %in% non.native) # accuracy only calculated in non.native regions
        TP <- gadm.2.ref.acc %>% dplyr::filter(occupied == 1 & pred.occ == 1) %>% nrow()
        TN <- gadm.2.ref.acc %>% dplyr::filter(occupied == 0 & pred.occ == 0) %>% nrow()
        FP <- gadm.2.ref.acc %>% dplyr::filter(occupied == 0 & pred.occ == 1) %>% nrow()
        FN <- gadm.2.ref.acc %>% dplyr::filter(occupied == 1 & pred.occ == 0) %>% nrow()
        
        accuracy <- 1 - (TP + TN) / (TP + FP + FN + TN)
        
      } else {# accuracy measurement with whole raster, better use this option:
        result.r.agg <- aggregate(result.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% # changing the resolution of the reference and output rasters can change the accuracy result. the ref.raster is aggregated and masked only once at beginning of the function loop
          terra::mask(vect(acc.vect))
        
        ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg) # the mean function in aggregate makes cell values which are between 0 and 1 so it is made to 0 and 1 here again.
        ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg)

        result.r.agg <- ifel(result.r.agg >= 0.5, 1, result.r.agg)
        result.r.agg <- ifel(result.r.agg < 0.5, 0, result.r.agg)
        
        ref.neg <- cells(ref.dist.r.agg, c(0))[[1]] # cell.IDs of reference absences
        ref.pos <- cells(ref.dist.r.agg, c(1))[[1]] # cell.IDs of reference presences
        pred.neg <- cells(result.r.agg, c(0))[[1]] # cell.IDs of predicted absences
        pred.pos <- cells(result.r.agg, c(1))[[1]] # cell.IDs of predicted presences
        
        TP <- sum(pred.pos %in% ref.pos) # uses cell.IDs to check which are correct/false
        TN <- sum(pred.neg %in% ref.neg) # uses cell.IDs to check which are correct/false
        FP <- sum(pred.pos %in% ref.neg) # uses cell.IDs to check which are correct/false
        FN <- sum(pred.neg %in% ref.pos) # uses cell.IDs to check which are correct/false
        
        accuracy <- 1 - (TP + TN) / (TP + FP + FN + TN) # is inverted because the optim algorithm needed it this way.
      }
    #
    #
    #
    
    
    gc()
    print(paste("time.step", t.s, "out of", time.steps, "done"))
    print(Sys.time() - a.int)
    if(t.s == 1){
      plot.stack <- result.r
      accuracy.list <- tibble(accuracy = accuracy, # all this is in the output of the function.
                              spread.val = spread.val,
                              thresh.disp.factor = thresh.disp.factor,
                              time.step = t.s,
                              initiation = initiation,
                              agg.acc.fact = agg.acc.fact, 
                              min.tr = min.tr, 
                              max.dist = max.dist)
    } else {
      plot.stack <- c(plot.stack, result.r)
      accuracy.list.update <- tibble(accuracy = accuracy, 
                              spread.val = spread.val,
                              thresh.disp.factor = thresh.disp.factor,
                              time.step = t.s,
                              initiation = initiation,
                              agg.acc.fact = agg.acc.fact, 
                              min.tr = min.tr, 
                              max.dist = max.dist)
      accuracy.list <- bind_rows(accuracy.list, 
                                 accuracy.list.update)
      }
  }
  if(plot.result == TRUE){
    return(
      list(
        result.r = result.r,
        updated.ua = updated.ua.i,
        accuracies = accuracy.list,
        plot.stack = plot.stack
      )
    )
  } else {
    return(accuracy) # maybe use only this one value for optimization, raster and nodes not necessary for optimization i guess... can be switched when accuracy is optimized to generate output raster
  }
}

# Tapinoma magnum with Seiferts data ====
## data preparation ====
# area of interest (the smaller the faster):
e <- ext(c(-11.8167128633999, 25.9369450325551, 35.5922736627615, 54.2843855437946))

### native and non.native areas ====
# Info from Manuelas Database:
native <- c("Corse", "Italy", "Sardegna", "Sicily", "Spain") # Sardegna = Sardinia, Corse = Corsica; "Morocco", "Tunisia","Algeria" left out because here only europe 
non.native <- c("Belgium", "France", "Germany", "Netherlands", "Slovenia", "Switzerland") # this line are known non.native countries
pot.countries <- c("Portugal", "United Kingdom", "Ireland", "Norway", "Sweden", "Finland", # these are potential countries
                   "Estonia", "Latvia", "Lithuania", "Poland", "Czechia", "Austria", "Croatia",
                   "Greece", "Slovenia", "San Marino", "Bosnia and Herzegovina", "Serbia", "Montenegro",
                   "Albania", "Bulgaria", "Romania", "Ukraine", "Moldova", "Belarus",
                   "Luxembourg", "Denmark", "Turkey", "Hungary", "Iceland", "Sicily","Slovakia", "Kosovo", "North Macedonia", "Montenegro") 
# "Azerbaijan" left out because no data available, majority in Europe

# filter with GADM instead of years => most records are from after 2010, even in the native region so I try here with regional separation, not temporal
gadm.0 <- st_read("simulation input data/world_gadm_410-levels.gpkg", 
                  layer = "ADM_0")
gadm.0.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% native)
gadm.0.n.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% non.native)
gadm.0.pot <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1 <- st_read("simulation input data/world_gadm_410-levels.gpkg", 
                  layer = "ADM_1")
gadm.1.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% native)
gadm.1.n.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% non.native)
gadm.1.pot <- gadm.1 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1.add <- gadm.1 %>% dplyr::filter(COUNTRY %in% c("North Macedonia" , "Moldova", "Montenegro", "Cyprus") | NAME_1 == "Kaliningrad") # are not in gadm.2 level.....

gadm.2 <- st_read("simulation input data/world_gadm_410-levels.gpkg", 
                  layer = "ADM_2") %>%
  dplyr::filter(COUNTRY %in% c(native, non.native, pot.countries))
gadm.2 <- dplyr::bind_rows(gadm.2, gadm.1.add)
gadm.2.native <- gadm.2 %>% 
  dplyr::filter(COUNTRY %in% native)
gadm.2.n.native <- gadm.2 %>% 
  dplyr::filter(COUNTRY %in% non.native)
gadm.2.ref <- gadm.2
gadm.2.pot <- gadm.2 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

### occurrence data from Seifert ====
# transform excel data into a shapefile.
tapintro <- readxl::read_excel("simulation input data/tapinoma/TAPINTRO.xlsx")
tapinoc <- readxl::read_excel("simulation input data/tapinoma/TAPINO_C.xlsx") %>% 
  dplyr::filter(HYP == "magn")

tap <- bind_rows(tapinoc, tapintro)

comp <- tibble(tap.mag = "tapinoma magnum",
               x = tap$LON,
               y = tap$LAT,
               year = tap$YEAR) %>% 
  na.omit()
comp$x <- as.numeric(comp$x)
comp$y <- as.numeric(comp$y)
comp <- st_as_sf(comp, coords = c("x", "y"), crs = 4326)

tap.mag.ini <- comp %>% # filter based on record day which occurrence records are used as initial distribution (i.e. starting points) in the simulation
  dplyr::filter(year <= 2015) %>%
  st_intersection(st_union(gadm.0.native, gadm.1.native)) 

st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
tap.mag.n.native <- comp %>% # create a subset of occurrence non.native records.
  st_as_sf() %>%
  st_intersection(gadm.0.n.native) %>%
  st_erase(gadm.0.native) %>%
  st_erase(gadm.1.native)

#
#
#
### read SDM map & make resistance raster for gridprocess::rawpsread() ====
# biomod2, gbm = boosted regression trees, variables are in the file name (merraclim database):
# BIO1: Annual Mean Temperature
# BIO10: Mean Temperature of Warmest Quarter
# BIO11: Mean Temperature of Coldest Quarter
# BIO13: Precipitation of Wettest Month
# BIO14: Precipitation of Driest Month
# LC = copernicus land-cover classes
# pop.density = population density raster (copernicus GHS)
gbm.r <- rast("simulation input data/tapinoma/biomod2_GBM.bio01.bio10.bio11.bio13.bio14.LC.pop.density.tif")
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1

gbm.r <- aggregate(gbm.r, 5,  # making the raster more coarse to avoid overly high resolution and reduce computing time; INCLUDE IN ASGRID() BELOW!!
                   fun = mean, 
                   na.rm = TRUE)

mask.thresh <- 0.35 
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>% # gets a spatraster that has only cells where the suitability is equal or below the mask.thresh, this will be used as unsuitability mask in the dispersal function
  cells() # gets the cell numbers of these cells; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r^2 # exponential conversion from habitat suitability to resistance instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 (again) irrespective of the transformation

gbm.r <- terra::mask(gbm.r,vect(gadm.0))
gbm.r.inv <- gbm.r*-1 + max(values(gbm.r[[1]]), na.rm = TRUE) # transform (i.e. invert) suitability raster to resistance raster.

spread.limit <- 9999 # set value for areas which cannot be crossed
gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below or it will crash, so replace with spread.limit to make these areas not crossable
                   from = NA, 
                   to = spread.limit)

gbm.r.inv <- asgrid(gbm.r.inv, # convert raster to grid format for spread function
                    xll = xmin(gbm.r.inf),
                    yll = ymin(gbm.r.inv),
                    cellsize = 5000) # update cellsize with the aggregate factor * 1000m
# creation of resistance layer for gridprocess::rawspread() done.
#
#
#

#
#
#
### plot to check if all data align and see what they look like. ====
plot(gbm.r)
plot(tap.mag.ini[,1], add = TRUE, col = "red", pch = 19)
plot(tap.mag.n.native[,1], add = TRUE , col = "yellow", pch = 19)
plot(st_union(st_buffer(tap.mag.n.native[,1], 120000)), add = TRUE , col = "yellow", pch = 19) # i chose this area as area inside which the accuracy is calculated to see whether sub-setting the geographic extend to an area which is within reach of the species makes the accuracy calculation more reasonable and trustworthy.
#
#
#

#
#
#
### read in traffic network ====
eu.links <- st_read("simulation input data/1031.1165ua.GHS.800km.gpkg", # this is a shapefile with all least-cost paths (i.e. open street map routes) between all pairs of urban areas
                    layer = "1031.1165ua.GHS.800km.exp.predict")
eu.links$length <- eu.links$original.dist
eu.links$link.id <- 1:nrow(eu.links) # paths need an ID for easier use later
eu.links <- eu.links %>% 
  rename(o.ID = start.ID)
eu.links <- eu.links %>% 
  rename(d.ID = dest.ID)

# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters:
eu.links <- eu.links[lengths(st_intersects(eu.links, 
                                           st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                        st_crs(eu.links)))) > 0, ]
eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links) # geometries not needed for use in the dispersal() function.
#
#
#

#
#
#
### read in urban areas ====
nodes <- st_read("simulation input data/eur.urban.areas.gpkg", 
                 layer = "GHS.29.10.2025.points")
nodes <- nodes[lengths(st_intersects(nodes, st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), st_crs(nodes)))) > 0, ]# crop the nodes to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters
#
#
#

#
#
#
### define input ====
empty.r <- gbm.r
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything, is used in dispersal function to make the result raster

native.dist.r <- empty.r %>% 
  terra::mask(vect(st_buffer(tap.mag.ini, 10000)),
              updatevalue = 1,
              inverse = TRUE)

ref.dist.r <- native.dist.r %>% # used as reference to calculate the accuracy of the dispersal()-output.
  terra::mask(vect(st_buffer(tap.mag.n.native, 10000)),
              updatevalue = 1,
              inverse = TRUE) 

ini.nodes <- nodes %>% 
  dplyr::filter(gc_ucn_mai_2025 %in% c(nodes[lengths(st_intersects(nodes, tap.mag.ini)) > 0,]$gc_ucn_mai_2025)) # nodes which overlap with the known initial distribution are used as ini.nodes (i.e. initial nodes)

add.nodes <- tibble(ID = nodes$ID, is.occupied = terra::extract(native.dist.r, vect(nodes), ID = FALSE)[[1]]) # nodes which lie within the known distributional range are added (number of added nodes depends on the buffer used above to make this raster)
add.nodes <- add.nodes %>% 
  dplyr::filter(is.occupied > 0)
ini.nodes <- unique(c(ini.nodes$ID, add.nodes$ID))

agg.acc.fact <-  10 # how much should the result and reference raster be aggregated to calculate the accuracy?
acc.vect <- st_union(st_buffer(tap.mag.n.native[,1], 120000)) #acc.vect = accuracy vector, inside this area, accuracy is calculated.
ext(ref.dist.r) == ext(empty.r) # check if data match.
ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% # aggregate the reference raster to avoid being overly fine-scaled
  terra::mask(vect(acc.vect))
ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg) # make binary (fun = "mean" in aggregate makes non-binary values)
ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg) # make binary (fun = "mean" in aggregate makes non-binary values)

spread.val <- 1 # budget for spread used in gridprocess::rawspread()
thresh.disp.factor <- 0.9 # threshold which has to be reached in a cell to be treated as occupied/presence
time.steps <- 2 # how many consecutive iterations?
min.tr <- quantile(eu.links$predicted, probs = c(0.95), na.rm = TRUE)[[1]] # probs defines which quantile of the traffic volumes is used as minimum traffic value that filter or paths which are used in the traffic network.
max.dist <- 350000 # paths in the network longer than this will not be used, in meter.
#
#
#

#
#
#
## run function ====
out <- dispersal(
  land.spread = TRUE,
  net.spread = TRUE,
  spread.val = spread.val,
  nodes.cut.off = nodes.cut.off,
  thresh.disp.factor = thresh.disp.factor, 
  time.steps = time.steps,
  dist.ini = tap.mag.ini,
  ini.nodes = ini.nodes,
  ref.raster = empty.r,
  result.r = empty.r,
  ref.dist.r = ref.dist.r,
  dist.red = FALSE,
  dist.red.param = 2,
  plot.result = TRUE, 
  sample.nodes.from.raster = TRUE,
  regional.approach = FALSE,
  regional.approach.acc = FALSE,
  unsuitability.mask = mask,
  acc.vect = acc.vect,
  min.tr = min.tr,
  max.dist = max.dist
)

#
#
#

#
#
#
## plot output ====
out[[1]] <- mask(out[[1]],vect(gadm.0))

par(mfrow = c(2,2))
plot(native.dist.r, 
     main = "initial distribution", 
     background = "darkgrey")
plot(out[[1]], 
     main = paste0(
       "s.v:", spread.val, "; ",
       "init: 1;",
       "\nmin.tr:", round(min.tr,2), "; ",
       "max.dist:", max.dist/1000, "km; ",
       #"\nn.c.o =", nodes.cut.off, ";",
       "t.s:", time.steps, ";",
       "\nt.d.f:", thresh.disp.factor, "; ",
       "transformation ²"), 
     background = "darkgrey")
plot(out[[1]], 
     main = paste(
       "s.v =", spread.val, ";",
       "init = 1;",
       #"\nn.c.o =", nodes.cut.off, ";",
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
out[[4]] <- mask(out[[4]], vect(gadm.0))
plot(out[[4]][[1:2]])
#
#
#


#
#
#
#
#
# Senecio inaequidens ====
# define inputs needed for function:
## data preparation ====
### geographic reference ====
ref.raster <- rast("simulation input data/senecio/dgm1000_utm32s.asc") # 1km² Digitales Geländemodell from https://gdz.bkg.bund.de/index.php/default/digitales-gelandemodell-gitterweite-1000-m-dgm1000.html ; 16.05.2025
ref.raster <- rast("simulation input data/senecio/biomod2_GBM.first.try.tif") %>%
  terra::aggregate(5)
names(ref.raster) <- "layer"
#terra::crs(ref.raster) <- "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # assign crs
#ref.raster <- terra::aggregate(ref.raster, 5) # do with all rasters that are used!

ger <- st_read("simulation input data/senecio/gadm41_DEU.gpkg", # to crop to germany
               layer = "ADM_ADM_0") %>% 
  st_transform(st_crs(ref.raster))


### occurrence data (CASPIAN & gbif) ====
sen.spread <- data.table::fread("simulation input data/senecio/SenecioSpread_HegerBoehmerCaspianSuppl.csv") # data from Hanno for the Caspian suppl map
sen.spread <- sen.spread[,-c(2,3)]
sen.gbif <- read.csv("simulation input data/senecio/sen.inaeq.gbif.csv", 
                     sep = "\t")
sen.gbif <- tibble(Long = sen.gbif$decimalLongitude,
                   Lat = sen.gbif$decimalLatitude,
                   year = sen.gbif$year) %>% 
  na.omit()
sen.spread <- bind_rows(sen.spread,
                   sen.gbif)

sen.spread <- st_as_sf(sen.spread, coords = c("Long","Lat"), crs = 4326) %>%
  st_transform(st_crs(ref.raster))

sen.spread$decade <- sen.spread$year
# convert to decades instead of years:
sen.spread$decade[sen.spread$decade < 1979] <- 1979 
sen.spread$decade[between(sen.spread$decade,1980,1989)] <- 1989 
sen.spread$decade[between(sen.spread$decade,1990,1999)] <- 1999 
sen.spread$decade[between(sen.spread$decade,2000,2009)] <- 2009
sen.spread$decade[between(sen.spread$decade,2010,2019)] <- 2019
sen.spread$decade[between(sen.spread$decade,2020,2025)] <- 2025

sen.ini <- sen.spread %>% # filter based on record day which occurrence records are used as initial distribution (i.e. starting points) in the simulation
  dplyr::filter(year <= 1979)

# final reference distribution:
# here try with regions invaded by 2004/2009 (i.e. time.step = 4)
ref.p <- sen.spread %>% 
  dplyr::filter(year <= 2009)
#
#
#
#
#

### read SDM map & make resistance raster for gridprocess::rawpsread() ====
gbm.r <- rast("simulation input data/senecio/biomod2_GBM.first.try.tif") %>%  # biomod SDM output from boosted regression trees (gbm)
  project(ref.raster)
names(gbm.r) <- "layer"
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r <- terra::mask(gbm.r,vect(ger))
gbm.r <- subst(gbm.r, NA, 0)

mask.thresh <- 0.4
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(ger)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY # deleted because if I do this, then the grid allows to cross borders over time
  cells()                               # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r^5 # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
gbm.r[mask] <- 0

gbm.r.inv <- gbm.r*-1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

limit <- max(values(gbm.r.inv), na.rm = TRUE)
#spread.limit <- limit*100 # set value for areas which are not crossable
gbm.r.inv <- subst(gbm.r.inv, # must not have NAs for the function below, so replace with spread.limit to make these areas not crossable
                   1, 
                   9999) # the final 
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network becaus of color-scale. reduce spread.limit to make it visible again.
gbm.r.inv <- asgrid(gbm.r.inv, # convert to grid for spread function
                    xll = xmin(gbm.r.inv),
                    yll = ymin(gbm.r.inv),
                    cellsize = 1000) # update cellsize with the aggregate factor * 1000m
#
#
#
#
#


#
#
#
### plot to check if all data align and see what they look like. ====
plot(gbm.r)
plot(ref.p[,1], add = TRUE , col = "yellow", pch = 19)
plot(sen.ini[,1], add = TRUE, col = "red", pch = 19)
#
#
#


### read in traffic network ====
eu.links <- st_read("simulation input data/1031.1165ua.GHS.800km.gpkg", # this is a shapefile with all least-cost paths (i.e. open street map routes) between all pairs of urban areas
                    layer = "1031.1165ua.GHS.800km.exp.predict")
eu.links$length <- eu.links$original.dist
eu.links$link.id <- 1:nrow(eu.links) # paths need an ID for easier use later
eu.links <- eu.links %>% 
  rename(o.ID = start.ID)
eu.links <- eu.links %>% 
  rename(d.ID = dest.ID)

# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters:
eu.links <- eu.links[lengths(st_intersects(eu.links, 
                                           ger)) > 0, ]
eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links) # geometries not needed for use in the dispersal() function.


#
#
#
### read in urban areas ====
nodes <- st_read("simulation input data/eur.urban.areas.gpkg", 
                 layer = "GHS.29.10.2025.points")
nodes <- nodes[lengths(st_intersects(nodes, ger)) > 0, ]# crop the nodes to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters
#
#
#


#
#
#
### define input ====
empty.r <- ref.raster
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything

ini.dist.r <- empty.r %>% 
  terra::mask(vect(st_buffer(sen.ini, 5000)),
              updatevalue = 1,
              inverse = TRUE)

ref.dist.r <- ini.dist.r %>% # used as reference to calculate the accuracy of the dispersal()-output.
  terra::mask(vect(st_buffer(ref.p, 5000)),
              updatevalue = 1,
              inverse = TRUE) 

nodes$occ <- terra::extract(ini.dist.r, vect(nodes))[,2]
ini.nodes <- nodes %>% dplyr::filter(occ == 1)
ini.nodes <- ini.nodes$ID
#
#
#

spread.val <- 1
nodes.cut.off <- 0.03
thresh.disp.factor <- 0.75
time.steps <- 4
agg.acc.fact <- 1
acc.vect <- ger
min.tr <- 0
max.dist <- 1000000

ext(ref.dist.r) == ext(empty.r)
ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% 
  terra::mask(vect(acc.vect))
ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg)
ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg)


out <- dispersal(
  land.spread = TRUE,
  net.spread = TRUE,
  spread.val = spread.val,
  nodes.cut.off = nodes.cut.off,
  thresh.disp.factor = thresh.disp.factor, 
  time.steps = time.steps,
  dist.ini = sen.ini,
  ini.nodes = ini.nodes,
  ref.raster = ref.raster,
  result.r = empty.r,
  ref.dist.r = ref.dist.r,
  dist.red = FALSE,
  dist.red.param = 2,
  plot.result = TRUE, 
  sample.nodes.from.raster = TRUE,
  regional.approach = FALSE, 
  regional.approach.acc = FALSE,
  unsuitability.mask = mask,
  acc.vect = acc.vect,
  min.tr = min.tr,
  max.dist = max.dist
)

par(mfrow = c(2,2))
plot(ini.dist.r, 
     main = "initial distribution", 
     background = "darkgrey")
plot(mask(out[[1]],vect(ger)), 
     main = paste(
       "s.v =", spread.val, ";",
       "\nn.c.o =", nodes.cut.off, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "transformation ^5"), 
     background = "darkgrey")
plot(mask(out[[1]],vect(ger)), 
     main = paste(
       "s.v =", spread.val, ";",
       "\nn.c.o =", nodes.cut.off, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "ini.points = green ;",
       "\nreached ua (sim.) = red"), 
     background = "darkgrey")
plot(sen.ini, add = TRUE, col = "darkgreen", pch = 19) # starting points
plot(dplyr::filter(nodes, ID %in% out[[2]])[,1], add = TRUE, col = "red", pch = 19) # reached ua
plot(ref.dist.r, 
     main = "final reference distribution", 
     background = "darkgrey")
out[[3]]
out[[4]] <- mask(out[[4]], vect(ger))
plot(out[[4]])

###
# end senecio
###
###
###
###
###




#
#
#
# example for optimization-by-hand approach ====
accuracy.list <- tibble(accuracy = numeric(), 
                        spread.val = numeric(),
                        thresh.disp.factor = numeric(),
                        time.step = integer(),
                        initiation = numeric(), 
                        agg.acc.fact = numeric(), 
                        min.tr = numeric(), 
                        max.dist = numeric())

parameters <- tidyr::crossing(spread.val = c(0.9, 1, 1.5, 2), #4 
         thresh.disp.factor = c(0.3,0.6,0.9), #3
         min.tr.quantile = c(0.1,0.5,0.9), #3
         max.dist = c(100000, 350000, 500000))


for(i.p in 1:nrow(parameters)){
  print(i.p)
  spread.val <- parameters[i.p,]$spread.val
  thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
  time.steps <- 27
  min.tr <- quantile(eu.links$predicted, probs = c(parameters[i.p,]$min.tr.quantile), na.rm = TRUE)[[1]]
  #min.tr <- 0.5 # derived from the reference distribution (see below)
  max.dist <- parameters[i.p,]$max.dist
  
  out <- dispersal(
    land.spread = TRUE,
    net.spread = TRUE,
    spread.val = spread.val,
    nodes.cut.off = nodes.cut.off,
    thresh.disp.factor = thresh.disp.factor, 
    time.steps = time.steps,
    dist.ini = tap.mag.ini,
    ini.nodes = ini.nodes,
    ref.raster = empty.r,
    result.r = empty.r,
    ref.dist.r = ref.dist.r,
    dist.red = FALSE,
    dist.red.param = 2,
    plot.result = TRUE, 
    sample.nodes.from.raster = TRUE,
    regional.approach = FALSE,
    regional.approach.acc = FALSE,
    unsuitability.mask = mask,
    acc.vect = acc.vect,
    min.tr = min.tr,
    max.dist = max.dist
  )
  
  if(i.p == 1){
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    optim.output <- out$accuracies
  } else {
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    optim.output <- bind_rows(optim.output, out$accuracies)
  }
  
}

save <- optim.output
# a hundred years later....
#
#
#


#
#
#
# better plots with tmap ====
library(tmap)
out.pol <- as.polygons(out[[4]][[23]]) %>%
  st_as_sf()
plot(st_as_sf(out.pol))

reference.pol <- as.polygons(ref.dist.r) %>%
  st_as_sf()

out.pol$layer <- as.character(out.pol$layer)
reference.pol$layer <- as.character(reference.pol$layer)

tm_shape(reference.pol) +
  tm_polygons(fill = "layer", 
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
