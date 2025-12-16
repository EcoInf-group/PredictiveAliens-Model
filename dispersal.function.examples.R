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
                      ini.nodes, # initial urban areas from which the species can spread through the traffic network and the landscape
                      #thresh.nodes, # threshold for which cells from spread are treated as origin in following step?
                      ref.raster, # reference raster (resolution etc) for the output
                      result.r, # output raster will be built from this empty raster
                      initiation = 1, # does not change relations (therefore no difference between plots when using different values), see nodes.cut.off
                      nodes.cut.off, # states how much traffic urban area "y" has to receive from urban area "x" to be treated as invaded. 
                      # each urban area gets an outgoing traffic budget of 1 and this budget is then distributed among all links going from this urban 
                      # area according to the traffic volume on this link. if the total traffic volume on all links going from urban area x is 1000 and the 
                      # link to urban area y has a volume of 10, then urban area y receives 1/1000*10 = 0.01 and it is checked whether or not this is above 
                      # the threshold for being treated as invaded. 
                      time.steps, # how many steps should the simulation run?
                      dist.red = FALSE, # in each step an output raster is generated with areas which are occupied. from these areas the border cells are used
                      # as starting points for further landscape spread in the following step. if dist.red = TRUE then only every xth cell is used as starting point
                      # with x = dist.red.param. greatly reduces computing time BUT MIGHT CHANGE OUTPUT PATTERNS!
                      dist.red.param, # define how strongly the initial cells should be reduced, numeric; 
                      # e.g. if 2, then every second cell/coordinate pair is chosen from the ini.update matrix
                      ref.dist.r, # reference raster with the final known distribution with cells being present (1) or absent (0), used for accuracy calculation
                      # (accuracy <- (True Positives + True Negatives) / (True Positives + False Positives + False Negatives + True Negatives))
                      # the result.r output raster will be compared with this and the values needed for accuracy calculation derived from the comparison.
                      plot.result = TRUE, 
                      node.accumulation = TRUE, # apparently deprecated. was formerly used to allow the urban areas to accumulate receiving traffic over time steps
                      # but now each urban area which is below the threshold at end of timestep is set to 0 (= not occupied)
                      sample.nodes.from.raster = TRUE, # if TRUE then urban areas which are within the occupied areas on the result raster will be treated as 
                      # occupied and can serve as starting point in the traffic network. this allows a species to enter the traffic network from the landscape.
                      regional.approach = FALSE, # not recommended. here calculations are based on gadm-0 levels, but gadm will have to be provided and it takes ages...
                      regional.approach.acc = FALSE, # calculates accuracy not on raster cells but based on gadm regions. however, i recommend to provided an aggregated
                      # raster to use for calculation since it is much faster and ecologically more meaningful than usein administrative boundaries
                      unsuitability.mask = NULL, # if provided (cell IDs) then every cell which is provided will be set to unoccupied at end of timestep. can be used
                      # to make sure that unsuitable cells will not become occupied. they can still be crossed though.
                      acc.vect,
                      agg.acc.fact = 1,
                      min.tr,
                      max.dist
) {
  a.int <- Sys.time()
  eu.links <- eu.links %>% # filter for all paths which are connected with the chosen start.node
    dplyr::filter(predicted >= min.tr) %>% # filter all paths that have at least the given traffic volume
    dplyr::filter(length <= max.dist) # filter all paths which are up to the max.dist length
  
  ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% # aggregate reference raster for accuracy measurement in the given area of interest (area of interest = acc.vect)
    terra::mask(vect(acc.vect))
  
  for(t.s in 1:time.steps){ # iterate over time-steps
    
    if(land.spread == TRUE){
      # start spread through landscape:
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
      # end spread through landscape
      
      if(t.s == 1){ # all reached raster-cells are marked as such in the final result.r. this is updated at the end of each landscape spread iteration.
        dist.r <- mask(empty.r, spread.m.r, updatevalue = 1, inverse = TRUE)
        result.r <- mask(result.r, spread.m.r, updatevalue = 1, inverse = TRUE)
      } else {
        dist.r <- mask(dist.r, spread.m.r, updatevalue = 1, inverse = TRUE)
        result.r <- mask(result.r, spread.m.r, updatevalue = 1, inverse = TRUE)
      }
      
      # reduce number of starting-points for consecutive time-step and update initial coordinates for raster spread:
      # IMPORTANT: I HAD CASES WHERE THE REDUCTION DID NOT CHANGE THE RESULT AND OTHERS WHERE THE OUTPUT WAS VERY DIFFERENT.
      dist.r <- result.r
      dist.r[is.na(dist.r)] <- 0 # sets the NA cells to 0 so that terra::boundaries can find the edges of the occupied patches below, otherwise it also identifies german border as edge (identifes all class differences between c(NA, 0, 1))
      b <- boundaries(dist.r, classes = TRUE, inner = FALSE)
      b.c <- cells(b, 1)[[1]] # identify boundary cells
      id.ini.update <- rowColFromCell(ref.raster, b.c) # in matrix: x,y ; lon, lat without resetting
      
      if(dist.red == TRUE){id.ini.update <- id.ini.update[c(seq(from = 0, 
                                                                to = nrow(id.ini.update), 
                                                                by = dist.red.param)), ]} else {} # takes only each nth cell (e.g. dist.red.param = 2 then each second cell)
      # end number reduction
    } else {}
    
    if(net.spread == TRUE) { # net.spread = spread through traffic network (i.e. IDs of urban areas and the traffic flows between each pair)
      # start network spread:
      if(t.s == 1){} else {
        ini.nodes <- updated.ua.i} # updated.ua.i are the ones whcih were reached in a previous time.step.
      
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
        # useless:
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
      # sample nodes from the distributional landscape raster which lie in invaded territory - optional to do this. allows the species to switch from the landscape to the traffic network. the opposite way is always active but can be limited with spread.value and thresh.disp.fact.
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
    

    
    if(!is.null(unsuitability.mask)){ # if there is an unsuitability mask provided all the cells which have a habitat suitability below a certain threshold are set to unoccupied at the end of each iteration. this way, a species might cross an unsuitable habitat but it can not establich there an each population will be deleted at the end of each time.step.
      result.r[mask] <- 0 # the unsuitability mask has IDs of each cells which were determined as being unsuitable (making an unsuitability mask is a preparatory step and not part of this function).
    } else {}
    
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
        
      } else {# accuracy measurement with whole raster :
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
                              #nodes.cut.off = nodes.cut.off, 
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
                              #nodes.cut.off = nodes.cut.off, 
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



# Tapinoma magnum ====
# define inputs needed for function:
#ref.raster <- rast("C:/Users/JLU-SU/Documents/PredictiveAliens/dgm1000_utm32s.asc") # 1km² Digitales Geländemodell from https://gdz.bkg.bund.de/index.php/default/digitales-gelandemodell-gitterweite-1000-m-dgm1000.html ; 16.05.2025
#terra::crs(ref.raster) <- "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # assign crs
#ref.raster <- terra::aggregate(ref.raster, 5) # do with all rasters that are used!
# Info from Manuelas Database:
native <- c("Corse", "Italy", "Sardegna", "Sicily", "Spain") # Sardegna = Sardinia, Corse = Corsica; "Morocco", "Tunisia","Algeria" left out because here only europe 
non.native <- c("Belgium", "France", "Germany", "Netherlands", "Slovenia", "Switzerland") # this line are known non.native countries
pot.countries <- c("Portugal", "United Kingdom", "Ireland", "Norway", "Sweden", "Finland", # these are potential countries
                   "Estonia", "Latvia", "Lithuania", "Poland", "Czechia", "Austria", "Croatia",
                   "Greece", "Slovenia", "San Marino", "Bosnia and Herzegovina", "Serbia", "Montenegro",
                   "Albania", "Bulgaria", "Romania", "Ukraine", "Moldova", "Belarus",
                   "Luxembourg", "Denmark", "Turkey", "Hungary", "Iceland", "Sicily","Slovakia", "Kosovo", "North Macedonia", "Montenegro") 
# "Azerbaijan" left out because no data available, majority in Europe

# filter with GADM? instead of years => most records are from after 2010, even in the native region so I try here with regional separation, not temporal
gadm.0 <- st_read("GADM/world_gadm_410-levels.gpkg", 
                  layer = "ADM_0") #%>% 
#st_transform(3035)
gadm.0.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% native)
gadm.0.n.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% non.native)
gadm.0.pot <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1 <- st_read("GADM/world_gadm_410-levels.gpkg", 
                  layer = "ADM_1") #%>% 
#st_transform(3035)
gadm.1.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% native)
gadm.1.n.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% non.native)
gadm.1.pot <- gadm.1 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1.add <- gadm.1 %>% dplyr::filter(COUNTRY %in% c("North Macedonia" , "Moldova", "Montenegro", "Cyprus") | NAME_1 == "Kaliningrad") # are not in gadm.2 level.....

gadm.2 <- st_read("GADM/world_gadm_410-levels.gpkg", 
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

e <- ext(c(-11.8167128633999, 25.9369450325551, 35.5922736627615, 54.2843855437946))

comp <- read.csv("Tapinoma magnum/gbif t.darioi t.ibericum.csv", sep = "\t")
comp <- tibble(spec = comp$species,
               x = comp$decimalLongitude,
               y = comp$decimalLatitude,
               year = comp$year) %>% 
  na.omit()
comp <- st_as_sf(comp, coords = c("x", "y"), crs = 4326)

#gbm.r <- rast("C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma magnum/biomod2_GBM.first.try.with.traffic.tif")# biomod SDM output from boosted regression trees (gbm)
gbm.r <- rast("Tapinoma magnum/biomod2_GBM.bio01.bio10.bio11.bio12.LC.tif")
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
#gbm.r <- terra::mask(gbm.r,vect(gadm.0))

gbm.r <- aggregate(gbm.r, 5, 
                   fun = mean, 
                   na.rm = TRUE) # making the raster more coarse; INCLUDE IN ASGRID() BELOW!!

mask.thresh <- 0.35 
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(gadm.0)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY
  cells()                               # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r^2 # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation

gbm.r <- terra::mask(gbm.r,vect(gadm.0))

gbm.r.inv <- gbm.r*-1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

#limit <- max(values(gbm.r.inv), na.rm = TRUE)
spread.limit <- 9999 #limit*100 # set value for areas which are not crossable
gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below, so replace with spread.limit to make these areas not crossable
                   from = NA, 
                   to = spread.limit) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network because of color-scale. reduce spread.limit to make it visible again.

gbm.r.inv <- asgrid(gbm.r.inv, # convert to grid for spread function
                    xll = xmin(gbm.r.inf),
                    yll = ymin(gbm.r.inv),
                    cellsize = 5000) # update cellsize with the aggregate factor * 1000m


tap.mag.spread <- readxl::read_xlsx("Tapinoma magnum/records tapinoma magnum.xlsx") # data from Benoit, AntMaps
tap.mag.spread$dec_lat <- as.numeric(tap.mag.spread$dec_lat)
tap.mag.spread$dec_lon <- as.numeric(tap.mag.spread$dec_lon)

tap.mag.spread$pres <- 1
tap.mag.spread <- tap.mag.spread %>% 
  dplyr::filter(!is.na(dec_lat) & !is.na(dec_lon))
tap.mag.spread <- st_as_sf(tap.mag.spread, coords = c("dec_lon","dec_lat"), crs = 4326) 

tap.mag.spread$suit <- terra::extract(gbm.r$layer, vect(tap.mag.spread))[,2]
tap.mag.spread <- tap.mag.spread %>% 
  dplyr::filter(suit > 0)


load("Traffic/ger.glm.rda")
eu.links <- st_read("Traffic/gravity.output/europe/shortest.paths.update.09.07.2025/661661ua.800km.gpkg", # includes all shortest paths with origin, destination and predicted traffic based on german GLM
                    "661661ua.800km")
eu.links$link.id <- 1:nrow(eu.links)
eu.links <- eu.links[lengths(st_intersects(eu.links, 
                                           st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                        st_crs(eu.links)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters


eu.links <- rename(eu.links, 
                   c(o.ID = start.ID,
                     d.ID = dest.ID))

eu.links$dist <- c(base::scale(log(eu.links$length)))
eu.links$o.gdp.per.cap <- c(base::scale(log(eu.links$s.gdp.2015 / eu.links$s.pop)))
eu.links$d.gdp.per.cap <- c(base::scale(log(eu.links$d.gdp.2015 / eu.links$d.pop)))

eu.links$o.poly.and.pts.n.poi.commercial.union.erase <- c(base::scale(log(eu.links$s.poly.and.pts.n.poi.commercial.union.erase+1)))
eu.links$d.poly.and.pts.n.poi.commercial.union.erase <- c(base::scale(log(eu.links$d.poly.and.pts.n.poi.commercial.union.erase+1)))

eu.links$o.poly.and.pts.n.poi.industrial.union.erase <- c(base::scale(log(eu.links$s.poly.and.pts.n.poi.industrial.union.erase+1)))
eu.links$d.poly.and.pts.n.poi.industrial.union.erase <- c(base::scale(log(eu.links$d.poly.and.pts.n.poi.industrial.union.erase+1)))

eu.links$predicted <- stats::predict(object = ger.glm, 
                                     newdata = eu.links,
                                     type = "response")
eu.links$predicted <- exp(eu.links$predicted)
eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links)


nodes <- st_read("Traffic/eur.urban.areas.gpkg", # urban areas with gdp, pop, n.poi
                 "updated.eu.centroids.09.07.2025.w.POIs")
nodes <- nodes[lengths(st_intersects(nodes, 
                                     st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                  st_crs(nodes)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters

sort(unique(tap.mag.spread$Country))
target <- c("Austria", "Belgium", "Croatia", "France", "Germany", "Greece", "Hungary", "Italy", "Netherlands", "Slovenia", "Spain", "Switzerland", "United Kingdom") # filtered for only european countries
tap.mag.spread <- tap.mag.spread %>% 
  dplyr::filter(Country %in% target)

tap.mag.gbif <- read.csv("Tapinoma magnum/gbif.01.07.2025.csv", sep = "\t") %>% 
  dplyr::filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  st_as_sf(coords = c("decimalLongitude","decimalLatitude"), crs = 4326)

tap.mag.spread %>% 
  count(`Year Collection`)
tap.mag.gbif %>% 
  count(year)

#tap.mag.spread.ini <- tap.mag.spread %>% 
#  dplyr::filter(`Year Collection` <= 2010)
#tap.mag.gbif.ini <- tap.mag.gbif %>% 
#  dplyr::filter(year <= 2010)

#st_write(tap.mag.gbif, 
#         dsn = "C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma magnum/dist.data.gpkg", 
#         layer = "gbif")
#st_write(tap.mag.spread, 
#         dsn = "C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma magnum/dist.data.gpkg", 
#         layer = "antmap")


tap.mag.ini <- tibble(species = "tapinoma.magnum", geom = c(tap.mag.spread$geometry, tap.mag.gbif$geometry), pres = 1, year = c(tap.mag.spread$`Year Collection`, tap.mag.gbif$year)) %>%
  st_as_sf() %>%
  st_intersection(st_union(gadm.0.native, gadm.1.native)) 

st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
tap.mag.n.native <- tibble(species = "tapinoma.magnum", geom = c(tap.mag.spread$geometry, tap.mag.gbif$geometry), pres = 1, year = c(tap.mag.spread$`Year Collection`, tap.mag.gbif$year)) %>%
  st_as_sf() %>%
  st_intersection(gadm.0.n.native) %>%
  st_erase(gadm.0.native) %>%
  st_erase(gadm.1.native)

# plot(gbm.r)
# plot(tap.mag.ini[,1], add = TRUE, col = "red", pch = 19)
# plot(tap.mag.n.native[,1], add = TRUE , col = "yellow", pch = 19)


empty.r <- gbm.r
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything

native.dist.r <- empty.r %>%
  terra::mask(vect(st_buffer(tap.mag.ini, 10000)),
              updatevalue = 1,
              inverse = TRUE)
# plot(native.dist.r)

ref.dist.r <- native.dist.r %>%
  terra::mask(vect(st_buffer(tap.mag.n.native, 10000)),
              updatevalue = 1,
              inverse = TRUE) 

# plot(ref.dist.r)

#gadm.2.ref$occupied <- "" # for reference and AUC calculation only the known non-native regions are used.
#gadm.2.ref$pred.occ <- ""
#for(i in 1:nrow(gadm.2.ref)){
#  ref.r <- terra::mask(ref.dist.r, 
#                       vect(gadm.2.ref[i,]), 
#                       inverse = FALSE, 
#                       touches = TRUE) 
#  ref.a <- as.polygons(ref.r) |> 
#    disagg() |> 
#    st_as_sf() |> 
#    st_intersection(gadm.2.ref[i,]) 
#  ref.a$area <- units::drop_units(st_area(ref.a))
#  
#  if(nrow(ref.a) > 0){
#    if(ref.a %>% 
#       dplyr::filter(layer == 1) %>% 
#       nrow() == 0) {gadm.2.ref[i,]$occupied <- 0} else {
#         if(100/sum(ref.a$area) * sum(ref.a[which(ref.a$layer == 1),]$area, na.rm = TRUE) >= 15){
#           # result.r <- terra::mask(result.r, vect(test.a), updatevalue = 1)
#           gadm.2.ref[i,]$occupied <- 1
#         } else {gadm.2.ref[i,]$occupied <- 0
#         }
#       }
#    
#  } else {}
#}
#gadm.2.ref$occupied <- as.numeric(gadm.2.ref$occupied)
#gadm.2.ref <- gadm.2.ref %>% 
#  dplyr::filter(!is.na(occupied)) # Islas Canarias removed
#gadm.2.ref$pred.occ <- as.numeric(gadm.2.ref$pred.occ)
#st_write(gadm.2.ref, 
#         "gadm.2.ref.gpkg")

#gadm.2.ref <- st_read("gadm.2.ref.gpkg")

#gadm.2.ref$accu.calc <- ""
#gadm.2.ref[lengths(st_intersects(gadm.2.ref, st_buffer(st_as_sf(bind_rows(tap.mag.n.native, tap.mag.ini)), 100000))) > 0, ]$accu.calc <- 1
#gadm.2.ref$accu.calc <- as.numeric(gadm.2.ref$accu.calc)  





ua <- st_read("Traffic/eur.urban.areas.gpkg", 
              layer = "urban.areas.update.09.07.2025.w.POIs")
ua <- ua[lengths(st_intersects(ua, st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), st_crs(ua)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters

#
#
#
#
#
# insert here: sample habitat suitability of ua's from gbmr and always remove the ones which are too low after each time.step
#ua$suitability <- ""
#ua$suitability <- as.numeric(ua$suitability)
#for(i.ua in 1 : nrow(ua)){
#  ua[i.ua,]$suitability <- terra::extract(x = gbm.r, y = vect(ua[i.ua,]), fun = mean, method = "simple", ID = FALSE, na.rm = TRUE)[[1]]
#}
#ua <- ua %>% dplyr::filter(suitability >= 0.2)
#
#
#
#
#

ini.nodes <- nodes %>% 
  dplyr::filter(Name %in% c(ua[lengths(st_intersects(ua, tap.mag.ini)) > 0,]$Name))
add.nodes <- tibble(ID = nodes$ID, is.occupied = terra::extract(native.dist.r, vect(nodes), ID = FALSE)[[1]])
add.nodes <- add.nodes %>% 
  dplyr::filter(is.occupied > 0)
ini.nodes <- unique(c(ini.nodes$ID, add.nodes$ID))

#gadm.2.n.native$area <- 
#  units::drop_units(st_area(gadm.2.n.native))
#gadm.2.n.native$occupied <- ""

agg.acc.fact <-  10
acc.vect <- gadm.0
ext(ref.dist.r) == ext(empty.r)
ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% 
  terra::mask(vect(acc.vect))
ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg)
ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg)


spread.val <- 1
thresh.disp.factor <- 0.9
time.steps <- 25
# min.tr <- quantile(eu.links$predicted, probs = c(0.3), na.rm = TRUE)
min.tr <- 50 # derived from the reference distribution (see below)
max.dist <- 350000



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


par(mfrow = c(2,2))
plot(native.dist.r, 
     main = "initial distribution", 
     background = "darkgrey")
plot(mask(out[[1]],vect(gadm.0)), 
     main = paste(
       "s.v =", spread.val, ";",
       "init = 1;",
       "\nn.c.o =", nodes.cut.off, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "transformation ²"), 
     background = "darkgrey")
plot(mask(out[[1]],vect(gadm.0)), 
     main = paste(
       "s.v =", spread.val, ";",
       "init = 1;",
       "\nn.c.o =", nodes.cut.off, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "ini.points = green ;",
       "\nreached ua (sim.) = red"), 
     background = "darkgrey")
plot(tap.mag.ini, add = TRUE, col = "green") # starting points
#plot(out[[2]], add = TRUE, col = "red") # reached ua
plot(dplyr::filter(nodes, ID %in% c(out[[2]])), add = TRUE, col = "red")
plot(ref.dist.r, 
     main = "final reference distribution", 
     background = "darkgrey")
# plot(out[[2]], add = TRUE, col = "red", pch = 19) # reached ua
#summary(out[[2]]$tot.prop.traff)
out[[3]]
#plot(mask(out[[4]], vect(gadm.0)))



## ID for tapinoma magnum invaded UAs ====
ref.pol <- as.polygons(ref.dist.r) %>%
  st_as_sf() %>%
  dplyr::filter(layer == 1) %>%
  st_cast("POLYGON")
target.ua <- st_intersection(ua, ref.pol)

target.links <- eu.links.geom %>% 
  dplyr::filter(d.ID %in% target.ua$ID & o.ID %in% target.ua$ID)
#> summary(target.links$predicted)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#10.49    129.83    355.90   2962.01   1425.52 132265.92 

x11()
plot(ref.dist.r)
plot(target.links["predicted"], add = TRUE)





# Tapinoma nigerrimu complex ====
# define inputs needed for function:
#ref.raster <- rast("C:/Users/JLU-SU/Documents/PredictiveAliens/dgm1000_utm32s.asc") # 1km² Digitales Geländemodell from https://gdz.bkg.bund.de/index.php/default/digitales-gelandemodell-gitterweite-1000-m-dgm1000.html ; 16.05.2025
#terra::crs(ref.raster) <- "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # assign crs
#ref.raster <- terra::aggregate(ref.raster, 5) # do with all rasters that are used!
# Info from Manuelas Database:
native <- c("Corse", "Italy", "Sardegna", "Sicily", "Spain") # Sardegna = Sardinia, Corse = Corsica; "Morocco", "Tunisia","Algeria" left out because here only europe 
non.native <- c("Belgium", "France", "Germany", "Netherlands", "Slovenia", "Switzerland") # this line are known non.native countries
pot.countries <- c("Portugal", "United Kingdom", "Ireland", "Norway", "Sweden", "Finland", # these are potential countries
                   "Estonia", "Latvia", "Lithuania", "Poland", "Czechia", "Austria", "Croatia",
                   "Greece", "Slovenia", "San Marino", "Bosnia and Herzegovina", "Serbia", "Montenegro",
                   "Albania", "Bulgaria", "Romania", "Ukraine", "Moldova", "Belarus",
                   "Luxembourg", "Denmark", "Turkey", "Hungary", "Iceland", "Sicily","Slovakia", "Kosovo", "North Macedonia", "Montenegro") 
# "Azerbaijan" left out because no data available, majority in Europe

# filter with GADM? instead of years => most records are from after 2010, even in the native region so I try here with regional separation, not temporal
gadm.0 <- st_read("GADM/world_gadm_410-levels.gpkg", 
                  layer = "ADM_0") #%>% 
#st_transform(3035)
gadm.0.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% native)
gadm.0.n.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% non.native)
gadm.0.pot <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1 <- st_read("GADM/world_gadm_410-levels.gpkg", 
                  layer = "ADM_1") #%>% 
#st_transform(3035)
gadm.1.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% native)
gadm.1.n.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% non.native)
gadm.1.pot <- gadm.1 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1.add <- gadm.1 %>% dplyr::filter(COUNTRY %in% c("North Macedonia" , "Moldova", "Montenegro", "Cyprus") | NAME_1 == "Kaliningrad") # are not in gadm.2 level.....

gadm.2 <- st_read("GADM/world_gadm_410-levels.gpkg", 
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

e <- ext(c(-11.8167128633999, 25.9369450325551, 35.5922736627615, 54.2843855437946))

comp <- readxl::read_xlsx("C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma magnum/records tapinoma magnum.xlsx")
comp <- tibble(tap.nig = "tapinoma magnum",
                   x = comp$dec_lon,
                   y = comp$dec_lat,
                   year = comp$`Year Collection`) %>% 
  na.omit()
comp$x <- as.numeric(comp$x)
comp$y <- as.numeric(comp$y)
comp[370,2] <- -3.36354
tap.mag <- comp

comp <- read.csv("tapinoma nigerrimum complex/nigerrimum complex.csv", sep = "\t")
comp <- tibble(species = comp$species,
                   pres = 1,
                   x = comp$decimalLongitude,
                   y = comp$decimalLatitude,
                   year = comp$year) %>% 
  na.omit()
comp$x <- as.numeric(comp$x)
comp$y <- as.numeric(comp$y)
comp <- bind_rows(tap.mag, comp)
comp <- st_as_sf(comp, coords = c("x", "y"), crs = 4326)






#gbm.r <- rast("C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma magnum/biomod2_GBM.first.try.with.traffic.tif")# biomod SDM output from boosted regression trees (gbm)
gbm.r <- rast("tapinoma nigerrimum complex/biomod2_GBM.bio01.bio10.bio11.bio13.bio14.LC.pop.density.tif")
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
#gbm.r <- terra::mask(gbm.r,vect(gadm.0))

gbm.r <- aggregate(gbm.r, 5, 
                   fun = mean, 
                   na.rm = TRUE) # making the raster more coarse; INCLUDE IN ASGRID() BELOW!!

mask.thresh <- 0.35 
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(gadm.0)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY
  cells()                               # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r^2 # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation

gbm.r <- terra::mask(gbm.r,vect(gadm.0))

gbm.r.inv <- gbm.r*-1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

#limit <- max(values(gbm.r.inv), na.rm = TRUE)
spread.limit <- 9999 #limit*100 # set value for areas which are not crossable
gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below, so replace with spread.limit to make these areas not crossable
                   from = NA, 
                   to = spread.limit) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network because of color-scale. reduce spread.limit to make it visible again.

gbm.r.inv <- asgrid(gbm.r.inv, # convert to grid for spread function
                    xll = xmin(gbm.r.inf),
                    yll = ymin(gbm.r.inv),
                    cellsize = 5000) # update cellsize with the aggregate factor * 1000m


load("Traffic/ger.glm.rda")
eu.links <- st_read("Traffic/gravity.output/europe/shortest.paths.update.09.07.2025/661661ua.800km.gpkg", # includes all shortest paths with origin, destination and predicted traffic based on german GLM
                    "661661ua.800km")
eu.links$link.id <- 1:nrow(eu.links)
eu.links <- eu.links[lengths(st_intersects(eu.links, 
                                           st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                        st_crs(eu.links)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters


eu.links <- rename(eu.links, 
                   c(o.ID = start.ID,
                     d.ID = dest.ID))

eu.links$dist <- c(base::scale(log(eu.links$length)))
eu.links$o.gdp.per.cap <- c(base::scale(log(eu.links$s.gdp.2015 / eu.links$s.pop)))
eu.links$d.gdp.per.cap <- c(base::scale(log(eu.links$d.gdp.2015 / eu.links$d.pop)))

eu.links$o.poly.and.pts.n.poi.commercial.union.erase <- c(base::scale(log(eu.links$s.poly.and.pts.n.poi.commercial.union.erase+1)))
eu.links$d.poly.and.pts.n.poi.commercial.union.erase <- c(base::scale(log(eu.links$d.poly.and.pts.n.poi.commercial.union.erase+1)))

eu.links$o.poly.and.pts.n.poi.industrial.union.erase <- c(base::scale(log(eu.links$s.poly.and.pts.n.poi.industrial.union.erase+1)))
eu.links$d.poly.and.pts.n.poi.industrial.union.erase <- c(base::scale(log(eu.links$d.poly.and.pts.n.poi.industrial.union.erase+1)))

eu.links$predicted <- stats::predict(object = ger.glm, 
                                     newdata = eu.links,
                                     type = "response")
eu.links$predicted <- exp(eu.links$predicted)
eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links)


nodes <- st_read("Traffic/eur.urban.areas.gpkg", # urban areas with gdp, pop, n.poi
                 "updated.eu.centroids.09.07.2025.w.POIs")
nodes <- nodes[lengths(st_intersects(nodes, 
                                     st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                  st_crs(nodes)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters


tap.mag.ini <- comp %>%
  dplyr::filter(year <= 1980) %>%
  st_intersection(st_union(gadm.0.native, gadm.1.native)) 

st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
tap.mag.n.native <- comp %>%
  st_as_sf() %>%
  st_intersection(gadm.0.n.native) %>%
  st_erase(gadm.0.native) %>%
  st_erase(gadm.1.native)

#plot(gbm.r)
#plot(tap.mag.ini[,1], add = TRUE, col = "red", pch = 19)
#plot(tap.mag.n.native[,1], add = TRUE , col = "yellow", pch = 19)


empty.r <- gbm.r
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything

native.dist.r <- empty.r %>%
  terra::mask(vect(st_buffer(tap.mag.ini, 10000)),
              updatevalue = 1,
              inverse = TRUE)
# plot(native.dist.r)

ref.dist.r <- native.dist.r %>%
  terra::mask(vect(st_buffer(tap.mag.n.native, 10000)),
              updatevalue = 1,
              inverse = TRUE) 

# plot(ref.dist.r)

#gadm.2.ref$occupied <- "" # for reference and AUC calculation only the known non-native regions are used.
#gadm.2.ref$pred.occ <- ""
#for(i in 1:nrow(gadm.2.ref)){
#  ref.r <- terra::mask(ref.dist.r, 
#                       vect(gadm.2.ref[i,]), 
#                       inverse = FALSE, 
#                       touches = TRUE) 
#  ref.a <- as.polygons(ref.r) |> 
#    disagg() |> 
#    st_as_sf() |> 
#    st_intersection(gadm.2.ref[i,]) 
#  ref.a$area <- units::drop_units(st_area(ref.a))
#  
#  if(nrow(ref.a) > 0){
#    if(ref.a %>% 
#       dplyr::filter(layer == 1) %>% 
#       nrow() == 0) {gadm.2.ref[i,]$occupied <- 0} else {
#         if(100/sum(ref.a$area) * sum(ref.a[which(ref.a$layer == 1),]$area, na.rm = TRUE) >= 15){
#           # result.r <- terra::mask(result.r, vect(test.a), updatevalue = 1)
#           gadm.2.ref[i,]$occupied <- 1
#         } else {gadm.2.ref[i,]$occupied <- 0
#         }
#       }
#    
#  } else {}
#}
#gadm.2.ref$occupied <- as.numeric(gadm.2.ref$occupied)
#gadm.2.ref <- gadm.2.ref %>% 
#  dplyr::filter(!is.na(occupied)) # Islas Canarias removed
#gadm.2.ref$pred.occ <- as.numeric(gadm.2.ref$pred.occ)
#st_write(gadm.2.ref, 
#         "gadm.2.ref.gpkg")

#gadm.2.ref <- st_read("gadm.2.ref.gpkg")

#gadm.2.ref$accu.calc <- ""
#gadm.2.ref[lengths(st_intersects(gadm.2.ref, st_buffer(st_as_sf(bind_rows(tap.mag.n.native, tap.mag.ini)), 100000))) > 0, ]$accu.calc <- 1
#gadm.2.ref$accu.calc <- as.numeric(gadm.2.ref$accu.calc)  





ua <- st_read("Traffic/eur.urban.areas.gpkg", 
              layer = "urban.areas.update.09.07.2025.w.POIs")
ua <- ua[lengths(st_intersects(ua, st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), st_crs(ua)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters

#
#
#
#
#
# insert here: sample habitat suitability of ua's from gbmr and always remove the ones which are too low after each time.step
#ua$suitability <- ""
#ua$suitability <- as.numeric(ua$suitability)
#for(i.ua in 1 : nrow(ua)){
#  ua[i.ua,]$suitability <- terra::extract(x = gbm.r, y = vect(ua[i.ua,]), fun = mean, method = "simple", ID = FALSE, na.rm = TRUE)[[1]]
#}
#ua <- ua %>% dplyr::filter(suitability >= 0.2)
#
#
#
#
#

ini.nodes <- nodes %>% 
  dplyr::filter(Name %in% c(ua[lengths(st_intersects(ua, tap.mag.ini)) > 0,]$Name))
add.nodes <- tibble(ID = nodes$ID, is.occupied = terra::extract(native.dist.r, vect(nodes), ID = FALSE)[[1]])
add.nodes <- add.nodes %>% 
  dplyr::filter(is.occupied > 0)
ini.nodes <- unique(c(ini.nodes$ID, add.nodes$ID))

#gadm.2.n.native$area <- 
#  units::drop_units(st_area(gadm.2.n.native))
#gadm.2.n.native$occupied <- ""

agg.acc.fact <-  10
acc.vect <- gadm.0
ext(ref.dist.r) == ext(empty.r)
ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% 
  terra::mask(vect(acc.vect))
ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg)
ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg)


spread.val <- 1
thresh.disp.factor <- 0.9
time.steps <- 16
# min.tr <- quantile(eu.links$predicted, probs = c(0.3), na.rm = TRUE)
min.tr <- 10 # derived from the reference distribution (see below)
max.dist <- 350000



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


par(mfrow = c(2,2))
plot(native.dist.r, 
     main = "initial distribution", 
     background = "darkgrey")
plot(mask(out[[1]],vect(gadm.0)), 
     main = paste(
       "s.v =", spread.val, ";",
       "init = 1;",
       "\nn.c.o =", nodes.cut.off, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "transformation ²"), 
     background = "darkgrey")
plot(mask(out[[1]],vect(gadm.0)), 
     main = paste(
       "s.v =", spread.val, ";",
       "init = 1;",
       "\nn.c.o =", nodes.cut.off, ";",
       "t.s =", time.steps, ";",
       "\nt.d.f =", thresh.disp.factor, ";",
       "ini.points = green ;",
       "\nreached ua (sim.) = red"), 
     background = "darkgrey")
plot(tap.mag.ini, add = TRUE, col = "green") # starting points
#plot(out[[2]], add = TRUE, col = "red") # reached ua
plot(dplyr::filter(nodes, ID %in% c(out[[2]])), add = TRUE, col = "red")
plot(ref.dist.r, 
     main = "final reference distribution", 
     background = "darkgrey")
# plot(out[[2]], add = TRUE, col = "red", pch = 19) # reached ua
#summary(out[[2]]$tot.prop.traff)
out[[3]]
out[[4]] <- mask(out[[4]], vect(gadm.0))
plot(out[[4]])



library(tmap)
out.pol <- as.polygons(out[[4]][[12]]) %>%
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

writeRaster(out[[4]], "C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma nigerrimum complex/10.11.2025.tiff")


### tapinoma seifert data ====
# define inputs needed for function:
#ref.raster <- rast("C:/Users/JLU-SU/Documents/PredictiveAliens/dgm1000_utm32s.asc") # 1km² Digitales Geländemodell from https://gdz.bkg.bund.de/index.php/default/digitales-gelandemodell-gitterweite-1000-m-dgm1000.html ; 16.05.2025
#terra::crs(ref.raster) <- "+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # assign crs
#ref.raster <- terra::aggregate(ref.raster, 5) # do with all rasters that are used!
# Info from Manuelas Database:
native <- c("Corse", "Italy", "Sardegna", "Sicily", "Spain") # Sardegna = Sardinia, Corse = Corsica; "Morocco", "Tunisia","Algeria" left out because here only europe 
non.native <- c("Belgium", "France", "Germany", "Netherlands", "Slovenia", "Switzerland") # this line are known non.native countries
pot.countries <- c("Portugal", "United Kingdom", "Ireland", "Norway", "Sweden", "Finland", # these are potential countries
                   "Estonia", "Latvia", "Lithuania", "Poland", "Czechia", "Austria", "Croatia",
                   "Greece", "Slovenia", "San Marino", "Bosnia and Herzegovina", "Serbia", "Montenegro",
                   "Albania", "Bulgaria", "Romania", "Ukraine", "Moldova", "Belarus",
                   "Luxembourg", "Denmark", "Turkey", "Hungary", "Iceland", "Sicily","Slovakia", "Kosovo", "North Macedonia", "Montenegro") 
# "Azerbaijan" left out because no data available, majority in Europe

# filter with GADM? instead of years => most records are from after 2010, even in the native region so I try here with regional separation, not temporal
gadm.0 <- st_read("GADM/world_gadm_410-levels.gpkg", 
                  layer = "ADM_0") #%>% 
#st_transform(3035)
gadm.0.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% native)
gadm.0.n.native <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% non.native)
gadm.0.pot <- gadm.0 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1 <- st_read("GADM/world_gadm_410-levels.gpkg", 
                  layer = "ADM_1") #%>% 
#st_transform(3035)
gadm.1.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% native)
gadm.1.n.native <- gadm.1 %>% 
  dplyr::filter(NAME_1 %in% non.native)
gadm.1.pot <- gadm.1 %>% 
  dplyr::filter(COUNTRY %in% pot.countries)

gadm.1.add <- gadm.1 %>% dplyr::filter(COUNTRY %in% c("North Macedonia" , "Moldova", "Montenegro", "Cyprus") | NAME_1 == "Kaliningrad") # are not in gadm.2 level.....

gadm.2 <- st_read("GADM/world_gadm_410-levels.gpkg", 
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

e <- ext(c(-11.8167128633999, 25.9369450325551, 35.5922736627615, 54.2843855437946))

tapintro <- readxl::read_excel("C:/Users/JLU-SU/Downloads/TAPINTRO.xlsx")
tapinoc <- readxl::read_excel("C:/Users/JLU-SU/Downloads/TAPINO_C.xlsx") %>% 
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

#gbm.r <- rast("C:/Users/JLU-SU/Documents/PredictiveAliens/tapinoma magnum/biomod2_GBM.first.try.with.traffic.tif")# biomod SDM output from boosted regression trees (gbm)
gbm.r <- rast("tapinoma nigerrimum complex/biomod2_GBM.bio01.bio10.bio11.bio13.bio14.LC.pop.density.tif")
names(gbm.r) <- "layer"
gbm.r <- crop(gbm.r, e)
gbm.r <- subst(gbm.r, NA, 0)

gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation
#gbm.r <- terra::mask(gbm.r,vect(gadm.0))

gbm.r <- aggregate(gbm.r, 5, 
                   fun = mean, 
                   na.rm = TRUE) # making the raster more coarse; INCLUDE IN ASGRID() BELOW!!

mask.thresh <- 0.35 
mask <- which.lyr(gbm.r[[1]] <= mask.thresh) %>%  # gets a spatraster that has only cells which are 0 in gbm.r (i.e. which are unsuitable)
  #terra::mask(vect(gadm.0)) %>%           # crops it to the area of interest -> CHECK IF NECESSARY
  cells()                               # gets the cell numbers; these are then set to 0 (i.e. unoccupied in the result.r in the dispersal() function)

gbm.r[gbm.r < mask.thresh] <- 0
gbm.r <- gbm.r^2 # exponential conversion instead of linear.
gbm.r <- 1 / max(values(gbm.r), na.rm = TRUE) * gbm.r # with this step it is set to a scale of 0 to 1 irrespective of the transformation

gbm.r <- terra::mask(gbm.r,vect(gadm.0))

gbm.r.inv <- gbm.r*-1 + max(values(gbm.r[[1]]), na.rm = TRUE) # invert raster for creation of resistance matrix

#limit <- max(values(gbm.r.inv), na.rm = TRUE)
spread.limit <- 9999 #limit*100 # set value for areas which are not crossable
gbm.r.inv <- subst(x = gbm.r.inv, # must not have NAs for the function below, so replace with spread.limit to make these areas not crossable
                   from = NA, 
                   to = spread.limit) # the final
# plot(gbm.r.inv) # no network visible in plot anymore if limit much higher than maximum value in network because of color-scale. reduce spread.limit to make it visible again.

gbm.r.inv <- asgrid(gbm.r.inv, # convert to grid for spread function
                    xll = xmin(gbm.r.inf),
                    yll = ymin(gbm.r.inv),
                    cellsize = 5000) # update cellsize with the aggregate factor * 1000m


eu.links <- st_read("C:/Users/JLU-SU/Documents/PredictiveAliens/Traffic/gravity.output/europe/shortest.paths.update.29.10.2025/1031.1165ua.GHS.800km.gpkg", 
                    layer = "1031.1165ua.GHS.800km.exp.predict")
eu.links$length <- eu.links$original.dist
eu.links$link.id <- 1:nrow(eu.links)
eu.links <- eu.links %>% 
  rename(o.ID = start.ID)
eu.links <- eu.links %>% 
  rename(d.ID = dest.ID)

eu.links <- eu.links[lengths(st_intersects(eu.links, 
                                           st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                        st_crs(eu.links)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters

eu.links.geom <- eu.links # to save it with geometries for later plotting
eu.links <- st_drop_geometry(eu.links)


nodes <- st_read("C:/Users/JLU-SU/Documents/PredictiveAliens/Traffic/eur.urban.areas.gpkg", 
                 layer = "GHS.29.10.2025.points")
nodes <- nodes[lengths(st_intersects(nodes, 
                                     st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), 
                                                  st_crs(nodes)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters


tap.mag.ini <- comp %>%
  dplyr::filter(year <= 2015) %>%
  st_intersection(st_union(gadm.0.native, gadm.1.native)) 

st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
tap.mag.n.native <- comp %>%
  st_as_sf() %>%
  st_intersection(gadm.0.n.native) %>%
  st_erase(gadm.0.native) %>%
  st_erase(gadm.1.native)

plot(gbm.r)
plot(tap.mag.ini[,1], add = TRUE, col = "red", pch = 19)
plot(tap.mag.n.native[,1], add = TRUE , col = "yellow", pch = 19)
plot(st_buffer(tap.mag.n.native[,1], 50000), add = TRUE , col = "yellow", pch = 19)


empty.r <- gbm.r
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything

native.dist.r <- empty.r %>%
  terra::mask(vect(st_buffer(tap.mag.ini, 10000)),
              updatevalue = 1,
              inverse = TRUE)
# plot(native.dist.r)

ref.dist.r <- native.dist.r %>%
  terra::mask(vect(st_buffer(tap.mag.n.native, 10000)),
              updatevalue = 1,
              inverse = TRUE) 

# plot(ref.dist.r)

#gadm.2.ref$occupied <- "" # for reference and AUC calculation only the known non-native regions are used.
#gadm.2.ref$pred.occ <- ""
#for(i in 1:nrow(gadm.2.ref)){
#  ref.r <- terra::mask(ref.dist.r, 
#                       vect(gadm.2.ref[i,]), 
#                       inverse = FALSE, 
#                       touches = TRUE) 
#  ref.a <- as.polygons(ref.r) |> 
#    disagg() |> 
#    st_as_sf() |> 
#    st_intersection(gadm.2.ref[i,]) 
#  ref.a$area <- units::drop_units(st_area(ref.a))
#  
#  if(nrow(ref.a) > 0){
#    if(ref.a %>% 
#       dplyr::filter(layer == 1) %>% 
#       nrow() == 0) {gadm.2.ref[i,]$occupied <- 0} else {
#         if(100/sum(ref.a$area) * sum(ref.a[which(ref.a$layer == 1),]$area, na.rm = TRUE) >= 15){
#           # result.r <- terra::mask(result.r, vect(test.a), updatevalue = 1)
#           gadm.2.ref[i,]$occupied <- 1
#         } else {gadm.2.ref[i,]$occupied <- 0
#         }
#       }
#    
#  } else {}
#}
#gadm.2.ref$occupied <- as.numeric(gadm.2.ref$occupied)
#gadm.2.ref <- gadm.2.ref %>% 
#  dplyr::filter(!is.na(occupied)) # Islas Canarias removed
#gadm.2.ref$pred.occ <- as.numeric(gadm.2.ref$pred.occ)
#st_write(gadm.2.ref, 
#         "gadm.2.ref.gpkg")

#gadm.2.ref <- st_read("gadm.2.ref.gpkg")

#gadm.2.ref$accu.calc <- ""
#gadm.2.ref[lengths(st_intersects(gadm.2.ref, st_buffer(st_as_sf(bind_rows(tap.mag.n.native, tap.mag.ini)), 100000))) > 0, ]$accu.calc <- 1
#gadm.2.ref$accu.calc <- as.numeric(gadm.2.ref$accu.calc)  





ua <- st_read("Traffic/eur.urban.areas.gpkg", 
              layer = "GHS.29.10.2025")
ua <- ua[lengths(st_intersects(ua, st_transform(st_as_sf(vect(e, crs = "EPSG:4326")), st_crs(ua)))) > 0, ]# crop the eu.links to the area of interest to reduce computing time and avoid errors when they lead out of the cropped rasters

ini.nodes <- nodes %>% 
  dplyr::filter(gc_ucn_mai_2025 %in% c(ua[lengths(st_intersects(ua, tap.mag.ini)) > 0,]$gc_ucn_mai_2025))
add.nodes <- tibble(ID = nodes$ID, is.occupied = terra::extract(native.dist.r, vect(nodes), ID = FALSE)[[1]])
add.nodes <- add.nodes %>% 
  dplyr::filter(is.occupied > 0)
ini.nodes <- unique(c(ini.nodes$ID, add.nodes$ID))

agg.acc.fact <-  10
#acc.vect <- gadm.0 # this is the mask for validation area, filter for non.native/known occupied areas
acc.vect <- st_union(st_buffer(tap.mag.n.native[,1], 120000))
ext(ref.dist.r) == ext(empty.r)
ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% 
  terra::mask(vect(acc.vect))
ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg)
ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg)



spread.val <- 1
thresh.disp.factor <- 0.9
time.steps <- 15
min.tr <- quantile(eu.links$predicted, probs = c(0.95), na.rm = TRUE)[[1]]
#min.tr <- 0.5 # derived from the reference distribution (see below)
max.dist <- 350000

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
plot(out[[4]][[17:30]])


# optimization trial ====
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

