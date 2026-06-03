dispersal <- function(land.spread = TRUE, # logical, is spread through the landscape allowed?
                      net.spread = FALSE, # logical, is spread through the network allowed?
                      
                      dist.ini, # initial distribution coordinates from where the species spreads through the landscape
                      spread.val = 1, # how far can the species spread in gridprocess::rawspread?
                      thresh.disp.factor = 0.9, # how much of the spread.val must a pixel receive to be treated as occupied?
                      ini.nodes, # initial urban areas from which the species can spread through the traffic network AND the landscape
                      ref.raster, # reference raster (resolution etc.) for the output
                      result.r, # empty raster that will be updated after each time.step and thus turned into the result raster.
                      gbm.r.inv, # provide resistance matrix
                      
                      initiation = 1, # the initial traffic budget that each used node gets per time.step. this will then be distributed among all outgoing paths relative to the traffic flow on each path.
                      # I tried to scale this with gdp of the respective node, but that did not improve the output. OPEN FOR DISCUSSION.
                      
                      time.steps = 2, # how many iterations should the simulation run?
                      
                      ref.dist.r, # reference raster with the final known distribution with cells being present (1) or absent (0), used for accuracy calculation
                      # the result.r output raster will be compared with this and the values needed for accuracy calculation derived from the comparison.
                      
                      sample.nodes.from.raster = TRUE, # if TRUE urban areas which are within occupied areas in the landscape raster will be treated as 
                      # occupied and can serve as starting point in the traffic network. this allows a species to enter the traffic network from the landscape.
                      
                      unsuitability.mask = NULL, # if provided then every cell which is provided will be set to unoccupied at end of each time step. can be used
                      # to make sure that unsuitable cells will not become occupied. they might still be crossed though.
                      
                      acc.vect = NULL, # accuracy vector; area inside which the accuracy will be calculated.
                      agg.acc.fact = 1, # integer, defines how much the output is be aggregated (terra::aggregate()) before calculation of accuracy (too avoid overly fine-scaled output and accuracy calculations).
                      min.tr = 0, # traffic network: minimum traffic volume that paths must have to be used at all.
                      max.dist = 500000, # traffic network: maximum length of paths in the traffic network (longer distances than this value will not be traveled).
                      
                      plot.result = TRUE # provides output raster if TRUE.
) {
  a.int <- Sys.time() # just to keep track of how much time the simulation needs.
  
  initiation = 1 # the initial traffic budget that each used node gets per time.step. this will then be distributed among all outgoing paths relative to the traffic flow on each path.
  # I tried to scale this with gdp of the respective node, but that did not improve the output. OPEN FOR DISCUSSION.
  
  eu.links <- eu.links %>% # filter for all paths which are connected with the chosen start.node
    dplyr::filter(predicted >= min.tr) %>% # filter all paths that have at least the given traffic volume
    dplyr::filter(length <= max.dist) # filter all paths which are up to the max.dist length
  
  if(agg.acc.fact != 1){
    ref.dist.r.agg <- aggregate(ref.dist.r, agg.acc.fact, fun = "mean", na.rm = TRUE) # aggregate the reference raster to avoid being overly fine-scaled
    ref.dist.r.agg <- ifel(ref.dist.r.agg >= 0.5, 1, ref.dist.r.agg) # make binary (fun = "mean" in aggregate makes non-binary values)
    ref.dist.r.agg <- ifel(ref.dist.r.agg < 0.5, 0, ref.dist.r.agg) # make binary (fun = "mean" in aggregate makes non-binary values)
  } else {
    ref.dist.r.agg <- ref.dist.r}
  
  if(!is.null(acc.vect)){
    ref.dist.r.agg <- terra::mask(ref.dist.r.agg, vect(acc.vect))
  } else{}
  
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
      
      # implement foreach here? ====
      for(i in 1:nrow(id.ini)){ # loop through all start-sites and determine where the species spreads to from here through the resistance layer.
        if(i == 1){spread.m <- rawspread(x = get(gbm.r.inv)[[1]],
                                         spread.value = spread.val,
                                         row = id.ini[i,1],
                                         col = id.ini[i,2]#,
                                         #sd = sd # sd = bandwidth, not used here. "In the standard Gaussian kernel, the “bandwidth” which controls the spread of the kernel is equal to one standard deviation and accounts for 39% of the kernel volume." from: doi:10.1007/s10980-018-0653-9 
        )}else{
          spread.m.c <- rawspread(
            x = get(gbm.r.inv)[[1]],
            # if it is not the first starting point then the output of this one will be added to the previous one, so that all reached cells are collected in one output raster. thus, the threshold can also be reached if a cell is reached just so from many starting locations.
            spread.value = spread.val,
            row = id.ini[i, 1],
            col = id.ini[i, 2]#,
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
      
    } else {} 
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
        
        #updated.ua.i <- dest.nodes$ID
        if(s.n == 1) {updated.ua.i <- dest.nodes$ID} else {updated.ua.i <- unique(c(updated.ua.i, dest.nodes$ID))}
        
        #
        #
        #
        # result generation for consecutive steps:
        # put into if statement above to use only if t.s > 1 to check if optim can better deal with it then
        id.ini.update <- rbind(id.ini.update, # these are used in consecutive landspread.
                               rowColFromCell(ref.raster, cellFromXY(
                                 ref.raster, st_coordinates(dplyr::filter(nodes, ID %in% updated.ua.i))
                               )))
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
      if (sample.nodes.from.raster == TRUE) {
        additional.nodes <- tibble(ID = nodes$ID,
                                   is.occupied = terra::extract(result.r, vect(nodes), ID = FALSE)[[1]])
        additional.nodes <- additional.nodes[which(additional.nodes[, 2] > 0), ]
        updated.ua.i <- sort(unique(c(updated.ua.i, additional.nodes$ID)))
        id.ini.update <- rbind(id.ini.update, rowColFromCell(ref.raster, cellFromXY(
          ref.raster, st_coordinates(dplyr::filter(nodes, ID %in% updated.ua.i))
        )))
        
      } else {
      }
      
    } else {updated.ua.i <- nodes[ini.nodes,]} # end of if net.spread == TRUE
    #
    #
    #
    #
    #
    
    
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
    
    # accuracy measurement with whole raster, better use this option
    if(!is.null(acc.vect)){
      result.r.agg <- aggregate(result.r, agg.acc.fact, fun = "mean", na.rm = TRUE) %>% # changing the resolution of the reference and output rasters can change the accuracy result. the ref.raster is aggregated and masked only once at beginning of the function loop
        terra::mask(vect(acc.vect))
    } else {
      result.r.agg <- aggregate(result.r, agg.acc.fact, fun = "mean", na.rm = TRUE) # changing the resolution of the reference and output rasters can change the accuracy result. the ref.raster is aggregated and masked only once at beginning of the function loop
    }
    
    result.r.agg <- ifel(result.r.agg >= 0.5, 1, result.r.agg)
    result.r.agg <- ifel(result.r.agg < 0.5, 0, result.r.agg)
    
    pred.neg <- cells(result.r.agg, c(0))[[1]] # cell.IDs of predicted absences
    pred.pos <- cells(result.r.agg, c(1))[[1]] # cell.IDs of predicted presences
    
    for(i.lyr in 1:nlyr(ref.dist.r)){
      ref.pos <- cells(ref.dist.r[[i.lyr]], c(1))[[1]] # cell.IDs of reference presences
      ref.neg <- cells(ref.dist.r[[i.lyr]], c(0))[[1]] # cell.IDs of reference absences
      
      TP <- sum(pred.pos %in% ref.pos) # uses cell.IDs to check which are correct/false
      TN <- sum(pred.neg %in% ref.neg) # uses cell.IDs to check which are correct/false
      FP <- sum(pred.pos %in% ref.neg) # uses cell.IDs to check which are correct/false
      FN <- sum(pred.neg %in% ref.pos) # uses cell.IDs to check which are correct/false
      
      precision <- TP / (TP + FP) #  how many of all positives are correctly classified as positive? range: 0 - 1
      sensitivity <- TP / (TP + FN) # of all predicted presences, how many are actual presences? range: 0 - 1
      F1 <- 2 * ((precision * sensitivity) / (precision + sensitivity)) # F-score, range: 0 - 1, especially suited for imbalanced datasets (e.g. were one class if overrepresented), right now seems to be the most adequate by visually comparing the outputs and the references
      assign(paste0("F.", names(ref.dist.r[[i.lyr]])), F1)
      
      if(i.lyr == 1){
        F.scores <- tibble(!!paste0("F.", names(ref.dist.r[[i.lyr]])) := F1)
      } else {
        F.scores.update <- tibble(!!paste0("F.", names(ref.dist.r[[i.lyr]])) := F1)
        F.scores <- cbind(F.scores, F.scores.update)
      }
    }
    
    #
    #
    #
    
    gc()
    print(paste("time.step", t.s, "out of", time.steps, "done"))
    print(Sys.time() - a.int)
    if(t.s == 1){
      plot.stack <- result.r
      accuracy.list <- tibble(F.scores,
                              spread.val = spread.val,
                              thresh.disp.factor = thresh.disp.factor,
                              time.step = t.s,
                              initiation = initiation,
                              agg.acc.fact = agg.acc.fact, 
                              min.tr = min.tr, 
                              max.dist = max.dist,
                              land.spread = land.spread,
                              net.spread = net.spread,
                              min.tr.quantile = min.tr,
                              gbm.r.inv = gbm.r.inv)
    } else {
      plot.stack <- c(plot.stack, result.r)
      accuracy.list.update <- tibble(F.scores,
                                     spread.val = spread.val,
                                     thresh.disp.factor = thresh.disp.factor,
                                     time.step = t.s,
                                     initiation = initiation,
                                     agg.acc.fact = agg.acc.fact, 
                                     min.tr = min.tr, 
                                     max.dist = max.dist,
                                     land.spread = land.spread,
                                     net.spread = net.spread,
                                     min.tr.quantile = min.tr,
                                     gbm.r.inv = gbm.r.inv)
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
