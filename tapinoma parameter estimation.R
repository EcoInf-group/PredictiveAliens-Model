### define input which is not varied -------------------------------------------

#
#
#
### define input ---------------------------------------------------------------
empty.r <- gbm.r
empty.r[!is.na(empty.r)] <- 0 # create empty raster with no connections or anything, is used in dispersal function to make the result raster

native.dist.r <- empty.r %>%
  terra::mask(vect(st_buffer(tap.mag.ini, 10000)), updatevalue = 1, inverse = TRUE)

ref.dist.r <- c(terra::mask(native.dist.r, vect(st_buffer(comp, 10000)), updatevalue = 1, inverse = TRUE),
                terra::mask(native.dist.r, vect(gadm.3[lengths(st_intersects(gadm.3, comp)) > 0, ]), updatevalue = 1, inverse = TRUE),
                terra::mask(native.dist.r, vect(gadm.2[lengths(st_intersects(gadm.2, comp)) > 0, ]), updatevalue = 1, inverse = TRUE)
)
names(ref.dist.r) <- c("points.10km", "gadm.3", "gadm.2")

acc.vect <- st_union(st_buffer(comp[, 1], 20 * 10 * 1000)) #acc.vect = accuracy vector, in meter, inside this area, accuracy is calculated. time.steps * resolution (or aggregation) * 1000 to make it in km
ref.dist.r.restr <- terra::mask(ref.dist.r, vect(acc.vect))

names(ref.dist.r.restr) <- c("points.10km.200km.radius", "gadm.3.200km.radius", "gadm.2.200km.radius")

ref.dist.r <- c(ref.dist.r, ref.dist.r.restr)

ini.nodes <- nodes %>%
  dplyr::filter(gc_ucn_mai_2025 %in% c(nodes[lengths(st_intersects(nodes, tap.mag.ini)) > 0, ]$gc_ucn_mai_2025)) # nodes which overlap with the known initial distribution are used as ini.nodes (i.e. initial nodes)

add.nodes <- tibble(
  ID = nodes$ID,
  is.occupied = terra::extract(native.dist.r, vect(nodes), ID = FALSE)[[1]]
) # nodes which lie within the known distributional range are added (number of added nodes depends on the buffer used above to make this raster)
add.nodes <- add.nodes %>%
  dplyr::filter(is.occupied > 0)
ini.nodes <- unique(c(ini.nodes$ID, add.nodes$ID))

ext(ref.dist.r) == ext(empty.r) # check if data match.

# acc.vect <- st_union(st_buffer(comp[, 1], 40 * 10 * 1000)) #acc.vect = accuracy vector, in meter, inside this area, accuracy is calculated. time.steps * resolution (or aggregation) * 1000 to make it in km

#
#
#

## parameter estimation --------------------------------------------------------
parameters <- tidyr::crossing(
  net.spread = c(TRUE), # with TRUE on other core
  spread.val = c(1), # , 1.5, 2
  thresh.disp.factor = c(.95), #0.5, 0.75, 1
  transformation = c("power.1"), # , "power.2", "power.3"
  min.tr.quantile = c(.5,.6,.7),
  max.dist.quantile = c(.7,.8,.9),
  time.steps = c(45)
  )
gc()

start.time <- Sys.time()
for(i.p in 1:nrow(parameters)){
  print(paste("i.p:", i.p))
  spread.val <- parameters[i.p,]$spread.val
  thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
  time.steps <- parameters[i.p,]$time.steps
  min.tr <- quantile(eu.links$predicted, probs = c(parameters[i.p,]$min.tr.quantile), na.rm = TRUE)[[1]]
  max.dist <- quantile(distances, probs = c(parameters[i.p,]$max.dist.quantile), na.rm = TRUE)[[1]]
  gbm.r.inv <- get(parameters[i.p,]$transformation)
  
  
  out <- par.dispersal(
    land.spread = TRUE,
    net.spread = parameters[i.p,]$net.spread,
    
    dist.ini = tap.mag.ini,
    spread.val = spread.val,
    thresh.disp.factor = thresh.disp.factor, 
    ini.nodes = ini.nodes,
    
    ref.raster = empty.r,
    result.r = empty.r,
    gbm.r.inv = gbm.r.inv,
    name.gbm.r.inv = parameters[i.p,]$transformation,
    
    time.steps = parameters[i.p, ]$time.steps,
    ref.dist.r = ref.dist.r,
    
    sample.nodes.from.raster = TRUE,
    
    unsuitability.mask = mask,
    #acc.vect = acc.vect,
    min.tr = min.tr,
    max.dist = max.dist,
    
    plot.result = TRUE,
    cores = 15
  )
  
  if(i.p == 1){
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    out$accuracies$thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
    out$accuracies$max.dist.quantile <- parameters[i.p,]$max.dist.quantile
    
    optim.output <- out$accuracies
    
  } else {
    out$accuracies$min.tr.quantile <- parameters[i.p,]$min.tr.quantile
    out$accuracies$transformation <- parameters[i.p,]$transformation
    out$accuracies$thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
    out$accuracies$max.dist.quantile <- parameters[i.p,]$max.dist.quantile
    
    optim.output <- bind_rows(optim.output, out$accuracies)
    
  }
}

finish <- Sys.time()
finish - start.time

write.csv(optim.output, 
          "C:/Users/JLU-SU/Nextcloud/Predictive Aliens/data/simulation output/tapinoma magnum/more.even.more.optim.output.net.spread.TRUE.csv",
          row.names = FALSE)
