### define input which is not varied -------------------------------------------

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


## parameter estimation --------------------------------------------------------
parameters <- tidyr::crossing(
  net.spread = c(TRUE, FALSE), 
  spread.val = c(1), # , 1.5, 2
  thresh.disp.factor = c(.05,.5,.95), #0.5, 0.75, 1
  transformation = c("power.1"), # , "power.2", "power.3"
  min.tr.quantile = c(.05,.5,.95),
  max.dist.quantile = c(.05,.5,.95),
  time.steps = c(15)
)
gc()

start.time <- Sys.time()
for(i.p in 1:nrow(parameters)){
  print(paste("i.p:", i.p))
  spread.val <- parameters[i.p,]$spread.val
  thresh.disp.factor <- parameters[i.p,]$thresh.disp.factor
  time.steps <- parameters[i.p,]$time.steps
  min.tr <- quantile(traffic, 
                     probs = c(parameters[i.p,]$min.tr.quantile), na.rm = TRUE)[[1]]
  max.dist <- quantile(distances, 
                       probs = c(parameters[i.p,]$max.dist.quantile), na.rm = TRUE)[[1]]
  gbm.r.inv <- get(parameters[i.p,]$transformation)
  
  
  out <- par.dispersal(
    land.spread = TRUE,
    net.spread = parameters[i.p,]$net.spread,
    
    dist.ini = myo.coy.xy.ini,
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
          "C:/Users/JLU-SU/Nextcloud/Predictive Aliens/data/simulation output/myocastor coypus/10.06.2026.parameter estimation.csv",
          row.names = FALSE)
