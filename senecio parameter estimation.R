### define input which is not varied -------------------------------------------


## parameter estimation --------------------------------------------------------
parameters <- tidyr::crossing(
  net.spread = c(TRUE), # with TRUE on other core
  spread.val = c(1), # , 1.5, 2
  thresh.disp.factor = c(.05,.5, .95), #0.5, 0.75, 1
  transformation = c("power.1", "power.1.5", "power.2"),
  min.tr.quantile = c(.05,.5,.95),
  max.dist.quantile = c(.05,.5,.95),
  time.steps = c(40)
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
    
    dist.ini = sen.ini,
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
    acc.vect = acc.vect,
    min.tr = min.tr,
    max.dist = max.dist,
    
    plot.result = TRUE,
    cores = 10
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
          "C:/Users/JLU-SU/Nextcloud/Predictive Aliens/data/simulation output/senecio inaequidens/more.optim.output.net.spread.TRUE.csv",
          row.names = FALSE)
