Migration_matrix <- function(working_folder, models_file, model_number, simulation_number, path_to_rasters,all_demes_coords, landscape_effect_file, dispersal_capacity, sp_information,species_assignment){
  
  library(raster)
  library(gtools)
  library(gdistance)
  library(scales)
  
  options(scipen=999)
  options(warn=-1)
  
  quiet <- function(x) {
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }
  
  scale_to_0_1 <- function(x){ 
    values <- getValues(x)
    x[which(values > as.numeric(quantile(x, 0.6)))] <- 1
    (x- cellStats(x,"min"))/(cellStats(x,"max") - cellStats(x,"min"))
  }
  
  scale_0_1 <- function(x){ 
    (x- cellStats(x,"min"))/(cellStats(x,"max") - cellStats(x,"min"))
  }
  
  scale <- function(x){
    (x-min(x))/(max(x)-min(x))
  }
  
  Models <- read.table(models_file,h=T)
  Models_sim <- Models[model_number,]
  layers_sim <- colnames(Models_sim[which(Models_sim == 1)])
  
  if(layers_sim != "Demographic_model"){
    
    layers_names <- list()
    
    setwd(file.path(path_to_rasters,"IBD"))
    n_files <- length(list.files(pattern=".tif"))
    layers_names <- c(layers_names, list(gsub(".tif","",list.files(pattern=".tif"))))
    
    IBD <- raster(list.files(pattern="tif"))
    
    colnames(all_demes_coords) <- c("lon","lat")
    all_demes_coords2 <- SpatialPointsDataFrame(coords = all_demes_coords, data = all_demes_coords,
                                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) 
    
    
    true_resistance_surface <- list()
    true_resistance_surface[[1]] <- IBD
    true_resistance_surface[[1]][] <- 0.1
    
    migration_matrix <- list()
    d.cap <- runif(1,dispersal_capacity[1],dispersal_capacity[2])
    tr <- transition(true_resistance_surface[[1]], function(x) 1/mean(x),8)
    resistance_model <- as.matrix(costDistance(tr, all_demes_coords2))
    resistance_model[which(resistance_model < 1)] <- 1
    resistance_dist <- 1/(resistance_model^3) * runif(1,dispersal_capacity[1],dispersal_capacity[2])
    diag(resistance_dist) <- 0
    migration_matrix[[1]] <- resistance_dist
    migration_matrix[[1]] <- round(migration_matrix[[1]], digits=6)
    
     if (layers_sim == "Demographic_Spatial_model"){
       sampling <- character()
       for (x in 1:length(species_assignment[,1])){
         set.seed(simulation_number)
         sampling <- c(sampling, sample(names(species_assignment),1,prob=as.numeric(species_assignment[x,])))
       }
       
       for (x in 1:length(names(species_assignment))){
         sampling[which(sampling == names(species_assignment)[x])] <- x
       }
       
       sampling <- as.numeric(sampling)
       
       m <- numeric()
       for (x in 1:(length(sampling)-1)){
         for (i in x:(length(sampling)-1)){
           if(sampling[x] == sampling[i+1]){
             m <- c(m,1)
           } else {
             m <- c(m,0)
           }
         }}
       
       makeSymm <- function(mt) {
         mt[upper.tri(mt)] <- t(mt)[upper.tri(mt)]
         return(mt)
       }
       
       index <- which(m == 0)
       sub <- migration_matrix[[1]][lower.tri(migration_matrix[[1]])]
       sub[index] <- sub[index] * runif(1,0.01,0.05)
       mt <- matrix(0,length(all_demes_coords2),length(all_demes_coords2))
       mt[lower.tri(mt,diag=F)] <- sub
       migration_matrix[[1]] <- makeSymm(mt)
     }
    
    migration_matrix[[1]] <- cbind(rep(0,length(all_demes_coords2)), migration_matrix[[1]])
    migration_matrix[[1]] <- rbind(rep(0,length(all_demes_coords2)+1), migration_matrix[[1]])
    
    matrix_zeros <- migration_matrix[[1]]
    matrix_zeros[] <- 0
    
    migration_matrix[[length(migration_matrix)+1]] <- matrix_zeros
    
    setwd(working_folder)
    
    if(!dir.exists("Migration_matrices")){
      dir.create("Migration_matrices")
    }
    
    setwd("Migration_matrices")
    
    save(migration_matrix, file=paste0("Simulation_",simulation_number,".RData"))
    options(warn=0)
  }}
