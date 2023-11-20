create_images <- function(working_folder, folder_save_simulations, model_number, simulation_number, pop_assignment, populations, summary_stats_folder, nsim, reps, n.snps, freqs, species_assignment, models_file){
  
  library(dplyr)
  library(raster)
  library(BEDASSLE)
  
  quiet <- function(x) {
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }
  
  setwd(working_folder)
  
  if (!dir.exists("Simulations")){
    dir.create("Simulations")
  }
  
  setwd(file.path(working_folder,"Simulations"))
  setwd(paste0("Simulation_",simulation_number,"/input"))
  
  Models <- read.table(models_file,h=T)
  Models_sim <- Models[model_number,]
  layers_sim <- colnames(Models_sim[which(Models_sim == 1)])
  
  if (layers_sim == "Demographic_model"){
    
    sampling <- character()
    for (x in 1:length(species_assignment[,1])){
      set.seed(simulation_number)
      sampling <- c(sampling, sample(names(species_assignment),1,prob=as.numeric(species_assignment[x,])))
    }
    
    pops <- table(populations)
    pops <- rep(sampling,table(populations))
    
    index <- numeric()
    for (x in unique(names(species_assignment))){
      index <- c(index,length(which(pops == x)))
    }
    
    index2 <- rep(unique(names(species_assignment)),index)
    
    for (x in unique(names(species_assignment))){
      index3<- sample(which(index2 == x),r=F)
      index4 <- sample(which(pops == x),r=F)
      pops[index4] <- index3
    }
    
    pops <- as.numeric(pops)
  }
  
  out <- readLines(list.files(pattern="arp"))
  lines <- which(grepl("^[[:digit:]]", out) == TRUE)
  SNPs <- out[lines]
  SNPs <- sub(".*\\s+(.*)", "\\1", SNPs)
  SNPs <- strsplit(SNPs, "")
  length.snps <- length(SNPs[[1]])
  n_seq <- length(lines)
  SNPs <- data.frame(matrix(unlist(SNPs), ncol = length.snps, nrow=n_seq, byrow = T))
  SNPs <- mutate_all(SNPs, function(x) as.numeric(as.character(x)))
  
  if (layers_sim == "Demographic_model"){
    SNPs <- SNPs[pops,]
  }
  
  freq <- numeric()
  for (r in 1:length.snps){
    freq <- c(freq, length(which(as.numeric(SNPs[,r])==1))/n_seq)
  }
  
  del <-which(freq < freqs[1] | freq > freqs[2])
  SNPs <- SNPs[,-del]
  
  keep <- sample(length(SNPs[1,]),n.snps,r=F)
  SNPs <- SNPs[,keep]
  
  freq <- numeric()
  for (r in 1:length(SNPs[1,])){
    freq <- c(freq, length(which(as.numeric(SNPs[,r])==1))/n_seq)
  }
  
  index <- sort(freq, index.return=TRUE, decreasing=T)
  SNPs <- SNPs[,index$ix]
  
  img <- raster(nrow=length(lines), ncol=n.snps)
  seq <- c(t(SNPs))
  img <- setValues(img,as.numeric(seq))
  
  setwd(working_folder)
  
  if (!dir.exists("Simulations")){
    dir.create("Simulations")
  }
  
  setwd(file.path(working_folder,"Simulations"))
  setwd(paste0("Simulation_",simulation_number,"/input"))
  
  setwd(folder_save_simulations)
  
  if (!dir.exists("Simulated_datasets")){
    dir.create("Simulated_datasets")
  }
  
  setwd("Simulated_datasets")
  
  if (!dir.exists("training_dataset")){
    dir.create("training_dataset")
  }
  
  if (!dir.exists("test_dataset")){
    dir.create("test_dataset")
  }
  
  if (simulation_number <= nsim*0.8){
    setwd(folder_save_simulations)
    setwd("Simulated_datasets")
    setwd("training_dataset")
    
    if (!dir.exists(paste0("Model",model_number))){
      dir.create(paste0("Model",model_number))
    }
    
    setwd(paste0("Model",model_number))
    
    png(filename=paste0("Image_",simulation_number,"_Model",model_number, ".png"),height=length(lines), width=n.snps)
    dat <- matrix(getValues(img), ncol = length(lines), nrow = n.snps)
    par(mai=c(0,0,0,0))
    image(dat, col=c("black", "white"), axes=F, frame.plot=T)
    quiet(dev.off())
    
  } else {
    
    setwd(folder_save_simulations)
    setwd("Simulated_datasets")
    setwd("test_dataset")
    
    if (!dir.exists(paste0("Model",model_number))){
      dir.create(paste0("Model",model_number))
    }
    
    setwd(paste0("Model",model_number))
    
    png(filename=paste0("Image_",simulation_number,"_Model",model_number, ".png"),height=length(lines), width=n.snps)
    dat <- matrix(getValues(img), ncol = length(lines), nrow = n.snps)
    par(mai=c(0,0,0,0))
    image(dat, col=c("black", "white"), axes=F, frame.plot=T)
    quiet(dev.off())
    
  }
  
}