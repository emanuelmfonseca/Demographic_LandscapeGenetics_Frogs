input_fastsimcoal_spatial <- function(working_folder,
                                      models_file,
                                      model_number,
                                      simulation_number,
                                      PopSizeAnc,
                                      TDIV,
                                      Number_ind_loci,
                                      number_of_demes,
                                      pop_size_range,
                                      number_of_individuals_sampled,
                                      path_to_migration_matrix,
                                      species_demography,
                                      species_assignment){
  setwd(working_folder)
  
  Models <- read.table(models_file,h=T)
  Models_sim <- Models[model_number,]
  layers_sim <- colnames(Models_sim[which(Models_sim == 1)])
  
  pop_size_range <- sample(pop_size_range[1]:pop_size_range[2],number_of_demes,r=T)
  
  if (layers_sim == "Demographic_Spatial_model"){
    species_demography <- read.table(species_demography,h=T)
    pop_demography <- character()
    for (x in 1:length(species_assignment[,1])){
      pop <- names(sample(species_assignment[x,],1,prob=as.numeric(species_assignment[x,])))
      index <- which(names(species_demography) %in% pop == TRUE)
      pop_demography <- c(pop_demography, as.vector(species_demography[as.numeric(index)][1,1]))
    }
    
    if(any(pop_demography %in% "Bottleneck")){
      index <- which(pop_demography %in% "Bottleneck" == TRUE)
      pop_size_range[index] <- round(runif(length(index),0.3,0.5) * pop_size_range[index])
    }}
  
  if (!dir.exists("Simulations")){
    dir.create("Simulations")
  }
  
  setwd("Simulations")
  
  if (!dir.exists(paste0("Simulation_",simulation_number))){
    dir.create(paste0("Simulation_",simulation_number))
  }
  
  setwd(paste0("Simulation_",simulation_number))
  
  Models <- read.table(models_file,h=T)
  Models_sim <- Models[model_number,]
  layers_sim <- colnames(Models_sim[which(Models_sim == 1)])
  
  PopSizeAnc_input <- ceiling(runif(1,PopSizeAnc[1], PopSizeAnc[2]))
  
  load(paste0(path_to_migration_matrix,"/Simulation_",simulation_number,".RData"))
  
  cat("", file="./input.par", append=FALSE, sep = "")
  
  part1 <- "//Parameters for the coalescence simulation program : simcoal.exe"
  cat(part1, file="./input.par", append=TRUE, sep = "\n")
  
  part2 <- paste(number_of_demes+1, "samples to simulate :")
  cat(part2, file="./input.par", append=TRUE, sep = "\n")
  
  part3 <- "//Population effective sizes (number of genes)"
  cat(part3, file="./input.par", append=TRUE, sep = "\n")
  
  for (x in 1:(number_of_demes+1)){
    if (x == 1){
      cat(PopSizeAnc_input, file="./input.par", append=TRUE, sep = "\n")
    } else {
      cat(pop_size_range[x-1], file="./input.par", append=TRUE, sep = "\n")
    }}
  
  part4 <- "//Samples sizes and samples age"
  cat(part4, file="./input.par", append=TRUE, sep = "\n")
  
  for (x in 1:(length(number_of_individuals_sampled)+1)){
    if (x == 1){
      cat("0", file="./input.par", append=TRUE, sep = "\n")
    } else {
      cat(number_of_individuals_sampled[x-1], file="./input.par", append=TRUE, sep = "\n")
    }}
  
  part5 <- "//Growth rates: negative growth implies population expansion"
  cat(part5, file="./input.par", append=TRUE, sep = "\n")
  
  cat("0", file="./input.par", append=TRUE, sep = "\n")
  for (x in 1:number_of_demes){
    cat("0", file="./input.par", append=TRUE, sep = "\n")
  }
  
  part6 <- "//Number of migration matrices : 0 implies no migration between demes"
  cat(part6, file="./input.par", append=TRUE, sep = "\n")
  cat(length(migration_matrix), file="./input.par", append=TRUE, sep = "\n")
  
  for (x in 1:length(migration_matrix)){
    cat("//Migration matrix", file="./input.par", append=TRUE, sep = "\n")
    write.table(migration_matrix[[x]], file="./input.par", append=TRUE,col.names = F, row.names = F)
  }
  
  part7 <- "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index"
  cat(part7, file="./input.par", append=TRUE, sep = "\n")
  
  if (length(migration_matrix) == 2){
    part8 <- paste((number_of_demes), "historical event")
    cat(part8, file="./input.par", append=TRUE, sep = "\n")
  } else {
    part8 <- paste((number_of_demes+length(migration_matrix)-2), "historical event")
    cat(part8, file="./input.par", append=TRUE, sep = "\n")
  }
  
  if (layers_sim == "Demographic_Spatial_model"){
    resize <- numeric()
    for (x in 1:length(pop_size_range)){
      if (pop_demography[x] == "Expansion"){
        anc_pop_size <- round(round(runif(1,0.05,0.2) * pop_size_range[x])/pop_size_range[x],digits = 5)
        resize <- c(resize, anc_pop_size)
      } else if (pop_demography[x] == "Bottleneck") {
        anc_pop_size <- round(round(runif(1,5,20) * pop_size_range[x])/pop_size_range[x],digits = 5)
        resize <- c(resize, anc_pop_size)
      } else if (pop_demography[x] == "Constant"){
        resize <- c(resize, 1)
      }}} else {
        resize <- rep(1,(length(TDIV)+1))
      }
  
  TDIV <- round(runif(1,TDIV[1],TDIV[2]))
  
  index <- 1
    for (x in 1:(number_of_demes)){
        cat(paste(TDIV, x ," 0 1 1 0 1"), file="./input.par", append=TRUE, sep = "\n")
      }
  
  
  part9 <- "//Number of independent loci [chromosome] "
  cat(part9, file="./input.par", append=TRUE, sep = "\n")
  
  cat(Number_ind_loci, file="./input.par", append=TRUE, sep = "\n")
  
  part10 <- "//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1"
  cat(part10, file="./input.par", append=TRUE, sep = "\n")
  
  part11 <- "//per Block:data type, number of loci, per gen recomb and mut rates"
  cat(part11, file="./input.par", append=TRUE, sep = "\n")
  
  part12 <- "SNP 1 0"
  cat(part12, file="./input.par", append=TRUE, sep = "\n")
  
}
