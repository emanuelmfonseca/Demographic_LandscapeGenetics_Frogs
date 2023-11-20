input_fastsimcoal_demographic <- function(working_folder,
                                      models_file,
                                      model_number,
                                      simulation_number,
                                      PopSize,
                                      TDIV,
                                      TEXP,
                                      Number_ind_loci,
                                      number_of_individuals_sampled,
                                      migration_rate,
                                      species_demography,
                                      species_assignment,
                                      pop){
  setwd(working_folder)
  
  pop_size_range <- sample(PopSize[1]:PopSize[2],length(unique(pop)),r=T)
  
  if (!dir.exists("Simulations")){
    dir.create("Simulations")
  }
  
  setwd("Simulations")
  
  if (!dir.exists(paste0("Simulation_",simulation_number))){
    dir.create(paste0("Simulation_",simulation_number))
  }
  
  setwd(paste0("Simulation_",simulation_number))
  
  PopSize_input <- ceiling(runif(1,PopSize[1], PopSize[2]))
  
  cat("", file="./input.par", append=FALSE, sep = "")
  
  part1 <- "//Parameters for the coalescence simulation program : simcoal.exe"
  cat(part1, file="./input.par", append=TRUE, sep = "\n")
  
  part2 <- paste(length(pop_size_range), "samples to simulate :")
  cat(part2, file="./input.par", append=TRUE, sep = "\n")
  
  part3 <- "//Population effective sizes (number of genes)"
  cat(part3, file="./input.par", append=TRUE, sep = "\n")
  
  for (x in 1:length(pop_size_range)){
      if(x != 3) {
        cat(pop_size_range[x], file="./input.par", append=TRUE, sep = "\n")
      } else {
        cat(round(pop_size_range[x]*runif(1,0.1,0.3)), file="./input.par", append=TRUE, sep = "\n")
      }}
  
  part4 <- "//Samples sizes and samples age"
  cat(part4, file="./input.par", append=TRUE, sep = "\n")
  
  sampling <- character()
  for (x in 1:length(species_assignment[,1])){
    set.seed(simulation_number)
    sampling <- c(sampling, sample(names(species_assignment),1,prob=as.numeric(species_assignment[x,])))
  }
  
  for (x in names(species_assignment)){
    index <- which(sampling == x)
    cat(sum(number_of_individuals_sampled[index]), file="./input.par", append=TRUE, sep = "\n")
  }
  
  part5 <- "//Growth rates: negative growth implies population expansion"
  cat(part5, file="./input.par", append=TRUE, sep = "\n")
  
  for (x in 1:length(pop_size_range)){
      cat("0", file="./input.par", append=TRUE, sep = "\n")
    }
  
  part6 <- "//Number of migration matrices : 0 implies no migration between demes"
  cat(part6, file="./input.par", append=TRUE, sep = "\n")
  cat(2, file="./input.par", append=TRUE, sep = "\n")
  
  mig_matrix <- matrix(0,(length(TDIV)/2+1),(length(TDIV)/2+1))
  
  index <- which(lower.tri(mig_matrix) == TRUE)
  
  mig_rate <- runif(length(index),migration_rate[1],migration_rate[2])
  mig_matrix[index] <- mig_rate
  
  makeSymm <- function(mt) {
    mt[upper.tri(mt)] <- t(mt)[upper.tri(mt)]
    return(mt)
  }
  
  mig_matrix <- makeSymm(mig_matrix)
  
  cat("//Migration matrix", file="./input.par", append=TRUE, sep = "\n")
  write.table(mig_matrix, file="./input.par", append=TRUE,col.names = F, row.names = F)
  
  mig_matrix <- matrix(0,(length(TDIV)/2+1),(length(TDIV)/2+1))
  cat("//Migration matrix", file="./input.par", append=TRUE, sep = "\n")
  write.table(mig_matrix, file="./input.par", append=TRUE,col.names = F, row.names = F)
  
  part7 <- "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index"
  cat(part7, file="./input.par", append=TRUE, sep = "\n")
  
  part8 <- paste((length(pop_size_range)+(length(TDIV)/2)), "historical event")
  cat(part8, file="./input.par", append=TRUE, sep = "\n")
  
  n_splits <- length(TDIV)/2
  index <- n_splits
  
  TVID_min <- seq(1,length(TDIV),by=2)
  TVID_max <- seq(2,length(TDIV),by=2)
  
    for (x in 1:n_splits){
    cat(paste(round(runif(1,TDIV[TVID_min[x]],TDIV[TVID_max[x]])), 0, index, "1 1 0 0"), file="./input.par", append=TRUE, sep = "\n")
      index <- index -1
  }
  
  for (x in 1:length(unique(pop))){
    cat(paste(round(runif(1,TEXP[1],TEXP[2])), (x-1), (x-1), "1", runif(1,0.01,0.1), "0 1"), file="./input.par", append=TRUE, sep = "\n")
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
