for (x in 1:3){
	if(!dir.exists('./Scripts_model_L_troglodytes')){
		dir.create('./Scripts_model_L_troglodytes')
	}

  setwd("./Scripts_model_L_troglodytes")
  
  if(!dir.exists(paste0("Model",x))){
    dir.create(paste0("Model",x))
  }
  
  setwd(paste0("Model",x))
  
  R <- paste0('library(parallel)
source("./Scripts/Landscape_effect.R")

Model_number <- ',x,'
n_sim <- 2500

if (Model_number != 1){
  landscape_effect(path_to_folder = paste0("./Models/Model",Model_number),
                   models_file = "./Models/Tested_models.txt",
                   model_number = Model_number,
                   nsim = n_sim)
}

leptodatylus_df <- read.table("./Species_information/Leptodactylus_troglodytes_species_information.txt",h=T)

Simulations <- function(sim, df){
  
  source("./Scripts/Migration_matrices.R")
  source("./Scripts/Input_fsc.R")
  source("./Scripts/Run_fsc.R")
  source("./Scripts/Create_images.R")
  source("./Scripts/Input_fsc_demographic.R")
  
  Model_number <- ',x,'
  
   if (sim == 1){
    cat("",file=paste0("./Simulations_ss/Summary_stats_Model_",Model_number,".txt"),append=F)
  }
  
  if (Model_number != 1){
    Migration_matrix(working_folder = paste0("./Model",Model_number),
                     models_file = "./Models/Tested_models.txt",
                     model_number = Model_number,
                     simulation_number = sim,
                     path_to_rasters = "./Rasters/Working_rasters",
                     all_demes_coords = df[,c("lon","lat")],
                     landscape_effect_file = paste0("./Models/Model",Model_number,"/Landscape_effect.txt"),
                     dispersal_capacity = c(.05,.1),
                     sp_information = df,
                     species_assignment = df[,c(4,5,6)])
  }
  
  if (Model_number != 1){ 
  input_fastsimcoal_spatial(working_folder = paste0("./Model",Model_number),
                    models_file = "./Models/Tested_models.txt",
                    model_number = Model_number,
                    simulation_number = sim,
                    PopSizeAnc = c(20000,100000),
                    TDIV = c(10000,50000),
                    Number_ind_loci = 10000,
                    number_of_demes = length(df[,1]),
                    pop_size_range = c(50,200),
                    number_of_individuals_sampled = df$n_seq*2,
                    path_to_migration_matrix = paste0("./Model",Model_number,"/Migration_matrices"),
                    species_demography = "./Species_information/Leptodactylus_troglodytes_demography.txt",
                    species_assignment = df[,c(4,5,6)])
  } else {
    input_fastsimcoal_demographic(working_folder = paste0("./Model",Model_number),
                                   models_file = "./Models/Tested_models.txt",
                                   model_number = Model_number,
                                   simulation_number = sim,
                                   PopSize = c(20000,100000),
                                   TDIV = c(500000,1000000,100000,250000),
                                   TEXP = c(10000,50000),
                                   Number_ind_loci = 10000,
                                   number_of_individuals_sampled = df$n_seq*2,
                                   migration_rate = c(1e-6,1e-5),
                                   species_demography = "./Species_information/Leptodactylus_troglodytes_demography.txt",
                                   species_assignment = df[,c(4,5,6)],
                                   pop = df$pop)
  }
  
  run_fsc(working_folder = paste0("./Model",Model_number),
          simulation_number = sim,
          path_to_fastsimcoal = "./fsc26")
  
  create_images(working_folder = paste0("./Model",Model_number),
                folder_save_simulations = ".",
                models_file = "./Models/Tested_models.txt",
                model_number = Model_number,
                simulation_number = sim,
                summary_stats_folder = ".",
                pop_assignment = df$pop,
                populations = rep(1:length(df$pop), df$n_seq*2),
                nsim = n_sim,
                n.snps = 782,
                freqs = c(0.01257862, 0.4968553),
                species_assignment = df[,c(4,5,6)])
  
  print(paste0(Model_number,"-",sim))
}

mclapply(1:n_sim, function(sim,df) Simulations(sim,leptodatylus_df), mc.cores=detectCores())

unlink(paste0("./Model",Model_number,"/Simulations"),recursive = T)
unlink(paste0("./Model",Model_number,"/Migration_matrices"),recursive = T)')
  
  PBS <- paste0('#!/usr/bin/bash
#
#PBS -l walltime=05:00:00
#PBS -l nodes=1:ppn=40
#PBS -N FonsecaModel
#PBS -j oe
#PBS -A PAA0202

cd $PBS_O_WORKDIR

cp -r ./Models/Model',x,' .

cd  ./Model',x,'

module load R

Rscript ./Scripts_model/Model',x,'/Model',x,'.R')
  
  cat(R,file=paste0("Model",x,".R"))
  cat(PBS,file="PBS_script.pbs")
  
}

