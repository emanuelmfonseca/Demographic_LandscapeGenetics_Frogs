source("./Scripts/Species_information.R")

R.granulosa_df <- species_information(species = "./Species_information/Rhinella_granulosa_data.txt",
                                        species_dem = "./Species_information/Rhinella_granulosa_demography.txt",
                                      pop_column = c(4,5))
rownames(R.granulosa_df) <- 1:length(rownames(R.granulosa_df))

write.table(R.granulosa_df, "./Species_information/Rhinella_granulosa_species_information.txt",row.names = F, quote = F)
