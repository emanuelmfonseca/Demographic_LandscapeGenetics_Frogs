source("./Scripts/Species_information.R")

l.troglodytes_df <- species_information(species = "./Species_information/Leptodactylus_troglodytes_data.txt",
                                        species_dem = "./Species_information/Leptodactylus_troglodytes_demography.txt",
                                        pop_column = c(4,5,6))
rownames(l.troglodytes_df) <- 1:length(rownames(l.troglodytes_df))

write.table(l.troglodytes_df, "./Species_information/Leptodactylus_troglodytes_species_information.txt",row.names = F, quote = F)
