species_information <- function(species, species_dem,pop_column){
  
  library(sp)
  library(raster)
  
  species <- read.table(species,h=T,sep="\t")
  species_dem <- read.table(species_dem,h=T,sep="\t")
  
  species_coords <- species[,c("lon","lat")]
  coordinates(species_coords) <- species_coords
  species_coords_unique <- remove.duplicates(species_coords)
  
  species_coords_unique <- as.data.frame(species_coords_unique@coords)
  
  for (x in 1:length(species_coords_unique[,1])){
    index1 <- which(species[,"lon"] == species_coords_unique[x,1])
    index2 <- which(species[,"lat"] == species_coords_unique[x,2])
    inters <- intersect(index1,index2)
    
    index_names <- which(names(species) %in% names(species_dem) == TRUE)
    populations <- species[,index_names]
    
    populations_reduced <- populations[inters,]
    populations_reduced <- colMeans(populations_reduced)
    
    if (x == 1){
      species_df <- data.frame(matrix(NA,length(species_coords_unique[,1]),3+length(names(populations_reduced))))
      colnames(species_df) <- c("lon","lat","n_seq",names(populations_reduced))
    }
    
    species_df[x, c(1:2)] <- species_coords_unique[x,]
    species_df[x, 3] <- length(inters)
    species_df[x,which(names(species_df) %in% names(populations_reduced) == TRUE)] <- populations_reduced
    
  }
  
  ml_coord_index <- which(species_df[,2] == max(species_df[,2]))
  coords_index <- 1:length(species_df[,2])
  coords_index <- coords_index[-ml_coord_index]
  
  dist <- pointDistance(species_df[ml_coord_index,c(1:2)], species_df[,c(1:2)], lonlat=F, allpairs=FALSE)
  dist_sorted <-sort(dist, index.return=TRUE)
  
  species_df <- species_df[dist_sorted$ix,]
  
  pops <- numeric()
  for (x in 1:length(species_df[,2])){
    pops <- c(pops, which(species_df[x,pop_column] == max(species_df[x,pop_column])))
  }
  
  pops_sorted <-sort(pops, index.return=TRUE)
  species_df2 <- species_df[pops_sorted$ix,]
  
  pop <- pops_sorted$x
  species_df2 <- cbind(species_df2,pop)
  
  return(species_df2)
}
