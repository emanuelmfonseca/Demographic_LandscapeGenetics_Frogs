library(raster)
library(BEDASSLE)
library(hierfstat)
source("./Scripts/Species_information.R")

quiet <- function(x) {
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
}

SNPs <- read.table("./Datasets/L_troglodytes.str")

del <- numeric()
for (x in 2:length(SNPs[1,])){
    if(any(which(SNPs[,x] == -9))){
        del <- c(del,x)
    }
}

SNPs <- SNPs[,-del]

l.troglodytes_df <- species_information(species = "./Species_information/Leptodactylus_troglodytes_data.txt",
                                        species_dem = "./Species_information/Leptodactylus_troglodytes_demography.txt")

l.troglodytes_info <- read.table("./Species_information/Leptodactylus_troglodytes_data.txt",h=T,sep="\t")


populations <- data.frame(matrix(NA,length(l.troglodytes_info[,1]),1))
colnames(populations) <- "Pops"
for (x in 1:length(l.troglodytes_df[,1])){
    coord <- l.troglodytes_df[x,c("lon","lat")]
    index1 <- which(l.troglodytes_info[,"lon"] %in% coord[1] == TRUE)
    index2 <- which(l.troglodytes_info[,"lat"] %in% coord[2] == TRUE)
    index <- intersect(index1, index2)
    
    vouchers <- l.troglodytes_info[index,1]
    populations[which(SNPs[,1] %in% vouchers == TRUE),1] <- x
}

index <- sort(populations[,1], index.return=TRUE, decreasing=F)
populations <- populations[index$ix,1]

SNPs <- SNPs[index$ix,]

SNPs <- SNPs[,-1]

for (x in 1:length(SNPs[1,])){
    len <- length(unique(SNPs[,x]))
    count <- data.frame(matrix(NA,1,len))
    for (i in 1:len){
        count[1,i] <- length(which(SNPs[,x] == unique(SNPs[,x])[i]))
        colnames(count)[i] <- unique(SNPs[,x])[i]
    }
    max <- colnames(count)[which(count == max(count))]
    SNPs[,x] <- ifelse(as.numeric(SNPs[,x]) == as.numeric(max), 0, 1)   
}

freq <- numeric()
for (r in 1:length(SNPs[1,])){
    freq <- c(freq, length(which(as.numeric(SNPs[,r])==1))/length(SNPs[,1]))
}

del <- which(freq < 0.05 | freq > 0.95)
SNPs <- SNPs[,-del]

freq <- numeric()
for (r in 1:length(SNPs[1,])){
    freq <- c(freq, length(which(as.numeric(SNPs[,r])==1))/length(SNPs[,1]))
}

index <- sort(freq, index.return=TRUE, decreasing=T)
SNPs <- SNPs[,index$ix]

img <- raster(nrow=length(SNPs[,1]), ncol=length(SNPs[1,]))
seq <- c(t(SNPs))
img <- setValues(img,as.numeric(seq))

if (!dir.exists("Observed/Leptodactylus_troglodytes")){
    dir.create("Observed/Leptodactylus_troglodytes", recursive = TRUE)
}

png(filename="./Observed/Leptodactylus_troglodytes/Observed_Leptodactylus_troglodytes.png",height=length(SNPs[,1]), width=length(SNPs[1,]))
dat <- matrix(getValues(img), ncol = length(SNPs[,1]), nrow = length(SNPs[1,]))
par(mai=c(0,0,0,0))
image(dat, col=c("black", "white"), axes=F, frame.plot=T)
dev.off()
