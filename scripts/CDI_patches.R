library(random)
library(stats)
library(plyr)
source("~/Google_Drive/functions.R")

# growth

setwd("")
cell_data<-read.csv('CDI_strips_normalised_3.csv',head=TRUE)

find_patches <- function(strip){
  patches <- c()
  cellType <- c()
  patch <- 1
  for(i in 1:length(strip)){
    if(i!= length(strip)){
      if(strip[i] == strip[i+1]){
        patch <- patch + 1
      }
      if(strip[i] != strip[i+1]){
        patches <- c(patches,patch)
        cellType <- c(cellType,strip[i])
        patch <- 1
      }
    }
    if(i== length(strip)){
      patches <- c(patches,patch)
      cellType <- c(cellType,strip[i])
    }
  }
  return(list('patch' = patches,'cellType' = cellType))
}


# group and rep
patch_sizes <- data.frame(growth_inhibited=double(),inhb_r=double(),density=integer(),rep=integer(),group=factor(),radius=integer(),patches=integer(),cellType = integer())

cell_data$sim <- interaction(cell_data$rep,cell_data$group,cell_data$radius)

ptm <- proc.time()
for(i in unique(cell_data$sim)){
  patches <- find_patches(cell_data[cell_data$sim == i,]$cellType)
  p<-data.frame(growth_inhibited=unique(cell_data[cell_data$sim == i,]$growth_inhibited),
                inhb_r=unique(cell_data[cell_data$sim == i,]$inhb_r),
             density = cell_data[cell_data$sim == i,]$density[1],
             rep = cell_data[cell_data$sim == i,]$rep[1],
             group = cell_data[cell_data$sim == i,]$group[1],
             radius = as.factor(cell_data[cell_data$sim == i,]$radius[1]),
             patches = patches$patch,
             cellType = as.factor(patches$cellType))
  patch_sizes <- rbind(patch_sizes,p)
}
proc.time() - ptm
write.csv(patch_sizes,'data/cdi_patches_normalised.csv',row.names = FALSE)

