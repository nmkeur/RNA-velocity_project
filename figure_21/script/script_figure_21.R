#load libraries
library(Scribe)
library(monocle)
library(cowplot)

#Load CDS object
load('~/Desktop/CDS_predicted_velocity_10_genes.RData')

#Change variable names for figures
names(pData(CDS_predicted_velocity))[names(pData(CDS_predicted_velocity)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted_velocity))[names(pData(CDS_predicted_velocity)) == 'stim'] <- 'Stimulation'
CDS_predicted_velocity$Batch <- as.character(CDS_predicted_velocity$Batch)
names(pData(CDS_predicted_velocity))[names(pData(CDS_predicted_velocity)) == 'sample'] <- 'Sample'

#This rescales the pseudotime to match the values used by monocle
ps <- seq(1,13,0.001148501)
CDS_predicted_velocity$Pseudotime <- ps[1:8707]

#Order CDS by Pseudotime
CDS_predicted_velocity <- CDS_predicted_velocity[, order(pData(CDS_predicted_velocity)$Pseudotime)]

a <- plot_cell_trajectory(CDS_predicted_velocity,
                          color_by = 'Pseudotime',
                          cell_size = 0.4)



b <- plot_genes_in_pseudotime(CDS_predicted_velocity['IRF1'],
                         color_by = 'Stimulation',
                         cell_size = 0.4)+
                         guides(colour = guide_legend(override.aes = list(size=3)))

png("~/RNA-velocity_project/figure_21/figure/figure_21.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a,b , labels = c("A","B"), nrow = 2)
dev.off()
