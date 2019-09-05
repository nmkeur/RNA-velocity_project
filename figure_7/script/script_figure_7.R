#Load required libraries
library(monocle)

#Load CDS file used in analysis
CDS_predicted <- readRDS('~/Desktop/CDS_predicted_fulldataset_INF.Rds')

#Change the variable names for figures
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'

#Plot the trajectory
plot_cell_trajectory(CDS_predicted)

png("~/RNA-velocity_project/figure_7/figure/figure_7.png", width = 1500, height = 1800, units='px', res = 200)
plot_cell_trajectory(CDS_predicted)
dev.off()
