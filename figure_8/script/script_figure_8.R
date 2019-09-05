#Load required libraries
library(monocle)
library(cowplot)

#Load CDS file used in analysis
CDS_predicted <- readRDS('~/Desktop/CDS_predicted_fulldataset_INF.Rds')

#Change the variable names for figures
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'

#Plot the trajectory
plot_cell_trajectory(CDS_predicted, color_by = 'State')
plot_cell_trajectory(CDS_predicted, color_by = 'Stimulation')
plot_cell_trajectory(CDS_predicted, color_by = 'Sample')
plot_cell_trajectory(CDS_predicted, color_by = 'Batch')

#Plot the trajectory with adjusted cell size and scaled legend size
a <- plot_cell_trajectory(CDS_predicted, color_by = 'State', cell_size = 0.4) +guides(colour = guide_legend(override.aes = list(size=3)))
b <- plot_cell_trajectory(CDS_predicted, color_by = 'Stimulation', cell_size = 0.4)+guides(colour = guide_legend(override.aes = list(size=3)))
c <- plot_cell_trajectory(CDS_predicted, color_by = 'Sample',cell_size = 0.4)+guides(colour = guide_legend(override.aes = list(size=3)))
d <- plot_cell_trajectory(CDS_predicted, color_by = 'Batch',cell_size = 0.4)+guides(colour = guide_legend(override.aes = list(size=3)))

#Write the figure to directory using 200 DPI for increased resolution
png("~/RNA-velocity_project/figure_8/figure/figure_8.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a,b,c,d, labels = c("A","B","C","D"))
dev.off()
