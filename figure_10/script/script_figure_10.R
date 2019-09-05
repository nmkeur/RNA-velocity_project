#Load required libraries
library(monocle)
library(cowplot)
library(RColorBrewer)
library(colorRamps)

#Load CDS file used in analysis
CDS_predicted <- readRDS('~/Desktop/CDS_predicted_fulldataset_INF.Rds')

#Change the variable names for figures
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'

#Add meanexpression colum to CDS data object
CDS_predicted$MeanExpression <- colMeans(exprs(CDS_predicted))

# Create color continuous color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             breaks=c(0,0.4,0.8,1.2,1.6,2.0,2.4),
                             limits=c(0.1, 2))
#Plot trajectory Pseudotime
a <- plot_cell_trajectory(CDS_predicted ,
                     color_by = 'Pseudotime',
                     cell_size = 0.5) +
  guides(colour = guide_legend(override.aes = list(size=2)))

#Plot mean expression trajectory with new color palette
b <- plot_cell_trajectory(CDS_predicted ,
                     color_by = 'MeanExpression',
                     cell_size = 0.5) +
                     guides(colour = guide_legend(override.aes = list(size=2))) +
                     sc

#Write the figure to directory using 200 DPI for increased resolution
png("~/RNA-velocity_project/figure_10/figure/figure_10.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a,b, labels=c("A","B"), nrow = 2)
dev.off()


