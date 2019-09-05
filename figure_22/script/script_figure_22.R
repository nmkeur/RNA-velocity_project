#load libraries
library(Scribe)
library(monocle)
library(cowplot)

#Load CDS object
load('~/Desktop/CDS_predicted_velocity_10_genes.RData')
load('~/Desktop/CDS_predicted_monocle_10_genes.RData')

#Change variable names for figures
names(pData(CDS_predicted_velocity))[names(pData(CDS_predicted_velocity)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted_velocity))[names(pData(CDS_predicted_velocity)) == 'stim'] <- 'Stimulation'
CDS_predicted_velocity$Batch <- as.character(CDS_predicted_velocity$Batch)
names(pData(CDS_predicted_velocity))[names(pData(CDS_predicted_velocity)) == 'sample'] <- 'Sample'

#Change variable names for figures
names(pData(CDS_predicted_monocle))[names(pData(CDS_predicted_monocle)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted_monocle))[names(pData(CDS_predicted_monocle)) == 'stim'] <- 'Stimulation'
CDS_predicted_monocle$Batch <- as.character(CDS_predicted_monocle$Batch)
names(pData(CDS_predicted_monocle))[names(pData(CDS_predicted_monocle)) == 'sample'] <- 'Sample'

a <- plot_genes_in_pseudotime(CDS_predicted_monocle['IRF3'],
                              color_by = 'Stimulation',
                              cell_size = 0.4)+
  guides(colour = guide_legend(override.aes = list(size=3)))

b <- plot_genes_in_pseudotime(CDS_predicted_velocity['IRF3'],
                              color_by = 'Stimulation',
                              cell_size = 0.4)+
  guides(colour = guide_legend(override.aes = list(size=3)))

png("~/RNA-velocity_project/figure_22/figure/figure_22.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a,b , labels = c("A","B"), nrow = 2)
dev.off()
