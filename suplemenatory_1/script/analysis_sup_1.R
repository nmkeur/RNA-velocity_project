library(monocle)
library(ggplot2)
require(gridExtra)
library(cowplot)
library(Scribe)

CDS_predicted <- readRDS('~/Desktop/github_data/suplemenatory_1/data/CDS_predicted_without_STAT1-IRF1.Rds')
#CDS2_predicted <- load('~/Desktop/CD4_INF_3000_PSEUDO_SCRIBE_RESULTS.Rdata')
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'batch_name'] <- 'Batch'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'stim'] <- 'Stimulation'
names(pData(CDS_predicted))[names(pData(CDS_predicted)) == 'sample'] <- 'Sample'
CDS_predicted$Batch <- as.character(CDS_predicted$Batch)

a <-plot_cell_trajectory(CDS_predicted, color_by = "State", cell_size = 0.4) + guides(colour = guide_legend(override.aes = list(size=2)))
b<- plot_cell_trajectory(CDS_predicted, color_by = "Sample",cell_size = 0.4)+ guides(colour = guide_legend(override.aes = list(size=2)))
c<- plot_cell_trajectory(CDS_predicted, color_by = "Stimulation",cell_size = 0.4)+ guides(colour = guide_legend(override.aes = list(size=2)))
d<- plot_cell_trajectory(CDS_predicted, color_by = "Batch",cell_size = 0.4)+ guides(colour = guide_legend(override.aes = list(size=2)))

png("~/Desktop/github_data/suplemenatory_1/figure/trajectory_without_stat1-irf1.png", width = 3500, height = 1800, units='px', res = 200)
plot_grid(a,b,c,d, labels = c("A","B","C","D"))
dev.off()
