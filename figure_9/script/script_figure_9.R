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



a <- plot_cell_trajectory(CDS_predicted[,CDS_predicted$Sample %in% '1_LLDeep_0253'],
                          cell_size = 0.2, 
                          color_by = 'Batch') +  
                          theme(legend.position = "none")+
                          theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank())

b <- plot_cell_trajectory(CDS_predicted[,CDS_predicted$Sample %in% '1_LLDeep_0358'],
                          cell_size = 0.2,
                          color_by = 'Batch') +
                          theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank())+
                          guides(colour = guide_legend(override.aes = list(size=3)))

c <- plot_cell_trajectory(CDS_predicted[,CDS_predicted$Sample %in% '1_LLDeep_0862'],
                          cell_size = 0.2, 
                          color_by = 'Batch') + theme(legend.position = "none") +
                          theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank())

d <-plot_cell_trajectory(CDS_predicted[,CDS_predicted$Sample %in% 'LLDeep_1603'],
                         cell_size = 0.2,
                         color_by = 'Batch') + 
                         theme(legend.position = "none")+
                         theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank())

e <- plot_cell_trajectory(CDS_predicted[,CDS_predicted$Sample %in% 'LLDeep_1622'],
                          cell_size = 0.2,
                          color_by = 'Batch') +
                          theme(legend.position = "none")+
                          theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank())

f <- plot_cell_trajectory(CDS_predicted[,CDS_predicted$Sample %in% 'LLDeep_1697'],
                          cell_size = 0.2,
                          color_by = 'Batch')  +
                          theme(legend.position = "none")+
                          theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank())

#Write the figure to directory using 200 DPI for increased resolution
png("~/RNA-velocity_project/figure_9/figure/figure_9.png", width = 2000, height = 1800, units='px', res = 200)
plot_grid(a,b,c,d,e,f, labels = c("A","B","C","D","E","F"), nrow=3)
dev.off()
