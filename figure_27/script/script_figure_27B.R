set.seed(60)
library(monocle)
library(VGAM)
library(Scribe)
library(plyr)
library(inflection)
library(igraph)
library(ggplot2)
library(GGally)
library(ggraph)
library(pracma)
library(lpSolveAPI)
library(cowplot)
load('~/Desktop/CDS_predicted_velocity_10_genes.RData')
ucrdi_list_velocity <- urdi_list_velocity$max_rdi_value
a <- pheatmap::pheatmap(ucrdi_list_velocity, cluster_rows = F, cluster_cols = F, annotation_names_col = T, border_color = NA,fontsize = 12)

#Filter network edges using CLR
ucrdi_list_velocity[ucrdi_list_velocity <= 0.05] <- 0
ucrdi_list_velocity[dn_sparsifer(ucrdi_list_velocity) ==1] <- 0
ucrdi_list <- clr(ucrdi_list_velocity)

load('~/Desktop/CDS_predicted_velocity_10_genes.RData')
ucrdi_list_velocity <- urdi_list_velocity$max_rdi_value
ucrdi_list_velocity[ucrdi_list <= 0] <- 0

b<- pheatmap::pheatmap(ucrdi_list_velocity, cluster_rows = F, cluster_cols = F, annotation_names_col = T, border_color = NA,fontsize = 12)

# Filter edges using delays
delay <- urdi_list_velocity$max_rdi_delays
d2 <- as.data.frame(as.table(urdi_list_velocity$max_rdi_delays))
for(i in 1:nrow(d2)) {
  row <- d2[i,c(1,2,3)]
  if(delay[row$Var1,row$Var2] > delay[row$Var2,row$Var1]){
    ucrdi_list_velocity[row$Var2,row$Var1] <- 0}
  else{ delay[row$Var1,row$Var2] <- 0}
}

#Create directed graph
net_v <- graph_from_adjacency_matrix(ucrdi_list_velocity, mode = 'directed',weighted = T)

# Create ggraph network plot
set.seed(60)
c <- ggraph(net_v, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name),
                     end_cap = label_rect(node2.name)), 
                 arrow = arrow(length = unit(4, 'mm'))) + 
  geom_node_text(aes(label = name))

# Write figure
png("~/RNA-velocity_project/figure_27/figure/figure_27B.png", width = 1300, height = 1800, units='px', res = 200)
plot(c)
dev.off()


