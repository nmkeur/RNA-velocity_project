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

load('~/Desktop/CDS_predicted_monocle_10_genes.RData')

ucrdi_list_monocle <- urdi_list_monocle$max_rdi_value
a <- pheatmap::pheatmap(ucrdi_list_monocle, cluster_rows = F, cluster_cols = F, annotation_names_col = T, border_color = NA,fontsize = 12)

# Filter using dn sparsifier and CLR
ucrdi_list_monocle[ucrdi_list_monocle <= 0.00] <- 0
ucrdi_list_monocle[dn_sparsifer(ucrdi_list_monocle) == 1] <- 0
ucrdi_list2 <- clr(ucrdi_list_monocle)

load('~/Desktop/CDS_predicted_monocle_10_genes.RData')
ucrdi_list_monocle <- urdi_list_monocle$max_rdi_value
ucrdi_list_monocle[ucrdi_list2 <= 0.0] <-0

b<- pheatmap::pheatmap(ucrdi_list_monocle, cluster_rows = F, cluster_cols = F, annotation_names_col = T, border_color = NA,fontsize = 12)

#Filter using delay
delay <- urdi_list_monocle$max_rdi_delays
d2 <- as.data.frame(as.table(urdi_list_monocle$max_rdi_delays))

for(i in 1:nrow(d2)) {
  row <- d2[i,c(1,2,3)]
  if(delay[row$Var1,row$Var2] > delay[row$Var2,row$Var1]){
    ucrdi_list_monocle[row$Var2,row$Var1] <- 0}
  else{ delay[row$Var1,row$Var2] <- 0}
}

net_m <- graph_from_adjacency_matrix(ucrdi_list_monocle, mode = 'directed',weighted = T)

# Create directed network graph
set.seed(60)
d <- ggraph(net_m, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name),
                     end_cap = label_rect(node2.name)), 
                 arrow = arrow(length = unit(4, 'mm'))) + 
  geom_node_text(aes(label = name))

# Write figure
png("~/RNA-velocity_project/figure_27/figure/figure_27A.png", width = 1500, height = 1800, units='px', res = 200)
plot(d)
dev.off()






#first_row = plot_grid(a[[4]],b[[4]], labels = c('A'))
#second_row = plot_grid(c, labels = c('B', 'C'), nrow = 1, rel_widths = c(0.5, 0.5))
#gg_all = plot_grid(first_row, second_row, labels=c('', ''), ncol=1)
