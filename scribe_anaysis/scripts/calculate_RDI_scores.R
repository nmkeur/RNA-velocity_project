# Load libraries
library(Scribe)
library(monocle)

# Read CDS object and order based on pseudotime
CDS_predicted <- readRDS('data/raw_data/CD4_expr/')
CDS_predicted <- CDS_predicted[, order(pData(CDS_predicted)$Pseudotime)]

# Select genes in the CDS object
BEAM <- c("IFIT3","HLA-A","HLA-E","IFIT1","IFI6","IRF2","XAF1","IRF1","IRF3","STAT1")
CDS_predicted <- CDS_predicted[CDS_predicted@featureData@data$gene_short_name %in% BEAM,]

# Create the network graph we want to estimate RDI (here we calculate all possible pair of gene regulation)
data <- t(exprs(CDS_predicted_velocity)[,])
tmp <- expand.grid(1:ncol(data), 1:ncol(data), stringsAsFactors = T)
super_graph <- tmp[tmp[, 1] != tmp[, 2], ] - 1 #
super_graph <- super_graph[, c(2, 1)]
super_graph[1:2,1:2]

# Calculate RDI scores
# super graph can be used to restrict calculations
# unifomrmization can be enables and delays can be provided
# different combination can be given such as; RDI uRDI cRDI ucRDI
urdi_list_velocity <- calculate_rdi(CDS_predicted,
                             delays = seq(400,8400,400),
                             super_graph = super_graph,
                             method=1,
                             log = FALSE,
                             uniformalize = FALSE) 

# Condition the score of top incomming nodes
ucrdi_list_velocity <- calculate_conditioned_rdi(CDS_predicted_velocity,
                                          rdi_list = urdi_list_velocity,
                                          super_graph = super_graph,
                                          top_incoming_k = 1,
                                          log = FALSE,
                                          uniformalize = FALSE) 

# Save objects
save(urdi_list_velocity, ucrdi_list_velocity, CDS_predicted_velocity, file='~/Desktop/CDS_predicted_velocity_10_genes.RData')
