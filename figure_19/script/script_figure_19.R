library(monocle)
library(Scribe)

#Read in csv files that have gamma values calculated using velocyto.
gammas_stim <- read.csv('~/Desktop/scribe_results/CD4_stim_gammas.csv', col.names = c("Genes",'Gamma'))
gammas_unstim <- read.csv('~/Desktop/scribe_results/CD4_unstim_gammas.csv', col.names = c("Genes",'Gamma'))

#Select only genes of interest
genes <- c("IFITM1","IFI6","ISG15","IFIT3","IRF1","STAT1","XAF1","IFIT1","RSAD2","IFIT2","MX1","BST2","ISG20","GBP2","PTPN1","HLA-E","HSP90AB1","HLA-A","HLA-B","LSM14A","YTHDF2","JAK1","IFNAR1","HLA-C","IRF2","CDC37","HLA-F","IRF3","ABCE1")
gammas_stim_sub <- gammas_stim[gammas_stim$Genes %in% genes,]
gammas_unstim_sub <-  gammas_unstim[gammas_unstim$Genes %in% genes,]

# Calculate difference in gamma using stimulation - unstimulation
gammas_stim_sub$DGamma <- gammas_stim_sub$Gamma - gammas_unstim_sub$Gamma

# Add ID to dataframe
gammas_stim_sub$ID <- seq.int(nrow(gammas_stim_sub))

# Create plot using difference in gamma's and ID's
png("~/RNA-velocity_project/figure_19/figure/figure_19.png", width = 2400, height = 1500, units='px', res = 200)
ggplot(gammas_stim_sub, aes(x=ID, y=DGamma)) +
  geom_point(size=2, aes(color=Genes) )+ xlab('Genes')+ ylab('Difference Gamma')+ggtitle('')+
  theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label=Genes),hjust=0, vjust=0,size=4)+
  theme(legend.position="none") +theme(text = element_text(size=10))
dev.off()

