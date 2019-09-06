# import libraries
import numpy as np
import pandas as pd
import scanpy as sc
import loompy
import velocyto as vlm

#Read in all seperate files
adata1 = sc.read_loom('data/raw_data/pilot4_lane_1.loom')
adata2 = sc.read_loom('data/raw_data/pilot4_lane_2.loom')
adata3 = sc.read_loom('data/raw_data/pilot4_lane_3.loom')

adata1.var_names_make_unique()
adata2.var_names_make_unique()
adata3.var_names_make_unique()

# Add batch ID
adata1.obs['batch_name'] = "1"
adata2.obs['batch_name'] = "2"
adata3.obs['batch_name'] = "3"

# Make the names matching between the different files.
adata1.obs.index= adata1.obs.index.to_series().astype(str).str.replace('pilot4_lane_1:', '')
adata1.obs.index= adata1.obs.index.to_series().astype(str).str.replace('x', '_lane1')
adata2.obs.index= adata2.obs.index.to_series().astype(str).str.replace('pilot4_lane_2:', '')
adata2.obs.index= adata2.obs.index.to_series().astype(str).str.replace('x', '_lane2')
adata3.obs.index= adata3.obs.index.to_series().astype(str).str.replace('pilot4_lane_3:', '')
adata3.obs.index= adata3.obs.index.to_series().astype(str).str.replace('x', '_lane3')

# Concatenate the seperate loom files.
adata = adata1.concatenate(adata2, adata3, join='inner', index_unique=None )

# Remove all duplicatedd data, we no longer need.
del adata.var['Accession-1']
del adata.var['Chromosome-1']
del adata.var['End-1']
del adata.var['Start-1']
del adata.var['Strand-1']

del adata.var['Accession-2']
del adata.var['Chromosome-2']
del adata.var['End-2']
del adata.var['Start-2']
del adata.var['Strand-2']

del adata.obs['_X']
del adata.obs['_Y']
del adata.obs['batch']

adata.var.columns = ['Accession', 'Chromosome', 'End', 'Start', 'Strand']

# Write loom file
adata.write_loom("data/raw_data/merged_pilot4.loom")


#Create index in adata file and join meta_data extracted from seurat object
adata.obs.index.name = 'name'
adata.obs = adata.obs.join(meta_data.set_index('name'))

# Select only the data that has no nan values.
adata.obs['sample'] = adata.obs['sample'].astype(str)
adata2 = adata[adata.obs['sample'] != 'nan']

#adata2= adata.obs.dropna()

#Check the number of dimensions in the matrix
adata.obs.shape # 15901, 7
adata2.obs.shape # (15085,7)
# adata[adata.obs == adata.obs.dropna()] #(15085,7)


# This can be used to filter genes from the anndata object


# genename_data = pd.read_csv('data/raw_data/genenames.csv')
# genename_data.keys()

#genlist = genename_data['genenames'].tolist()
#adata2.var.index = adata.var['Accession']
#adata3 = adata2[:,adata2.var.index.isin(genlist)]
#adata3.var.index = adata3.var['Gene']

#Write new loom file
adata3.write_loom("data/raw_data/filtered_scaled_pilot4v03.loom")
