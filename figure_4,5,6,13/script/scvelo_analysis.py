import scvelo as scv
import scanpy as sc


# read adata object
adata = sc.read_loom('path_to_file')

# Normalize and calculate the moments
scv.pp.normalize_per_cell(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#Calculate the velocities
scv.tl.velocity(adata)


# Here you could first filter cells.
# So they dont get included in the velocity graph
scv.tl.velocity_graph(adata)

#Command in scvelo and scanpy both work on the same object
# so we can switch between the two modules

# calculate tsne, with optional arguments
sc.tl.tsne(adata)
scv.tl.tsne(adata,perplexity=50, early_exaggeration=100,learning_rate=100,random_state=108, n_jobs=12)

# write tsne to file.
scv.pl.tsne(adata, color='Stimulation', save='path_to_save_figure')
scv.pl.tsne(adata, color='Batch', save='path_to_save_figure')

#calculate the velocity embedding by project velocities on the tsne embedding
scv.tl.velocity_embedding(adata, basis='tsne')
scv.pl.velocity_embedding(adata, basis='tsne', save='path_to_save_figure')

# The monocle embedding can be added to the
# obsm layer and then projected using the follow command
scv.tl.velocity_embedding(adata, basis='name_of_embedding')
scv.pl.tsne(adata, basis='name_of_embedding', color='Stimulation', save='path_to_save_figure')
