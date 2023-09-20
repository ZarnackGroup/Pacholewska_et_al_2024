# run Umap

# libs
library(umap)

# files
smoothMM_par1 = readRDS("./data/smoothMM_par1.rds")
smoothMM_par2 = readRDS("./data/smoothMM_par2.rds")
smoothMM_par3 = readRDS("./data/smoothMM_par3.rds")

# umap settings
set.seed(1234)
custom.settings = umap.defaults
custom.settings$n_epochs = 5000
custom.settings$n_components = 2
custom.settings$min_dist = 0.01
custom.settings$n_neighbors = 5

# run umaps
print("UMAP 1")
umapDf_par1 = umap(smoothMM_par1, config = custom.settings)
print("UMAP 2")
umapDf_par2 = umap(smoothMM_par2, config = custom.settings)
print("UMAP 3")
umapDf_par3 = umap(smoothMM_par3, config = custom.settings)

# save results
saveRDS(umapDf_par1, file = "./data/umapDf_par1.rds")
saveRDS(umapDf_par2, file = "./data/umapDf_par2.rds")
saveRDS(umapDf_par3, file = "./data/umapDf_par3.rds")