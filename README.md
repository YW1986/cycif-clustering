# cycif-clustering
Cell state calling docker app for cycif data
# Dockerized app for cell clustering and top marker identification from Cycif data. 
## Preprocessing and normalization
log2 transformation => robust quantile normalization

## Dimension reduction
UMAP (Uniform Manifold Approximation and Projection)
https://github.com/lmcinnes/umap

## Clustering and differential expression analysis
KMeans and HDBSCAN
Differential expression based on fold change against overall population (self vs non-self).
