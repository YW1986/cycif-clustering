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

## Usage
### Install
Build the docker image by `docker build -t drclust .`

### Use
Run the app by `docker-compose run --rm drclust`

### Paremeters
#### Positional arguments: 

`file_path` points to the name of the segmented cycif pseudo-single-cell data in the `/input` folder.

#### Optional arguments:

`-h` or `--help`: Show this help message and exit.

`-algo`: Clustering algorithm, select from `['KMeans','HDBSCAN']`. KMeans by default.

`-mcs`: Minimal cluster size for HDBSCAN clustering, only for HDBSCAN.

`-nc`: Number of clusters for KMeans clustering.

# Example
#### Data projection on 2D UMAP transformed data.
![alt text](/output/Clustering_on_2D.png)

#### Marker expression in the transformed 2D data.
![alt text](/output/raw%20expr%20on%202D.png)
