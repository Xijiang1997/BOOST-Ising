# BOOST-Potts

BOOST-Potts is a method to detect genes with spatial expression patterns (SE genes) in spatial transcriptomics datasets. It is a Bayesian modeling framework for the analysis of single-gene count data on a large-scale lattice defined by a number of array spots. For each gene, expression counts are clustered into two groups - high-expression level and low-expression level, and spatial pattern is defined by the interaction between these two groups in Potts model. 

# How to use BOOST-Potts functions

We use `MouseOB dataset` (Spatial Transcriptomics assay of a slice of Mouse Olfactory Bulb) as an example. This dataset can be found in data file.

Firstly, we need to load data and functions

```r
load("data/olfactory_bulb_11.Rdata")
source("functions/Boost_HP_function.R")
```

`MouseOB dataset` includes two parts: `count data` and `location data`. In count data, each column is the expression counts for a gene. Location data is the coordinates to indecate which locations of the tissue slice has been sampled.

Before detecting SE genes, we need to filter the dataset, which can remove sample locations and genes with few expression points. 

```r
filter_result <- filter_count(count, loc, min_total = 10, min_percentage = 0.1)
loc_f <- filter_result[[1]]
count_f <- filter_result[[2]]
```
In the above function, `min_total` is the minimum total counts, and locations are selected if the total counts for all genes in this location is not less than it. `min_percentage` is the minimum percentage of non-zero counts for genes. If a gene has so many zero counts that the percentage of non-zero count is less than this threshold, this gene will be removed. 

After filteration, we can run the main detection function `Boost_HP`. 
```r
detect_result <- Boost_HP(count_f, loc_f, normalization = 2, clustermethod = 'Mclust')
```
In this function, we need to determine which normalization method is used. If `normalization = 1`, counts data are devided by the summation of total counts for each location. If  `normalization = 2`, we use the normalization method mentioned in [SpatialDE](https://www.nature.com/articles/nmeth.4636), which is set as default. For clustering method, model-based clustering method is applied. We can choose K-means by setting ` clustermethod = 'Kmeans'`.

To obtain detected SE genes, we can check the adjusted p-values. 

```r
SE_gene <- colnames(detect_result)[which(detect_result$pvalues_neg_ad< 0.05)
```

