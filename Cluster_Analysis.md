---
title: "Cluster Analysis"
author: "Quynh Nhu Nguyen"
date: "2023-04-05"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# **CLUSTER ANALYSIS IN R ASSIGNMENT**

```{r}
# Loaded library 
library(tidyverse)
library(utils)
library(factoextra)
library(NbClust)
library(cluster)
```

Download the BetaMatrix.tsv data from: <https://drive.google.com/file/d/1tOdeLpEzhEcsDPU6Vz_dIQsV0UZy0bz> 0/view?usp=share_link. Base on value of 200 CpG sites, do the following requests:

```{r}
# Specify URL where file is stored
url <- "https://drive.google.com/file/d/1tOdeLpEzhEcsDPU6Vz_dIQsV0UZy0bz0&export=download"

file <- basename(url)

# Specify destination where file should be saved
destfile <- "C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data/COVID_assignment.tar.gz"

# Apply download.file function in R
download.file(url, destfile)

untar(destfile, list = TRUE)

.exdir <- "C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data"
dir.create(.exdir)
tarfile <- file.path(.exdir, "COVID_assignment.tar")

untar(tarfile, list = TRUE)

# Extract the two TSV files from the tar archive
tar_files <- tar(tarfile, files = c("SampleSheet.tsv.tsv", "BetaMatrix.tsv"))

# Read the first TSV file into R as a data frame
file1 <- read.table(tar_files[1], header = TRUE, sep = "\t")

# Read the second TSV file into R as a data frame
file2 <- read.table(tar_files[2], header = TRUE, sep = "\t")

```

```{r}
# Load data 
tsv_betamatrix <- read_tsv("C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data/COVID/BetaMatrix.tsv") 

# Convert the tibble into a data frame
data_tsv <- as.data.frame(tsv_betamatrix)

# Check data
head(data_tsv)
dim(data_tsv)
str(data_tsv)
```

```{r}
# Check and remove NA 
colSums(is.na(data_tsv))
data_rm <- na.omit(data_tsv)
head(data_tsv)
```

```{r}
# Standardization
dt_scale <- as.data.frame(scale(data_rm[,2:301]))
dt_scale
```

## 1. Draw dendrogram.

```{r}
# Distance calculation
dist.ma <- dist(dt_scale, method = "euclidean")
dist.ma
```

```{r}
#Perform hierarchical clustering
hc <- hclust(dist.ma, method="complete")
fviz_dend(hc, k=3,
          k_color=c("darkorange", "purple", "cyan3"),
          color_labels_by_k = TRUE,
          type="rectangle",
          cex = 0.15,
          main = "Cluster Dendrogram\nRectangle",
          xlab = "",
          ylab = "Distance")

fviz_dend(hc, k=3,
          k_color=c("darkorange", "purple", "cyan3"),
          type="circular",
          cex = 0.4,
          main = "Cluster Dendrogram\nCircular")

fviz_dend(hc, k=3,
          k_color=c("darkorange", "purple", "cyan3"),
          type="phylogenic",
          cex = 0.5,
          repel = TRUE,
          max.overlaps=Inf,
          main = "Cluster Dendrogram\nPhylogenic")
```

## 2. How many optimal k clusters does each method (Elbow, average silhouette, gap statistic method) tell you?

```{r}
# K estimate

## Elbow method
fviz_nbclust(dt_scale, kmeans, method="wss") + labs(title= "Elbow method")

## Silhouette method
fviz_nbclust(dt_scale, kmeans, method="silhouette") + labs(title= "Average silhouette method")

## Gap statistic method
fviz_nbclust(dt_scale, kmeans, method="gap_stat") + labs(title= "Gap statistic method")
```

According to the three graphs above, the optimal k for each method is as follows: 
- Elbow method: 2 
- Silhouette method: 2 
- Gap statistic method: 7

## 3. Draw clustering results from K-means clustering algorithm.

```{r}
# Visualiztion of k-means clusters (k=2)
km = kmeans(dt_scale, center=2, nstart=20)
fviz_cluster(km,
             dt_scale,
             ellipse.type="convex",
             repel = TRUE,
             max.overlaps=Inf,
             show.clust.cent = TRUE,
             geom = "point", #Show points only
             palette = "Set2", ggtheme = theme_minimal(),## Change the color palette and theme
             main = "Cluster plot with 2 clusters k"
)
```

```{r}
# Visualiztion of k-means clusters (k=7)
km = kmeans(dt_scale, center=7, nstart=20)
fviz_cluster(km,
             dt_scale,
             ellipse.type="convex",
             repel = TRUE,
             max.overlaps=Inf,
             show.clust.cent = TRUE,
             geom = "point", #Show points only
             palette = "Set2", ggtheme = theme_minimal(),## Change the color palette and theme
             main = "Cluster plot with 7 clusters k"
)
```
