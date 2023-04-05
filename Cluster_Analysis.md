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
library(data.table)
library(factoextra) # Package for calculation PCA
library(ComplexHeatmap)
library(circlize) # Package of Color for heatmap visualization
library(reshape) # reshapes a data frame to 'wide' format
```

Download the BetaMatrix.tsv data from: <https://drive.google.com/file/d/1tOdeLpEzhEcsDPU6Vz_dIQsV0UZy0bz> 0/view?usp=share_link Base on value of 200 CpG sites, do the following requests: 
```{r}
# Specify URL where file is stored
url <- "https://drive.google.com/file/d/1tOdeLpEzhEcsDPU6Vz_dIQsV0UZy0bz0/view"

file <- basename(url)

# Specify destination where file should be saved
destfile <- "C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data/COVID_assignment.tar.gz"

# Apply download.file function in R
download.file(url, destfile)

untar(destfile, list = TRUE)
```
