---
title: "ggplot2 in R Assignment"
author: "Quynh Nhu Nguyen"
date: "2023-03-26"
output: html_document
---
# **GGPLOT2 IN R ASSIGNMENT**

```{r}
# Loaded library 
library(tidyverse)
library(data.table)
library(factoextra) # Package for calculation PCA
library(ComplexHeatmap)
library(circlize) # Package of Color for heatmap visualization
library(reshape) # reshapes a data frame to 'wide' format
```

## 1. DATA "DNA METHYLATION"

### Prepare data

#### Method 1:

```{r}
# Set the working directory to the folder containing the TSV files
setwd("C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data/Basic R5_ ggplot2 package/DNAmethylation")
getwd()

# Load data 
tsv_files <- list.files(pattern = "\\.tsv$") # List all the TSV files in the folder
combined_data <- data.frame() # Create an empty data frame to store the combined data

# Loop through each TSV file, read it in as a data frame, and add it to the combined data
for (file in tsv_files) {
  data <- read_tsv(file)
  data_combine <- bind_rows(combined_data, data)
}

write_tsv(data_combine, "combined_data.tsv") # Write the combined data to a new TSV file
```

#### Method 2:

```{r}
# List all the TSV files in the folder
filename1 <- list.files("C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data/Basic R5_ ggplot2 package/DNAmethylation", pattern = '\\.tsv$', full.names = TRUE)
# Read all files
data_f <- lapply(filename1, function(x) read.table(x,header=T))
# merge all files using do.call
data_combine <- do.call("rbind", data_f )
```

```{r}
# Check data
data_combine <- data_combine[, c(14, 1:13)] # Move the "subtype" column to 1st column
head(data_combine)
dim(data_combine)
str(data_combine)
```

```{r}
# Check and remove NA 
colSums(is.na(data_combine))
data_rm <- na.omit(data_combine)
head(data_rm)
```

### 1.1 Plot PCA following subtypes.

```{r}
# Numeric data 
number_data <- data_rm[,2:14]
head(number_data)
str(number_data)
```

```{r}
### Normalise data 
colMeans(number_data)
data_normalized <- scale(number_data)
colMeans(data_normalized)
```

```{r}
# Caculate PCA 
Xpca <- prcomp(data_normalized)
fviz_eig(Xpca,addlabels=TRUE, ncp=20)
fviz_cos2(Xpca, choice = "var", axes = 1:2)
```

```{r}
# Draw PCA 
data_PCA <- as.data.frame(Xpca$x)
head(data_PCA)
data_PCA$subtype <- data_rm$subtype
ggplot(data_PCA, aes(x = PC1, y = PC2, col=subtype)) + geom_point(size=4, alpha=2) +
  theme(legend.position="bottom") + theme(text = element_text(size = 15))+
  labs(x="PC1 (35.9%) ", y="PC1 (17.6%)", title="DNA methylation")
```

### 1.2 Plot heatmap between four subtypes FA, FC, NIFTP and fv.

```{r}
# Convert the "subtype" column to a factor
data_rm$subtype <- as.factor(data_rm$subtype)
```

The `aggregate()` function is used to calculate the mean of each column in the data frame for each level of the "x" factor variable. The `by` argument is set to `List(df[, 1])` to indicate that the means should be calculated for each level of the factor variable, which is the first column of the data frame. The `FUN` argument is set to mean to indicate that the mean should be calculated for each column. Note that `aggregate()` only works with numeric columns, so any non-numeric columns will be ignored. If your data frame contains missing values (NAs), you may want to use the `na.rm` argument of `mean()` to exclude missing values from the calculation.

```{r}
# Calculate the mean of each column for each level of the factor variable
means <- aggregate(data_rm[, 2:14], by = list(data_rm[, 1]), FUN = mean)
head(means)
str(means)
```

The `t()` function is used to transpose the data frame, so that the rows become columns and the columns become rows. The `[, -1]` indexing syntax is used to exclude the first column, which is assumed to be the row labels.

```{r}
# Transpose the data frame from rows become columns
means1 <- t(means[, -1])
colnames(means1) <- means[, 1] # Names of "4 subtypes" are names of columns in new data-frame
```

```{r}
# Convert a data frame to a matrix
class(means1)
mat <- as.matrix(means1)
class(mat)
```

```{r}
# Draw basic Heatmap
Heatmap(mat)
```

```{r}
# Create annotation color 
col_fun = colorRamp2(c(0, 0.4, 0.8), c("purple", "white", "red"))

# create annotation for column
column_ha = HeatmapAnnotation(
    dist1 = anno_barplot(
    colSums(mat), 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "#2a0134"), 
    border = FALSE,
    axis_param = list(at = c(0, 2, 4, 6, 8), labels = c("0", "2", "4", "6", "8")),
    height = unit(2, "cm")), 
    show_annotation_name = FALSE)

# create annotation for row 
row_ha = rowAnnotation(
     dist2 = anno_barplot(
    rowSums(mat), 
    bar_width = 1, 
    gp = gpar(col = "white", fill = "#2a0134"), 
    border = FALSE,
    axis_param = list(at = c(0, 1, 2, 3), labels = c("0", "1", "2", "3")),
    width = unit(2, "cm")), 
    show_annotation_name = FALSE)

# deal with column variables 
sub_type = colnames(mat)
ha_column = HeatmapAnnotation(
  subtype = anno_text(sub_type, rot = 0, location = unit(1, "npc"), just = "top"))

# Draw plot
ht_list = Heatmap(mat, 
        name = "values", # Legend
        col = col_fun, # Generate the color mapping function
        cluster_columns = FALSE, #  Turn off row clustering
        show_row_dend = FALSE, # Hide column dendrogram
        rect_gp = gpar(col = "white", lwd = 0.1), # Set cell borders
        show_column_names = FALSE,
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 10),
        column_title = 'DNA methylation between four subtype FA, FC, NIFTP and fvPTC',
        column_title_gp = gpar(fontsize = 11, fontface = "bold", hjust = 0.5),
        top_annotation = column_ha, 
        bottom_annotation = ha_column, 
        heatmap_legend_param = list(at = c(0, 0.2, 0.4, 0.6, 0.8), 
                  labels = c("0", "0.2", "0.4", "0.6", "0.8"))) + row_ha
draw(ht_list, ht_gap = unit(1, "mm"))

```

## 2. DATA "STUDENTS' PERFORMANCE"

### Prepare data.

```{r}
# Load data 
csv_files <- read_csv("C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/NGS3-W9-19-03-2023/Data/Basic R5_ ggplot2 package/StudentsPerformance.csv") 

# Convert the tibble into a data frame
data_csv <- as.data.frame(csv_files)

# Check data
head(data_csv)
dim(data_csv)
str(data_csv)
```

```{r}
# Check and remove NA 
colSums(is.na(data_csv)) 
data_rm_data2 <- na.omit(data_csv)
head(data_csv)
```

### 2.1 Plot PCA.

```{r}
# Numeric data 
number_data2 <- data_rm_data2[,6:8] 
head(number_data2) 
str(number_data2)
```

```{r}
# Caculate PCA 
Xpca_data2 <- prcomp(number_data2) 
fviz_eig(Xpca_data2, addlabels=TRUE, ncp=20)
fviz_cos2(Xpca_data2, choice = "var", axes = 1:2)
```

```{r}
# Draw PCA 
data_PCA_data2 <- as.data.frame(Xpca_data2$x) 
head(data_PCA_data2) 
```

```{r}
## Plot PCA following gender:
data_PCA_data2$gender <- data_rm_data2$gender
ggplot(data_PCA_data2, aes(x = PC1, y = PC2, col=gender)) + geom_point(size=4, alpha=2) + theme(legend.position="bottom") + theme(text = element_text(size = 15))+ labs(x="PC1 (90.5%)", y="PC1 (8%)", title="Students Performance in Exams following gender") +
  theme(plot.title = element_text(face="bold",size=12, hjust =0.5))
```

```{r}
## Plot PCA following race/ethnicity:
data_PCA_data2$race <- data_rm_data2$`race/ethnicity`
ggplot(data_PCA_data2, aes(x = PC1, y = PC2, col=race)) + geom_point(size=4, alpha=2) + theme(legend.position="bottom") + theme(text = element_text(size = 15))+ labs(x="PC1 (90.5%)", y="PC1 (8%)", title="Students Performance in Exams following race/ethnicity") +
  theme(plot.title = element_text(face="bold",size=12, hjust =0.5, colour = "purple"))
```

### 2.2 Plot heat.

```{r}
# Plot heat following gender
data_melt <- melt(data_rm_data2)                                           
head(data_melt)
ggp <- ggplot(data_melt, aes(gender, variable)) +# Create heatmap with ggplot2
  geom_tile(aes(fill = value)) +
  labs(title="Students Performance in Exams following gender", 
       y="Coures", 
       x="Gender",
       fill="Score") + # Name of legend
  theme_set(theme_classic()) + 
  theme(plot.title = element_text(face="bold",size=13, hjust =0.5, colour = "purple"))

ggp + scale_fill_gradient(low = "violet", high = "black")
```

```{r}
# Plot heat following race/ethnicity
ggp <- ggplot(data_melt, aes(`race/ethnicity`, variable)) +# Create heatmap with ggplot2
  geom_tile(aes(fill = value)) +
  labs(title="Students Performance in Exams following race/ethnicity", 
       y="Coures", 
       x="Race / Ethnicity",
       fill="Score") + # Name of legend
  theme_set(theme_classic()) +
  theme(plot.title = element_text(face="bold",size=13, hjust =0.5, colour = "darkgreen"))
ggp + scale_fill_gradient(low = "green", high = "black")
```
