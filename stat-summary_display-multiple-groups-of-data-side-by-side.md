---
title: "Reshape and plot"
author: "Quynh Nhu Nguyen"
date: "2023-04-08"
output: html_document
---

# **RESHAPE AND PLOT**

```{r}
# Load library 
library(tidyverse)
library(reshape) # reshape a data frame to 'wide' format
```

```{r setup, include=FALSE}
# Set the root directory for notebook chunks
knitr::opts_knit$set(root.dir = "C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/exam_data/exam_data")
```

## 1. Read file "SNPs_before_filter.csv" into R

```{r}
# Load data 
csv_file <- read_csv("C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/exam_data/exam_data/SNPs_before_filter.csv") 

# Convert the tibble into a data frame
data_csv <- as.data.frame(csv_file)

# Check data
head(data_csv)
dim(data_csv)
str(data_csv)
```

```{r}
# Check and remove NA 
colSums(is.na(data_csv))
data_rm <- na.omit(data_csv)
head(data_csv)
```

```{r}
# Rename the first column
names(data_csv)[1] <- "Cases"
```

```{r}
# Reshape the data frame into a "long" format with separate columns for the tool and its values
data_csv_long <- data_csv %>%
  pivot_longer(cols=c(FreeBayes, Mutect2, Strelka), names_to="Tools", values_to="SNPs") %>%
  na.omit()
print(data_csv_long)

# Check summary statistics of the data
summary(data_csv_long)
```

## 2. Compute and plot the mean and median of number of SNP of FreeBayes, Mutect2, Strelka

### Compute the mean and median of number of SNP of FreeBayes, Mutect2, Strelka

#### Method 1:

```{r}
for (i in 2:4) {
  mean_snps <- round(mean(data_csv[,i]),3)
  median_snps <- round(median(data_csv[,i]),3)
  # Get the name of each columns
  column_name<- colnames(data_csv)[i]
  # Note that the print() function is usually used to print objects in R, not individual values.
  print(paste("The mean", mean_snps, "median", median_snps, "of number of SNP of", column_name))
}
```

#### Method 2:

```{r}
mean <- apply(data_csv[-1], 2, mean)
print("The mean are")
mean

median <- apply(data_csv[-1], 2, median)
print("The median are")
median
```

`data_csv[-1]` is a command in R that returns a new data frame that is identical to the original `data_csv` data frame, except that the first column has been removed. The `[-1]` notation is used to remove the first column, as the minus sign indicates that the first column should be excluded.

In R, data frames are represented as two-dimensional objects, with rows and columns. The `[i, j]` notation can be used to select individual elements of a data frame, where i is the row index and j is the column index. The `[, j]` notation can be used to select an entire column of a data frame, where j is the column index. By using the `[-1]` notation instead of the `[1]` notation, we are indicating that we want to exclude the first column rather than include it.

#### Method 3:

```{r}
data_summary <- data_csv_long %>%
  group_by(Tools) %>%
  summarize(mean = mean(SNPs), median = median(SNPs)) %>%
  as.data.frame()
data_summary
```

### Plot the mean and median of number of SNP of FreeBayes, Mutect2, Strelka

```{r}
# Reshape a data frame to 'wide' format
data_melt <- melt(data_summary)
data_melt
```

```{r}
# Create a boxplot with the mean and median values for each tools
data_melt %>% ggplot(aes(x=Tools, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Mean and Median number of SNPs identified by different Tools", x="Tools", y="Number of SNPs") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(25000, 100000, 200000, 700000)) +
  theme_classic()
```

## 3. Plot barplot with X axis is Cases and Y axis is number of SNPs of FreeBayes, Mutect2, Strelka

```{r}
# Create a barplot with X axis is Cases and Y axis is number of SNP
data_csv_long %>% 
  mutate(Number = as.numeric(gsub("Case", "", Cases))) %>%
  mutate(Cases = reorder(Cases, Number)) %>% # reorder
  ggplot(aes(x=Cases, y=SNPs, fill=Tools)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Number of SNPs identified by different Tools", x="Cases", y="Number of SNPs") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_classic()
```
In this example, we are creating a bar chart using the `SNPs_before_filter` dataset included in the `ggplot2` package. We are mapping the `Cases` variable to the x-axis, the `SNPs` variable to the y-axis, and the `Tools` variable to the fill. We are using `geom_bar` to create the bars, and setting the `stat` parameter to `"identity"` to plot the actual values rather than a summary statistic. Finally, we are setting the position parameter to `position_dodge()` to place the bars side by side.

In `ggplot2`, `position_dodge` is a positioning option that is commonly used to **display multiple groups of data side by side**, particularly in bar charts. When `position_dodge` is applied, the bars for each group are placed adjacent to each other, rather than overlapping.

## 4. Plot boxplot with X axis is FreeBayes, Mutect2, Strelka and Y axis is number of SNP

```{r}
# Create a boxplot for each tools
data_csv_long %>%
  mutate(Tools = reorder(Tools, SNPs, fun = "median")) %>% # reorder
  ggplot(aes(x=Tools, y=SNPs)) +
  geom_boxplot(fill="lightblue", alpha=0.5, na.rm = TRUE) +
  stat_summary(fun.y="mean", geom="point", shape=18, size=4, color="blue", na.rm = TRUE) +
  stat_summary(fun.y="median", geom="point", shape=18, size=4, color="red", na.rm = TRUE) +
  scale_y_continuous(breaks = c(25000, 100000, 200000, 700000)) +
  labs(title="Number of SNPs identified by different Tools", x="Tools", y="Number of SNPs") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_bw()
```
