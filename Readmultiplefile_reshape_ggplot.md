---
title: "Read multiple files, reshape and plot"
author: "Quynh Nhu Nguyen"
date: "2023-04-09"
output: html_document
---
# **READ MULTIPLE FILES, RESHAPE AND PLOT**
```{r}
# Load library 
library(tidyverse)
library(data.table)
library(reshape) # reshape a data frame to 'wide' format
```

```{r setup, include=FALSE}
# Set the root directory for notebook chunks
knitr::opts_knit$set(root.dir = "C:/Users/lucia/Downloads/Lecturer from Dr. Loi/NGS3/exam_data/exam_data")
```

## 1. Read file “tools_overlap_filter_1.csv” to “tools_overlap_filter_11.csv” into R and merge them into a single data frame
```{r}
# Create an empty list to store the data frames
tools_overlap_filter_list <- list()

# Loop over the file names and read them into the list
for (i in 1:11) {
  file_name <- paste0("tools_overlap_filter_", i, ".csv")
  df <- read.csv(file_name)
  
  # Add a case column with the current iteration number
  df$Cases <- paste0("Case_", i)
  
  # Add the data frame to the list
  tools_overlap_filter_list[[i]] <- df
}

# Merge the data frames into a single data frame
tools_overlap_filter_combine <- do.call(rbind, tools_overlap_filter_list)
```
In R, the double square bracket notation `[[` is used to extract or replace elements of a list by their index or name. It is different from the single square bracket notation `[` in that it returns the actual element of the list, rather than a sub-list containing that element.

In the line `tools_overlap_filter_list[[i]] <- df`, the double square brackets are used to set the`i-th` element of the `tools_overlap_filter_list` list to be the data frame `df`. The double brackets are used instead of single brackets because we want to set the actual element of the list to be the data frame, rather than creating a sub-list containing the data frame.

```{r}
# Check data
head(tools_overlap_filter_combine)
dim(tools_overlap_filter_combine)
str(tools_overlap_filter_combine)
```

```{r}
# Check and remove NA 
colSums(is.na(tools_overlap_filter_combine))
data_rm <- na.omit(tools_overlap_filter_combine)
head(tools_overlap_filter_combine)
```

```{r}
# Move the "subtype" column to 1st column
tools_overlap_filter_combine <- tools_overlap_filter_combine[, c(5, 1:4)]

# Reshape the data frame into a "long" format with separate columns for the tool and its values
data_melt_tools_overlap_filter<- reshape2::melt(tools_overlap_filter_combine)

# rename the first column
names(data_melt_tools_overlap_filter)[2] <- "Tools"
names(data_melt_tools_overlap_filter)[3] <- "SNPs"

# View the melt data frame
data_melt_tools_overlap_filter
```
## 2. Compute and plot the mean and median of number of overlapped SNP of Mutect2/FreeBayes, MuTect2/Strelka, FreeBayes/Mutect2/Strelka
```{r}
data_summary_tools_overlap_filter <- data_melt_tools_overlap_filter %>%
  group_by(Tools) %>%
  summarize(mean = mean(SNPs), median = median(SNPs)) %>%
  as.data.frame()
data_summary_tools_overlap_filter
```

```{r}
# Reshape a data frame to 'wide' format
data_melt_summary_tools_overlap_filter <- reshape2::melt(data_summary_tools_overlap_filter)

# View the melt data frame
data_melt_summary_tools_overlap_filter
```
# Create a boxplot with the mean and median values for each tools
```{r}
data_melt_summary_tools_overlap_filter %>% 
  ggplot(aes(x=Tools, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Mean and Median number of SNPs identified by different Tools", x="Tools", y="Number of overlapped SNPs") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(25000, 100000, 200000, 700000)) +
  theme_classic()
```
## 3. Plot barplot with X axis is Cases and Y axis is number of overlapped SNP of FreeBayes/Mutect2/Strelka
```{r}
# Create a barplot with X axis is Cases and Y axis is number of SNP
data_melt_tools_overlap_filter %>%
  filter(Tools == "FreeBayes.Mutect2.Strelka") %>%
  mutate(Number = as.numeric(gsub("Case_", "", Cases))) %>%
  mutate(Cases = reorder(Cases, Number)) %>% # reorder
  ggplot(aes(x=Cases, y=SNPs, fill=Cases)) +
  geom_bar(stat="identity") +
  labs(title="Number of SNPs identified by FreeBayes.Mutect2.Strelka", x="Cases", y="Number of SNPs") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_classic()
```
The `stat` parameter specifies the statistical transformation to apply to the data before plotting. In the case of `geom_bar()`, `stat="identity`" means that the height of each bar in the chart represents the actual count or value of the data, rather than any statistical summary such as the mean or median.

## 4. Plot boxplot with X axis is Mutect2/FreeBayes, Mutect2/Strelka,  FreeBayes/Mutect2/Strelka and Y axis is number of overlapped SNP
```{r}
# Create a boxplot for each tools
data_melt_tools_overlap_filter %>%
  filter(Tools != "FreeBayes.Strelka") %>%
  mutate(Tools = reorder(Tools, SNPs, fun = "median")) %>% # reorder
  ggplot(aes(x=Tools, y=SNPs, fill=Tools)) +
  geom_boxplot(alpha=1, na.rm = TRUE) +
  labs(title="Number of overlapped SNPs identified by different Tools", x="Tools", y="Number of overlapped SNPs") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_bw()
```
