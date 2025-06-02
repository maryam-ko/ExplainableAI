#!/bin/Rscript

# Setup

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.19")

BiocManager::install('limma')
library(limma)

install.packages('dplyr')
library(dplyr)




# set current working directory
setwd('/data/home/bt24990/maryam-ko-QMUL-MSc-Project')

# load matrix with column names as header (note, this matrix has column 1 as phosphosite_IDs)
matrix <- read.csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/RawMatrix_NoOutliers.csv', header = TRUE)
matrix <- as.data.frame(matrix)

phosphosite_IDs <- matrix$DatasetName

# save matrix without phosphosite_ID column (all columns need to be numeric for normalizeBetweenArrays function)
matrix_numeric <- select(matrix, -DatasetName)

# Transpose matrix so samples are rows, features are columns
matrix_numeric_t <- t(matrix_numeric)

# perform normalizeBetweenArrays function
quantile_matrix <- normalizeBetweenArrays(matrix_numeric_t, method = 'quantile')

# Transpose back: phosphosites are rows
quantile_matrix_t <- t(quantile_matrix)

# append the phosphosite_ID column to the normalised matrix
quantile_normalised <- data.frame(DatasetName = phosphosite_IDs, quantile_matrix_t)

write.csv(quantile_normalised, '/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NBA-Matrix_Quantile.csv', row.names = F)
```

