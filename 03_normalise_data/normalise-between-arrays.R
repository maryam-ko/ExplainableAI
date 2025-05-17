#!/bin/Rscript

# Set up personal library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

# Setup
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.r-project.org")
BiocManager::install(version = "3.19", ask = FALSE)
BiocManager::install('limma', ask = FALSE)

install.packages('dplyr', repos = "http://cran.r-project.org")
library(limma)
library(dplyr)




# set current working directory
setwd('/data/home/bt24990/maryam-ko-QMUL-MSc-Project')

# load matrix with column names as header (note, this matrix has column 1 as phosphosite_IDs)
matrix <- read.csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/RawMatrix_NoOutliers.csv', header = T)
head(matrix)

# save matrix without phosphosite_ID column (all columns need to be numeric for normalizeBetweenArrays function)
matrix_numeric <- select(matrix, -DatasetName)



# perform normalizeBetweenArrays function
quantile_matrix <- normalizeBetweenArrays(matrix_numeric, method = 'quantile')


# append the phosphosite_ID column to the normalised matrix
quantile_normalised <- cbind(quantile_matrix, matrix$DatasetName)

quantile_normalised <- quantile_normalised[, c(782, 1:781)] 
colnames(quantile_normalised)[colnames(quantile_normalised) == ''] <- 'DatasetName'

write.csv(quantile_normalised, '/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NBA-Matrix_Quantile.csv', row.names = F)
```

