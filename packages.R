if(!require("BiocManager",quietly=T)){install.packages("BiocManager")}
if(!require("xcms",  quietly=T)){BiocManager::install("xcms")}
if(!require("CAMERA",quietly=T)){BiocManager::install("CAMERA")}

library(ggplot2)
library(pals)
library(reshape2)
library(pheatmap)
library(randomForest)
library(cluster)
library(umap)
library(CAMERA)
library(data.table)
library(ggplot2)
library(ggraph)
library(igraph)

source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")