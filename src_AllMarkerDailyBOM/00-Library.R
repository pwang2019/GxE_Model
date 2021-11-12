######################################################################################################################################
#Set working directory to Honours.R folder

source_dir = "D:\\Works_2020\\GxE_Models\\GxE_Models_Oct2021\\"

setwd(source_dir)

plot_path = paste(source_dir, "outputs/plots", sep="")

options("experssion" = 5e5)

######################################################################################################################################
##Library
library(BGLR)
library(tidyverse)
library(tidyquant)
library(readxl)
library(data.table)
library(slga)
library(BiocManager)
library(snpReady)
library(caret)
library(factoextra)
library(gdata)
library(geosphere)
library(cowplot)
library(doParallel)
library(randomForestSRC)
library(VIM)

######################################################################################################################################
