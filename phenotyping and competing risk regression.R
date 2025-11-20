rm(list = ls())
setwd("...\\duke")
library(SNFtool)
# 140 IBSI features
#tumor core
dk_ce_r1<-read.csv("feature_data_cluster_mask1.csv",header = TRUE,stringsAsFactors = FALSE)
# leading edge
dk_ce_r2<-read.csv("feature_data_cluster_mask2.csv",header = TRUE,stringsAsFactors = FALSE)
# normalization
normalized_df1 <- standardNormalization(as.matrix(df1))
normalized_df2 <- standardNormalization(as.matrix(df2))
# calculating similarity matrix
W1 <- affinityMatrix(SNFtool::dist2(as.matrix(normalized_df1), normalized_df1), K = 30, sigma = 0.8)
W2 <- affinityMatrix(SNFtool::dist2(as.matrix(normalized_df2), normalized_df2), K = 30, sigma = 0.8)
# similarity network fusion
fused_similarity_matrix <- SNF(list(W1, W2), K = 30, t = 20)
# spectral clustering
clusters <- spectralClustering(fused_similarity_matrix, K = 3)
# print
table(clusters)

# Fine-Gray competing risk regression analyses 
library(tidycmprsk)
data<-read.csv(".../duke_sur_data_stage2_3_complete.csv", header = TRUE,stringsAsFactors = FALSE)
crr_mod <- tidycmprsk::crr(Surv(dfs,dfs_e) ~ phenotype, failcode = 1, data=data) # interest event=1
crr_mod
