library(missMDA)
library(sva)
library(FactoMineR)
library(factoextra)

#loading datasets
mbc.ht.to.impute <- load(mbc_ht_merged_with_gaps.RData)
mbc.ht.native <- load(mbc_ht_merged_native_only.Rdata)
batch.mbc <- load(batch_for_combat.RData)

#imputing missing values
imputed_data <- imputePCA(mbc.ht.to.impute, method = "EM", ncp = 10)

#batch correcting dataset with imputed values 
imputed.data.combat <- ComBat(imputed_data$completeObs, batch = batch.mbc$batch)

#batch correcting dataset without imputed values for future comparison
native.data.combat <- ComBat(mbc.ht.native, batch = batch.mbc$batch)

#running PCA on imputed dataset
res.pca.imputed <- PCA(imputed.data.combat)

#creating PCA plot on the imputed dataset selected for top 20 genes by contribution
pca.contrib.top20.imputed <- fviz_pca_ind(res.pca.imputed, 
                                  col.ind="contrib", 
                                  select.ind = list(contrib = 20),
                                  legend.title = "Confidence",
                                  repel = TRUE,
                                  title = 'PCA Plot of Imputed Data') +
  scale_color_gradient2(low="forestgreen", mid="goldenrod",
                        high="firebrick4", midpoint=4)
                        
#running PCA on native dataset
res.pca.native <- PCA(native.data.combat)

#creating PCA plot on the native dataset selected for top 20 genes by contribution
pca.contrib.top20.native <- fviz_pca_ind(res.pca.native, 
                                  col.ind="contrib", 
                                  select.ind = list(contrib = 20),
                                  legend.title = "Confidence",
                                  repel = TRUE,
                                  title = 'PCA Plot of Native Data') +
  scale_color_gradient2(low="forestgreen", mid="goldenrod",
                        high="firebrick4", midpoint=4)
