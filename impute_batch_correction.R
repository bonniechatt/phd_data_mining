library(missMDA)
library(sva)

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

