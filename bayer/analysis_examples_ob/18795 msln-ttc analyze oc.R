source('bayer/analysis_examples_ob/func_ana_simres.r')
#library(crmPack)
library(openxlsx)

res <- vector(mode="list",length=15)
sim_names_ob <- paste(c('safe.', 'late.', 'early.', 'toxic.', 'peak.'), 18795, sep = "")

truth_names <- c('safe', 'late', 'early', 'toxic', 'peak')

truth_all <- mget(truth_names)
target <- 0.3
true.mtd <- c(15, 4.5, 2.5, 1, 3.5)
sim_all_ob <- mget(sim_names_ob)

for (i in 1:5){
  res[[i]] <-simsum(sim_all_ob[[i]], i, truth_all[[i]], target, true.mtd[i])
}

r <- do.call("cbind", res)
colnames(r) <- c(truth_names)
r
