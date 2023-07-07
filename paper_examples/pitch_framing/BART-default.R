# BART simulation that drops constant columns
library(BART)
library(dbarts)

############
# Change for different folds and seasons!
fold_num <- 1
year <- 2013
############

load(paste0("data/bart_data_", year, ".RData"))

train_index <- folds[[fold_num]][["train_index"]]
test_index <- folds[[fold_num]][["test_index"]]

tmp_mm <- dbarts::makeModelMatrixFromDataFrame(X_df_all, drop = FALSE)
Y_train <- Y_all[train_index]
Y_test <- Y-all[test_index]
mm_train <- tmp_mm[train_index,]
mm_test <- tmp_mm[test_index,]

timing <-
  system.time(
    fit <- 
      BART::pbart(x.train = mm_train,
                  y.train = Y_train,
                  x.test = mm_test,
                  sparse = FALSE, rm.const = TRUE,
                  ntree = 200, ndpost = 1000, nskip = 1000, nkeeptreedraws = 0,
                  keepevery = 1, printevery = 2001))
misclass <- c(train = compute_misclass(fit$prob.train.mean, Y_train),
              test = compute_misclass(fit$prob.test.mean, Y_test))
brier <- c(train = compute_brier(fit$prob.train.mean, Y_train),
           test = compute_brier(fit$prob.test.mean, Y_test))
logloss <- c(train = compute_logloss(fit$prob.train.mean, Y_train),
             test = compute_logloss(fit$prob.test.mean, Y_test))


assign(paste0("BART_default_", year, "_", fold_num),
       list(year = year, fold_num = fold_num,
            timing = timing,
            misclass = misclass,
            brier = brier,
            logloss = logloss))
save(list = paste0("BART_default_", year, "_", fold_num),
     file = paste0("BART_default_", year, "_", fold_num, ".RData"))

