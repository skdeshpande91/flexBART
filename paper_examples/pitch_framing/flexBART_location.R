# flexBART but only using location
load("data/sim_settings.RData")


############
# Change for different folds and seasons!
fold_num <- 1
year <- 2013
############

##############
load(paste0("data/bart_data_", year, ".RData"))


train_index <- folds[[fold_num]][["train_index"]]
test_index <- folds[[fold_num]][["test_index"]]

X_cont_train <- X_cont_all[train_index,]
X_cont_test <- X_cont_all[test_index,]

X_cat_train <- X_cat_all[train_index,]
X_cat_test <- X_cat_all[test_index,]

Y_train <- Y_all[train_index]
Y_test <- Y_all[test_index]

#########
# To use only location, suffices to not pass X_cat's
timing <- 
  system.time(
    fit <- flexBART::probit_flexBART(Y_train = Y_train,
                                     X_cont_train = X_cont_train,
                                     X_cont_test = X_cont_test,
                                     unif_cuts = unif_cuts, cutpoints_list = cutpoints_list,
                                     sparse = FALSE, verbose = FALSE,
                                     save_samples = FALSE, save_trees = FALSE))

misclass <- c(train = compute_misclass(fit$prob.train.mean, Y_train),
              test = compute_misclass(fit$prob.test.mean, Y_test))
brier <- c(train = compute_brier(fit$prob.train.mean, Y_train),
           test = compute_brier(fit$prob.test.mean, Y_test))
logloss <- c(train = compute_logloss(fit$prob.train.mean, Y_train),
             test = compute_logloss(fit$prob.test.mean, Y_test))

assign(paste0("flexBART_location_", year, "_", fold_num),
       list(year = year, fold_num = fold_num,
            timing = timing,
            misclass = misclass,
            brier = brier,
            logloss = logloss))
save(list = paste0("flexBART_location_", year, "_", fold_num),
     file = paste0("flexBART_location_", year, "_", fold_num, ".RData"))
