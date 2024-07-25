load("sim_settings.RData")
source("targetBART.R")
args <- commandArgs(TRUE)
job_id <- as.numeric(args[1])

year <- sim_settings[job_id, "year"]
method <- sim_settings[job_id, "method"]
sim_number <- sim_settings[job_id, "sim_number"]

load(paste0("data/pitchFraming_", year, ".RData"))
year_ix <- which(2013:2019 == year)
set.seed(718*year_ix + sim_number)
n_all <- length(Y_all)
train_ix <- sort(sample(1:n_all, size = floor(0.9*n_all), replace = FALSE))
test_ix <- sort( (1:n_all)[-train_ix])

Y_train <- Y_all[train_ix]
Y_test <- Y_all[test_ix]

X_cont_train <- X_cont_all[train_ix,]
X_cont_test <- X_cont_all[test_ix,]

X_cat_train <- X_cat_all[train_ix,]
X_cat_test <- X_cat_all[test_ix,]

bart_df_train <- X_df_all[train_ix,]
bart_df_test <- X_df_all[test_ix,]

nd <- 1000
burn <- 1000

p_cont <- ncol(X_cont_all)
p_cat <- ncol(X_cat_all)

if(method == "flexBART"){
  timing <- 
    system.time(
      fit <- 
        flexBART::probit_flexBART(Y_train = Y_train,
                                  X_cont_train = X_cont_train,
                                  X_cat_train = X_cat_train,
                                  X_cont_test = X_cont_test,
                                  X_cat_test = X_cat_test,
                                  unif_cuts = unif_cuts, cutpoints_list = cutpoints_list,
                                  cat_levels_list = cat_levels_list,
                                  sparse = FALSE, verbose = FALSE,
                                  nd = nd, burn = burn, 
                                  save_samples = FALSE, save_trees = FALSE))[["elapsed"]]
} else if(method == "flexBART_loc"){
  timing <- 
    system.time(
      fit <- 
        flexBART::probit_flexBART(Y_train = Y_train,
                                  X_cont_train = X_cont_train,
                                  X_cont_test = X_cont_test,
                                  unif_cuts = unif_cuts, cutpoints_list = cutpoints_list,
                                  cat_levels_list = NULL,
                                  sparse = FALSE, verbose = FALSE,
                                  nd = nd, burn = burn, 
                                  save_samples = FALSE, save_trees = FALSE))[["elapsed"]]
} else if(method == "BART"){
  timing <-
    system.time(
      fit <-
        BART::pbart(x.train = bart_df_train,
                    y.train = Y_train,
                    x.test = bart_df_test,
                    ndpost = nd, nskip = burn, 
                    printevery = nd+burn+1))[["elapsed"]]
  
  
} else if(method == "targetBART"){
  timing <-
    system.time(
      fit <-
        targetBART(Y_train = Y_train,
                   df_train = bart_df_train,
                   df_test = bart_df_test,
                   p_cont = p_cont, p_cat = p_cat,
                   cat_levels_list = cat_levels_list,
                   nd = nd, burn = burn))[["elapsed"]]
} else{
  stop("bad method")
}

misclass <- c(train = compute_misclass(fit$prob.train.mean, Y_train),
              test = compute_misclass(fit$prob.test.mean, Y_test))
brier <- c(train = compute_brier(fit$prob.train.mean, Y_train),
           test = compute_brier(fit$prob.test.mean, Y_test))
logloss <- c(train = compute_logloss(fit$prob.train.mean, Y_train),
             test = compute_logloss(fit$prob.test.mean, Y_test))
assign(paste0("pitchFraming_", job_id),
       list(job_id = job_id,
            year = year,
            sim_number = sim_number,
            misclass = misclass,
            brier = brier,
            logloss = logloss,
            timing = timing))
save(list = paste0("pitchFraming_", job_id),
     file = paste0("pitchFraming_", job_id, ".RData"))
