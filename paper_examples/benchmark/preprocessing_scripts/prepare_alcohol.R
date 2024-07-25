source("preprocessing_scripts/preprocess.R")
orig_data <- read.table("raw_data/ktdata.dat", fill = TRUE)

n <- 2467
# 12 continuous preds
# 6 cat preds (3,3,3,4,4,6)


# race (3 black, white, other)
# marital (4 married, widow, divsep, single)
# employment (3, employed, unemployed, retired?)
# region (4, northest, midwest, south, other??)
# insurance source: (3, medicaid, champus, other)
# Age (6, < 30, 30-40, 40-50, 50-60, 60-70, > 70)

# regmed
# dri


# Continuous: editinc, educ, diab, heartcond, stroke, hlthins
# major daily act, some daily act, A,  regmed, dri,

raw_data <- data.frame(D = rep(NA, times = n),
                       A = rep(NA, times = n),
                       EDITINC = rep(NA, times = n),
                       AGE = rep(NA, times = n),
                       EDUC = rep(NA, times = n),
                       RACE = rep(NA, times = n),
                       MARITAL = rep(NA, times = n),
                       EMPLOY = rep(NA, times = n),
                       REGION = rep(NA, times = n),
                       MEDICARE = rep(NA, times = n),
                       INS_SOURCE = rep(NA, times = n),
                       INSURED = rep(NA, times = n),
                       REGMED = rep(NA, times = n),
                       DRI = rep(NA, times = n),
                       MAJOR_LIM = rep(NA, times = n),
                       SOME_LIM = rep(NA, times = n),
                       DIAB = rep(NA, times = n),
                       HEART = rep(NA, times = n),
                       STROKE = rep(NA, times = n))

for(i in 1:n){
  
  # need to recode lots of values
  start_index <- 1 + 4 * (i-1)
  end_index <- 4*i
  
  raw_data[i,"D"] <- orig_data[start_index,1]
  raw_data[i,"A"] <- orig_data[start_index,2]
  raw_data[i,"EDITINC"] <- orig_data[start_index,3]
  
  

  # figure out age
  tmp_age <- orig_data[start_index, 4:8]
  if(!sum(tmp_age) %in% c(0,1)){
    print(paste("i = ", i))
    print(tmp_age)
    stop("two_age_categories")
  } else{
    if(all(tmp_age == 0)){
      raw_data[i,"AGE"] <- "age_25" # less than 30
    } else{
      if(tmp_age[1] == 1) raw_data[i,"AGE"] <- "age_35"# 30-40
      if(tmp_age[2] == 1) raw_data[i, "AGE"] <- "age_45" # 40-50
      if(tmp_age[3] == 1) raw_data[i, "AGE"] <- "age_55" # 50-60
      if(tmp_age[4] == 1) raw_data[i, "AGE"] <- "age_65" # 60-70
      if(tmp_age[5] == 1) raw_data[i, "AGE"] <- "age_75" # over 70
    }
  }
  
  # on line start_index+1 we have
  # EDUC, 2 entries for race, 3 entries for marital, 2 for employ/unemploy
  raw_data[i, "EDUC"] <- orig_data[1 + start_index,1]
  
  tmp_race <- orig_data[1 + start_index,2:3]
  if(!sum(tmp_race) %in% c(0,1)){
    print(paste("i = ", i))
    print(tmp_age)
    stop("Invalid coding of race")
  } else{
    if(all(tmp_race == 0)) raw_data[i,"RACE"] <- "white"
    else{
      if(tmp_race[1] == 1) raw_data[i, "RACE"] <- "black"
      if(tmp_race[2] == 1) raw_data[i, "RACE"] <- "other"
    }
  }
  
  tmp_marital <- orig_data[1 + start_index,4:6]
  if(!sum(tmp_marital) %in% c(0,1)){
    print(paste("i = ", i))
    print(tmp_marital)
    stop("Invalid coding of marital status")
  } else{
    if(all(tmp_marital == 0)) raw_data[i, "MARITAL"] <- "single"
    if(tmp_marital[1] == 1) raw_data[i, "MARITAL"] <- "married"
    if(tmp_marital[2] == 1) raw_data[i, "MARITAL"] <- "widowed"
    if(tmp_marital[3] == 1) raw_data[i,"MARITAL"] <- "divorced"
  }
  
  tmp_employ <- orig_data[1 + start_index, 7:8]
  if(!sum(tmp_employ) %in% c(0,1)){
    print(paste("i = ", i))
    print(tmp_employ)
    stop("Invalid coding of employment")
  } else{
    if(all(tmp_employ == 0)) raw_data[i,"EMPLOY"] <- "other_employ"
    if(tmp_employ[1] == 1) raw_data[i,"EMPLOY"] <- "employed"
    if(tmp_employ[2] == 1) raw_data[i, "EMPLOY"] <- "unemployed"
  }
  # line start_index+2:
  # 3 entries for region, medicare, medicaid, champus, htlhins, regmed
  tmp_region <- orig_data[2+start_index,1:3]
  if(!sum(tmp_region) %in% c(0,1)){
    print(paste("i = ", i))
    print(tmp_region)
    stop("Invalid coding of region")
  } else{
    if(all(tmp_region == 0)) raw_data[i,"REGION"] <- "other"
    if(tmp_region[1] == 1) raw_data[i,"REGION"] <- "Northeast"
    if(tmp_region[2] == 1) raw_data[i,"REGION"] <- "Midwest"
    if(tmp_region[3] == 1) raw_data[i,"REGION"] <- "South"
  }
  raw_data[i, "MEDICARE"] <- orig_data[2 + start_index,4]
  
  tmp_ins <- orig_data[2 + start_index,5:6]
  if(!sum(tmp_ins) %in% c(0,1)){
    print(paste("i = ", i))
    print(tmp_ins)
    if(tmp_ins[1] == 1 & tmp_ins[2] == 1) raw_data[i,"INS_SOURCE"] <- "Military"
    #stop("Invalid coding of insurance source")
  } else{
    if(all(tmp_ins == 0)) raw_data[i,"INS_SOURCE"] <- "private"
    if(tmp_ins[1] == 1) raw_data[i,"INS_SOURCE"] <- "Medicaid"
    if(tmp_ins[2] == 1) raw_data[i,"INS_SOURCE"] <- "Military"
  }
  raw_data[i,"INSURED"] <- orig_data[2+start_index,7]
  raw_data[i,"REGMED"] <- orig_data[2+start_index,8]
  
  
  # line start_index + 3: 
  # dri, majorlim, somelim, diab, heart, stroke
  raw_data[i,"DRI"] <- orig_data[3 + start_index,1]
  raw_data[i,"MAJOR_LIM"] <- orig_data[3 + start_index,2]
  raw_data[i,"SOME_LIM"] <- orig_data[3 + start_index,3]
  raw_data[i,"DIAB"] <- orig_data[3 + start_index,4]
  raw_data[i,"HEART"] <- orig_data[3 + start_index,5]
  raw_data[i,"STROKE"] <- orig_data[3 + start_index,6]
}

raw_data[,"AGE"] <- factor(raw_data[,"AGE"])
raw_data[,"RACE"] <- factor(raw_data[,"RACE"])
raw_data[,"MARITAL"] <- factor(raw_data[,"MARITAL"])
raw_data[,"EMPLOY"] <- factor(raw_data[,"EMPLOY"])
raw_data[,"REGION"] <- factor(raw_data[,"REGION"])
raw_data[,"INS_SOURCE"] <- factor(raw_data[,"INS_SOURCE"])


out_name <- c("D")
cont_names <- c("A", "EDITINC", "EDUC", "MEDICARE", "INSURED", "REGMED",
                "DRI", "MAJOR_LIM", "SOME_LIM", "DIAB", "HEART", "STROKE")
cat_names <- c("AGE", "RACE", "MARITAL", "EMPLOY", "REGION", "INS_SOURCE")
use_cut_names <- c("A", "EDUC", "MEDICARE", "INSURED", "REGMED", "DRI", "MAJOR_LIM",
                   "SOME_LIM", "DIAB", "HEART", "STROKE")

Alcohol_data <- preprocess(raw_data, cont_names, cat_names, out_name, use_cut_names, log_y = FALSE)

save(Alcohol_data, file = "data/Alcohol.RData")
