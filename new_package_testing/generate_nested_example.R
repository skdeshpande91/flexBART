###############################
# Example 1: unit of observation is a single student
# students nested in classrooms, classrooms nested in school, school nested in district
# students also nested in blkgroup, blkgroup nested in tract
# students from same blkgroup could go to different schools and can also go to different districts
set.seed(125)
n <- 500

train_data <- 
  data.frame(X1 = runif(n, min = -1, max = 1),
             X2 = sample(seq(-1,1, by = 0.1), size = n, replace = TRUE))

train_data[,"Classroom"] <-  sample(c(0:119,sample(0:119, size = n-120, replace = TRUE)))
train_data[,"School"] <- floor(train_data$Classroom/5) # 5 classrooms per school
train_data[,"District"] <- floor(train_data$School/6) # 6 schools per district

train_data[,"Blkgrp"] <- sample(0:49, size = n, replace = TRUE)
train_data[,"Tract"] <- floor(train_data$Blkgrp/10) # 10 blkgrps per tract
train_data[,"Race"] <- sample(0:2, size = n, replace = TRUE)

train_data$Classroom <- factor(train_data$Classroom, levels = 0:119)
train_data$School <- factor(train_data$School, levels = 0:23)
train_data$District <- factor(train_data$District, levels = 0:3)

train_data$Blkgrp <- factor(train_data$Blkgrp, levels = 0:49)
train_data$Tract <- factor(train_data$Tract, levels = 0:4)
train_data$Race <- factor(train_data$Race, levels = 0:2)

train_data$Y <- rnorm(n, mean = 0, sd = 1)
#cov_ensm <- 
#  matrix(0, nrow = ncol(train_data), ncol = 4,
#         dimnames = list(colnames(train_data), rep(NA, times = 4)))
#cov_ensm[c("Classroom", "School", "District"),1] <- 1
#cov_ensm[c("Blkgrp", "Tract"),2] <- 1
#cov_ensm[c("X1", "X2", "Classroom", "District", "Blkgrp"), 3] <- 1
#cov_ensm[,4] <- 1

#train_data$Y <- rep(0, times = n)
