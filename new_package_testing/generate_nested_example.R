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

train_data$Classroom <-  sample(c(0:119,sample(0:119, size = n-120, replace = TRUE)))
train_data$School <- floor(train_data$Classroom/5) # 5 classrooms per school
train_data$District <- floor(train_data$School/6) # 6 schools per district

train_data$Blkgrp <- sample(0:49, size = n, replace = TRUE)
train_data$Tract <- floor(train_data$Blkgrp/10) # 10 blkgrps per tract
train_data$Race <- sample(0:2, size = n, replace = TRUE)

train_data$Classroom <- factor(train_data$Classroom, levels = 0:119)
train_data$School <- factor(train_data$School, levels = 0:23)
train_data$District <- factor(train_data$District, levels = 0:3)

train_data$Blkgrp <- factor(train_data$Blkgrp, levels = 0:49)
train_data$Tract <- factor(train_data$Tract, levels = 0:4)
train_data$Race <- factor(train_data$Race, levels = 0:2)

train_data$Y <- rnorm(n, mean = 0, sd = 1)

#####################
# Example 2: 
# student nested in classrooms, classrooms nested in school, school nested in district
# student also nested in blkgrp, blkgrp in tract
# # schools also nested in census tracts

set.seed(125)
train_data2 <- 
  data.frame(X1 = runif(n, min = -1, max = 1),
             X2 = sample(seq(-1,1, by = 0.1), size = n, replace = TRUE))

train_data2$Classroom <-  sample(c(0:119,sample(0:119, size = n-120, replace = TRUE)))
train_data2$School <- floor(train_data2$Classroom/5) # 5 classrooms per school
train_data2$District <- floor(train_data2$School/6) # 6 schools per district

train_data2$Tract <- floor(train_data2$School/5) # at most 5 schools per tract
train_data2$Blkgrp <- rep(NA, times = n)
train_data2$Blkgrp[which(train_data2$Tract == 0)] <- 
  sample(0:9, size = sum(train_data2$Tract == 0), replace = TRUE)
train_data2$Blkgrp[which(train_data2$Tract==1)] <-
  sample(10:19, size = sum(train_data2$Tract==1), replace = TRUE)
train_data2$Blkgrp[which(train_data2$Tract==2)] <-
  sample(20:29, size = sum(train_data2$Tract==2), replace = TRUE)
train_data2$Blkgrp[which(train_data2$Tract==3)] <-
  sample(30:39, size = sum(train_data2$Tract==3), replace = TRUE)
train_data2$Blkgrp[which(train_data2$Tract==4)] <-
  sample(40:49, size = sum(train_data2$Tract==4), replace = TRUE)

train_data2$Race <- sample(0:2, size = n, replace = TRUE)
  
train_data2$Classroom <- factor(train_data2$Classroom, levels = 0:119)
train_data2$School <- factor(train_data2$School, levels = 0:23)
train_data2$District <- factor(train_data2$District, levels = 0:3)

train_data2$Blkgrp <- factor(train_data2$Blkgrp, levels = 0:49)
train_data2$Tract <- factor(train_data2$Tract, levels = 0:4)
train_data2$Race <- factor(train_data2$Race, levels = 0:2)
train_data2$Y <- rnorm(n = n, mean = 0, sd = 1)

train_data2 <- train_data2[,colnames(train_data)]
