#This file contains R code for analyzing the NTP dataset using the approaches 
#in the article entitled "Doubly robust estimation and sensitivity analysis with
#outcomes truncated by death in multi-arm clinical trials".

library(geex)
source("SampleCodeFun.R")

#Data pre-processing
load("./ntp.rdata")
df <- ntp
#Reorder the treatment values such that higher level indicates higher survival probability
df$Z <- 5-df$Z 
df$C <- as.factor(df$C)
df$C2 <- as.factor(df$C2)

d234 <- point_est(df,g=2,z=3,zp=4)
d323 <- point_est(df,g=3,z=2,zp=3)
d324 <- point_est(df,g=3,z=2,zp=4)
d334 <- point_est(df,g=3,z=3,zp=4)
d412 <- point_est(df,g=4,z=1,zp=2)
d413 <- point_est(df,g=4,z=1,zp=3)
d414 <- point_est(df,g=4,z=1,zp=4)
d423 <- point_est(df,g=4,z=2,zp=3)
d424 <- point_est(df,g=4,z=2,zp=4)
d434 <- point_est(df,g=4,z=3,zp=4)
final_re <- rbind(d234,d323,d324,d334,d412,d413,d414,d423,d424,d434)
colnames(final_re) <- c('PSW','PSWSE','OR','ORSE','DR','DRSE')
rownames(final_re) <- c(1:nrow(final_re))