###########################SA under violation of monotonicity#############################

library(geex)
library(parallel)
source("SampleCodeFunMono.R")

#Data pre-processing
load("./ntp.rdata")
df <- ntp
#Reorder the treatment values such that higher level indicates higher survival probability
df$Z <- 5-df$Z 
df$C <- as.factor(df$C)
df$C2 <- as.factor(df$C2)

rho <- seq(0, 0.52, by = 0.02)
rho1011 <- rho
rho0101 <- rho
rho0010 <- rho
d412 <- d413 <- d414 <- d423 <- d424 <- d434 <-  c()

# Define a function that performs the computation for a given index i
compute_point_est <- function(i) {
  point_est(df, g = 4, z = 1, zp = 2, rho[i], rho[i], rho[i])
}

# Use mclapply to perform the computations in parallel
results <- mclapply(1:length(rho1011), function(i) {
  list(
    d412 = compute_point_est(i),
    d413 = compute_point_est(i),
    d414 = compute_point_est(i),
    d423 = compute_point_est(i),
    d424 = compute_point_est(i),
    d434 = compute_point_est(i)
  )
}, mc.cores = 6)

# Initialize the result data frames
d412 <- d413 <- d414 <- d423 <- d424 <- d434 <- data.frame()

# Combine the results
for (res in results) {
  d412 <- rbind(d412, res$d412)
  d413 <- rbind(d413, res$d413)
  d414 <- rbind(d414, res$d414)
  d423 <- rbind(d423, res$d423)
  d424 <- rbind(d424, res$d424)
  d434 <- rbind(d434, res$d434)
}

# Optionally set column names
colnames(d412) <- colnames(d413) <- colnames(d414) <- colnames(d423) <- colnames(d424) <- colnames(d434) <- c("rho", "rho", "rho", "PSW", "PSWSE", "OR", "ORSE", "DR", "DRSE")

d412 <- d412[,-c(1,2)]
d413 <- d413[,-c(1,2)]
d414 <- d414[,-c(1,2)]
d423 <- d423[,-c(1,2)]
d424 <- d424[,-c(1,2)]
d434 <- d434[,-c(1,2)]