###########################SA under violation of PI#############################

library(geex)
library(latex2exp)
source("SampleCodeFunPI.R")

#Data pre-processing
load("./ntp.rdata")
df <- ntp
#Reorder the treatment values such that higher level indicates higher survival probability
df$Z <- 5-df$Z 
df$C <- as.factor(df$C)
df$C2 <- as.factor(df$C2)

######################################fix delta1=1 ########################################
#setup for the coordinates for the plot when fixing delta1 and varying delta2 and delta3
x_seq <- seq(0.51, 2.00, length = 150) #delta2
y_seq <- seq(0.51, 2.00, length = 150) #delta3

par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
#panel PSW delta1=1 
plot_sa_pi_psw12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=2, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_psw12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=3, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_psw13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_psw14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_psw23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=2, zp=3, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_psw23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=2, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_psw24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=3, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_psw34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Principal Score Weighting", side = 3, line = -3, cex = 1.5, outer = TRUE)


#panel OR delta1=1 
par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
plot_sa_pi_or12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=2, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_or12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=3, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_or13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_or14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_or23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=2, zp=3, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_or23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=2, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_or24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=3, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_or34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Outcome Regression", side = 3, line = -3, cex = 1.5, outer = TRUE)


#panel DR delta1=1 
par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
plot_sa_pi_dr12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=2, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_dr12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=3, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_dr13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_dr14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_dr23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=2, zp=3, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_dr23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=2, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_dr24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=3, zp=4, 1, x, y)))
contour(x_seq, y_seq, plot_sa_pi_dr34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_2$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Doubly Robust", side = 3, line = -3, cex = 1.5, outer = TRUE)



######################################fix delta2=1 ########################################
#setup for the coordinates for the plot when fixing delta2 and varying delta1 and delta3
x_seq <- seq(0.51, 2.00, length = 150) #delta1
y_seq <- seq(0.51, 2.00, length = 150) #delta3

par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
#panel PSW delta2=1 
plot_sa_pi_psw12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=2, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_psw12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=3, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_psw13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_psw14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_psw23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=2, zp=3, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_psw23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=2, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_psw24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=3, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_psw34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Principal Score Weighting", side = 3, line = -3, cex = 1.5, outer = TRUE)


#panel OR delta2=1 
par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
plot_sa_pi_or12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=2, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_or12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=3, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_or13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_or14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_or23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=2, zp=3, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_or23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=2, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_or24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=3, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_or34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Outcome Regression", side = 3, line = -3, cex = 1.5, outer = TRUE)


#panel DR delta2=1 
par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
plot_sa_pi_dr12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=2, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_dr12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=3, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_dr13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_dr14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_dr23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=2, zp=3, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_dr23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=2, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_dr24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=3, zp=4, x, 1, y)))
contour(x_seq, y_seq, plot_sa_pi_dr34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_3$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Doubly Robust", side = 3, line = -3, cex = 1.5, outer = TRUE)



######################################fix delta3=1 ########################################
#setup for the coordinates for the plot when fixing delta3 and varying delta1 and delta2
x_seq <- seq(0.51, 2.00, length = 150) #delta1
y_seq <- seq(0.51, 2.00, length = 150) #delta2

par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
#panel PSW delta3=1 
plot_sa_pi_psw12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=2, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_psw12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=3, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_psw13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=1, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_psw14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_psw23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=2, zp=3, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_psw23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=2, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_psw24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_psw34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_psw_sa_pi(df, g=4, z=3, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_psw34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Principal Score Weighting", side = 3, line = -3, cex = 1.5, outer = TRUE)


#panel OR delta3=1 
par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
plot_sa_pi_or12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=2, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_or12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=3, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_or13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=1, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_or14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_or23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=2, zp=3, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_or23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=2, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_or24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_or34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_or_sa_pi(df, g=4, z=3, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_or34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Outcome Regression", side = 3, line = -3, cex = 1.5, outer = TRUE)


#panel DR delta3=1 
par(mar = c(5, 5, 7, 5), mfrow = c(2, 3))
plot_sa_pi_dr12 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=2, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_dr12, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,2)$"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr13 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=3, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_dr13, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr14 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=1, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_dr14, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(1,4)"), side = 3, line = 1, cex = 1.5)


plot_sa_pi_dr23 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=2, zp=3, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_dr23, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,3)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr24 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=2, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_dr24, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(2,4)"), side = 3, line = 1, cex = 1.5)

plot_sa_pi_dr34 <- outer(x_seq, y_seq, Vectorize(function(x, y)point_dr_sa_pi(df, g=4, z=3, zp=4, x, y, 1)))
contour(x_seq, y_seq, plot_sa_pi_dr34, labcex = 1, cex.lab = 1.5, cex.axis = 1.5, col = "blue", lwd = 2)
mtext(TeX("$\\delta_1$"), side = 1, line = 3, cex = 1.5)
mtext(TeX("$\\delta_2$"), side = 2, line = 3, cex = 1.5, las = 1)
mtext(TeX("$\\Delta_4(3,4)"), side = 3, line = 1, cex = 1.5)
mtext("Doubly Robust", side = 3, line = -3, cex = 1.5, outer = TRUE)


