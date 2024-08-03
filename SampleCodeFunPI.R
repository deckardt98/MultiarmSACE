###########################SA under violation of PI#############################

#Function calculating point estimates if PI is violated
point_est_sa_pi <- function(df, g, z, zp, delta1, delta2, delta3){
  #Input:
  #df: data set analyzed
  #g: principal strata variable
  #z: z, treatment comparison candidate 1  
  #zp: z^\prime, treatment comparison candidate 2
  #delta1, delta2, delta3: sensitivity parameters
  
  delta <- c(delta1,delta2,delta3,1)
  J <- 4
  n <- nrow(df)
  Z <- df$Z
  S <- df$S
  Y <- df$Y
  
  #Estimate principal scores
  prinscore0 <-  rep(0, n)
  ps_m1 <- glm(S ~ A+C, data = df, subset = (Z==1), family = binomial)
  ps_m2 <- glm(S ~ A+C, data = df, subset = (Z==2), family = binomial)
  ps_m3 <- glm(S ~ A+C, data = df, subset = (Z==3), family = binomial)
  ps_m4 <- glm(S ~ A+C, data = df, subset = (Z==4), family = binomial)
  prinscore1 <-  predict(ps_m1, newdata=df, type="response")
  prinscore2 <-  predict(ps_m2, newdata=df, type="response")
  prinscore3 <-  predict(ps_m3, newdata=df, type="response")
  prinscore4 <-  predict(ps_m4, newdata=df, type="response")
  ps_model <- list(ps_m1,ps_m2,ps_m3,ps_m4)
  prinscore5 <-  rep(1, n)
  prinscore <-  list(prinscore0, prinscore1, prinscore2, prinscore3, prinscore4, prinscore5)
  
  #Treatment probabilities
  proscore <-  list(1/4, 1/4, 1/4, 1/4)
  
  #Omega_{z,g}(X)
  omegadenominator <- as.matrix(cbind(delta[1]*(prinscore[[J+1]]-prinscore[[J]]),
                                      delta[2]*(prinscore[[J]]-prinscore[[J-1]]),
                                      delta[3]*(prinscore[[J-1]]-prinscore[[J-2]]),
                                      delta[4]*(prinscore[[J-2]]-prinscore[[J-3]])))
  omegaxz <- delta[g]*prinscore[[z+1]]/rowSums(as.matrix(omegadenominator[,((J+1-z):4)]))
  omegaxzp <- delta[g]*prinscore[[zp+1]]/rowSums(as.matrix(omegadenominator[,((J+1-zp):4)]))
  
  #Estimate psi_s, z
  psi_s0 <-  rep(0, n)
  psi_s1 <-  I(Z==1)*(S-prinscore[[2]])/proscore[[1]]+prinscore[[2]]
  psi_s2 <-  I(Z==2)*(S-prinscore[[3]])/proscore[[2]]+prinscore[[3]]
  psi_s3 <-  I(Z==3)*(S-prinscore[[4]])/proscore[[3]]+prinscore[[4]]
  psi_s4 <-  I(Z==4)*(S-prinscore[[5]])/proscore[[4]]+prinscore[[5]]
  psi_s5 <-  rep(1, n)
  psi_s <-  list(psi_s0, psi_s1, psi_s2, psi_s3, psi_s4, psi_s5)
  
  #Estimate outcome mean
  om1 <- glm(Y~A+C, data=df, subset = (Z==1&S==1), family = gaussian)
  om2 <- glm(Y~A+C, data=df, subset = (Z==2&S==1),family = gaussian)
  om3 <- glm(Y~A+C, data=df, subset = (Z==3&S==1), family = gaussian)
  om4 <- glm(Y~A+C, data=df, subset = (Z==4&S==1), family = gaussian)
  om_model <- list(om1,om2,om3,om4)
  omf <-  predict(om_model[[z]], newdata=df, type="response")
  omfp <-  predict(om_model[[zp]], newdata=df, type="response")
  
  #Calculate principal score weighting estimators
  if (g==J){
    w1g <- (prinscore[[2]])/prinscore[[z+1]]*mean(I(Z==z)*S)/mean(S*I(Z==1))*omegaxz
    w1gp <- (prinscore[[2]])/prinscore[[zp+1]]*mean(I(Z==zp)*S)/mean(S*I(Z==1))*omegaxzp
  } else {
    w1g <- (prinscore[[J-g+2]]-prinscore[[J-g+1]])/prinscore[[z+1]]*mean(I(Z==z)*S)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))*omegaxz
    w1gp <- (prinscore[[J-g+2]]-prinscore[[J-g+1]])/prinscore[[zp+1]]*mean(I(Z==zp)*S)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))*omegaxzp
  }
  est1 <-  mean(I(Z==z)*S*w1g*Y)/mean(I(Z==z)*S)
  est1p <-  mean(I(Z==zp)*S*w1gp*Y)/mean(I(Z==zp)*S)
  #est1c: point estimate for contrast using principal score weighting
  est1c <- est1-est1p

  #Calculate outcome regression estimators
  factor0 <-  prinscore0
  factor1 <-  S*I(Z==1)/proscore[[1]]
  factor2 <-  S*I(Z==2)/proscore[[2]]
  factor3 <-  S*I(Z==3)/proscore[[3]]
  factor4 <-  S*I(Z==4)/proscore[[4]]
  factor5 <-  prinscore5
  factor <-  list(factor0, factor1, factor2, factor3, factor4, factor5)
  if (g==J){
    est2 <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omf*omegaxz)/(mean(S*I(Z==1)))
    est2p <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omfp*omegaxzp)/(mean(S*I(Z==1)))
    #est2c: point estimate for contrast using outcome regression weighting
    est2c <- est2 - est2p
  } else {
    est2 <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omf*omegaxz)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))
    est2p <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omfp*omegaxzp)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))
    est2c <- est2 - est2p
  }
  
  #Calculate doubly robust estimator
  deltapsiprod <- as.matrix(cbind(delta[1]*(psi_s[[J+1]]-psi_s[[J]]),
                                  delta[2]*(psi_s[[J]]-psi_s[[J-1]]),
                                  delta[3]*(psi_s[[J-1]]-psi_s[[J-2]]),
                                  delta[4]*(psi_s[[J-2]]-psi_s[[J-3]])))
  psiysz <- 4*I(Z==z)*(Y*S-omf*prinscore[[z+1]])+omf*prinscore[[z+1]]
  dr <- mean((prinscore[[J-g+2]]-prinscore[[J-g+1]])*omegaxz/prinscore[[z+1]]*(psiysz-omegaxz/delta[g]*omf*rowSums(as.matrix(deltapsiprod[,((J+1-z):4)])))+omegaxz*omf*(psi_s[[J-g+2]]-psi_s[[J-g+1]]))/mean(psi_s[[J-g+2]]-psi_s[[J-g+1]])
  psiyszp <- 4*I(Z==zp)*(Y*S-omfp*prinscore[[zp+1]])+omfp*prinscore[[zp+1]]
  drp <- mean((prinscore[[J-g+2]]-prinscore[[J-g+1]])*omegaxzp/prinscore[[zp+1]]*(psiyszp-omegaxzp/delta[g]*omfp*rowSums(as.matrix(deltapsiprod[,((J+1-zp):4)])))+omegaxzp*omfp*(psi_s[[J-g+2]]-psi_s[[J-g+1]]))/mean(psi_s[[J-g+2]]-psi_s[[J-g+1]])
  #drc: point estimate for contrast using doubly robust weighting
  drc <- dr - drp
  
  return(c(est1c,est2c,drc))
}
point_psw_sa_pi <- function(df, g, z, zp, delta1, delta2, delta3){
  est <- point_est_sa_pi(df, g, z, zp, delta1, delta2, delta3)
  return(est[1])
}
point_or_sa_pi <- function(df, g, z, zp, delta1, delta2, delta3){
  est <- point_est_sa_pi(df, g, z, zp, delta1, delta2, delta3)
  return(est[2])
}
point_dr_sa_pi <- function(df, g, z, zp, delta1, delta2, delta3){
  est <- point_est_sa_pi(df, g, z, zp, delta1, delta2, delta3)
  return(est[3])
}