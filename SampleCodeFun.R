#This file contains R functions for analyzing the NTP dataset using the approaches 
#in the article entitled "Doubly robust estimation and sensitivity analysis with
#outcomes truncated by death in multi-arm clinical trials".

#Please consult R documentations for geex package for more information and examples 
#for defining the estimating equations

#Joint estimating equations for the principal score weighting estimator
#using simple proportions for marginal principal score
eq_psw <- function(data, models, g, z, zp, J){
  #Input:
  #data: data set analyzed
  #models: fitted models for the nuisance regression functions
  #g: principal strata variable
  #z: z, treatment comparison candidate 1  
  #zp: z^\prime, treatment comparison candidate 2
  #J: number of treatment arms
  
  Z <- data$Z
  S <- data$S
  Y <- data$Y
  
  #ind: indicator for classifying different cases
  ind <- length(unique(c(J-g,J-g+1,z)))
  
  if (ind==3){
    if (g==J){
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))
      Xz <- grab_design_matrix(data = data,
                               rhs_formula = grab_fixed_formula(models$pz_model))
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      
      pJg1_pos <- 1:ncol(XJg1)
      pz_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xz))
      pzp_pos <- (max(pz_pos) + 1):(max(pz_pos) + ncol(Xzp))
      
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pz_scores <- grab_psiFUN(models$pz_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
    } else {
      XJg <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pJg_model))
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))
      Xz <- grab_design_matrix(data = data,
                               rhs_formula = grab_fixed_formula(models$pz_model))
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      
      pJg_pos <- 1:ncol(XJg)
      pJg1_pos <- (max(pJg_pos) + 1):(max(pJg_pos) + ncol(XJg1))
      pz_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xz))
      pzp_pos <- (max(pz_pos) + 1):(max(pz_pos) + ncol(Xzp))
      
      pJg_scores <- grab_psiFUN(models$pJg_model, data)
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pz_scores <- grab_psiFUN(models$pz_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
    }
  } else if (ind==2) {
    if (g==J){
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      
      pJg1_pos <- 1:ncol(XJg1)
      pzp_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xzp))
      
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
    } else {
      XJg <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pJg_model))
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      
      pJg_pos <- 1:ncol(XJg)
      pJg1_pos <- (max(pJg_pos) + 1):(max(pJg_pos) + ncol(XJg1))
      pzp_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xzp))
      
      pJg_scores <- grab_psiFUN(models$pJg_model, data)
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
    }
  }
  
  function(theta, g, z, zp, J){
    
    ind <- length(unique(c(J-g,J-g+1,z)))
    
    if (ind==3){
      if (g==J){
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pz <- plogis(Xz %*% theta[pz_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        c(pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pz_scores(theta[pz_pos])*I(Z==z),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          4*S*I(Z==1)-theta[max(pzp_pos)+1],
          4*I(Z==z)*S*(pJg1)/pz/(theta[max(pzp_pos)+1])*Y-4*I(Z==zp)*S*(pJg1)/pzp/(theta[max(pzp_pos)+1])*Y-theta[max(pzp_pos)+2]) 
      } else {
        pJg <- plogis(XJg %*% theta[pJg_pos])
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pz <- plogis(Xz %*% theta[pz_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        c(pJg_scores(theta[pJg_pos])*I(Z==J-g),
          pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pz_scores(theta[pz_pos])*I(Z==z),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          4*S*I(Z==J-g)-theta[max(pzp_pos)+1],
          4*S*I(Z==J-g+1)-theta[max(pzp_pos)+2],
          4*I(Z==z)*S*(pJg1-pJg)/pz/(theta[max(pzp_pos)+2]-theta[max(pzp_pos)+1])*Y-4*I(Z==zp)*S*(pJg1-pJg)/pzp/(theta[max(pzp_pos)+2]-theta[max(pzp_pos)+1])*Y-theta[max(pzp_pos)+3]) 
      }
    } else {
      if (g==J){
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        c(pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          4*S*I(Z==1)-theta[max(pzp_pos)+1],
          4*I(Z==z)*S*Y/theta[max(pzp_pos)+1]-4*I(Z==zp)*S*(pJg1)/pzp/(theta[max(pzp_pos)+1])*Y-theta[max(pzp_pos)+2])
      } else {
        pJg <- plogis(XJg %*% theta[pJg_pos])
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        c(pJg_scores(theta[pJg_pos])*I(Z==J-g),
          pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          4*S*I(Z==J-g)-theta[max(pzp_pos)+1],
          4*S*I(Z==J-g+1)-theta[max(pzp_pos)+2],
          4*I(Z==z)*S*(pJg1-pJg)/pJg1/(theta[max(pzp_pos)+2]-theta[max(pzp_pos)+1])*Y-4*I(Z==zp)*S*(pJg1-pJg)/pzp/(theta[max(pzp_pos)+2]-theta[max(pzp_pos)+1])*Y-theta[max(pzp_pos)+3]) 
      }
    }
  }
}

#Joint estimating equations for the outcome regression estimator
eq_or <- function(data, models, g, z, zp, J){
  #Input:
  #data: data set analyzed
  #models: fitted models for the nuisance regression functions
  #g: principal strata variable
  #z: z, treatment comparison candidate 1  
  #zp: z^\prime, treatment comparison candidate 2
  #J: number of treatment arms
  
  Z <- data$Z
  S <- data$S
  Y <- data$Y
  
  Xom <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$om_model))  
  Xomp <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$omp_model))  
  om_pos <- 1:ncol(Xom)
  omp_pos <- (max(om_pos)+1):(max(om_pos)+ncol(Xomp))
  
  om_scores <- grab_psiFUN(models$om_model, data)
  omp_scores <- grab_psiFUN(models$omp_model, data)
  
  function(theta, g, z, zp, J){
    if (g==J){
      om <- Xom %*% theta[om_pos]
      omp <- Xomp %*% theta[omp_pos]
      c(om_scores(theta[om_pos])*I(Z==z)*I(S==1),
        omp_scores(theta[omp_pos])*I(Z==zp)*I(S==1),
        4*S*I(Z==1)-theta[max(omp_pos)+1],
        (S*I(Z==1)*4)/(theta[max(omp_pos)+1])*(om-omp)-theta[max(omp_pos)+2])
    } else {
      om <- Xom %*% theta[om_pos]
      omp <- Xomp %*% theta[omp_pos]
      c(om_scores(theta[om_pos])*I(Z==z)*I(S==1),
        omp_scores(theta[omp_pos])*I(Z==zp)*I(S==1),
        4*S*I(Z==J-g)-theta[max(omp_pos)+1],
        4*S*I(Z==J-g+1)-theta[max(omp_pos)+2],
        (S*I(Z==J-g+1)*4-S*I(Z==J-g)*4)/(theta[max(omp_pos)+2]-theta[max(omp_pos)+1])*(om-omp)-theta[max(omp_pos)+3])
    }
  }
}

#Joint estimating equations for the doubly robust estimator
eq_dr <- function(data, models, g, z, zp, J){
  #Input:
  #data: data set analyzed
  #models: fitted models for the nuisance regression functions
  #g: principal strata variable
  #z: z, treatment comparison candidate 1  
  #zp: z^\prime, treatment comparison candidate 2
  #J: number of treatment arms
  
  Z <- data$Z
  S <- data$S
  Y <- data$Y
  
  ind <- length(unique(c(J-g,J-g+1,z)))
  
  if (ind == 3){
    if (g==J){
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))  
      Xz <- grab_design_matrix(data = data,
                               rhs_formula = grab_fixed_formula(models$pz_model))
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      Xom <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$om_model))  
      Xomp <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$omp_model))  
      
      pJg1_pos <- 1: ncol(XJg1)
      pz_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xz))
      pzp_pos <- (max(pz_pos) + 1):(max(pz_pos) + ncol(Xzp))
      om_pos <- (max(pzp_pos) + 1):(max(pzp_pos) + ncol(Xom)) 
      omp_pos <- (max(om_pos) + 1):(max(om_pos) + ncol(Xomp)) 
      
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pz_scores <- grab_psiFUN(models$pz_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
      om_scores <- grab_psiFUN(models$om_model, data)
      omp_scores <- grab_psiFUN(models$omp_model, data)
    } else {
      XJg <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pJg_model))
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))  
      Xz <- grab_design_matrix(data = data,
                               rhs_formula = grab_fixed_formula(models$pz_model))
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      Xom <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$om_model))  
      Xomp <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$omp_model))  
      
      pJg_pos <- 1:ncol(XJg)
      pJg1_pos <-(max(pJg_pos) + 1):(max(pJg_pos) + ncol(XJg1))
      pz_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xz))
      pzp_pos <- (max(pz_pos) + 1):(max(pz_pos) + ncol(Xzp))
      om_pos <- (max(pzp_pos) + 1):(max(pzp_pos) + ncol(Xom)) 
      omp_pos <- (max(om_pos) + 1):(max(om_pos) + ncol(Xomp)) 
      
      pJg_scores <- grab_psiFUN(models$pJg_model, data)
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pz_scores <- grab_psiFUN(models$pz_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
      om_scores <- grab_psiFUN(models$om_model, data)
      omp_scores <- grab_psiFUN(models$omp_model, data)
      
    }
  } else {
    if (g==J) {
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))  
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      Xom <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$om_model))  
      Xomp <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$omp_model))  
      
      pJg1_pos <- 1: ncol(XJg1)
      pzp_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xzp))
      om_pos <- (max(pzp_pos) + 1):(max(pzp_pos) + ncol(Xom)) 
      omp_pos <- (max(om_pos) + 1):(max(om_pos) + ncol(Xomp)) 
      
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
      om_scores <- grab_psiFUN(models$om_model, data)
      omp_scores <- grab_psiFUN(models$omp_model, data)
    } else {
      XJg <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pJg_model))
      XJg1 <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$pJg1_model))  
      Xzp <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$pzp_model))
      Xom <- grab_design_matrix(data = data,
                                rhs_formula = grab_fixed_formula(models$om_model))  
      Xomp <- grab_design_matrix(data = data,
                                 rhs_formula = grab_fixed_formula(models$omp_model))  
      
      pJg_pos <- 1:ncol(XJg)
      pJg1_pos <-(max(pJg_pos) + 1):(max(pJg_pos) + ncol(XJg1))
      pzp_pos <- (max(pJg1_pos) + 1):(max(pJg1_pos) + ncol(Xzp))
      om_pos <- (max(pzp_pos) + 1):(max(pzp_pos) + ncol(Xom)) 
      omp_pos <- (max(om_pos) + 1):(max(om_pos) + ncol(Xomp)) 
      
      pJg_scores <- grab_psiFUN(models$pJg_model, data)
      pJg1_scores <- grab_psiFUN(models$pJg1_model, data)
      pzp_scores <- grab_psiFUN(models$pzp_model, data)
      om_scores <- grab_psiFUN(models$om_model, data)
      omp_scores <- grab_psiFUN(models$omp_model, data)
    }
  }
  
  function(theta, g, z, zp, J){
    
    if (ind==3){
      if (g==J){
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pz <- plogis(Xz %*% theta[pz_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        om <- Xom %*% theta[om_pos]
        omp <- Xomp %*% theta[omp_pos]
        psiJg1 <- I(Z==(J-g+1))*(S-pJg1)*4+pJg1
        c(pJg1_scores(theta[pJg1_pos])*I(Z==(J-g+1)),
          pz_scores(theta[pz_pos])*I(Z==z),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          om_scores(theta[om_pos])*I(Z==z)*I(S==1),
          omp_scores(theta[omp_pos])*I(Z==zp)*I(S==1),
          4*I(Z==(J-g+1))*(S-pJg1)+pJg1-theta[max(omp_pos)+1],
          (pJg1)*S*I(Z==z)/(pz*1/4*(theta[max(omp_pos)+1]))*(Y-om)+om*(psiJg1)/(theta[max(omp_pos)+1])-
            (pJg1)*S*I(Z==zp)/(pzp*1/4*(theta[max(omp_pos)+1]))*(Y-omp)-omp*(psiJg1)/(theta[max(omp_pos)+1]) - theta[max(omp_pos)+2])
      } else {
        pJg <- plogis(XJg %*% theta[pJg_pos])
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pz <- plogis(Xz %*% theta[pz_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        om <- Xom %*% theta[om_pos]
        omp <- Xomp %*% theta[omp_pos]
        psiJg <- I(Z==J-g)*(S-pJg)*4+pJg
        psiJg1 <- I(Z==J-g+1)*(S-pJg1)*4+pJg1
        c(pJg_scores(theta[pJg_pos])*I(Z==J-g),
          pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pz_scores(theta[pz_pos])*I(Z==z),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          om_scores(theta[om_pos])*I(Z==z)*I(S==1),
          omp_scores(theta[omp_pos])*I(Z==zp)*I(S==1),
          4*I(Z==(J-g))*(S-pJg)+pJg-theta[max(omp_pos)+1],
          4*I(Z==(J-g+1))*(S-pJg1)+pJg1-theta[max(omp_pos)+2],
          (pJg1-pJg)*S*I(Z==z)/(pz*1/4*(theta[max(omp_pos)+2]-theta[max(omp_pos)+1]))*(Y-om)+om*(psiJg1-psiJg)/(theta[max(omp_pos)+2]-theta[max(omp_pos)+1])-
            (pJg1-pJg)*S*I(Z==zp)/(pzp*1/4*(theta[max(omp_pos)+2]-theta[max(omp_pos)+1]))*(Y-omp)-omp*(psiJg1-psiJg)/(theta[max(omp_pos)+2]-theta[max(omp_pos)+1])-theta[max(omp_pos)+3])
      }
    } else if (ind==2){
      if (g==J){
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        om <- Xom %*% theta[om_pos]
        omp <- Xomp %*% theta[omp_pos]
        psiJg1 <- I(Z==J-g+1)*(S-pJg1)*4+pJg1
        c(pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          om_scores(theta[om_pos])*I(Z==z)*I(S==1),
          omp_scores(theta[omp_pos])*I(Z==zp)*I(S==1),
          4*I(Z==(J-g+1))*(S-pJg1)+pJg1-theta[max(omp_pos)+1],
          S*I(Z==z)/(1/4*(theta[max(omp_pos)+1]))*(Y-om)+om*(psiJg1)/(theta[max(omp_pos)+1])-
            pJg1*S*I(Z==zp)/(pzp*1/4*(theta[max(omp_pos)+1]))*(Y-omp)-omp*(psiJg1)/(theta[max(omp_pos)+1])-theta[max(omp_pos)+2]) 
      } else {
        pJg <- plogis(XJg %*% theta[pJg_pos])
        pJg1 <- plogis(XJg1 %*% theta[pJg1_pos])
        pzp <- plogis(Xzp %*% theta[pzp_pos])
        om <- Xom %*% theta[om_pos]
        omp <- Xomp %*% theta[omp_pos]
        psiJg <- I(Z==J-g)*(S-pJg)*4+pJg
        psiJg1 <- I(Z==J-g+1)*(S-pJg1)*4+pJg1
        c(pJg_scores(theta[pJg_pos])*I(Z==J-g),
          pJg1_scores(theta[pJg1_pos])*I(Z==J-g+1),
          pzp_scores(theta[pzp_pos])*I(Z==zp),
          om_scores(theta[om_pos])*I(Z==z)*I(S==1),
          omp_scores(theta[omp_pos])*I(Z==zp)*I(S==1),
          4*I(Z==(J-g))*(S-pJg)+pJg-theta[max(omp_pos)+1],
          4*I(Z==(J-g+1))*(S-pJg1)+pJg1-theta[max(omp_pos)+2],
          (pJg1-pJg)*S*I(Z==z)/(pJg1*1/4*(theta[max(omp_pos)+2]-theta[max(omp_pos)+1]))*(Y-om)+om*(psiJg1-psiJg)/(theta[max(omp_pos)+2]-theta[max(omp_pos)+1])-
            (pJg1-pJg)*S*I(Z==zp)/(pzp*1/4*(theta[max(omp_pos)+2]-theta[max(omp_pos)+1]))*(Y-omp)-omp*(psiJg1-psiJg)/(theta[max(omp_pos)+2]-theta[max(omp_pos)+1])-theta[max(omp_pos)+3])
      }
    }
  }
}

#Function calculating point estimates and interval estimates
point_est <- function(df, g, z, zp){
  #Input:
  #df: data set analyzed
  #g: principal strata variable
  #z: z, treatment comparison candidate 1  
  #zp: z^\prime, treatment comparison candidate 2
  
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
  
  #calculate principal score weighting estimator
  if (g==J){
    w1g <- (prinscore[[2]])/prinscore[[z+1]]*mean(I(Z==z)*S)/mean(S*I(Z==1))
    w1gp <- (prinscore[[2]])/prinscore[[zp+1]]*mean(I(Z==zp)*S)/mean(S*I(Z==1))
  } else {
    w1g <- (prinscore[[J-g+2]]-prinscore[[J-g+1]])/prinscore[[z+1]]*mean(I(Z==z)*S)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))
    w1gp <- (prinscore[[J-g+2]]-prinscore[[J-g+1]])/prinscore[[zp+1]]*mean(I(Z==zp)*S)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))
  }
  est1 <-  mean(I(Z==z)*S*w1g*Y)/mean(I(Z==z)*S)
  est1p <-  mean(I(Z==zp)*S*w1gp*Y)/mean(I(Z==zp)*S)
  #est1c: point estimate for contrast using principal score weighting
  est1c <- est1-est1p
  
  #calculate outcome regression estimator
  factor0 <-  prinscore0
  factor1 <-  S*I(Z==1)/proscore[[1]]
  factor2 <-  S*I(Z==2)/proscore[[2]]
  factor3 <-  S*I(Z==3)/proscore[[3]]
  factor4 <-  S*I(Z==4)/proscore[[4]]
  factor5 <-  prinscore5
  factor <-  list(factor0, factor1, factor2, factor3, factor4, factor5)
  if (g==J){
    est2 <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omf)/(mean(S*I(Z==1)))
    est2p <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omfp)/(mean(S*I(Z==1)))
    #est2c: point estimate for contrast using outcome regression weighting
    est2c <- est2 - est2p
  } else {
    est2 <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omf)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))
    est2p <-  1/4*mean((factor[[J-g+2]]-factor[[J-g+1]])*omfp)/(mean(S*I(Z==J-g+1))-mean(S*I(Z==J-g)))
    est2c <- est2 - est2p
  }
  
  #calculate doubly robust estimator
  dr <- mean((prinscore[[J-g+2]]-prinscore[[J-g+1]])*S*I(Z==z)/prinscore[[z+1]]*4*(Y-omf)+omf*(psi_s[[J-g+2]]-psi_s[[J-g+1]]))/mean(psi_s[[J-g+2]]-psi_s[[J-g+1]])
  drp <- mean((prinscore[[J-g+2]]-prinscore[[J-g+1]])*S*I(Z==zp)/prinscore[[zp+1]]*4*(Y-omfp)+omfp*(psi_s[[J-g+2]]-psi_s[[J-g+1]]))/mean(psi_s[[J-g+2]]-psi_s[[J-g+1]])
  #drc: point estimate for contrast using doubly robust weighting
  drc <- dr - drp
  
  ##################################################variance estimates###############################################################
  #For this part, we use geex package to implement the calculations for sandwich variance estimator. Please consult the documentations
  #for geex package for more information and examples.
  
  #ind: indicator for classifying difference scenarios of estimating equations
  ind <- length(unique(c(J-g,J-g+1,z)))
  theta_ps <- 
    if (ind==3){
      if (g==J){
        c(as.vector(coef(ps_model[[1]])),
          as.vector(coef(ps_model[[z]])),
          as.vector(coef(ps_model[[zp]])),
          4*mean(I(Z==1)*S),
          est1c) 
      } else {
        c(as.vector(coef(ps_model[[J-g]])),
          as.vector(coef(ps_model[[J-g+1]])),
          as.vector(coef(ps_model[[z]])),
          as.vector(coef(ps_model[[zp]])),
          4*mean(I(Z==J-g)*S),
          4*mean(I(Z==J-g+1)*S),
          est1c) 
      }
    } else {
      if (g==J){
        c(as.vector(coef(ps_model[[1]])),
          as.vector(coef(ps_model[[zp]])),
          4*mean(I(Z==1)*S),
          est1c) 
      } else {
        c(as.vector(coef(ps_model[[J-g]])),
          as.vector(coef(ps_model[[J-g+1]])),
          as.vector(coef(ps_model[[zp]])),
          4*mean(I(Z==J-g)*S),
          4*mean(I(Z==J-g+1)*S),
          est1c) 
      }
    }

  theta_or <- 
    if (g==J){
      c(as.vector(coef(om_model[[z]])),
        as.vector(coef(om_model[[zp]])),
        4*mean(I(Z==1)*S),
        est2c)
    } else {
      c(as.vector(coef(om_model[[z]])),
        as.vector(coef(om_model[[zp]])),
        4*mean(I(Z==J-g)*S),
        4*mean(I(Z==J-g+1)*S),
        est2c) 
    }
  
  theta_dr <-
    if (ind==3){
      if (g==J){
        c(as.vector(coef(ps_model[[J-g+1]])),
          as.vector(coef(ps_model[[z]])),
          as.vector(coef(ps_model[[zp]])),
          as.vector(coef(om_model[[z]])),
          as.vector(coef(om_model[[zp]])),
          mean(psi_s[[J-g+2]]),
          drc) 
      } else {
        c(as.vector(coef(ps_model[[J-g]])),
          as.vector(coef(ps_model[[J-g+1]])),
          as.vector(coef(ps_model[[z]])),
          as.vector(coef(ps_model[[zp]])),
          as.vector(coef(om_model[[z]])),
          as.vector(coef(om_model[[zp]])),
          mean(psi_s[[J-g+1]]),
          mean(psi_s[[J-g+2]]),
          drc) 
      }
    } else {
      if (g==J){
        c(as.vector(coef(ps_model[[J-g+1]])),
          as.vector(coef(ps_model[[zp]])),
          as.vector(coef(om_model[[z]])),
          as.vector(coef(om_model[[zp]])),
          mean(psi_s[[J-g+2]]),
          drc)
      } else {
        c(as.vector(coef(ps_model[[J-g]])),
          as.vector(coef(ps_model[[J-g+1]])),
          as.vector(coef(ps_model[[zp]])),
          as.vector(coef(om_model[[z]])),
          as.vector(coef(om_model[[zp]])),
          mean(psi_s[[J-g+1]]),
          mean(psi_s[[J-g+2]]),
          drc)
      }
    }
  
  estimate_psw <- function(data, g, z, zp, J){
    pJg1_model <- glm(S~A+C, subset = (Z==J-g+1), data=data, family = binomial) 
    if (g<J){
      pJg_model <- glm(S~A+C, subset = (Z==J-g), data=data, family = binomial)
    } else {
      pJg_model <- pJg1_model
    }
    pz_model <- glm(S~A+C, subset = (Z==z), data=data, family = binomial) 
    pzp_model <- glm(S~A+C, subset = (Z==zp), data=data, family = binomial) 
    models <- list(pJg_model=pJg_model, pJg1_model=pJg1_model, pz_model=pz_model,
                   pzp_model=pzp_model)
    m_estimate(estFUN = eq_psw, data = data,
               outer_args = list(models = models, g=g, z=z, zp=zp, J=J),
               compute_roots = FALSE,
               roots = theta_ps,
               inner_args = list(g = g, z = z, zp=zp, J = J)
    )
  }
  
  estimate_or <- function(data, g, z, zp, J){
    om_model <- glm(Y~A+C, subset = (Z==z&S==1), data=data, family = gaussian)
    omp_model <- glm(Y~A+C, subset = (Z==zp&S==1), data=data, family = gaussian)
    models <- list(om_model=om_model,omp_model=omp_model)
    m_estimate(estFUN = eq_or, data = data,
               outer_args = list(models = models, g=g, z=z, zp=zp, J=J),
               compute_roots = FALSE,
               roots = theta_or,
               inner_args = list(g = g, z = z, zp = zp, J = J)
    )
  }
  
  estimate_dr <- function(data, g, z, zp, J){
    om_model <- glm(Y~A+C, subset = (Z==z&S==1), data=data, family = gaussian)
    omp_model <- glm(Y~A+C, subset = (Z==zp&S==1), data=data, family = gaussian)
    pJg1_model <- glm(S~A+C, subset = (Z==J-g+1), data=data, family = binomial) 
    if (g<J){
      pJg_model <- glm(S~A+C, subset = (Z==J-g), data=data, family = binomial)
    } else {
      pJg_model <- pJg1_model
    }
    pz_model <- glm(S~A+C, subset = (Z==z), data=data, family = binomial) 
    pzp_model <- glm(S~A+C, subset = (Z==zp), data=data, family = binomial) 
    models <- list(pJg_model=pJg_model, pJg1_model=pJg1_model, pz_model=pz_model, pzp_model=pzp_model,
                   om_model=om_model, omp_model=omp_model)
    m_estimate(estFUN = eq_dr, data = data,
               outer_args = list(models = models, g=g, z=z, zp=zp, J=J),
               compute_roots = FALSE,
               roots = theta_dr,
               inner_args = list(g = g, z = z, zp=zp, J = J)
    )
  }
  
  cov_matrix1 <- vcov(estimate_psw(data=df, g = g, z = z, zp=zp, J = 4))
  cov_matrix2 <- vcov(estimate_or(data=df, g = g, z = z, zp=zp, J = 4))
  cov_matrix_dr <- vcov(estimate_dr(data=df, g = g, z = z, zp=zp, J = 4))
  
  se1 <- sqrt(cov_matrix1[ncol(cov_matrix1),ncol(cov_matrix1)]) #se1: standard error estimate for the principal score weighting estimator
  se2 <- sqrt(cov_matrix2[ncol(cov_matrix2),ncol(cov_matrix2)])  #se2: standard error estimate for the outcome regression estimator
  se_dr <- sqrt(cov_matrix_dr[ncol(cov_matrix_dr),ncol(cov_matrix_dr)])  #se_dr: standard error estimate for the doubly robust estimator
  
  return(c(est1c,se1,est2c,se2,drc,se_dr))
}