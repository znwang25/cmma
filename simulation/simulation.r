library(parallel)
library(MASS)
library(tidyverse)
library(stats4)
library(knitr)
library(kableExtra)
library(lfe)
library(aod)
library(sandwich)
library(ggthemes)
library(stats)
set.seed(42)


simulate_data <- function(y_model, samplesize_per_trial, no_trials,  sd.m00 = 3,
                          sd.y000 = 3, cor.m00.y000 = 0.95,
                          num_treatment_group = 3,
                          mu.tau = c(0.5,1,2.5),
                          mu.gamma = c(0,1.5,3),
                          violate_A2 = FALSE,
                          violate_A3 = FALSE,
                          violate_A4 = FALSE
){
  # individual error structure
  mu <- c(m00=0, y000=0, noise_beta_correlated = 0) # means
  
  sd <- c(sd.m00, sd.y000, 0.5) # standard deviations
  cormat <- rbind(c(1,cor.m00.y000, 0.8),c(cor.m00.y000,1, 0.8), c(0.8, 0.8, 1))
  # var-cov matrix
  covmat <- diag(sd) %*% cormat %*% diag(sd) 
  
  df <- mvrnorm(samplesize_per_trial*no_trials,mu,covmat)%>%
    as_tibble()
  
  # Set model parameters
  phi <- runif(no_trials,-2,2)
  theta <- runif(no_trials,-2,2)
  trial_tau <- runif(no_trials,-3, 3)
  trial_gamma <- rep(0,no_trials)
  tau_innovation <- rnorm(samplesize_per_trial*no_trials,0,0.5)
  gamma_innovation <- rnorm(samplesize_per_trial*no_trials,0,0.5)
  beta_innovation <- rnorm(samplesize_per_trial*no_trials,0,0.5)
  
  if (violate_A2){
    beta_innovation <- df$noise_beta_correlated
  }
  if (violate_A3){
    complement <- function(y, rho, x) {
      # Generate a new vector correlates with y
      if (missing(x)) x <- runif(length(y),-1, 1) # Optional: supply a default if `x` is not given
      y.perp <- residuals(lm(x ~ y))
      rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
    }
    trial_gamma <- complement(trial_tau, 0.8)
  }
  if (violate_A4){
    trial_tau <- rep(0,no_trials)
  }
  df <- 
    df %>% mutate(    trial_index=row_number()%%no_trials+1,
                      treatment_group=trial_index%%num_treatment_group+1,
                      Ti = rbernoulli(samplesize_per_trial*no_trials),
                      phi_s = phi[as.vector(trial_index)], 
                      phi_s = phi[as.vector(trial_index)], 
                      phi_s = phi[as.vector(trial_index)], 
                      tau_s = mu.tau[as.vector(treatment_group)]+ trial_tau[as.vector(trial_index)],
                      theta_s = theta[as.vector(trial_index)],
                      gamma_s = mu.gamma[as.vector(treatment_group)]+trial_gamma[as.vector(trial_index)],
                      noise_tau=tau_innovation,
                      noise_gamma=gamma_innovation,
                      noise_beta=beta_innovation,
                      Mi = phi_s + tau_s*Ti+ m00 + noise_tau*Ti,
                      Yi = y_model(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta)
    )%>%
    mutate(treatment_group=factor(treatment_group),
           trial_index=factor(trial_index))

  return(df)
}


regress_one_trial <- function(data){
  data <- data%>%mutate(Mi2=Mi^2,Mi3=Mi^3)
  result <- lm(Yi~1+Ti, data=data)%>%summary
  delta_y <- result$coefficients['TiTRUE',1]
  delta_y_ste <- result$coefficients['TiTRUE',2]
  result <- lm(Mi~1+Ti, data=data)%>%summary
  delta_m <- result$coefficients['TiTRUE',1]
  delta_m_ste <- result$coefficients['TiTRUE',2]
  result <- lm(Mi2~1+Ti, data=data)%>%summary
  delta_m2 <- result$coefficients['TiTRUE',1]
  delta_m2_ste <- result$coefficients['TiTRUE',2]
  result <- lm(Mi3~1+Ti, data=data)%>%summary
  delta_m3 <- result$coefficients['TiTRUE',1]
  delta_m3_ste <- result$coefficients['TiTRUE',2]
  out <- tibble(delta_y,delta_y_ste,delta_m,delta_m_ste,delta_m2,delta_m2_ste,delta_m3,delta_m3_ste)
  return(out)
}


mc_model_main <- function(y_model, samplesize_per_trial,no_trials, no_sims=100, seed=42, proposed_estimator_only=FALSE, DRF_degree=1, wald_test=FALSE,violate_A2 = FALSE,
                          violate_A3 = FALSE, violate_A4 = FALSE){
  set.seed(seed)
  
  MC_trial<- function(index){
    if(index %%10==0){
      print(paste("Prcessing simulation no.",index))
    }
    sim_data <- simulate_data(y_model,samplesize_per_trial,no_trials,violate_A2 = violate_A2,
                              violate_A3 = violate_A3, violate_A4 = violate_A4)
    sim_data <- sim_data%>%mutate(Mi2=Mi^2,Mi3=Mi^3)
    sim_data <- sim_data%>%bind_cols(
      as_tibble(model.matrix(~0+treatment_group,data = sim_data)*sim_data$Ti)%>%
        rename_at(vars(matches("treatment_group")),list(~paste(.,"Ti",sep =""))))
    
    meta_data <- sim_data%>%
      group_by(trial_index,treatment_group)%>%
      do(regress_one_trial(.))%>%
      ungroup
    
    out <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(out) <- c("valname", "value", "model")
    out <- as_tibble(out)
    
    # only M
    if (DRF_degree==1) {

    ## Proposed method
    regres <- lm(delta_y~1+delta_m+treatment_group,data=meta_data)%>%summary
    # Point estimates
    beta2.est <- regres$coefficients['delta_m',1]
    # Standard Error
    beta2.est.std <- regres$coefficients['delta_m',2]
    out <- out %>% 
      add_row(valname = "beta2.est", value=beta2.est, model ='proposed')%>%
      add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed')
    
      if (proposed_estimator_only == FALSE) {
        ## Simple pooled OLS (OLS SI)
        regres <- felm(Yi~1+Mi|factor(Ti):trial_index +trial_index|0|trial_index,data=sim_data)%>%summary
        # Point estimates from pooled ols
        beta2.est <- regres$coefficients['Mi',1]
        # Standard Error from pooled ols
        beta2.est.std <- regres$coefficients['Mi',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='pooled ols')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='pooled ols')
        
        ## IV method (small)
        regres <- felm(Yi~1+Ti|trial_index|(Mi~trial_index:Ti)|trial_index,data=sim_data)%>%summary
        # Point estimates from pooled ols
        beta2.est <- regres$coefficients['`Mi(fit)`',1]
        # Standard Error from pooled ols
        beta2.est.std <- regres$coefficients['`Mi(fit)`',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='IV (small)')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='IV (small)')
        
        ## proposed IV equivalent
        iv_form <- as.formula(paste("Yi~1+", 
                                    paste(grep("treatment_group.*Ti", names(sim_data), value=TRUE), collapse = " + "),
                                    "|trial_index|(Mi~trial_index:factor(Ti))|0") )
        regres <- felm(iv_form,data=sim_data)%>%summary    
        # Point estimates from pooled ols
        beta2.est <- regres$coefficients['`Mi(fit)`',1]
        # Standard Error from pooled ols
        beta2.est.std <- regres$coefficients['`Mi(fit)`',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='proposed IV equivalent')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed IV equivalent')
        
        ## Liml
        regres <- felm(iv_form,data=sim_data,kclass='liml')%>%summary    
        # Point estimates from pooled ols
        beta2.est <- regres$coefficients['Mi',1]
        # Standard Error from pooled ols
        beta2.est.std <- regres$coefficients['Mi',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='Liml equivalent')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='Liml equivalent')
      }
      
    }
    else if (DRF_degree==2) {
      ## Proposed method
      regres <- lm(delta_y~1+delta_m+delta_m2+treatment_group,data=meta_data)%>%summary
      # Point estimates
      beta2.est <- regres$coefficients['delta_m',1]
      beta3.est <- regres$coefficients['delta_m2',1]
      # Standard Error
      beta2.est.std <- regres$coefficients['delta_m',2]
      beta3.est.std <- regres$coefficients['delta_m2',2]
      out <- out %>% 
        add_row(valname = "beta2.est", value=beta2.est, model ='proposed with m2')%>%
        add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed with m2')%>%
        add_row(valname = "beta3.est", value=beta3.est, model ='proposed with m2')%>%
        add_row(valname = "beta3.est.std", value=beta3.est.std, model ='proposed with m2')
      
      if (proposed_estimator_only == FALSE) {
        ## Simple pooled OLS (OLS SI)
        regres <- felm(Yi~1+Mi+Mi2|factor(Ti):trial_index +trial_index|0|trial_index,data=sim_data)%>%summary
        # Point estimates
        beta2.est <- regres$coefficients['Mi',1]
        beta3.est <- regres$coefficients['Mi2',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['Mi',2]
        beta3.est.std <- regres$coefficients['Mi2',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='ols with m2')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='ols with m2')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='ols with m2')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='ols with m2')
        
        ## IV method (small)
        regres <- felm(Yi~1+Ti|trial_index|(Mi|Mi2~trial_index:Ti)|trial_index,data=sim_data)%>%summary
        # Point estimates
        beta2.est <- regres$coefficients['`Mi(fit)`',1]
        beta3.est <- regres$coefficients['`Mi2(fit)`',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['`Mi(fit)`',2]
        beta3.est.std <- regres$coefficients['`Mi2(fit)`',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='IV (small) with m2')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='IV (small) with m2')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='IV (small) with m2')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='IV (small) with m2')
        
        ## proposed IV equivalent
        iv_form <- as.formula(paste("Yi~1+", 
                                    paste(grep("treatment_group.*Ti", names(sim_data), value=TRUE), collapse = " + "),
                                    "|trial_index|(Mi|Mi2~trial_index:factor(Ti))|0") )
        
        regres <- felm(iv_form,data=sim_data)%>%summary
        # Point estimates
        beta2.est <- regres$coefficients['`Mi(fit)`',1]
        beta3.est <- regres$coefficients['`Mi2(fit)`',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['`Mi(fit)`',2]
        beta3.est.std <- regres$coefficients['`Mi2(fit)`',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='proposed IV equivalent with m2')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed IV equivalent with m2')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='proposed IV equivalent with m2')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='proposed IV equivalent with m2')
        
        ## Liml
        regres <- felm(iv_form,data=sim_data,kclass='liml')%>%summary    
        # Point estimates
        beta2.est <- regres$coefficients['Mi',1]
        beta3.est <- regres$coefficients['Mi2',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['Mi',2]
        beta3.est.std <- regres$coefficients['Mi2',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='Liml equivalent with m2')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='Liml equivalent with m2')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='Liml equivalent with m2')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='Liml equivalent with m2')
      }
    }
    else if (DRF_degree==3) {
      ## Proposed Method
      regres <- lm(delta_y~1+delta_m+ delta_m2+ delta_m3 + treatment_group,data=meta_data)%>%summary
      # Point estimates
      beta2.est <- regres$coefficients['delta_m',1]
      beta3.est <- regres$coefficients['delta_m2',1]
      beta4.est <- regres$coefficients['delta_m3',1]
      # Standard Error
      beta2.est.std <- regres$coefficients['delta_m',2]
      beta3.est.std <- regres$coefficients['delta_m2',2]
      beta4.est.std <- regres$coefficients['delta_m3',2]
      out <- out %>% 
        add_row(valname = "beta2.est", value=beta2.est, model ='proposed with m2 m3')%>%
        add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed with m2 m3')%>%
        add_row(valname = "beta3.est", value=beta3.est, model ='proposed with m2 m3')%>%
        add_row(valname = "beta3.est.std", value=beta3.est.std, model ='proposed with m2 m3')%>%
        add_row(valname = "beta4.est", value=beta4.est, model ='proposed with m2 m3')%>%
        add_row(valname = "beta4.est.std", value=beta4.est.std, model ='proposed with m2 m3')
      
      if (proposed_estimator_only == FALSE) {
        ## Simple pooled OLS (OLS SI)
        # M and M^2 and M^3
        regres <- felm(Yi~1+Mi+Mi2+Mi3|factor(Ti):trial_index+trial_index|0|trial_index,data=sim_data)%>%summary
        # Point estimates
        beta2.est <- regres$coefficients['Mi',1]
        beta3.est <- regres$coefficients['Mi2',1]
        beta4.est <- regres$coefficients['Mi3',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['Mi',2]
        beta3.est.std <- regres$coefficients['Mi2',2]
        beta4.est.std <- regres$coefficients['Mi3',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='ols with m2 m3')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='ols with m2 m3')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='ols with m2 m3')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='ols with m2 m3')%>%
          add_row(valname = "beta4.est", value=beta4.est, model ='ols with m2 m3')%>%
          add_row(valname = "beta4.est.std", value=beta4.est.std, model ='ols with m2 m3')
        
        ## IV method (small)
        regres <- felm(Yi~1+Ti|trial_index|(Mi|Mi2|Mi3~trial_index:Ti)|trial_index,data=sim_data)%>%summary
        # Point estimates
        beta2.est <- regres$coefficients['`Mi(fit)`',1]
        beta3.est <- regres$coefficients['`Mi2(fit)`',1]
        beta4.est <- regres$coefficients['`Mi3(fit)`',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['`Mi(fit)`',2]
        beta3.est.std <- regres$coefficients['`Mi2(fit)`',2]
        beta4.est.std <- regres$coefficients['`Mi3(fit)`',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='IV (small) with m2 m3')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='IV (small) with m2 m3')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='IV (small) with m2 m3')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='IV (small) with m2 m3')%>%
          add_row(valname = "beta4.est", value=beta4.est, model ='IV (small) with m2 m3')%>%
          add_row(valname = "beta4.est.std", value=beta4.est.std, model ='IV (small) with m2 m3')
        
        ## proposed IV equivalent
        iv_form <- as.formula(paste("Yi~1+", 
                                    paste(grep("treatment_group.*Ti", names(sim_data), value=TRUE), collapse = " + "),
                                    "|trial_index|(Mi|Mi2|Mi3~trial_index:factor(Ti))|0") )
        
        regres <- felm(iv_form,data=sim_data)%>%summary
        # Point estimates
        beta2.est <- regres$coefficients['`Mi(fit)`',1]
        beta3.est <- regres$coefficients['`Mi2(fit)`',1]
        beta4.est <- regres$coefficients['`Mi3(fit)`',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['`Mi(fit)`',2]
        beta3.est.std <- regres$coefficients['`Mi2(fit)`',2]
        beta4.est.std <- regres$coefficients['`Mi3(fit)`',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='proposed IV equivalent with m2 m3')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed IV equivalent with m2 m3')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='proposed IV equivalent with m2 m3')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='proposed IV equivalent with m2 m3')%>%
          add_row(valname = "beta4.est", value=beta4.est, model ='proposed IV equivalent with m2 m3')%>%
          add_row(valname = "beta4.est.std", value=beta4.est.std, model ='proposed IV equivalent with m2 m3')
        
        ## Liml
        regres <- felm(iv_form,data=sim_data,kclass='liml')%>%summary    
        # Point estimates
        beta2.est <- regres$coefficients['Mi',1]
        beta3.est <- regres$coefficients['Mi2',1]
        beta4.est <- regres$coefficients['Mi3',1]
        # Standard Error
        beta2.est.std <- regres$coefficients['Mi',2]
        beta3.est.std <- regres$coefficients['Mi2',2]
        beta4.est.std <- regres$coefficients['Mi3',2]
        out <- out %>% 
          add_row(valname = "beta2.est", value=beta2.est, model ='Liml equivalent with m2 m3')%>%
          add_row(valname = "beta2.est.std", value=beta2.est.std, model ='Liml equivalent with m2 m3')%>%
          add_row(valname = "beta3.est", value=beta3.est, model ='Liml equivalent with m2 m3')%>%
          add_row(valname = "beta3.est.std", value=beta3.est.std, model ='Liml equivalent with m2 m3')%>%
          add_row(valname = "beta4.est", value=beta4.est, model ='Liml equivalent with m2 m3')%>%
          add_row(valname = "beta4.est.std", value=beta4.est.std, model ='Liml equivalent with m2 m3')
      }
    }
    else {
    }
    if (wald_test){
      ## Proposed method
      regres <- lm(delta_y~1+delta_m+treatment_group,data=meta_data)%>%summary
      # Point estimates
      beta2.est <- regres$coefficients['delta_m',1]
      # Standard Error
      beta2.est.std <- regres$coefficients['delta_m',2]
      out <- out %>% 
        add_row(valname = "beta2.est", value=beta2.est, model ='proposed')%>%
        add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed')
      
      ## Proposed method
      regres <- lm(delta_y~1+delta_m+delta_m2+treatment_group,data=meta_data)%>%summary
      # Point estimates
      beta2.est <- regres$coefficients['delta_m',1]
      beta3.est <- regres$coefficients['delta_m2',1]
      # Standard Error
      beta2.est.std <- regres$coefficients['delta_m',2]
      beta3.est.std <- regres$coefficients['delta_m2',2]
      out <- out %>% 
        add_row(valname = "beta2.est", value=beta2.est, model ='proposed with m2')%>%
        add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed with m2')%>%
        add_row(valname = "beta3.est", value=beta3.est, model ='proposed with m2')%>%
        add_row(valname = "beta3.est.std", value=beta3.est.std, model ='proposed with m2')
      
      
      # M and M^2 and M^3
      regres_wald <- lm(delta_y~1+delta_m+ delta_m2+ delta_m3 + treatment_group,data=meta_data)
      regres <- regres_wald%>%summary
      # Point estimates
      beta2.est <- regres$coefficients['delta_m',1]
      beta3.est <- regres$coefficients['delta_m2',1]
      beta4.est <- regres$coefficients['delta_m3',1]
      # Standard Error
      beta2.est.std <- regres$coefficients['delta_m',2]
      beta3.est.std <- regres$coefficients['delta_m2',2]
      beta4.est.std <- regres$coefficients['delta_m3',2]
      out <- out %>% 
        add_row(valname = "beta2.est", value=beta2.est, model ='proposed with m2 m3')%>%
        add_row(valname = "beta2.est.std", value=beta2.est.std, model ='proposed with m2 m3')%>%
        add_row(valname = "beta3.est", value=beta3.est, model ='proposed with m2 m3')%>%
        add_row(valname = "beta3.est.std", value=beta3.est.std, model ='proposed with m2 m3')%>%
        add_row(valname = "beta4.est", value=beta4.est, model ='proposed with m2 m3')%>%
        add_row(valname = "beta4.est.std", value=beta4.est.std, model ='proposed with m2 m3')
      
      # Wald Test
      # regress cubic
      wald <-  wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = 2)
      out <- out %>% 
        add_row(valname = "p value, m", value=wald$result$chi2['P'], model ='wald test')
      wald <-  wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = 3)
      out <- out %>% 
        add_row(valname = "p value, m2", value=wald$result$chi2['P'], model ='wald test')
      wald <-  wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = 4)
      out <- out %>% 
        add_row(valname = "p value, m3", value=wald$result$chi2['P'], model ='wald test')
      wald <- wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = 2:3)
      out <- out %>% 
        add_row(valname = "p value, m m2", value=wald$result$chi2['P'], model ='wald test')
      wald <- wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = c(2,4))
      out <- out %>% 
        add_row(valname = "p value, m m3", value=wald$result$chi2['P'], model ='wald test')
      wald <- wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = 3:4)
      out <- out %>% 
        add_row(valname = "p value, m2 m3", value=wald$result$chi2['P'], model ='wald test')
      wald <- wald.test(b = coef(regres_wald), Sigma = vcov(regres_wald), Terms = 2:4)
      out <- out %>% 
        add_row(valname = "p value, m m2 m3", value=wald$result$chi2['P'], model ='wald test')
    }
    out <- out%>%mutate(simulation_index=index)
    return(out)
  }
  # parameters for monte carlo 
  outcomes <- mclapply(1:no_sims, MC_trial)
  outcomes <- outcomes %>% bind_rows%>%
    mutate("samplesize_per_trial"=samplesize_per_trial,"no_trials"=no_trials)
  return(outcomes)
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
# Compute coverage probabilities for the 90% and 95% CI and other stats from simulation results
stats_summary <- function(data,model_name,est_valname,std_valname, true_value,level=0.95){
  df <- data%>%filter(model==model_name)%>%
    spread(valname,value)%>%
    dplyr::select(est_valname,std_valname,samplesize_per_trial,no_trials)
  
  df <- df%>%rename(estimate=1, std=2)%>% 
    mutate_at(vars(estimate,std), remove_outliers)%>%
    mutate(CI.up = estimate + qnorm(0.5+level/2)*std,
           CI.lo = estimate - qnorm(0.5+level/2)*std,
           in_ci=CI.lo < true_value & true_value < CI.up)%>%
    summarise(ci_covprob = mean(in_ci,na.rm=T), 
              point.est = mean(estimate,na.rm=T),
              mse=mean((estimate-true_value)^2,na.rm=T),
              samplesize_per_trial=mean(samplesize_per_trial),
              no_trials=mean(no_trials),
              n= sum(!is.na(estimate)))%>%
    mutate(bias = point.est-true_value,
           rmse=sqrt(mse),
           Estimator=model_name)
  return(df)
}

simulation_summary <- function(data, models,est_valname,std_valname, true_value, level=0.95){
  lapply(models,
         stats_summary,
         data=data,est_valname=est_valname,std_valname=std_valname,true_value=true_value)%>%
  bind_rows()
}
get_avg_simulated_estimate <- function(data, model_name,est_valnames){
  data%>%filter(model==model_name)%>%
              filter(valname %in% est_valnames)%>%
              spread(valname,value)%>%
              dplyr::select(!!est_valnames)%>%
              summarise_all(mean)
}
# Get wald test result
wald_test_summary <- function(data, level=0.95){
  data%>%filter(model=='wald test')%>%
    spread(valname,value)%>%
    dplyr::select(-model,-samplesize_per_trial,-no_trials)%>%
    summarise_all(function(x) mean(x<(1-level)))
}

# First model linear
# Y_i =  beta*Mi + theta_s + gamma_s*Ti + y000 


simu_model_0 <- function(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta){
  beta <- 4
  Yi =  beta*Mi + theta_s + gamma_s*Ti + y000
  return(Yi)
}

simu_model_1 <- function(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta){
  beta <- 4
  Yi =  beta*Mi + theta_s + gamma_s*Ti + y000 + noise_gamma*Ti + noise_beta*Mi
  return(Yi)
}

# Homogeneous direct treatment effects
sim_data <- simulate_data(simu_model_1,1000,10,
                          num_treatment_group = 1,
                          mu.gamma = c(1.5))

#  OLS SI
felm(Yi~1+Mi|factor(Ti):trial_index +trial_index|0|trial_index,data=sim_data) %>% summary

#  IV
felm(Yi~1|trial_index|(Mi ~ factor(Ti):trial_index)|trial_index,data=sim_data) %>% summary

# IV (small)
felm(Yi~1+Ti|trial_index|(Mi~trial_index:factor(Ti))|trial_index,data=sim_data)%>%summary    

# Liml
felm(Yi~1+Ti|trial_index|(Mi~trial_index:factor(Ti)),data=sim_data, kclass='liml')%>%summary    


# Heterogeneous direct treatment effects
sim_data <- simulate_data(simu_model_1,1000,10, violate_A4 = T)
sim_data <- sim_data%>%bind_cols(
  as_tibble(model.matrix(~0+treatment_group,data = sim_data)*sim_data$Ti)%>%
    rename_at(vars(matches("treatment_group")),list(~paste(.,"Ti",sep =""))))

#  OLS SI
felm(Yi~1+Mi|factor(Ti):trial_index +trial_index|0|trial_index,data=sim_data) %>% summary

#  IV
felm(Yi~1|trial_index|(Mi ~ factor(Ti):trial_index)|trial_index,data=sim_data) %>% summary

# IV (small)
felm(Yi~1+Ti|trial_index|(Mi~trial_index:factor(Ti))|trial_index,data=sim_data)%>%summary    

# IV equivalent
iv_form <- as.formula(paste("Yi~1+", 
                            paste(grep("treatment_group.*Ti", names(sim_data), value=TRUE), collapse = " + "),
                            "|trial_index|(Mi~trial_index:factor(Ti))|trial_index") )
felm(iv_form,data=sim_data)%>%summary 

# Liml
felm(iv_form,data=sim_data, kclass='liml')%>%summary    
meta_data <- sim_data%>%
       group_by(trial_index,treatment_group)%>%
       do(regress_one_trial(.))%>%
       ungroup

regres <- lm(delta_y~ delta_m+treatment_group,data=meta_data)%>%summary

# 

# Second model second degree polynomial
# Y_i = beta[1]*Mi + beta[2]*Mi^2 + theta_s + gamma_s*Ti + y000

simu_model_2 <- function(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta){
  beta <- c(4, 2)
  Yi =  beta[1]*Mi + beta[2]*Mi^2 + theta_s + gamma_s*Ti + y000 + noise_gamma*Ti + noise_beta*Mi
  return(Yi)
}

sim_data <- simulate_data(simu_model_2,1000,10) %>% mutate(Mi2 = Mi^2)
sim_data <- sim_data%>%bind_cols(
  as_tibble(model.matrix(~0+treatment_group,data = sim_data)*sim_data$Ti)%>%
    rename_at(vars(matches("treatment_group")),list(~paste(.,"Ti",sep =""))))



#  OLS SI
felm(Yi~1+Mi + Mi2|factor(Ti):trial_index +trial_index|0|trial_index,data=sim_data) %>% summary

#  IV
felm(Yi~1|trial_index|(Mi + Mi2 ~ factor(Ti):trial_index)|trial_index,data=sim_data) %>% summary

# IV (small)
felm(Yi~1+Ti|trial_index|(Mi+Mi2~trial_index:factor(Ti))|trial_index,data=sim_data)%>%summary    

# IV equivalent
iv_form <- as.formula(paste("Yi~1+", 
                            paste(grep("treatment_group.*Ti", names(sim_data), value=TRUE), collapse = " + "),
                            "|trial_index|(Mi+Mi2~trial_index:factor(Ti))|trial_index") )
felm(iv_form,data=sim_data)%>%summary 

# Liml
felm(iv_form,data=sim_data, kclass='liml')%>%summary    


# Third model third degree polynomial
# Y_i = beta[1]*Mi + beta[2]*Mi^2 + beta[3]*Mi^3  + theta_s + gamma_s*Ti + y000 
simu_model_3 <- function(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta){
  beta <- c(4, 0, 5)
  Yi =  beta[1]*Mi + beta[2]*Mi^2 + beta[3]*Mi^3  + theta_s + gamma_s*Ti + y000 + noise_gamma*Ti
  return(Yi)
}
sim_data <- simulate_data(simu_model_3,1000,10) %>% mutate(Mi2 = Mi^2, Mi3 = Mi^3)
sim_data <- sim_data%>%bind_cols(
  as_tibble(model.matrix(~0+treatment_group,data = sim_data)*sim_data$Ti)%>%
    rename_at(vars(matches("treatment_group")),list(~paste(.,"Ti",sep =""))))

# IV (small)
felm(Yi~1+Ti|trial_index|(Mi+Mi2+Mi3~trial_index:factor(Ti))|trial_index,data=sim_data)%>%summary    

# IV equivalent
iv_form <- as.formula(paste("Yi~1+", 
                            paste(grep("treatment_group.*Ti", names(sim_data), value=TRUE), collapse = " + "),
                            "|trial_index|(Mi+Mi2+Mi3~trial_index:factor(Ti))|trial_index") )
felm(iv_form,data=sim_data)%>%summary 

# Liml
felm(iv_form,data=sim_data, kclass='liml')%>%summary    


out = mc_model_main(simu_model_1,200,10,100)
out = mc_model_main(simu_model_1,200,10,100, violate_A2 =T)
out = mc_model_main(simu_model_1,200,10,100, violate_A3 =T)
out = mc_model_main(simu_model_1,200,10,100, violate_A4 =T)

out%>%
  simulation_summary(models=c('proposed','proposed IV equivalent','pooled ols','IV (small)', 'Liml equivalent'),est_valname='beta2.est',std_valname='beta2.est.std', 4)


# Table output
# ## Increase samplesize per trial
out2 <- lapply(c(200,500,1000), function(n) mc_model_main(simu_model_1,n,50))
# out2%>%
#   map(simulation_summary,
#       models=c('proposed','proposed IV equivalent','pooled ols','IV (small)', 'Liml equivalent'),
#       est_valname='beta2.est',std_valname='beta2.est.std', 4)%>%
#   bind_rows()%>%
#   dplyr::select(samplesize_per_trial,no_trials,Estimator,bias,ci_covprob)%>%
#   group_by(Estimator)%>%
#   group_split()%>%
#   rev()%>%
#   bind_cols()%>%
#   dplyr::select(-samplesize_per_trial1,-no_trials1,-samplesize_per_trial2,-no_trials2, -Estimator,-Estimator1,-Estimator2)%>%
#   kable(format="latex",digits = 2,booktabs = T,
#         col.names = c("Sample size per trial`","No. of trials",
#                       "Bias","95% CI coverage", 
#                       "Bias","95% CI coverage",
#                       "Bias","95% CI coverage"))%>%
#   add_header_above(c(" " = 2, "Proposed estimator" = 2, "Pooled OLS estimator" = 2,"IV" = 2))

out2%>%
  map(simulation_summary,
      models=c('proposed','proposed IV equivalent','pooled ols','IV (small)', 'Liml equivalent'),
      est_valname='beta2.est',std_valname='beta2.est.std', 4)%>%
  bind_rows()%>%
  mutate(ci_covprob = paste(round(ci_covprob*100,1),'\\%',sep = ''))%>%
  dplyr::select(samplesize_per_trial,no_trials,Estimator,bias,ci_covprob)%>%
  group_by(samplesize_per_trial)%>%
  group_split()%>%
  bind_cols()%>%
  dplyr::select(-matches('Estimator.'),
                -starts_with('no_trials'),
                -starts_with('samplesize_per_trial'))%>%
  mutate(Estimator = factor(Estimator, 
                            labels = c('LIML','CMMA','Full Sample 2SLS','\\citet{Sobel2008IdentificationVariables}' ,"LSEM \\citep{Baron1986}"),
                            levels=c('Liml equivalent','proposed', 'proposed IV equivalent', 'IV (small)','pooled ols')))%>%
  arrange(Estimator)%>%
  kable(format="latex",digits = 3,booktabs = T,
        col.names = NULL,escape = F)%>%
  add_header_above(escape = F,line_sep = 0, c("Estimator" = 1, 
                                "Bias"=1,"\\\\makecell{ 95\\\\% CI \\\\\\\\coverage }"=1, 
                                "Bias"=1,"\\\\makecell{ 95\\\\% CI \\\\\\\\coverage }"=1, 
                                "Bias"=1,"\\\\makecell{ 95\\\\% CI \\\\\\\\coverage }"=1))%>%
  add_header_above(escape = F,c(" " = 1, "$N_{per}=200$" = 2,"$N_{per}=500$" = 2,"$N_{per}=1000$"=2))


# Increase number of trials
out1 <- lapply(c(10,20, 50, 100), function(n) mc_model_main(simu_model_1,500,n))
out1%>%  
  map(simulation_summary,
      c('proposed','pooled ols'),
             est_valname='beta2.est',std_valname='beta2.est.std', 4)%>%
  bind_rows()%>%
  group_by(Estimator)%>%
  group_split()%>%
  rev()%>%
  bind_cols()%>%
  dplyr::select(samplesize_per_trial,no_trials,bias,RMSE,ci_covprob, bias1,RMSE1,ci_covprob1)%>%
  kable(format="latex",digits = 2,booktabs = T,
        col.names = c("Sample size per trial","No. of trials",
                      "Bias","RMSE","\makecell{ 95\% CI \\coverage }", 
                      "Bias","RMSE","\makecell{ 95\% CI \\coverage }"))%>%
  add_header_above(c(" " = 2, "Proposed estimator" = 3, "Pooled OLS estimator" = 3))

# Table output
# ## Assumption Violation
out3 <- lapply(c(200,500,1000), 
               function(n) mc_model_main(simu_model_1,n,50, proposed_estimator_only = T))%>%
  map(mutate,violation='None')
out3_A2 <-  lapply(c(200,500,1000), 
                   function(n) mc_model_main(simu_model_1,n,50, proposed_estimator_only = T,  violate_A2 =T))%>%
  map(mutate,violation='A2')
out3_A3 <-  lapply(c(200,500,1000), 
                   function(n) mc_model_main(simu_model_1,n,50, proposed_estimator_only = T,  violate_A3 =T))%>%
  map(mutate,violation='A3')
out3_A4 <-  lapply(c(200,500,1000), 
                   function(n) mc_model_main(simu_model_1,n,50, proposed_estimator_only = T,  violate_A4 =T))%>%
  map(mutate,violation='A4')


lapply(list(out3, out3_A2, out3_A3, out3_A4), function(data){
  violation_str <- data[[1]]%>%select(violation)%>%pluck(1,1)
  data%>%map(simulation_summary,
             models=c('proposed'),
             est_valname='beta2.est',std_valname='beta2.est.std', 4)%>%
    bind_rows()%>%
    mutate(ci_covprob = paste(round(ci_covprob*100,1),'%',sep = ''),
           violation=!!violation_str)
}) %>% 
  bind_rows()%>%
  dplyr::select(samplesize_per_trial,no_trials,violation,bias,ci_covprob)%>%
  group_by(samplesize_per_trial)%>%
  group_split()%>%
  bind_cols()%>%
  dplyr::select(-matches('violation.'),
                -starts_with('no_trials'),
                -starts_with('samplesize_per_trial'))%>%
  mutate(violation = factor(violation,
                            levels=c('None','A2', 'A3', 'A4')))%>%
  arrange(violation)%>%
  kable(format="latex",digits = 3,booktabs = T,
        col.names = c("Estimator",
                      "Bias","\\makecell{ 95% CI \\\\coverage }", 
                      "Bias","\\makecell{ 95% CI \\\\coverage }",
                      "Bias","\\makecell{ 95% CI \\\\coverage }"))%>%
  add_header_above(escape = F,c(" " = 1, "$N_{per}=200$" = 2,"$N_{per}=500$" = 2,"$N_{per}=1000$"=2))


# Wald test table
## mu(m) = beta2*m
out_model1 <-  mc_model_main(simu_model_1,1000,100, DRF_degree=0, wald_test= TRUE)
## mu(m) = beta2*m + beta3*m^2
out_model2 <-  mc_model_main(simu_model_2,1000,100,DRF_degree=0, wald_test= TRUE)
## mu(m) = beta2*m + 0 + beta4*m^3
out_model3 <-  mc_model_main(simu_model_3,1000,100, DRF_degree=0, wald_test= TRUE)

bind_rows(
  out_model1%>%
  wald_test_summary%>%
  dplyr::select(`p value, m3`,`p value, m2 m3`,`p value, m m2 m3`)%>%
  mutate(`$\\mu(m)$`='$4m$',
         degree = 1,
         estimate=out_model1 %>%
           get_avg_simulated_estimate('proposed',c('beta2.est'))%>%
           as_vector()%>%round(3)%>%paste(collapse=', ')),
  out_model2%>%wald_test_summary%>%
  dplyr::select(`p value, m3`,`p value, m2 m3`,`p value, m m2 m3`)%>%
  mutate(`$\\mu(m)$`='$4m + 2m^2$',
         degree = 2,
         estimate=out_model2 %>%
           get_avg_simulated_estimate('proposed with m2',c('beta2.est','beta3.est'))%>%
           as_vector()%>%round(3)%>%paste(collapse=', ')),
  out_model3%>%wald_test_summary%>%
  dplyr::select(`p value, m3`,`p value, m2 m3`, `p value, m m2 m3`)%>%
  mutate(`$\\mu(m)$`='$4m + 5m^3$',
         degree = 3,
         estimate=out_model3 %>%
           get_avg_simulated_estimate('proposed with m2 m3',c('beta2.est','beta3.est', 'beta4.est'))%>%
           as_vector()%>%round(3)%>%paste(collapse=', ')))%>%
  dplyr::select(`$\\mu(m)$`,everything())%>%
  kable(format="latex",digits = 2,booktabs = T,escape=FALSE,
        col.names = c("$\\mu(m)$","$H_0: \\beta_3 = 0$","$H_0: \\beta_2 = \\beta_3 = 0$",
                      "$H_0: \\beta_1=\\beta_2=\\beta_3= 0$", "Highest degree used",
                      "Results"))%>%
  add_header_above(c(" "=1, "Percentage of wald test rejected null hypothesis (95% confidence)" = 3,"Estimation"=2))


# Bias vs sample size per trial
grid <- seq(2,4,0.25)
out2 <- lapply(as.integer(10^(grid)), function(n) mc_model_main(m_model, y_model,n,50))
df <- out2%>%bind_rows()
df%>%write_rds('bias_vs_sample_size.rds')
ggplot(data=df,aes(x=log(samplesize_per_trial),y=bias,color=Estimator)) +
  geom_point()+geom_line()


# Draw graph
## Simple case, single treatment group
sim_data <- simulate_data(simu_model_1, 400,20, 
                          sd.y000 = 6,
                          num_treatment_group = 1,
                          mu.gamma = c(3))
true_dose <- function(m) 4*m
meta_data <- sim_data%>%
  group_by(trial_index,treatment_group)%>%
  do(regress_one_trial(.))%>%
  ungroup
samp_df <- sim_data%>%sample_n(100)
draw_df <- bind_rows(samp_df%>%filter(Ti==F)%>%
                       mutate(id= row_number())%>%
                       mutate(Mi=range(samp_df%>%filter(Ti==F)%>%dplyr::select(Mi))[1],
                              Yi = simu_model_1(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta)),
                     samp_df%>%filter(Ti==F)%>%
                       mutate(id= row_number())%>%
                       mutate(Mi=range(samp_df%>%filter(Ti==F)%>%dplyr::select(Mi))[2],
                              Yi = simu_model_1(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta)))

ggplot(data = samp_df%>%filter(Ti==F), aes(x=Mi, y=Yi)) +
  geom_smooth(data =draw_df, 
              aes(x=Mi, y=Yi,group=id), 
              method='lm',se = F,
              color='grey')+
  geom_point()+
  geom_point(data=meta_data,aes(x=delta_m,y=delta_y),color='red')+  
  geom_smooth(data=meta_data,aes(x=delta_m,y=delta_y),color='red',method='lm',se = F)+
  geom_smooth(method='lm',se = F)+
  stat_function(fun=true_dose)+
  labs(x='M',y='Y')

ggplot(data = samp_df%>%filter(Ti==F), aes(x=Mi, y=Yi)) +
  geom_smooth(data =draw_df, 
              aes(x=Mi, y=Yi,group=id), 
              method='lm',se = F,
              color='grey',
              size =0.5)+
  geom_point(size =0.5)+
  geom_smooth(method='lm',se = F, color = 'black', size =0.5)+
  stat_function(fun=true_dose, color='blue', size=0.7)+
  labs(x='M',y='Y')+
  theme_few()

draw_df <- bind_rows(samp_df%>%
                       mutate(id= row_number())%>%
                       mutate(Mi=range(samp_df%>%dplyr::select(Mi))[1],
                              Yi = simu_model_1(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta)),
                     samp_df%>%
                       mutate(id= row_number())%>%
                       mutate(Mi=range(samp_df%>%dplyr::select(Mi))[2],
                              Yi = simu_model_1(theta_s, gamma_s, Ti, Mi, y000, noise_gamma, noise_beta)))

fit <- lm(Yi ~ 1+ Mi + Ti, data=samp_df)
scale = 0.3
g <- ggplot(data = samp_df, aes(x=Mi, y=Yi)) +
  geom_smooth(data =draw_df, 
              aes(x=Mi, y=Yi,group=id), 
              method='lm',se = F,
              color='grey',
              alpha = 0.5,
              size =1.2* scale)+
  geom_point(size =1.8 * scale)+
  # geom_smooth(method='lm',formula = str(fit$call), se = F, color = 'black', size =0.2)+
  geom_smooth(method='loess', se = F, color = 'black', size =1.8* scale)+
  stat_function(fun=true_dose, color='blue', size=2*scale)+
  labs(x='M',y='Y')+
  theme_few()
ggsave('potential_outcomes.png',g,scale=scale)


## quadrtic, single treatment group
sim_data <- simulate_data(simu_model_2, 400,20, 
                          sd.m00 = 1,
                          sd.y000 = 6,
                          cor.m00.y000 = 0.95,
                          num_treatment_group = 1,
                          mu.gamma = c(3)
                        )
true_dose <- function(m) 4*m+2*m^2
meta_data <- sim_data%>%
  group_by(trial_index,treatment_group)%>%
  do(regress_one_trial(.))%>%
  ungroup

samp_df <- sim_data%>%sample_n(100)

m_range <- samp_df%>%filter(Ti==F)%>%dplyr::select(Mi)%>%range()
m_seq <- tibble(t = rep(1:20,200),
       m = rep(seq(m_range[1],m_range[2],length.out = 20),200))%>%
  group_by(t)%>%
  mutate(id= row_number())%>%
  ungroup

draw_df <- samp_df%>%filter(Ti==F)%>%
  mutate(id= row_number())%>%
  left_join(m_seq)%>%
  mutate(Yi = simu_model_2(theta_s, gamma_s, Ti, m, y000, noise_gamma, noise_beta))%>%
  dplyr::select(-Mi,Mi=m)

ggplot(data = samp_df%>%filter(Ti==F), aes(x=Mi, y=Yi)) +
  geom_smooth(data =draw_df, 
              aes(x=Mi, y=Yi,group=id), 
              method='loess',se = F,
              color='grey')+
  geom_point()+
  geom_smooth(method='loess',se = F)+
  stat_function(fun=true_dose)+
  labs(x='M',y='Y')

# Interaction between mediation and ov
m_model <- function(alpha1,Ti,eta){
  alpha0 <- 1
  Mi = alpha0 + alpha1*Ti + eta
  return(Mi)
}

y_model <- function(beta1, Ti, Mi,epsilon){
  beta0 <- 0.5
  beta2 <- 4
  Yi = beta0 + beta1*Ti + beta2*Mi + beta2*Mi*epsilon + epsilon
  return(Yi)
}
sim_data <- simulate_data(m_model,y_model, 10,20, sd.epsilon = 6)
true_dose <- function(m) 4*m
draw_df <- bind_rows(sim_data%>%filter(Ti==F)%>%
                       mutate(id= row_number())%>%
                       mutate(Mi=range(sim_data%>%filter(Ti==F)%>%dplyr::select(Mi))[1],
                              Yi = y_model(beta1, Ti, Mi,epsilon)),
                     sim_data%>%filter(Ti==F)%>%
                       mutate(id= row_number())%>%
                       mutate(Mi=range(sim_data%>%filter(Ti==F)%>%dplyr::select(Mi))[2],
                              Yi = y_model(beta1, Ti, Mi,epsilon)))

ggplot(data = sim_data%>%filter(Ti==F), aes(x=Mi, y=Yi)) +
  geom_smooth(data =draw_df, 
              aes(x=Mi, y=Yi,group=id), 
              method='lm',se = F,
              color='grey')+
  geom_point()+
  geom_smooth(method='lm',se = F)+
  stat_function(fun=true_dose)