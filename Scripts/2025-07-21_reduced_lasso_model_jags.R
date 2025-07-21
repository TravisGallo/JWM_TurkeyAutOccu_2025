model{
  
  for(site in 1:nsite){
   
    # First Season
    
    # latent state model
    logit(psi[site,1]) <- psi0 + inprod(psi_beta, X[site,])
    z[site,1] ~ dbern(psi[site,1])
    
    # data model
    for(occasion in 1:noccasion){
    logit(p[site,occasion,1]) <- a0 + a1 * A_a[site,occasion,1]
    y[site, occasion,1] ~ dbern(z[site,1] * p[site,occasion,1])
    }
    
    # Remaining seasons
    for(season in 2:nseason){
      
      # latent state model
      logit(psi[site,season]) <- psi0 + inprod(psi_beta, X[site,]) +
        theta * z[site,season-1]
      z[site,season] ~ dbern(psi[site,season])
      
      # data model
      for(occasion in 1:noccasion){
        logit(p[site,occasion,season]) <- a0 + a1 * A_a[site,occasion,season]
        y[site,occasion,season] ~ dbern(z[site,season] * p[site,occasion,season])
      }
    }
  }
  
  ## Priors
  
  # Intercepts
  psi0 ~ dlogis(0,1)
  a0 ~ dlogis(0,1)
  
  # # Slope term for occupancy
  # # Use spike and slab for variable selection
  for(covar in 1:ncovar){
    #ridge[covar] ~ dnorm(0, lambda)
    psi_beta[covar] ~ ddexp(0, lambda) # lasso regression
  }
  
  # scaling parameter for double exponential priors
  # smaller mean more mass at zero
  lambda ~ dgamma(0.001,.001)
  
  # slope term for detection
  a1 ~ dlogis(0,1)
  
  # First-order autoregressive term
  theta ~ dlogis(0,1)
}
