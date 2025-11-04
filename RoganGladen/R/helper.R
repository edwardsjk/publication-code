## helper functions for RG ##

# helper functions ----

#' m-estimator engine
#' @param ef an estimating function
#' @param init initial parameter values

mestimator <- function(ef, init, ...){
  
  require(rootSolve)
  
  # get estimated coefficients
  ef_colsums <- function(...) {
    return(colSums(ef(...)))
  }
  fit <- multiroot(ef_colsums, start = init,
                   rtol = 1e-6, atol = 1e-8, ctol = 1e-8, ...  )
  betahat <- fit$root
  
  # bread
  efhat <- ef(betahat, ...)
  n <- nrow(efhat)
  # derivative of estimating function at betahat
  pd <- gradient(ef_colsums, betahat, ...) 
  pd <- as.matrix(pd)
  bread <- -pd/n
  
  # meat1
  meat1 <- (t(efhat) %*% efhat)/n
  
  # sandwich
  sandwich1 <- (solve(bread) %*% meat1 %*% t(solve(bread)))/n
  se_sandwich1 <- sqrt(diag(sandwich1))
  
  results <- list(data.frame(estimates = betahat, stderr = se_sandwich1), covariance = sandwich1)
  return(results)
  
}


##' estimating function to account for misclassification

ef <- function(theta, ystar_, v_){
  alpha <- theta[1] 
  beta <- theta[2] 
  rhostar <- theta[3] 
  rho <- theta[4] 
  
  #efs
  ef1 <- as.numeric(v_ == 1) * (ystar_ - alpha)
  ef2 <- as.numeric(v_ == 2) * ((1 - ystar_) - beta)
  ef3 <- as.numeric(v_ == 0) * (ystar_ - rhostar)
  ef4 <- rho * (alpha + beta - 1) - (rhostar + beta - 1)
  
  stacked <- cbind(ef1, ef2, ef3, ef4)
  return(stacked)
  
}

##' estimating function to account for misclassification and missing data 
##' 
ef_m <- function(theta, ystar_, v_, xmat){
  alpha <- theta[1]
  beta <- theta[2] 
  rhostar <- theta[3]
  rho <- theta[4] 
  gamma <- theta[5:(5+ncol(xmat))]
  
  r <- as.numeric(!is.na(ystar_))
  ystar_ <- replace_na(ystar_, 0)
  x <- as.matrix(rbind(cbind(rep(1, nrow(xmat)), xmat), 
                       matrix(0, ncol = ncol(xmat) + 1, nrow = sum(v_ > 0))) 
                 # add matrix of 0s for rows in validation data
  )
  
  pi <- as.vector(plogis(x %*% gamma)) 
  
  #efs
  ef1 <- as.numeric(v_ == 1) * (ystar_ - alpha)
  ef2 <- as.numeric(v_ == 2) * ((1 - ystar_) - beta)
  ef3 <- as.numeric(v_ == 0) * (r) * (ystar_ - rhostar) * 1/pi
  ef4 <- rho * (alpha + beta - 1) - (rhostar + beta - 1)
  ef5 <- as.numeric(v_ == 0) * x * (r - pi)
  
  stacked <- cbind(ef1, ef2, ef3, ef4, ef5)
  return(stacked)
  
}


##' a function to implement the RG estimator
rgfunc <- function(pstar, se, sp){
  p <- (pstar + sp - 1)/(se+sp-1)
  return(p)
}

##' a function to compute sensitivity and specificity
getsesp <- function(v){
  se <- mean(v[v$y==1 ,]$ystar)
  sp <- 1 - mean(v[v$y==0 ,]$ystar)
  return(c(se, sp))
}

##' a function to account for misclassification and missing data
rgall <- function(dat, vd){
  sesp <- getsesp(vd)
  se <- sesp[1]
  sp <- sesp[2]
  # fit ipw
  den <- glm(r ~  factor(edu) + drinking + selfreportsw + age_g30 + sti, data = dat, 
             family = "binomial"(link = "logit"))
  
  # compute weights
  cc <- dat |> filter(!is.na(ystar))
  cc$weight <- 1/predict(den, newdata = cc, type = "response")
  
  # plug in to RG
  pstar_w <- weighted.mean(cc$ystar, cc$weight)
  p <- rgfunc(pstar_w, se, sp)
  return(p)
}



