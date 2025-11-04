### Code to replicate example ####
source("R/helper.R") 
require(fastDummies)
require(dplyr)
require(tidyr)

# read in data ----
dat <- read.csv(file = "data/rg_dat.csv")

## validation data ----
# from https://www.sciencedirect.com/science/article/abs/pii/S1386653212000893

valdat <- data.frame(ystar = c(1, 0, 1, 0),
                     y = c(1, 0, 0, 1),
                     n = c(208, 31, 1, 13)) |> 
          uncount(n)

# complete case analysis ----
main <- dat %>% filter(!is.na(ystar))

## compute se and sp ----
se <- mean(valdat[valdat$y==1 ,]$ystar)
sp <- 1 - mean(valdat[valdat$y==0 ,]$ystar)


## the RG ----
pstar <- mean(main$ystar, na.rm = T)
sehat <- sqrt(pstar*(1 - pstar)/(nrow(main)))

adjprev <- (pstar + sp - 1)/(se+sp-1)
adjprev

## the var ---
rgvar <- function(pstar, p, alpha, beta, n, r, m){
  var <- (pstar * (1 - pstar))/(n * (alpha + beta - 1)^2) +
    (alpha*(1 - alpha)*p^2)/(m * ((alpha + beta - 1)^2)) +
    (beta * (1 - beta) * (1 - p)^2)/(r * ((alpha + beta - 1)^2))
}

var <- rgvar(pstar, adjprev, se, sp, nrow(main), sum(valdat$y==0), sum(valdat$y==1))
var
stderr <- sqrt(var)
stderr

## repeat using m-estimation ----

ystar_main <- main$ystar
ystar_1 <- valdat[valdat$y == 1, ]$ystar
ystar_2 <- valdat[valdat$y == 0, ]$ystar
ystar <- c(ystar_main, ystar_1, ystar_2)

v <- c(rep(0, nrow(main)), 
       rep(1, length(ystar_1)), 
       rep(2, length(ystar_2)))
init_theta <- c(0.9, 0.9, 0.5, 0.5)
m_est <- mestimator(ef, init = init_theta, ystar_ = ystar, v_ = v)

phat <- m_est[[1]][4,]
phat


# accounting for missing data and misclassification ----

# fit ipw
den <- glm(r ~  factor(edu) + drinking + selfreportsw + age_g30 + sti, data = dat, 
           family = "binomial"(link = "logit"))

# compute weights
main$weight <- 1/predict(den, newdata = main, type = "response")

# plug in to RG
pstar_w <- weighted.mean(main$ystar, main$weight)
adjprev <- (pstar_w + sp - 1)/(se+sp-1)
adjprev


## compute variance via m-estimation----

ystar_main <- dat$ystar
ystar_1 <- valdat[valdat$y == 1, ]$ystar
ystar_2 <- valdat[valdat$y == 0, ]$ystar
ystar <- c(ystar_main, ystar_1, ystar_2)

v <- c(rep(0, length(ystar_main)), 
       rep(1, length(ystar_1)), 
       rep(2, length(ystar_2)))

dat2 <- fastDummies::dummy_cols(dat, select_columns = "edu")
covs <- as.matrix(dat2 %>% 
                    select(age_g30, selfreportsw, drinking, edu_2, edu_3, edu_4, edu_5, edu_6, sti)) 

init_theta <- c(0.9, 0.9, 0.05, 0.05, rep(0.5, ncol(covs)+1))
m_est_withmissing <- mestimator(ef_m, init = init_theta, ystar_ = ystar, v_ = v, xmat = covs)
m_est_withmissing[[1]]
phat2 <- m_est_withmissing[[1]][4,]
phat2

