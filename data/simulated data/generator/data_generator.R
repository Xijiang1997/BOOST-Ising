library(ggplot2)
library(reshape2)
library(gridExtra)

simu_data_discrete <- function(loc, pattern, sigma, beta1, beta0, phi_mean = 10, p = 100, p_gamma = 10, prop_zero = 0.4, seed = NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  n <- nrow(loc)
  s <- exp(rnorm(n, mean = 0, sd = 0.2))
  
  # Generate simulated data
  gamma <- rep(FALSE, p)
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)
  phi <- rexp(p, 1 / phi_mean)
  
  # low region
  loglambda <- beta0 + epsilon
  count_null <- matrix(rnbinom(n * p, 
                               mu = s * exp(loglambda), 
                               size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                       nrow = n, ncol = p)
  gamma[sample(1:p, p_gamma)] <- TRUE
  
  # high region
  loglambda[pattern, gamma] <- loglambda[pattern, gamma] + (beta1 - beta0)
  count <- matrix(rnbinom(n * p, 
                          mu = s * exp(loglambda), 
                          size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                  nrow = n, ncol = p)
  
  # zero index
  H <- matrix(0, nrow = n, ncol = p)
  if (prop_zero > 0) {
    H <- matrix(rbinom(n * p, 1, prop_zero), nrow = n, ncol = p)
    count_null[which(H == 1)] <- 0
    count[which(H == 1)] <- 0
  }
  
  return(list(count = count, loc = loc, gamma = gamma, count_null = count_null, zero_index = H, phi = phi))
}

simu_data_continous <- function(loc, pattern, sigma, beta1, beta0, phi_mean = 10, p = 100, p_gamma = 10, prop_zero = 0.4, seed = NA) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  n <- nrow(loc)
  s <- exp(rnorm(n, mean = 0, sd = 0.2))
  
  # Generate simulated data
  gamma <- rep(FALSE, p)
  epsilon <- matrix(rnorm(n * p, mean = 0, sd = sigma), nrow = n, ncol = p)
  phi <- rexp(p, 1 / phi_mean)
  
  # low region
  loglambda <- beta0 + epsilon
  count_null <- matrix(rnbinom(n * p, 
                               mu = s * exp(loglambda), 
                               size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                       nrow = n, ncol = p)
  gamma[sample(1:p, p_gamma)] <- TRUE
  
  # high region
  if (pattern == 'spot') {
    c.x <- mean(loc$x)
    c.y <- mean(loc$y)
    r <- min(diff(range(loc$x)), diff(range(loc$y))) / 3
    # add pattern
    d <- sqrt((loc$x - c.x)^2 + (loc$y - c.y)^2)
  } else if(pattern == 'linear') {
    c.x <- min(loc$x)
    c.y <- min(loc$y)
    r <- max(diff(range(loc$x)), diff(range(loc$y)))
    # add pattern
    d <- loc$x - c.x + loc$y - c.y
  }
  beta.diff <- (beta0 - beta1) / r * d + beta1 - beta0
  loglambda[, gamma] <- loglambda[, gamma] + sapply(beta.diff, max, 0)
  count <- matrix(rnbinom(n * p, 
                          mu = s * exp(loglambda), 
                          size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), 
                  nrow = n, ncol = p)
  
  # zero index
  H <- matrix(0, nrow = n, ncol = p)
  if (prop_zero > 0) {
    H <- matrix(rbinom(n * p, 1, prop_zero), nrow = n, ncol = p)
    count_null[which(H == 1)] <- 0
    count[which(H == 1)] <- 0
  }
  
  return(list(count = count, loc = loc, gamma = gamma, count_null = count_null, zero_index = H, phi = phi, s = s))
}


output_path <- "data/simulated data/";

# generate data with mob_i_pattern
pattern_name <- 'mob_i_pattern'
load("data/simulated data/generator/mob_patterns.Rdata")

n <- nrow(pattern_loc)
pattern <- patterns[, 2]

# generate data with discrete patterns
rpls <- 30
p <- 100
p_gamma <- 15
beta0 <- 2
phi_mean <- 10
zero_p <- c(0, 0.1, 0.3, 0.5)
sigma <- 0.3
beta1 <- beta0 + log(3)
  for (z.i in zero_p){
    for(r in 1:rpls){
      tmp <- simu_data_discrete(loc = pattern_loc, pattern = pattern, sigma = sigma, p = p, p_gamma = p_gamma, 
                                beta1 = beta1, beta0 = beta0, prop_zero = z.i, phi_mean = phi_mean, seed = 123 + r)
      count <- tmp$count
      loc <- as.matrix(tmp$loc)
      gamma <- tmp$gamma 
      s <- tmp$s
      parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean)     
      save(count, loc, gamma, s, parameters, file = paste0(output_path, pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
    }
  }


# generate data with mob_ii_pattern
pattern_name <- 'mob_ii_pattern'
pattern <- patterns[, 1]

rpls <- 30
p <- 100
p_gamma <- 15
beta0 <- 2
phi_mean <- 10
zero_p <- c(0, 0.1, 0.3, 0.5)
sigma <- 0.3
beta1 <- beta0 + log(3)

for (z.i in zero_p) {
  for(r in 1:rpls) {
    tmp <- simu_data_discrete(loc = pattern_loc, pattern = pattern, sigma = sigma, p = p, p_gamma = p_gamma, 
                              beta1 = beta1, beta0 = beta0, prop_zero = z.i, phi_mean = phi_mean, seed = 123 + r)
    count <- tmp$count
    loc <- as.matrix(tmp$loc)
    gamma <- tmp$gamma 
    s <- tmp$s
    parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean)     
    save(count, loc, gamma, s, parameters, file = paste0(output_path, pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
  }
}

# generate data with bc pattern
pattern_name <- 'bc_pattern'
load("data/simulated data/generator/bc_patterns.Rdata")

n <- nrow(pattern_loc)

# 30 replicates 
rpls <- 30
p <- 100
p_gamma <- 15
beta0 <- 2
phi_mean <- 10
zero_p <- c(0, 0.1, 0.3, 0.5)
sigma <- 0.3
beta1 <- beta0 + log(3)

for (z.i in zero_p) {
  for(r in 1:rpls) {
    tmp <- simu_data_discrete(loc = pattern_loc, pattern = pattern, sigma = sigma, p = p, p_gamma = p_gamma, 
                              beta1 = beta1, beta0 = beta0, prop_zero = z.i, phi_mean = phi_mean, seed = 123 + r)
    count <- tmp$count
    loc <- as.matrix(tmp$loc)
    gamma <- tmp$gamma 
    s <- tmp$s
    parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean)     
    save(count, loc, gamma, s, parameters, file = paste0(output_path, pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
  }
}

# generate data with linear pattern
pattern_name <- 'linear_pattern'
pattern_loc <- data.frame(x = rep(1:16, times = 16), y = rep(1:16, each = 16))

# 30 replicates 
rpls <- 30
p <- 100
p_gamma <- 15
beta0 <- 2
beta1 <- beta0 + log(6)
phi_mean <- 10
zero_p <- c(0, 0.1, 0.3, 0.5)
zero_names <- c('none','low', 'med', 'high')
sigma <- 0.3

for (z.i in zero_p) {
  for(r in 1:rpls) {
    tmp <- simu_data_continous(loc = pattern_loc, pattern = 'linear', sigma = sigma, p = p, p_gamma = p_gamma, 
                               beta1 = beta1, beta0 = beta0, prop_zero = z.i, phi_mean = phi_mean, seed = 123 + r)
    count <- tmp$count
    loc <- as.matrix(tmp$loc)
    gamma <- tmp$gamma 
    s <- tmp$s
    parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean)     
    save(count, loc, gamma, s, parameters, file = paste0(output_path, pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
  }
}

# generate data with spot pattern
pattern_name <- 'spot_pattern'
pattern_loc <- data.frame(x = rep(1:16, times = 16), y = rep(1:16, each = 16))

# 30 replicates 
rpls <- 30
p <- 100
p_gamma <- 15
beta0 <- 2
beta1 <- beta0 + log(6)
phi_mean <- 10
zero_p <- c(0, 0.1, 0.3, 0.5)
zero_names <- c('none','low', 'med', 'high')
sigma <- 0.3

for (z.i in zero_p) {
  for(r in 1:rpls) {
    tmp <- simu_data_continous(loc = pattern_loc, pattern = 'spot', sigma = sigma, p = p, p_gamma = p_gamma, 
                               beta1 = beta1, beta0 = beta0, prop_zero = z.i, phi_mean = phi_mean, seed = 123 + r)
    count <- tmp$count
    loc <- as.matrix(tmp$loc)
    gamma <- tmp$gamma 
    s <- tmp$s
    parameters <- list(sigma = sigma, p=p, p_gamma = p_gamma, beta1 = beta1, beta0=beta0, prop_zero = z.i, phi_mean=phi_mean)     
    save(count, loc, gamma, s, parameters, file = paste0(output_path, pattern_name, '_zero_', z.i*100, '_replicate_', r, '.RData'))
  }
}
