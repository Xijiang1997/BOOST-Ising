# load packages

library(mclust)
library(lattice)

source("src/functions_potts.R");
Rcpp::sourceCpp('src/functions_potts_omega.cpp');

# functions
normalization2 <- function(counts) {
  varx = apply(counts, 2, var)
  meanx = apply(counts, 2, mean)
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
  
  ## regress out log total counts
  norm_counts <- log(counts + 1/(2 * phi))
  total_counts <- apply(counts, 1, sum)
  log_total_counts <- log(total_counts)

  coeff <- apply(norm_counts, 2, function(x){ x - log_total_counts} )
 
  return(coeff)
}

CombinePValues <- function(pvalues, weights=NULL){
  if(!is.matrix(pvalues)){pvalues <- as.matrix(pvalues)}
  ## to avoid extremely values
  pvalues[which(pvalues==0)] <- 5.55e-17
  pvalues[which((1-pvalues)<1e-3)] <- 0.99
 
  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  if(is.null(weights)){
    weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  }# end fi
  if( (nrow(weights) != num_gene) || (ncol(weights) != num_pval)){
    stop("the dimensions of weights does not match that of combined pvalues")
  }# end fi
 
  Cstat <- tan((0.5 - pvalues)*pi)
 
  wCstat <- weights*Cstat
  Cbar <- apply(wCstat, 1, sum)
  #combined_pval <- 1.0/2.0 - atan(Cbar)/pi
  combined_pval <- 1.0 - pcauchy(Cbar) 
  combined_pval[which(combined_pval <= 0)] <- 5.55e-17
  return(combined_pval)
}# end func


transform_data <- function(count, loc){
  array_y <- array(0, dim = c((max(loc[,1])-min(loc[,1])+1),(max(loc[,2])-min(loc[,2])+1),ncol(count)))
  loc1 <- min(loc[,1])-1
  loc2 <- min(loc[,2])-1
  for (i in 1:ncol(count)){
    for (j in 1:nrow(loc)){
      array_y[loc[j,1]-loc1,loc[j,2]-loc2, i] <- count[j,i]
    }
  }
  return(array_y)
}


# filter data
filter_count <- function(count, sample_info, min_total = 10,min_percentage = 0.1){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  sample_f <- sample_info[-which(rowSums(count)<min_total),]
  count_f <- count[-which(rowSums(count)<min_total),-(which(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num))]
  
  return(list(sample_f,count_f))
}

# main function

Boost_HP_4chain <- function(count,sample_info, normalization = 1, clustermethod = 'Kmeans', dout = 3,  sigma = 1, prior = 'N', prior_omega = 'noninfo')
{
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
# normalization
  
  if (normalization == 1)
  {
    count_nor <- matrix(0, nrow = sample_num, ncol = gene_num)
    for (i in 1:sample_num){
      count_nor[i,] <- 10000*count[i,]/count_rowsum[i]
    }

  }
  
  else if (normalization == 2)
  {
    count_nor <- normalization2(count)

  }
  
  else if (normalization == 3)
  {
    count_nor <- matrix(0, nrow = sample_num, ncol = gene_num)
    for (i in 1:sample_num){
      count_q <- quantile(count[i,],0.75)
      count_nor[i,] <- count[i,]/count_q
    }
  }
  
  # clustering
  
  if (clustermethod == 'MGC'){
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
    mIQR <- IQR(count_nor[,i])
    list_large <- which(count_nor[,i] > median(count_nor[,i]) + dout*mIQR)
    list_small <- NULL

    if (mIQR == 0){list_large <- NULL}
    list_remain <- setdiff(1:sample_num, list_large)
      if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
      {
        count_binary[list_remain,i] <- 1
      }
      else
      {
      k <- Mclust(count_nor[list_remain,i], G=2)
      if (k$parameters$mean[1]>k$parameters$mean[2])
      {count_binary[list_remain,i] <- 3 - k$classification}
      else {count_binary[list_remain,i] <- k$classification}
    }
      count_binary[list_large,i] <- 2
      
      # If cannot split into two groups, we assume that there is no spatial pattern
      if (sum(count_binary[,i]) == 2*sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
      if (sum(count_binary[,i]) == sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
    }
  }
    
  # Kmeans clustering
  
  if (clustermethod == 'Kmeans'){
   
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
      mIQR <- IQR(count_nor[,i])
    list_large <- which(count_nor[,i] > median(count_nor[,i]) + dout*mIQR)
    list_small <- NULL

    if (mIQR == 0){list_large <- NULL}
    list_remain <- setdiff(1:sample_num, list_large)
    
      if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
      {
        count_binary[list_remain,i] <- 1
      }
      else{
      k <- kmeans(count_nor[list_remain,i], 2, nstart = 25)
      if (k$centers[1]>k$centers[2])
      {count_binary[list_remain,i] <- 3 - k$cluster}
      else {count_binary[list_remain,i] <- k$cluster}}
      count_binary[list_large,i] <- 2

      count_binary[list_small,i] <- 1
    }

    if (max(count_binary[,i]) == 1 )
    {
      count_binary[sample(1:sample_num, 100), i] <- 2
    }
    else if (min(count_binary[,i]) == 2 )
    {
      count_binary[sample(1:sample_num, 100), i] <- 1
    }
  }
  
  # transform to array
  P_nor <- transform_data(count_binary,sample_info)
  
  # estimate parameters in Potts model
  theta_est <- array(0,dim=c(5000, gene_num, 4))
  omega_est <- array(0,dim=c(5000, gene_num, 4))
  
  theta_CI_low <- numeric(gene_num)
  theta_CI_high <- numeric(gene_num)
  theta_mean <- numeric(gene_num)
  omega_mean <- numeric(gene_num)
  omega_CI_low <- numeric(gene_num)
  omega_CI_high <- numeric(gene_num)

  accept_theta <- numeric(gene_num)
  accept_omega <- numeric(gene_num)
  
  potential_R_theta <- numeric(gene_num) + 2
  potential_R_omega <- numeric(gene_num) 

  iter <- numeric(gene_num)

  for (i in 1:gene_num){
    if (prior_omega == 'noninfo')
    {
      mean1 = 1
      sigma1 = 2.5
      is_normal = 1
    }
    else if (prior_omega == 'info') {
      propor <- sum(P_nor[,,i] == 1)/sum(P_nor[,,i] != 0)
      mean1 = log((propor/(1-propor))*exp(1))
      sigma1 = 0.1
      is_normal = 1
    }
    else if (prior_omega == 'uniform')
    {
      is_normal = 0
      mean1 = 1
      sigma1 = 2.5
    }

    resu <- TRUE
  
    while (resu){

    if (prior == 'N'){
    res1 <- potts_2_omega(P_nor[,,i], sigma, mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000);
    res2 <- potts_2_omega(P_nor[,,i], sigma, mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000);
    res3 <- potts_2_omega(P_nor[,,i], sigma, mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000);
    res4 <- potts_2_omega(P_nor[,,i], sigma, mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000);
  }
    else if (prior == 'L'){
      res1 <- potts_2_omega_laplace(P_nor[,,i], sigma/sqrt(2), mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000)
      res2 <- potts_2_omega_laplace(P_nor[,,i], sigma/sqrt(2), mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000)
      res3 <- potts_2_omega_laplace(P_nor[,,i], sigma/sqrt(2), mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000)
      res4 <- potts_2_omega_laplace(P_nor[,,i], sigma/sqrt(2), mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 10000)
    }
    
    if (dim(res1$theta)[2]+dim(res2$theta)[2] + dim(res3$theta)[2] + dim(res4$theta)[2] == 4 & dim(res1$omega)[2] != 0 & dim(res2$omega)[2] != 0 & dim(res3$omega)[2] != 0 & dim(res4$omega)[2] != 0  ){
    theta_est[,i, 1] <- res1$theta[5001:10000,1]
    omega_est[,i, 1] <- res1$omega[5001:10000,1]

    theta_est[,i, 2] <- res2$theta[5001:10000,1]
    omega_est[,i, 2] <- res2$omega[5001:10000,1]

    theta_est[,i, 3] <- res3$theta[5001:10000,1]
    omega_est[,i, 3] <- res3$omega[5001:10000,1]

    theta_est[,i, 4] <- res4$theta[5001:10000,1]
    omega_est[,i, 4] <- res4$omega[5001:10000,1]
  }
    
    #theta_CI_low[i] <- quantile(theta_est[,i, 1],0.025)
    #theta_CI_high[i] <- quantile(theta_est[,i,1],0.975)
    
    #omega_CI_low[i] <- quantile(omega_est[,i,1],0.025)
    #omega_CI_high[i] <- quantile(omega_est[,i,1],0.975)
    accept_theta[i] <- 0.25*(res1$accept + res2$accept + res3$accept + res4$accept)
    accept_omega[i] <- 0.25*(res1$accept_omega + res2$accept_omega + res3$accept_omega + res4$accept_omega)

    w <- 0.25*(var(theta_est[,i, 1]) + var(theta_est[,i, 2]) + var(theta_est[,i, 3]) + var(theta_est[,i, 4])) + 0.0000000001
    b <- 5000/3 * ((mean(theta_est[,i, 1]) - mean(theta_est[,i,]))^2 + (mean(theta_est[,i, 2]) - mean(theta_est[,i,]))^2 + (mean(theta_est[,i, 3]) - mean(theta_est[,i,]))^2 + (mean(theta_est[,i, 4]) - mean(theta_est[,i,]))^2)
    potential_R_theta[i] <- sqrt((4999/5000 * w + b/5000)/w)

    w2 <- 0.25*(var(omega_est[,i, 1]) + var(omega_est[,i, 2]) + var(omega_est[,i, 3]) + var(omega_est[,i, 4])) + 0.0000000001
    b2 <- 5000/3 * ((mean(omega_est[,i, 1]) - mean(omega_est[,i,]))^2 + (mean(omega_est[,i, 2]) - mean(omega_est[,i,]))^2 + (mean(omega_est[,i, 3]) - mean(omega_est[,i,]))^2 + (mean(omega_est[,i, 4]) - mean(omega_est[,i,]))^2)
    potential_R_omega[i] <- sqrt((4999/5000 * w2 + b2/5000)/w2)

    iter[i] <- iter[i] + 1

    resu <- potential_R_theta[i] > 1.1 && iter[i] < 3

    if (is.logical(resu) == FALSE && iter[i] < 3)
    {resu <- TRUE}
    else if (is.logical(resu) == FALSE && iter[i] >= 3) {resu <- FALSE}
  }}
  
  theta_mean <- 0.25*(apply(theta_est[,,1], 2, mean) + apply(theta_est[,,2], 2, mean) + apply(theta_est[,,3], 2, mean) + apply(theta_est[,,4], 2, mean))
  omega_mean <- 0.25*(apply(omega_est[,,1], 2, mean) + apply(omega_est[,,2], 2, mean) + apply(omega_est[,,3], 2, mean) + apply(omega_est[,,4], 2, mean))
  
  # detect SV genes
  pvalues_neg <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_neg[i] <- sum(theta_est[,i,]>=0)/20000
  }
  pvalues_neg_ad <- p.adjust(pvalues_neg, "BH")
  
  pvalues_cauchy <- numeric(gene_num)
  pp <- matrix(0, nrow = gene_num, ncol = 4)
   for (i in 1:gene_num){
    pp[i, 1] <- sum(theta_est[,i,1]>=0)/5000
    pp[i, 2] <- sum(theta_est[,i,2]>=0)/5000
    pp[i,3] <- sum(theta_est[,i,3]>=0)/5000
    pp[i,4]<- sum(theta_est[,i,4]>=0)/5000
  }
 pvalues_cauchy <- CombinePValues(pp)

  pvalues_pos <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_pos[i] <- sum(theta_est[,i,]<=0)/20000
  }
  pvalues_pos_ad <- p.adjust(pvalues_pos, "BH")

  pvalues_double <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_double[i] <- 2 * min(sum(theta_est[,i,]<=0),sum(theta_est[,i,]>=0)) /20000
  }
  pvalues_double_ad <- p.adjust(pvalues_double, "BH")

  
  results <- data.frame(theta_mean, theta_CI_low, theta_CI_high, omega_mean,omega_CI_low,omega_CI_high, pvalues_neg,pvalues_cauchy, pvalues_neg_ad, pvalues_pos, pvalues_pos_ad, pvalues_double, pvalues_double_ad, accept_theta,accept_omega, potential_R_theta, potential_R_omega, iter)
  rownames(results) <- colnames(count)
  return(results)
}  


Boost_HP <- function(count,sample_info, normalization = 1, clustermethod = 'Kmeans', dout = 3,  sigma = 1, prior = 'N', prior_omega = 'noninfo')
{
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
# normalization
  
  if (normalization == 1)
  {
    count_nor <- matrix(0, nrow = sample_num, ncol = gene_num)
    for (i in 1:sample_num){
      count_nor[i,] <- 10000*count[i,]/count_rowsum[i]
    }

  }
  
  else if (normalization == 2)
  {
    count_nor <- normalization2(count)

  }
  
  else if (normalization == 3)
  {
    count_nor <- matrix(0, nrow = sample_num, ncol = gene_num)
    for (i in 1:sample_num){
      count_q <- quantile(count[i,],0.75)
      count_nor[i,] <- count[i,]/count_q
    }
  }
  
  # clustering
  
  if (clustermethod == 'MGC'){
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
    mIQR <- IQR(count_nor[,i])
    list_large <- which(count_nor[,i] > median(count_nor[,i]) + dout*mIQR)
    list_small <- NULL

    if (mIQR == 0){list_large <- NULL}
    list_remain <- setdiff(1:sample_num, list_large)
      if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
      {
        count_binary[list_remain,i] <- 1
      }
      else
      {
      k <- Mclust(count_nor[list_remain,i], G=2)
      if (k$parameters$mean[1]>k$parameters$mean[2])
      {count_binary[list_remain,i] <- 3 - k$classification}
      else {count_binary[list_remain,i] <- k$classification}
    }
      count_binary[list_large,i] <- 2
      
      # If cannot split into two groups, we assume that there is no spatial pattern
      if (sum(count_binary[,i]) == 2*sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
      if (sum(count_binary[,i]) == sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
    }
  }
    
  # Kmeans clustering
  
  if (clustermethod == 'Kmeans'){
   
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
      mIQR <- IQR(count_nor[,i])
    list_large <- which(count_nor[,i] > median(count_nor[,i]) + dout*mIQR)
    list_small <- NULL

    if (mIQR == 0){list_large <- NULL}
    list_remain <- setdiff(1:sample_num, list_large)
    
      if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
      {
        count_binary[list_remain,i] <- 1
      }
      else{
      k <- kmeans(count_nor[list_remain,i], 2, nstart = 25)
      if (k$centers[1]>k$centers[2])
      {count_binary[list_remain,i] <- 3 - k$cluster}
      else {count_binary[list_remain,i] <- k$cluster}}
      count_binary[list_large,i] <- 2

      count_binary[list_small,i] <- 1
    }

    if (max(count_binary[,i]) == 1 )
    {
      count_binary[sample(1:sample_num, 100), i] <- 2
    }
    else if (min(count_binary[,i]) == 2 )
    {
      count_binary[sample(1:sample_num, 100), i] <- 1
    }
  }
  
  # transform to array
  P_nor <- transform_data(count_binary,sample_info)
  
  # estimate parameters in Potts model
  theta_est <- array(0,dim=c(10000, gene_num))
  omega_est <- array(0,dim=c(10000, gene_num))
  
  theta_CI_low <- numeric(gene_num)
  theta_CI_high <- numeric(gene_num)
  
  theta_mean <- numeric(gene_num)
  omega_mean <- numeric(gene_num)
  
  omega_CI_low <- numeric(gene_num)
  omega_CI_high <- numeric(gene_num)

  accept_theta <- numeric(gene_num)
  accept_omega <- numeric(gene_num)
  
  iter <- numeric(gene_num)

  for (i in 1:gene_num){
    if (prior_omega == 'noninfo')
    {
      mean1 = 1
      sigma1 = 2.5
      is_normal = 1
    }
    else if (prior_omega == 'info') {
      propor <- sum(P_nor[,,i] == 1)/sum(P_nor[,,i] != 0)
      mean1 = log((propor/(1-propor))*exp(1))
      sigma1 = 0.1
      is_normal = 1
    }
    else if (prior_omega == 'uniform')
    {
      is_normal = 0
      mean1 = 1
      sigma1 = 2.5
    }

    if (prior == 'N'){
    res <- potts_2_omega(P_nor[,,i], sigma, mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 20000);
    }
    else if (prior == 'L'){
      res <- potts_2_omega_laplace(P_nor[,,i], sigma/sqrt(2), mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5), 20000)
    }
    
    
    theta_est[,i] <- res$theta[10001:20000,1]
    omega_est[,i] <- res$omega[10001:20000,1]
    
    theta_CI_low[i] <- quantile(theta_est[,i],0.025)
    theta_CI_high[i] <- quantile(theta_est[,i],0.975)
    
    omega_CI_low[i] <- quantile(omega_est[,i],0.025)
    omega_CI_high[i] <- quantile(omega_est[,i],0.975)
    accept_theta[i] <- res$accept 
    accept_omega[i] <- res$accept_omega 

    theta_mean[i] <- mean(res$theta[10001:20000,1])
    omega_mean[i] <- mean(res$omega[10001:20000,1])
    }
  
  
  # detect SV genes
  pvalues_neg <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_neg[i] <- sum(theta_est[,i]>=0)/10000
  }
  pvalues_neg_ad <- p.adjust(pvalues_neg, "BH")

  pvalues_pos <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_pos[i] <- sum(theta_est[,i]<=0)/20000
  }
  pvalues_pos_ad <- p.adjust(pvalues_pos, "BH")
  
  results <- data.frame(theta_mean, theta_CI_low, theta_CI_high, omega_mean,omega_CI_low,omega_CI_high, pvalues_neg,pvalues_neg_ad, pvalues_pos, pvalues_pos_ad, accept_theta,accept_omega)
  rownames(results) <- colnames(count)
  return(results)
}  



