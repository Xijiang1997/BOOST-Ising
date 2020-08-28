# load packages

library(mclust)
library(lattice)

source("functions/functions_potts.R");
Rcpp::sourceCpp('functions/functions_potts_omega.cpp');

# functions

normalization2 <- function(counts) {
  varx = apply(counts, 2, var)
  meanx = apply(counts, 2, mean)
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
  
  ## regress out log total counts
  norm_counts <- log(counts + 1/(2 * phi))
  total_counts <- apply(counts, 1, sum)
  log_total_counts <- log(total_counts)
  coeff <- apply(norm_counts, 2, function(x){coef(lm(x ~ log_total_counts))} )
  prod <- log_total_counts %*% t(as.matrix(coeff[2,]))
  res_norm_counts <- norm_counts - prod
  
  return(res_norm_counts)
}



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

Boost_HP <- function(count,sample_info, normalization = 2, clustermethod = 'Mclust', dout = 1.5)
{
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  sample_info <- round(sample_info)
  
  # normalization
  
  if (normalization == 1)
  {
    count_rowsum <- rowSums(count)
    count_nor <- matrix(0, nrow = sample_num, ncol = gene_num)
    for (i in 1:sample_num){
      count_nor[i,] <- 10000*count[i,]/count_rowsum[i]
    }
    mIQR <- IQR(count_nor[,i])
    if (mIQR != 0)
      {
      list_large <- which(count_nor[,i] > median(count_nor[,i]) + dout*mIQR)
      list_small <- NULL
      }
    else
      {
      list_small <- NULL
      list_large <- NULL
    }
  }
  
  else if (normalization == 2)
  {
    count_nor <- normalization2(count)
    mIQR <- IQR(count_nor[,i])
    list_large <- which(count_nor[,i] > median(count_nor[,i]) + dout*mIQR)
    list_small <- which(count_nor[,i] < median(count_nor[,i]) - dout*mIQR)
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
  
  if (clustermethod == 'Mclust'){
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
      list_remain <- setdiff(1:sample_num, union(list_large, list_small))
      k <- Mclust(count_nor[list_remain,i], G=2)
      if (k$parameters$mean[1]>k$parameters$mean[2])
      {count_binary[list_remain,i] <- 3 - k$classification}
      else {count_binary[list_remain,i] <- k$classification}
      count_binary[list_large,i] <- 2
      count_binary[list_small, i] <- 1
      
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
      list_remain <- setdiff(1:sample_num, union(list_large, list_small))
      k <- kmeans(count_nor[list_remain,i], 2, nstart = 25)
      if (k$centers[1]>k$centers[2])
      {count_binary[list_remain,i] <- 3 - k$cluster}
      else {count_binary[list_remain,i] <- k$cluster}
      count_binary[list_large,i] <- 2
      count_binary[list_small, i] <- 1
    }
    
  }
  
  # transform to array
  P_nor <- transform_data(count_binary,sample_info)
  
  # estimate parameters in Potts model
  theta_est <- matrix(0,ncol =gene_num, nrow = 5000)
  omega_est <- matrix(0,ncol =gene_num, nrow = 5000)
  
  theta_CI_low <- numeric(gene_num)
  theta_CI_high <- numeric(gene_num)
  
  omega_CI_low <- numeric(gene_num)
  omega_CI_high <- numeric(gene_num)
  for (i in 1:gene_num){
    res <- potts_2_omega(P_nor[,,i]);
    theta_est[,i] <- res$theta[5001:10000,1]
    omega_est[,i] <- res$omega[5001:10000,1]
    
    theta_CI_low[i] <- quantile(theta_est[,i],0.025)
    theta_CI_high[i] <- quantile(theta_est[,i],0.975)
    
    omega_CI_low[i] <- quantile(omega_est[,i],0.025)
    omega_CI_high[i] <- quantile(omega_est[,i],0.975)
  }
  
  theta_mean <- apply(theta_est, 2, mean)
  omega_mean <- apply(omega_est, 2, mean)
  
  # detect SV genes
  pvalues_neg <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_neg[i] <- sum(theta_est[,i]>=0)/5000
  }
  pvalues_neg_ad <- p.adjust(pvalues_neg, "BH")
  
  pvalues_pos <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_pos[i] <- sum(theta_est[,i]<=0)/5000
  }
  pvalues_pos_ad <- p.adjust(pvalues_pos, "BH")
  
  results <- data.frame(theta_mean, theta_CI_low, theta_CI_high, omega_mean,omega_CI_low,omega_CI_high, pvalues_neg,pvalues_neg_ad, pvalues_pos, pvalues_pos_ad)
  rownames(results) <- colnames(count)
  return(results)
  
}  


