
data_name <- "breast";

if (data_name == 'breast')
{
pattern_loc <- read.csv('pattern_loc_breast.csv')[,2:3]
patterns <- read.csv('pattern_breast.csv')[,2:5]
}  else if (data_name == 'bulb'){
  pattern_loc <- read.csv('pattern_loc_bulb.csv')[,2:3]
  patterns <- read.csv('pattern_bulb.csv')[,2:5]
}
# loc <- pattern_loc;
# n <- dim(loc)[1];

# Load size factor
if (data_name == "bulb") {
  count = read.csv("count_bulb.csv")[,2:16219];
  loc <- read.csv("loc_bulb.csv")[,2:3]
} else if (data_name == "breast") {
  count = read.csv("count_breast.csv")[,2:14790];
  loc <- read.csv("loc_breast.csv")[,2:3]
}

loc <- as.matrix(loc)
count <- as.matrix(count)
black_index <- which(as.character(rowSums(loc)) %in% setdiff(as.character(rowSums(loc)), as.character(rowSums(pattern_loc))));
loc <- loc[-black_index,];
count <- count[-black_index,];

n <- dim(loc)[1];
s <- exp(rnorm(n, mean = 0, sd = 0.2))
# data <- data.frame(lambda = s, x = loc[, 1], y = loc[, 2]);
# ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = lambda), size = 4) + coord_fixed(ratio = 1) + scale_color_distiller(palette = "Spectral") + theme_classic() + labs(color = "")


# Load settings
Rep <- 10;
p <- 100; # **********************************
p_gamma <- 15; # **********************************
output_path <- paste0("simulated_data/", data_name, "_p_", p, "_");
beta0 <- 2;
phi_mean <- 10;
pis <- c(0, 0.3);
names(pis) <- c("none", "low");
sigmas <- c(0.3); # the paper claims to use c(sqrt(0.2), sqrt(0.35), sqrt(0.6))
beta1s <- c( beta0 + log(3));
names(beta1s) <- c( 'med');
for (model in c("nb")) {
  for (pattern_id in 1:4) {
    for (sigma in sigmas) {
      for (beta1 in beta1s) {
        for (pi in pis) {
          for (i in 1:Rep) {
            set.seed(i);
            # Generate simulated data
            epsilon <- matrix(rnorm(n*p, mean = 0, sd = sigma), nrow = n, ncol = p);
            phi <- rexp(p, 1/phi_mean);
            logLambda <- beta0 + epsilon;
            if (model == "poi") {
              set.seed(i);
              count_null <- matrix(rpois(n*p, s*exp(logLambda)), nrow = n, ncol = p);
            } else if (model == "nb") {
              set.seed(i);
              count_null <- matrix(rnbinom(n*p, mu = s*exp(logLambda), size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), nrow = n, ncol = p);
            }
            gamma <- rep(FALSE, p);
            set.seed(10+i);
            gamma[sample(1:p, p_gamma)] <- TRUE;
            logLambda[patterns[, pattern_id], gamma] <- logLambda[patterns[, pattern_id], gamma] + (beta1 - beta0);
            if (model == "poi") {
              set.seed(i);
              count <- matrix(rpois(n*p, s*exp(logLambda)), nrow = n, ncol = p);
            } else if (model == "nb") {
              set.seed(i);
              count <- matrix(rnbinom(n*p, mu = s*exp(logLambda), size = matrix(phi, nrow = n, ncol = p, byrow = TRUE)), nrow = n, ncol = p);
            }
            H <- matrix(0, nrow = n, ncol = p);
            if (pi != 0) {
              set.seed(i);
              pi_temp <- pi;
              # pi_temp <- 1 - 1/(1 + exp(-10*(logLambda - quantile(logLambda, pi))));
              H <- matrix(rbinom(n*p, 1, pi_temp), nrow = n, ncol = p);
              count_null[which(H == 1)] <- 0;
              count[which(H == 1)] <- 0
            }
            save(gamma, phi, count, count_null, loc, logLambda, s, H, file = paste0(output_path, "model=", model, "_pattern=", pattern_id, "_zero=", names(pis)[which(pis == pi)], "_rep=", i, ".Rdata"))
          }

        }
      }
    }
  }
}




# Check data
# source('code/functions.R')
# library(ggplot2)
# load("simulated_data/simulated_poi_pattern=3_noise=low_fc=high_rep=1.Rdata")
# y <- anscombe_transformer(count);
# j <- sample(which(gamma), 1);
# data <- data.frame(lambda = y[, j], x = loc[, 1], y = loc[, 2]);
# data <- data.frame(lambda = logLambda[, j], x = loc[, 1], y = loc[, 2]);
# ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = lambda), size = 4) + coord_fixed(ratio = 1) + scale_color_distiller(palette = "Spectral") + theme_classic() + labs(color = "")
