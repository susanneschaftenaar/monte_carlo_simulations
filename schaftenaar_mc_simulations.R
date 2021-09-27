######################################################################
######################################################################
########################## ASSIGNMENT 4 ##############################
####################### SUSANNE SCHAFTENAAR ##########################
############################# DPCR, UU ###############################
######################################################################

### clean workspace
rm(list=ls(all=TRUE)) # Clear the workspace

### set working directory
setwd("~/Documents/PHD/Course work/aqm/ass_4")

### Load the required libraries
library(tidyverse)
library(lme4)
library(stargazer)
library(AER)
library(ivpack)
library(MASS)

### data set-up
## set seed
set.seed(53572)

##number of obs
n=5000

## true values b1 and b2
b1 <- -2 # true value X1(A, B, C)
b2 <- 1.5 # true value W (covariates)

## Generate X* and W as correlated normal random variables
pred_vars <- mvrnorm(n, c(-2, -1.5), matrix(c(1, 0.5, 0.5, 1), 2, 2)) # correlation of 0.5 between x_star and W
x_star <- pred_vars[, 1]  
W <- pred_vars[, 2]

df1 <- cbind(W, x_star)
cor(df1)

df2 <- data.frame(id = rep(1:100,50), 
                  a = rep(runif(100,-1,1),50), # 100 intercepts by country, 50 years p/ country
                  Z1 = rep(rnorm(100, 0, 1),50) # instrument only on country-level, 50 years p/ country
)

sort_df <- arrange(df2, id)
sort_df # looks correct

df <- cbind(sort_df, df1) # combine sets
cor(df) # assess current correlation

### set-up three different hypothetical situations to test
##1 create a correlation between Z1 (country-level instrument without violating exclusion criterion) and X1_A
X1_A <- x_star + df$Z1 + rnorm(n, 0, 1)
##2 create an instrument with correlation to the intercept/ clustered data (i.e a violation of the exclusion criterion)
Z2 <- df$Z1 + df$a 
## create a correlation between Z2 (country_level instrument violating exclusion criterion through being correlated with the intercept) and X1
X1_B <- x_star + Z2 + rnorm(n, 0, 1)
##3 create an instrument that is not grouped (i.e. varies per row, emulates country-year variation) and does not violate exclustion criterion
Z3 <- rnorm(n) 
## create a correlation between Z3 (country-year instrument does not violate exclusion criterion) and X1_C
X1_C <- x_star + Z3 + + rnorm(n, 0, 1) 

### set-up three true DGPs based on the above three instruments and correlated main independent variables
Y1 <- df$a + b1*X1_A + b2*W + rnorm(n, 0, 1) # true DGP N(0,1) error, Z1 and X1_A
Y2 <- df$a + b1*X1_B + b2*W + rnorm(n, 0, 1) # true DGP N(0,1) error, Z2 and X1_B
Y3 <- df$a + b1*X1_C + b2*W + rnorm(n, 0, 1) # true DGP N(0,1) error, Z3 and X1_C

df <- cbind(df, X1_A, X1_B, X1_C, Z2, Z3, Y1, Y2, Y3) # put all the above together in a dataframe
cor(df) # assess if set-up works: i.e. right correlations? yes.

df <- dplyr::select(df, c(id, a, Z1, X1_A, Z2, X1_B, Z3, X1_C, W, Y1, Y2, Y3)) # select only those variables to be used in analysis
cor(df) # double check set-up


######################################################################
########################### OLS REGRESSIONS ##########################
######################################################################
## below i draw realizations to assess the conditions under which OLS works with the three set-ups
## 1) a bivariate OLS regression not controlling for W, not accounting for potential clustered data/ fixed intercepts
## 2) a multivariate OLS regression controlling for W, not accounting for potential clustered data/ fixed intercepts
## 3) a bivariate OLS regression not controlling for W, but accounting for potential clustered data/ fixed intercepts
## 4) a multivariate OLS regression controlling for W and accounting for potential clustered data/ fixed intercepts

### estimating first DGP Y1 <- df$a + b1*X1_A + b2*W + rnorm(n, 0, 1)
## 1) OLS, without W
m1_A <- lm(Y1 ~ X1_A)
summary(m1_A)
## 2) OLS, with W
m2_A <- lm(Y1 ~ X1_A + W)
summary(m2_A)
## 3) fixed-effects OLS, without W
m1_fixed_A <- lm(Y1 ~ X1_A + factor(id), data=df)
summary(m1_fixed_A)
## 4) fixed-effects OLS, with W
m2_fixed_A <- lm(Y1 ~ X1_A + W + factor(id), data=df)
summary(m2_fixed_A)

stargazer(m1_A, m1_fixed_A, m2_A, m2_fixed_A, omit = "id", out = "OLS_X1A.tex")

### estimating first DGP Y2 <- df$a + b1*X1_B + b2*W + rnorm(n, 0, 1)
## 1) OLS, without W
m1_B <- lm(Y2 ~ X1_B)
summary(m1_B)
## 2) OLS, with W
m2_B <- lm(Y2 ~ X1_B + W)
summary(m2_B)
## 3) fixed-effects OLS, without W
m1_fixed_B <- lm(Y2 ~ X1_B + factor(id), data=df)
summary(m1_fixed_B)
## 4) fixed-effects OLS, with W
m2_fixed_B <- lm(Y2 ~ X1_B + W + factor(id), data=df)
summary(m2_fixed_B)
stargazer(m1_B, m1_fixed_B, m2_B, m2_fixed_B, omit = "id", out = "OLS_X1B.tex")

### estimating first DGP Y3 <- df$a + b1*X1_C + b2*W + rnorm(n, 0, 1) 
## 1) OLS, without W
m1_C <- lm(Y3 ~ X1_C)
summary(m1_C)
## 2) OLS, with W
m2_C <- lm(Y3 ~ X1_C + W)
summary(m2_C)
## 3) fixed-effects OLS, without W
m1_fixed_C <- lm(Y3 ~ X1_C + factor(id), data=df)
summary(m1_fixed_C)
## 4) fixed-effects OLS, with W
m2_fixed_C <- lm(Y3 ~ X1_C + W + factor(id), data=df)
summary(m2_fixed_C)
stargazer(m1_C, m1_fixed_C, m2_C, m2_fixed_C, omit = "id", out = "OLS_X1C.tex")


######################################################################
################################ IV ################################## 
############################# ANALYSES ###############################
######################################################################

### 1. DGP: Y1 <- df$a + b1*X1_A + b2*W + rnorm(n, 0, 1): low correlation Z1 (-0.12) with id due to country-level Z1 
## first stage
first_stage_A <- lm(X1_A ~ Z1, data=df)
summary(first_stage_A)
## 2SLS
iv1 <- ivreg(formula = Y1  ~ X1_A | Z1, data=df)
summary(iv1, vcov = sandwich, diagnostics = TRUE)
## with clustered standard errors by id (can't use fixed effects due to country-level instrument)
iv1_clust <- cluster.robust.se(iv1, df$id)
stargazer( iv1, iv1_clust, out = "IV1.tex")

### 2. DGP: Y2 <- df$a + b1*X1_B + b2*W + rnorm(n, 0, 1). Z2 is country-level AND induced correlation (close to 0.4) with DGP intercept 
## first stage
first_stage_B <- lm(X1_B ~ Z2, data = df)
summary(first_stage_B)
## 2SLS
iv2 <- ivreg(formula = Y2  ~ X1_B | Z2, data=df)
summary(iv2, vcov = sandwich, diagnostics = TRUE)
## with clustered standard errors by id (can't use fixed effects due to country-level instrument)
iv2_clust <- cluster.robust.se(iv2, df$id)
stargazer( iv2, iv2_clust, out = "IV2.tex")

### 3. DGP: Y3 <- df$a + b1*X1_C + b2*W + rnorm(n, 0, 1). Z3 is country-year. No correlation to intercept (-0.02).
## first stage
first_stage_C <- lm(X1_C ~ Z3, data = df)
summary(first_stage_C)
## 2SLS
iv3 <- ivreg(formula = Y3  ~ X1_C | Z3, data=df)
summary(iv3, vcov = sandwich, diagnostics = TRUE)
## with clustered standard errors by id 
iv3_clust <- cluster.robust.se(iv3, df$id)
stargazer( iv3, iv3_clust, out = "IV3.tex")

stargazer(iv1, iv1_clust, iv2, iv2_clust, iv3, iv3_clust, out = "IVs.tex")

######################################################################
############################ BOOTSTRAPPED ############################
################################ IV ################################## 
############################# ANALYSES ###############################
######################################################################

### 1. bootstrapped version of iv-estimation with Z1, X1_A

## set seed
set.seed(1234)

b1A_est <- summary(iv1)$coef[2,1] # store B1
b1A_est
b1A_se <- summary(iv1)$coef[2,2] # store se B1
b1A_se
b1A_clus_se <- iv1_clust[2,2] # store clustered se B1
b1A_clus_se
b1A_sim <- rnorm(1000,b1A_est,b1A_se) # create 1000 B1s based on B1_A_est + unadjusted se
b1A_sim_clus <- rnorm(1000,b1A_est, b1A_clus_se) # create 1000 B1s based on B1_A_est + clustered se

## 1000 bootstraps
boots = 1000

## Setup bootstrap loop
b1A_bs <- rep(NA,boots) # empty vector to add bootstrapped B1s
for(i in 1:boots){ # for loop number of times sim to run 
  bs_index <- sample(1:nrow(df),boots, replace=T) #then steps from slide: index 1 to no rows in df
  df_bs <- df[bs_index,] # original data, then rows of bootstrap index
  m1_bs <- ivreg(Y1 ~ X1_A | Z1 , data=df_bs) # run the model, with df_bs
  b1A_bs[i] <- m1_bs$coefficients[2] # extract coefficient for b1, store in b1_bs that was the empty vector. 1000 times, 1 for each bootstrap
}

## create density plot: comparison B1_A with unadjusted se, with clustered se, and bootstrapped
pdf("b1A_bootstrapped.pdf")
df_density_A <- data.frame(regular = b1A_sim, # reg similated b1s using normal distr
                         regular_clustered = b1A_sim_clus,
                         bootstrap = b1A_bs) %>% gather(key="type",value="value") # call bootstrap. Throw in gather, otherwise doesnt work in ggplot

ggplot(df_density_A,aes(x=value,col=type))+ # corresponds to df_density vars
  geom_density()
dev.off()

## parametric confidence intervals
df_density_A %>% group_by(type) %>% summarize(mean=mean(value),
                                            sd = sd(value))
## nonparametric confidence intervals
quantile(b1A_bs, 0.025)
quantile(b1A_bs, 0.975)



### 2. bootstrapped version of iv-estimation with Z2, X1_B

## set seed
set.seed(1234)


b1B_est <- summary(iv2)$coef[2,1] # store B1
b1B_est
b1B_se <- summary(iv2)$coef[2,2] # store B1 se 
b1B_se
b1B_clus_se <- iv2_clust[2,2] # store clustered se B1
b1B_clus_se
b1B_sim <- rnorm(1000,b1B_est,b1B_se) # create 1000 B1s based on B1_A_est + unadjusted se
b1B_sim_clus <- rnorm(1000,b1B_est, b1B_clus_se) # create 1000 B1s based on B1_A_est + clustered se


## Setup bootstrap loop
b1B_bs <- rep(NA,boots) # empty vector to adds bootstrap
for(i in 1:boots){ # for loop number of times sim to run (1000)
  bs_index <- sample(1:nrow(df),boots, replace=T) #then steps from slide: index 1 to no rows in df
  df_bs <- df[bs_index, ] # original data, then rows of bootstrap index
  m1_bs <- ivreg(Y2 ~ X1_B | Z2 , data=df_bs) # run the model, with df_bs
  b1B_bs[i] <- m1_bs$coefficients[2] # extract coefficient for b1, store in b1B_bs that was the empty vector. 1000 times, 1 for each bootstrap
}

## create density plot: comparison B1_B with unadjusted se, with clustered se, and bootstrapped
pdf("b1B_bootstrapped.pdf")
df_density_B <- data.frame(regular = b1B_sim, # reg similated b3s using normal distr
                         regular_clustered = b1B_sim_clus,
                         bootstrap = b1B_bs) %>% gather(key="type",value="value") # call bootstrap. Throw in gather, otherwise doesnt work in ggplot

ggplot(df_density_B,aes(x=value,col=type))+ # corresponds to df_density vars
  geom_density()
dev.off()

## parametric confidence intervals
df_density_B %>% group_by(type) %>% summarize(mean=mean(value),
                                            sd = sd(value))
## nonparametric confidence intervals
quantile(b1B_bs, 0.025)
quantile(b1B_bs, 0.975)



### 3. bootstrapped version of iv-estimation with Z2, X1_B

b1C_est <- summary(iv3)$coef[2,1] # store B1
b1C_est
b1C_se <- summary(iv3)$coef[2,2] # store se B1
b1C_se
b1C_clus_se <- iv3_clust[2,2] # store clustered se B1
b1C_clus_se
b1C_sim <- rnorm(1000,b1C_est,b1C_se) # create 1000 B1s based on B1_C_est + unadjusted se
b1C_sim_clus <- rnorm(1000,b1C_est, b1C_clus_se) # create 1000 B1s based on B1_C_est + clustered se


## setup bootstrap loop
b1C_bs <- rep(NA,boots) # empty vector to adds bootstrap
for(i in 1:boots){ # for loop number of times sim to run (1000)
  bs_index <- sample(1:nrow(df),boots, replace=T) #then steps from slide: index 1 to no rows in df
  df_bs <- df[bs_index, ] # original data, then rows of bootstrap index
  m1_bs <- ivreg(Y3 ~ X1_C | Z3 , data=df_bs) # run the model, with df_bs
  b1C_bs[i] <- m1_bs$coefficients[2] # extract coefficient for b1, store in b1C_bs that was the empty vector. 1000 times, 1 for each bootstrap
}

## create density plot: comparison B1_B with unadjusted se, with clustered se, and bootstrapped
pdf("b1C_bootstrapped.pdf")
df_density_C <- data.frame(regular = b1C_sim, # reg similated b3s using normal distr
                         regular_clustered = b1C_sim_clus,
                         bootstrap = b1C_bs) %>% gather(key="type",value="value") # call bootstrap. Throw in gather, otherwise doesnt work in ggplot

ggplot(df_density_C,aes(x=value,col=type))+ # corresponds to df_density vars
  geom_density()
dev.off()

# parametric confidence intervals
df_density_C %>% group_by(type) %>% summarize(mean=mean(value),
                                            sd = sd(value))
## nonparametric confidence intervals
quantile(b1C_bs, 0.025)
quantile(b1C_bs, 0.975)


######################################################################
################################ MC ##################################
############################ SIMULATION ##############################
################################ IV ################################## 
############################# ANALYSES ###############################
######################################################################

reps = 1000
par.est.iv <- matrix(NA, nrow=reps, ncol=15) # matrix to store estimates in (B1s + se)
par_est <- matrix(NA, nrow=reps, ncol=6)# matrix to store estimates in (only B1s)

for(i in 1:reps){ # Start the loop
  # set-up data
  pred_vars <- mvrnorm(n, c(-2, -1.5), matrix(c(1, 0.5, 0.5, 1), 2, 2))
  x_star <- pred_vars[, 1]  
  W <- pred_vars[, 2]
  
  df1 <- cbind(W, x_star)

  df2 <- data.frame(id = rep(1:100,50), 
                    a = rep(runif(100,-1,1),50), # 100 intercepts by country, 50 years p/ country
                    Z1 = rep(rnorm(100, 0, 1),50) # instrument only on country-level, 50 years p/ country
  )
  
  sort_df <- arrange(df2, id)
  df <- cbind(sort_df, df1) # combine sets
  # create an instrument with correlation to the clustered data (i.e a violation of the exclusion criterion)
  Z2 <- df$Z1 + df$a 
  # create a correlation between Z1 (instrument without violating exclusion criterion) and X1
  X1_A <- x_star + df$Z1 + rnorm(n, 0, 1)
  # create a correlation between Z2 (instrument violating exclusion criterion) and X1
  X1_B <- x_star + Z2 + rnorm(n, 0, 1)
  # create an instrument that is not grouped
  Z3 <- rnorm(n) 
  X1_C <- x_star + Z3 + + rnorm(n, 0, 1) 
  
  Y1 <- df$a + b1*X1_A + b2*W + rnorm(n, 0, 1) # true DGP N(0,1) error
  Y2 <- df$a + b1*X1_B + b2*W + rnorm(n, 0, 1) # true DGP N(0,1) error
  Y3 <- df$a + b1*X1_C + b2*W + rnorm(n, 0, 1) # true DGP N(0,1) error
  
  df <- cbind(df, X1_A, X1_B, X1_C, Z2, Z3, Y1, Y2, Y3)
  df <- dplyr::select(df, c(id, a, Z1, X1_A, Z2, X1_B, Z3, X1_C, W, Y1, Y2, Y3))

  
  ## run iv model 1: Z1 grouped by id, X1_A multiple observations per id
  iv1 <- ivreg(formula = Y1  ~ X1_A | Z1, data=df)
  iv1_clust <- cluster.robust.se(iv1, df$id)
  
  par.est.iv[i, 1] <- summary(iv1)$coef[2,1] # store B1 for ivmodel 1
  par.est.iv[i, 2] <- summary(iv1)$coef[2,2] # store se ivmodel 1
  par.est.iv[i, 3] <- iv1_clust[2,2] # store se ivmodel 1 clustered
  
  
  ## run iv model 2: Z2 grouped by id and covaries with intercept, X1_B multiple observations per id
  iv2 <- ivreg(formula = Y2  ~ X1_B | Z2, data=df)
  iv2_clust <- cluster.robust.se(iv2, df$id)
  
  par.est.iv[i, 4] <- summary(iv2)$coef[2,1] # store B1 for ivmodel 2
  par.est.iv[i, 5] <- summary(iv2)$coef[2,2] # store se ivmodel 2
  par.est.iv[i, 6] <- iv2_clust[2,2] # store se ivmodel 2 clustered
  
  
  ## run iv model 3: Z3 not grouped (multiple obs per id), X1_C multiple observations per id
  iv3 <- ivreg(formula = Y3  ~ X1_C | Z3, data=df)
  iv3_clust <- cluster.robust.se(iv3, df$id)
  
  par.est.iv[i, 7] <- summary(iv3)$coef[2,1] # store B1 for ivmodel 3
  par.est.iv[i, 8] <- summary(iv3)$coef[2,2] # store se ivmodel 3
  par.est.iv[i, 9] <- iv3_clust[2,2] # store se ivmodel 3 clustered
  
  par_est[i, 1] <- summary(iv1)$coef[2,1] # store B1 for ivmodel 1
  par_est[i, 2] <- summary(iv2)$coef[2,1] # store B1 for ivmodel 2
  par_est[i, 3] <- summary(iv3)$coef[2,1] # store B1 for ivmodel 3
  
  
  ## bootstrapping iv1, 2, 3
  boots_sim = 1000 # no of bootstraps
  
  ## Setup bootstrap loop
  b1A_bs_sim <- rep(NA,boots_sim) # empty vector to adds bootstrap
  for(j in 1:boots_sim){ # for loop number of times sim to run (1000)
    bs_index <- sample(1:nrow(df),boots_sim, replace=T) #then steps from slide: index 1 to no rows in df
    df_bs <- df[bs_index,] # original data, then rows of bootstrap index
    m1_bs_sim <- ivreg(Y1 ~ X1_A | Z1 , data=df_bs) # run the model, with df_bs
    b1A_bs_sim[j] <- m1_bs_sim$coefficients[2] # extract coefficient for b1, store in b1_bs that was the empty vector. 1000 times, 1 for each bootstrap
  }
  
  par.est.iv[i, 10] <- mean(b1A_bs_sim) # store B1 bootstrapped ivmodel 1 
  par.est.iv[i, 11] <- sd(b1A_bs_sim) # store B1 se bootstrapped ivmodel 1 
  
  par_est[i, 4] <- mean(b1A_bs_sim) # store B1 for bootstrapped ivmodel 1
  
  b1B_bs_sim <- rep(NA,boots_sim) # empty vector to adds bootstrap
  for(j in 1:boots_sim){ # for loop number of times sim to run (1000)
    bs_index <- sample(1:nrow(df),boots_sim, replace=T) #then steps from slide: index 1 to no rows in df
    df_bs <- df[bs_index,] # original data, then rows of bootstrap index
    m2_bs_sim <- ivreg(Y2 ~ X1_B | Z2 , data=df_bs) # run the model, with df_bs
    b1B_bs_sim[j] <- m2_bs_sim$coefficients[2] # extract coefficient for b1, store in b1_bs that was the empty vector. 1000 times, 1 for each bootstrap
  }
  
  par.est.iv[i, 12] <- mean(b1B_bs_sim) # store coef bootstrapped iv2
  par.est.iv[i, 13] <- sd(b1B_bs_sim) # store B1 se bootstrapped ivmodel 2 
  
  par_est[i, 5] <- mean(b1B_bs_sim) # store coef bootstrapped iv2
  

  b1C_bs_sim <- rep(NA,boots_sim) # empty vector to adds bootstrap
  for(j in 1:boots_sim){ # for loop number of times sim to run (1000)
    bs_index <- sample(1:nrow(df),boots_sim, replace=T) #then steps from slide: index 1 to no rows in df
    df_bs <- df[bs_index,] # original data, then rows of bootstrap index
    m3_bs_sim <- ivreg(Y3 ~ X1_C | Z3 , data=df_bs) # run the model, with df_bs
    b1C_bs_sim[j] <- m3_bs_sim$coefficients[2] # extract coefficient for b1, store in b1_bs that was the empty vector. 1000 times, 1 for each bootstrap
  }
  
  par.est.iv[i, 14] <- mean(b1C_bs_sim) # store B1 boostrapped ivmodel 3 
  par.est.iv[i, 15] <- sd(b1C_bs_sim) # store B1 se bootstrapped ivmodel 3 
  
  par_est[i, 6] <- mean(b1C_bs_sim) # store B1 for bootstrapped ivmodel 1
  
  {
    cat(i,"\n")
  }
} # End the loop

colnames(par.est.iv) <- c("IV1", "IV1 se", 
                          "IV1 se cl", 
                          "IV2", "IV2 se",
                          "IV2 se cl",
                          "IV3", "IV3 se",
                          "IV3 se cl",
                          "BS IV1", "BS IV1 se",
                          "BS IV2", "BS IV2 se",
                          "BS IV3", "BS IV3 se")

### overview stored B1s, se, clustered se
mean(par.est.iv[, 1]) # B1
mean(par.est.iv[, 2]) # se
mean(par.est.iv[, 3]) # se cl
mean(par.est.iv[, 4]) # B1
mean(par.est.iv[, 5]) # se
mean(par.est.iv[, 6]) # se cl
mean(par.est.iv[, 7]) # B1 
mean(par.est.iv[, 8]) # se
mean(par.est.iv[, 9]) # se cl
mean(par.est.iv[, 10]) # B1 bs IV1
mean(par.est.iv[, 11]) # B1 bs IV1 se
mean(par.est.iv[, 12]) # B1 bs IV2
mean(par.est.iv[, 13]) # B1 bs IV2 se
mean(par.est.iv[, 14]) # B1 bs IV3
mean(par.est.iv[, 15]) # B1 bs IV3 se

### are the estimands biased: assess MSE
mse_b1_A <- mean((par.est.iv[, 1]-b1)^2)
mse_b1_B <- mean((par.est.iv[, 4]-b1)^2)
mse_b1_C <- mean((par.est.iv[, 7]-b1)^2)

mse_b1_A_bstr <- mean((par.est.iv[, 10]-b1)^2)
mse_b1_B_bstr <- mean((par.est.iv[, 12]-b1)^2)
mse_b1_C_bstr <- mean((par.est.iv[, 14]-b1)^2)


### are the estimands biased: visualise by plotting B1s for different models

pdf("IV1_MC.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.iv[ , 1]), lty = 1, xlim = c(-2.5, -1.5),
     ylim = c(0, 10), lwd = 3, xlab = "", ylab = "", main = "", axes = TRUE)
lines(density(par.est.iv[ , 10]), lwd = 3, lty = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.7, 7, expression("True"~beta[1]~"= -2"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("2SLS IV1"), expression("2 SLS bootstrapped IV1")),
       lty = c(1, 2), lwd = 3, cex = 1)
dev.off()

pdf("IV2_MC.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.iv[ , 4]), lty = 1, xlim = c(-2.5, -1.5),
     ylim = c(0, 10), lwd = 3, xlab = "", ylab = "", main = "", axes = TRUE)
lines(density(par.est.iv[ , 12]), lwd = 3, lty = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.7, 7, expression("True"~beta[1]~"= -2"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("2SLS IV2"), expression("2SLS bootstrapped IV2")),
       lty = c(1, 2), lwd = 3, cex = 1)
dev.off()

pdf("IV3_MC.pdf")
par(mar = c(5, 5.25, .5, .5))
plot(density(par.est.iv[ , 7]), lty = 1, xlim = c(-2.5, -1.5),
     ylim = c(0, 15), lwd = 3, xlab = "", ylab = "", main = "", axes = TRUE)
lines(density(par.est.iv[ , 14]), lwd = 3, lty = 2)
title(xlab = expression(hat(beta[1])), cex.lab = 1.5)
title(ylab = expression("Density"), line = 3.75, cex.lab = 1.5)
abline(v = b1, lwd = 2)
text(.7, 7, expression("True"~beta[1]~"= -2"), cex = 1.5)
box()
legend("topright", bty = "n", c(expression("2SLS"), expression("2SLS bootstrapped IV3")),
       lty = c(1, 2), lwd = 3, cex = 1)
dev.off()

### what about the standard errors: coverage

# coverage function
# CP Function
coverage <- function(b, se, true, level = .95, df = Inf){ # Estimate, 
  # standard error,
  # true parameter, 
  # confidence level, 
  # and df  
  qtile <- level + (1 - level)/2 # Compute the proper quantile
  lower.bound <- b - qt(qtile, df = df)*se # Lower bound
  upper.bound <- b + qt(qtile, df = df)*se # Upper bound 
  # Is the true parameter in the confidence interval? (yes = 1)
  true.in.ci <- ifelse(true >= lower.bound & true <= upper.bound, 1, 0)
  cp <- mean(true.in.ci) # The coverage probability
  mc.lower.bound <- cp - 1.96*sqrt((cp*(1 - cp))/length(b)) # Monte Carlo error  
  mc.upper.bound <- cp + 1.96*sqrt((cp*(1 - cp))/length(b))  
  return(list(coverage.probability = cp, # Return results
              true.in.ci = true.in.ci,
              ci = cbind(lower.bound, upper.bound),
              mc.eb = c(mc.lower.bound, mc.upper.bound)))
}


iv_1 <- coverage(par.est.iv[ , 1], par.est.iv[ , 2], b1,
                   df = n - iv1$rank)
iv_1_clust <- coverage(par.est.iv[ , 1], par.est.iv[ , 3], b1,
                  df = n - iv1$rank)
iv_1_bs <- coverage(par.est.iv[ , 10], par.est.iv[ , 11], b1,
                   df = n - m1_bs_sim$rank)


iv_2 <- coverage(par.est.iv[ , 4], par.est.iv[ , 5], b1,
                 df = n - iv2$rank)
iv_2_clust <- coverage(par.est.iv[ , 4], par.est.iv[ , 6], b1,
                       df = n - iv2$rank)
iv_2_bs <- coverage(par.est.iv[ , 12], par.est.iv[ , 13], b1,
                    df = n - m2_bs_sim$rank)

iv_3 <- coverage(par.est.iv[ , 7], par.est.iv[ , 8], b1,
                 df = n - iv3$rank)
iv_3_clust <- coverage(par.est.iv[ , 7], par.est.iv[ , 9], b1,
                       df = n - iv3$rank)
iv_3_bs <- coverage(par.est.iv[ , 14], par.est.iv[ , 15], b1,
                    df = n - m3_bs_sim$rank)


# coverage plot IV1
pdf("IV1_cp.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(1, iv_1$coverage.probability, pch = 19, xlim = c(0, 8),
     ylim = c(0, 1), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
segments(1, iv_1$mc.eb[1], 1, iv_1$mc.eb[2], lwd = 2)
points(3, iv_1_clust$coverage.probability, pch = 1, lwd = 3)
segments(3, iv_1_clust$mc.eb[1], 3, iv_1_clust$mc.eb[2], lwd = 2) 
points(5, iv_1_bs$coverage.probability, pch = 19, lwd = 3) 
segments(5, iv_1_bs$mc.eb[1], 5, iv_1_bs$mc.eb[2], lwd = 2) 
axis(1, at = c(1, 3, 5), labels = c(expression("IV1"), expression("IV1 clust se"),
                                       expression("IV1 bs")), cex.axis = 1.25)
axis(2, at = seq(0, 1, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Estimator"), cex.lab = 1.5)
title(ylab = expression("Coverage Probability"), line = 3.75, cex.lab = 1.5)
abline(h = .95, lwd = 2, lty = 2)
box()
dev.off()


# coverage plot IV2
pdf("IV2_cp.pdf")

par(mar = c(5, 5.25, .5, .5))
plot(1, iv_2$coverage.probability, pch = 19, xlim = c(0, 8),
     ylim = c(0, 1), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
segments(1, iv_2$mc.eb[1], 1, iv_2$mc.eb[2], lwd = 2)
points(3, iv_2_clust$coverage.probability, pch = 1, lwd = 3)
segments(3, iv_2_clust$mc.eb[1], 3, iv_2_clust$mc.eb[2], lwd = 2) 
points(5, iv_2_bs$coverage.probability, pch = 19, lwd = 3) 
segments(5, iv_2_bs$mc.eb[1], 5, iv_2_bs$mc.eb[2], lwd = 2) 
axis(1, at = c(1, 3, 5), labels = c(expression("IV2"), expression("IV2 clust se"),
                                    expression("IV2 bs")), cex.axis = 1.25)
axis(2, at = seq(0, 1, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Estimator"), cex.lab = 1.5)
title(ylab = expression("Coverage Probability"), line = 3.75, cex.lab = 1.5)
abline(h = .95, lwd = 2, lty = 2)
box()
dev.off()

# coverage plot IV3
pdf("IV3_cp.pdf")
par(mar = c(5, 5.25, .5, .5))
plot(1, iv_3$coverage.probability, pch = 19, xlim = c(0, 8),
     ylim = c(0, 1), lwd = 3, xlab = "", ylab = "", main = "", axes = FALSE)
segments(1, iv_3$mc.eb[1], 1, iv_3$mc.eb[2], lwd = 2)
points(3, iv_3_clust$coverage.probability, pch = 1, lwd = 3)
segments(3, iv_3_clust$mc.eb[1], 3, iv_3_clust$mc.eb[2], lwd = 2) 
points(5, iv_3_bs$coverage.probability, pch = 19, lwd = 3) 
segments(5, iv_3_bs$mc.eb[1], 5, iv_3_bs$mc.eb[2], lwd = 2) 
axis(1, at = c(1, 3, 5), labels = c(expression("IV3"), expression("IV3 clust se"),
                                    expression("IV3 bs")), cex.axis = 1.25)
axis(2, at = seq(0, 1, .05), cex.axis = 1.25, las = 2)
title(xlab = expression("Estimator"), cex.lab = 1.5)
title(ylab = expression("Coverage Probability"), line = 3.75, cex.lab = 1.5)
abline(h = .95, lwd = 2, lty = 2)
box()
dev.off()
