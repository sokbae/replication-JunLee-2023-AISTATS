### updated on 22 Feb 2023
### Monte Carlo Experiments for Average Adjusted Association
rm(list=ls(all=T))
gc()

start_time = Sys.time()

# need to set the working directory 
#setwd("TO BE ADDED")

library("fastDummies")
library("splines")
library("MASS")
library("Matrix")
library("glmnet")


source("src/AAA_DML_pro.R")
source("src/AAA_DML_retro.R")
source("src/Lasso_pro.R")
source("src/Lasso_retro.R")
source("src/Logit_pro.R")
source("src/Logit_retro.R")

seed = 432624
set.seed(seed) 

ACS_case = read.table("data/ACS_cases_HD.csv",header=T,sep=",")
colnames(ACS_case) = c("age", "married", "ind", "baplus", "topincome")

ACS_control = read.table("data/ACS_controls_HD.csv",header=T,sep=",")
colnames(ACS_control) = c("age", "married", "ind", "baplus", "topincome")

dt = rbind(ACS_case,ACS_control)
 n = nrow(dt)
dt = dt[sample.int(n),] # shuffle the data 
 
attach(dt)

   age_bs_HD = bs(age, df=20)              # high-dimensional b-splines for age
       
     ind_var = dummy_cols(ind)  
 ind_dummies = ind_var[,3:ncol(ind_var)]  
 ind_dummies_matrix = as.matrix(ind_dummies)

      T_data = baplus
      Y_data = topincome
           X = age
        HD_X = cbind(age_bs_HD,ind_dummies_matrix) 

prop_topincome = mean(Y_data)

### Simulation Input ###

age_quad = cbind(age, age^2)
fit_treat = glm(formula = baplus ~ age_quad, 
               family = "binomial")
fit_treat_coef = fit_treat$coefficients
Pr_treat = predict(fit_treat, newx = age_quad, type = "response")

fit_outcome = glm(formula = topincome ~ baplus + age_quad, 
                family = "binomial")
fit_outcome_coef = fit_outcome$coefficients
Pr_outcome = predict(fit_outcome, newx = age_quad, type = "response")

true_coefficient = fit_outcome$coefficients["baplus"]

alpha = 0.10
cv = qnorm(1-alpha/2)

n_sim = 500
m = 5000    # sample size in the Monte Carlo experiments

results = {}
for (sim_i in 1:n_sim){
  
# print out sim_i to check the progress
  print(sim_i)

si = sample.int(n,m)
HD_X_s = HD_X[si,]
X_s = cbind(rep(1,m),age_quad[si,])
fit_index = X_s%*%fit_treat_coef
Pr_treat_s = exp(fit_index)/(1+exp(fit_index)) 
T_s = as.integer(runif(m) < Pr_treat_s)
TX_s = cbind(rep(1,m),T_s,age_quad[si,])
fit_index = TX_s%*%fit_outcome_coef
Pr_outcome_s = exp(fit_index)/(1+exp(fit_index)) 
Y_s = as.integer(runif(m) < Pr_outcome_s)

### Lasso logit estimation ###

pro_est_lasso = Lasso_pro(outcome=Y_s,treat=T_s,X=HD_X_s)
retro_est_lasso = Lasso_retro(outcome=Y_s,treat=T_s,X=HD_X_s)

### DML estimation ###

    pro_est_results_HD = AAA_DML_pro(outcome=Y_s,treat=T_s,X=HD_X_s,K=5)
            pro_est_HD = pro_est_results_HD$est
             pro_se_HD = pro_est_results_HD$se       

             
  retro_est_results_HD = AAA_DML_retro(outcome=Y_s,treat=T_s,X=HD_X_s,K=5)
          retro_est_HD = retro_est_results_HD$est
           retro_se_HD = retro_est_results_HD$se 


### Save estimation results ###
           
estimates = c(pro_est_lasso, retro_est_lasso,
              pro_est_HD,retro_est_HD,
              pro_se_HD,retro_se_HD)

results = rbind(results, estimates)

}

bias = slam::col_means(results[,1:4]-true_coefficient)
std = sd(results[,1:4])

pro_lb_HD = results[,3] - cv*results[,5]
pro_ub_HD = results[,3] + cv*results[,5]
retro_lb_HD = results[,4] - cv*results[,6]
retro_ub_HD = results[,4] + cv*results[,6]
  pro_size = (pro_lb_HD <= true_coefficient)*(pro_ub_HD >= true_coefficient)
retro_size = (retro_lb_HD <= true_coefficient)*(retro_ub_HD >= true_coefficient)

    results_table = matrix(NA,nrow=4,ncol=2)
results_table[,1] = bias
results_table[,2] = std

colnames(results_table) = c("Mean Bias", "Std. Dev.") 
rownames(results_table) = c("Pro. Lasso", "Retro. Lasso",
                            "Pro. DML", "Retro. DML")

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
results_xtable = xtable(results_table, digits=3)

### Save inference results ###

results_table_CI = matrix(NA,nrow=1,ncol=2)
results_table_CI[1,] = c(mean(pro_size),mean(retro_size))

colnames(results_table_CI) = c("Prospective DML", "Retrospective DML") 
rownames(results_table_CI) = c("Size")

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
results_xtable = xtable(results_table, digits=2)
results_xtable_CI = xtable(results_table_CI, digits=2)


output_file_name <- paste("results/MC",".Rdata", sep = "")
save(true_coefficient, results, file = output_file_name)

        
sink(file = "results/MC_AAA_2022.txt", append = FALSE)
options(digits=3)
    
print("Original Sample Size")
print(n)
    
print("Proportion of Top-Coded")
print(prop_topincome)
    
print(results_xtable)

print(results_xtable_CI)

time_taken = Sys.time() - start_time 
print("Time Taken")
print(time_taken)
sink()

detach(dt)


