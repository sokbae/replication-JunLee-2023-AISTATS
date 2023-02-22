### updated on 22 Feb 2023
### DML estimators of Average Adjusted Association
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
library("FNN")

source("src/AAA_DML_pro.R")
source("src/AAA_DML_retro.R")

seed = 345469
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

### DML estimation ###

alpha = 0.05
cv = qnorm(1-alpha/2)

    pro_est_results_HD = AAA_DML_pro(outcome=Y_data,treat=T_data,X=HD_X,K=10)
            pro_est_HD = pro_est_results_HD$est
             pro_se_HD = pro_est_results_HD$se       
             pro_lb_HD = pro_est_HD - cv*pro_se_HD
             pro_ub_HD = pro_est_HD + cv*pro_se_HD
             
  retro_est_results_HD = AAA_DML_retro(outcome=Y_data,treat=T_data,X=HD_X,K=10)
          retro_est_HD = retro_est_results_HD$est
           retro_se_HD = retro_est_results_HD$se 
           retro_lb_HD = retro_est_HD - cv*retro_se_HD
           retro_ub_HD = retro_est_HD + cv*retro_se_HD 

### Save estimation results ###
    
    results_table = matrix(NA,nrow=2,ncol=2)
results_table[1,] = c(pro_est_HD,retro_est_HD)
results_table[2,] = c(pro_se_HD,retro_se_HD)

colnames(results_table) = c("Prospective DML", "Retrospective DML") 
rownames(results_table) = c("Estimate", "Standard Error")

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
results_xtable = xtable(results_table, digits=3)

### Save inference results ###

results_table_CI = matrix(NA,nrow=3,ncol=2)
results_table_CI[1,] = c(exp(pro_est_HD),exp(retro_est_HD))
results_table_CI[2,] = c(exp(pro_lb_HD),exp(retro_lb_HD))
results_table_CI[3,] = c(exp(pro_ub_HD),exp(retro_ub_HD))

colnames(results_table_CI) = c("Prospective DML", "Retrospective DML") 
rownames(results_table_CI) = c("Estimate", "CI (lower)", "CI (upper)")

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
results_xtable = xtable(results_table, digits=2)
results_xtable_CI = xtable(results_table_CI, digits=2)

        
sink(file = "results/ACS2018_AAA_w_DML.txt", append = FALSE)
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


