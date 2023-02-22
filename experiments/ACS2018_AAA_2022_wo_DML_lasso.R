### updated on 9 October 2022
### Plug-in estimators of Average Adjusted Association without Regularization
### Lasso based estimation
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

source("src/Lasso_pro.R")
source("src/Lasso_retro.R")

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

### Lasso logit estimation ###

   pro_est_logit = Lasso_pro(outcome=Y_data,treat=T_data,X=HD_X)
 retro_est_logit = Lasso_retro(outcome=Y_data,treat=T_data,X=HD_X)

 ### Save estimation results ### 
 
 sink(file = "results/ACS2018_AAA_2022_wo_DML_lasso.txt", append = FALSE)
 options(digits=3)
 
 print("Sample Size")
 print(n)
 
 print("Proportion of Top-Coded")
 print(prop_topincome) 
 
 print("Prospective Estimate: theta --- exp(theta)")
 print(c(pro_est_logit, exp(pro_est_logit)))
 
 print("Retrospective Estimate: theta --- exp(theta)")
 print(c(retro_est_logit, exp(retro_est_logit)))

time_taken = Sys.time() - start_time 
print("Time Taken")
print(time_taken)
sink()

detach(dt)


