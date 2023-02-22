########################################################################## 
####  High-dimensional retrospective DML estimation of                ####
####  average adjusted association under random sampling              #### 
####  Inputs: (Y_i,T_i,X_i): i=1,...,n and K (no. of folds)           ####
####  Outputs: estimate and standard error                            ####  
##########################################################################
AAA_DML_retro = function(outcome=outcome,treat=treat,X=X,K=K){

    h_hat = mean(outcome)

  ### Construct k-fold partition using random splitting ###

  unif_tmp = runif(length(outcome), min = 0, max = K)      
kfoldindex = ceiling(unif_tmp)

  n_kfold = {}
 
 for(k_j in 1:K){
   
   n_main = length(outcome[kfoldindex==k_j])
  n_kfold = c(n_kfold, n_main)
 }  
   n_kmax = max(n_kfold)
 
 est_matrix = matrix(0,nrow=n_kmax,ncol=K)

 for(k_j in 1:K){
   
  Y_main = outcome[kfoldindex==k_j]
   Y_est = outcome[kfoldindex!=k_j]
  T_main = treat[kfoldindex==k_j]
   T_est = treat[kfoldindex!=k_j] 
  X_main = X[kfoldindex==k_j,]
   X_est = X[kfoldindex!=k_j,]    
   
  ### Estimation of Pr(treat=1|X=x,outcome=1) using glmnet package ###   
    T_case = T_est[Y_est==1]
    X_case = X_est[Y_est==1,]
    
cvfit_case = cv.glmnet(X_case, T_case, family = "binomial", type.measure = "deviance")
 Pr_T_case = predict(cvfit_case, newx = X_main, s = "lambda.min", type = "response")       
     
 ### Estimation of Pr(treat=1|X=x,outcome=0) using glmnet package ### 
    T_control = T_est[Y_est==0]
    X_control = X_est[Y_est==0,]
cvfit_control = cv.glmnet(X_control, T_control, family = "binomial", type.measure = "deviance")
 Pr_T_control = predict(cvfit_control, newx = X_main, s = "lambda.min", type = "response") 

 ### Estimation of Pr(outcome=1|X=x) using glmnet package ###   
 
 cvfit_Y = cv.glmnet(X_est, Y_est, family = "binomial", type.measure = "deviance")
    Pr_Y = predict(cvfit_Y, newx = X_main, s = "lambda.min", type = "response")     
  
 ### Conditional odds ratio ###   

      OR = (Pr_T_case*(1-Pr_T_control))/((1-Pr_T_case)*Pr_T_control)
   term1 = log(OR)
   term2 = (Y_main/Pr_Y)*(T_main - Pr_T_case)/(Pr_T_case*(1-Pr_T_case))
   term3 = -((1-Y_main)/(1-Pr_Y))*(T_main - Pr_T_control)/(Pr_T_control*(1-Pr_T_control))
     est = term1  + term2 + term3

  est_matrix[1:nrow(est),k_j] = est

} 

      est = mean(colSums(est_matrix)/n_kfold) 
     avar = mean(colSums((est_matrix-est)^2)/n_kfold) 
       se = sqrt(avar/length(outcome))

output_list <- list("est"=est,"se"=se,"est_matrix"=est_matrix)
return(output_list) 
}
