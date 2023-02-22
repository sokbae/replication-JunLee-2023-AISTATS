########################################################################## 
####  High-dimensional prospective DML estimation of                  ####
####  average adjusted association under random sampling              #### 
####  Inputs: (Y_i,T_i,X_i): i=1,...,n and K (no. of folds)           ####
####  Outputs: estimate and standard error                            ####  
##########################################################################
AAA_DML_pro = function(outcome=outcome,treat=treat,X=X,K=K){

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
   
  ### Estimation of Pr(outcome=1|X=x,treat=1) using glmnet package ###   
    Y_case = Y_est[T_est==1]
    X_case = X_est[T_est==1,]
    
cvfit_case = cv.glmnet(X_case, Y_case, family = "binomial", type.measure = "deviance")
 Pr_Y_case = predict(cvfit_case, newx = X_main, s = "lambda.min", type = "response")       
     
 ### Estimation of Pr(outcome=1|X=x,treat=0) using glmnet package ### 
    Y_control = Y_est[T_est==0]
    X_control = X_est[T_est==0,]
cvfit_control = cv.glmnet(X_control, Y_control, family = "binomial", type.measure = "deviance")
 Pr_Y_control = predict(cvfit_control, newx = X_main, s = "lambda.min", type = "response") 

 ### Estimation of Pr(treat=1|X=x) using glmnet package ###   
 
 cvfit_T = cv.glmnet(X_est, T_est, family = "binomial", type.measure = "deviance")
    Pr_T = predict(cvfit_T, newx = X_main, s = "lambda.min", type = "response")     
  
 ### Conditional odds ratio ###   

      OR = (Pr_Y_case*(1-Pr_Y_control))/((1-Pr_Y_case)*Pr_Y_control)
   term1 = log(OR)
   term2 = (T_main/Pr_T)*(Y_main - Pr_Y_case)/(Pr_Y_case*(1-Pr_Y_case))
   term3 = -((1-T_main)/(1-Pr_T))*(Y_main - Pr_Y_control)/(Pr_Y_control*(1-Pr_Y_control))
     est = term1  + term2 + term3

  est_matrix[1:nrow(est),k_j] = est

} 

      est = mean(colSums(est_matrix)/n_kfold) 
     avar = mean(colSums((est_matrix-est)^2)/n_kfold) 
       se = sqrt(avar/length(outcome))

output_list <- list("est"=est,"se"=se,"est_matrix"=est_matrix)
return(output_list) 
}
