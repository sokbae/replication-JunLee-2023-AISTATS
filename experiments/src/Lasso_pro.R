########################################################################## 
####  Plug-in prospective lasso logit estimation of                   ####
####  average adjusted association under random sampling              #### 
####  Inputs: (Y_i,T_i,X_i): i=1,...,n                                ####
####  Outputs: estimate                                               ####  
##########################################################################
Lasso_pro = function(outcome=outcome,treat=treat,X=X){
  
  ### Estimation of Pr(outcome=1|X=x,treat=1) using glmnet package ###   
  Y_case = outcome[treat==1]
  X_case = X[treat==1,]
  cvfit_case = cv.glmnet(X_case, Y_case, family = "binomial", type.measure = "deviance")
  Pr_Y_case = predict(cvfit_case, newx = X, s = "lambda.min", type = "response")       
  
  ### Estimation of Pr(outcome=1|X=x,treat=0) using glmnet package ### 
  Y_control = outcome[treat==0]
  X_control = X[treat==0,]
  cvfit_control = cv.glmnet(X_control, Y_control, family = "binomial", type.measure = "deviance")
  Pr_Y_control = predict(cvfit_control, newx = X, s = "lambda.min", type = "response") 
  
  ### Conditional odds ratio ###   
  OR = (Pr_Y_case*(1-Pr_Y_control))/((1-Pr_Y_case)*Pr_Y_control)
  term1 = log(OR)
  est = mean(term1)

return(est) 
}
