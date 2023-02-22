########################################################################## 
####  Plug-in retrospective lasso logit estimation of                 ####
####  average adjusted association under random sampling              #### 
####  Inputs: (Y_i,T_i,X_i): i=1,...,n                                ####
####  Outputs: estimate                                               ####  
##########################################################################
Lasso_retro = function(outcome=outcome,treat=treat,X=X){
  
  X = scale(X)
  X_matrix = cbind(rep(1,nrow(X)),X)

  ### Estimation of Pr(treat=1|X=x,outcome=1) using glmnet package ###   
  T_case = treat[outcome==1]
  X_case = X[outcome==1,]
  cvfit_case = cv.glmnet(X_case, T_case, family = "binomial", type.measure = "deviance")
  Pr_T_case = predict(cvfit_case, newx = X, s = "lambda.min", type = "response")       
  
  ### Estimation of Pr(treat=1|X=x,outcome=0) using glmnet package ### 
  T_control = treat[outcome==0]
  X_control = X[outcome==0,]
  cvfit_control = cv.glmnet(X_control, T_control, family = "binomial", type.measure = "deviance")
  Pr_T_control = predict(cvfit_control, newx = X, s = "lambda.min", type = "response") 
 
  ### Conditional odds ratio ###   
  OR = (Pr_T_case*(1-Pr_T_control))/((1-Pr_T_case)*Pr_T_control)
  term1 = log(OR)
  est = mean(term1)

return(est) 
}
