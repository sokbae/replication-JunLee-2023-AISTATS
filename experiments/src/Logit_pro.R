########################################################################## 
####  Plug-in prospective logit estimation of                         ####
####  average adjusted association under random sampling              #### 
####  Inputs: (Y_i,T_i,X_i): i=1,...,n                                ####
####  Outputs: estimate                                               ####  
##########################################################################
Logit_DML_pro = function(outcome=outcome,treat=treat,X=X){
  
  X = scale(X)
  X_matrix = cbind(rep(1,nrow(X)),X)
  epsilon = 1e-06  

  ### Estimation of Pr(outcome=1|X=x,treat=1) using logit ###   
  fit_case = glm(formula = outcome[treat==1] ~ X[treat==1,], 
                 family = "binomial")
  fit_coef = fit_case$coefficients
  fit_coef[is.na(fit_coef)] = 0    # replace NA with 0      
  case_index = X_matrix%*%fit_coef
  Pr_Y_case = exp(case_index)/(1+exp(case_index)) 
  Pr_Y_case[Pr_Y_case <= epsilon] = epsilon # lower bound by epsilon 
  
  ### Estimation of Pr(outcome=1|X=x,treat=0) using logit ###     
  fit_control = glm(formula = outcome[treat==0] ~ X[treat==0,], 
                    family = "binomial")
  control_coef = fit_control$coefficients
  control_coef[is.na(control_coef)] = 0    # replace NA with 0 
  control_index = X_matrix%*%control_coef
  Pr_Y_control = exp(control_index)/(1+exp(control_index)) 
  Pr_Y_control[Pr_Y_control <= epsilon] = epsilon  # lower bound by epsilon  
  
  ### Conditional odds ratio ###   
  OR = (Pr_Y_case*(1-Pr_Y_control))/((1-Pr_Y_case)*Pr_Y_control)
  term1 = log(OR)
  est = mean(term1)

return(est) 
}
