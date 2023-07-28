##############################################################################
#                                                                            #
#                                                                            #
#  Self-defined functions used for "Enhanced Inference for Finite Population #
#  Sampling-Based Prevalence Estimation with Misclassification Errors        #
#                                                                            #
#                                                                            #
##############################################################################

## A function to simulate individual level imperfect testing data based on 
## pre-specified parameters

library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)

simu_sym_RS = function(N,N_true,p.RS,Se,Sp){
  
  ## Parameters
  ## N: total population size
  ## N_true: true disease people (determined by the disease prevalence)
  ## p.RS: the sampling rate for random samples
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # calculate the sample size
  n.RS = round(N*p.RS)
  
  # generate the individual testing data by sampling from the total population
  test = rep(0,N)
  id_test = sort(sample(N,n.RS))
  test[id_test] = 1
  n_test = length(id_test)
  
  # generate the individual imperfect test results by the given Se and Sp
  testpos = rep(0,N)
  testpos[id_test] = N_true[id_test]
  id_testpos = which(testpos == 1)
  n_testpos =length(id_testpos)
  testpos[which(test == 1 & testpos == 0)] = rbinom(n_test-n_testpos,1,1-Sp)
  testpos[id_testpos] = rbinom(n_testpos,1,Se)
  
  return(list(test = test,testpos = testpos))
}


## A function for the proposed Bayesian Credible interval approach

RS_BC = function(N,n_pos.RS,n.RS,Se,Sp){
  
  ## Parameters
  ## N: total population size (N=1, output the BC for prevalence estimator;
  ##    N=Ntot, output the BC for case count estimator)
  ## n_pos.RS: number of positive test results from the imperfect test
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # set the number of posterior samples for BC interval
  m = 1000
  
  # calculate the estimate of test positive frequency ("pi" in paper) and the 
  # bias-corrected disease prevalence ("pi_c" in paper) with threshold [0, 1]
  p_star_RS = min(max(n_pos.RS/n.RS,0.0001), 0.9999)     # 1>= p_star_RS >= 0
  p_RS = min(max((p_star_RS+Sp-1)/(Se+Sp-1),0.0001), 0.9999)  # 1>= p_RS >= 0
  
  # V1(pi) in eqn.(1)
  V_p_star_RS = p_star_RS*(1-p_star_RS)/n.RS
  
  # V3(pi) in eqn.(12)
  fpc = n.RS*(N-n.RS)/(N*(n.RS-1)) 
  V_p_RS = fpc*V_p_star_RS+(p_RS*Se*(1-Se) + (1-p_RS)*Sp*(1-Sp))/N
  
  # calculate the scale and shift parameter a and b
  a=sqrt(V_p_RS/V_p_star_RS)
  b=p_star_RS*(1-a)
  
  # generate beta posterior samples and adjust by a and b
  p_star_post = rbeta(m,n_pos.RS+1/2,n.RS-n_pos.RS+1/2)
  p_true_post = a*p_star_post + b
  p_true_post = (p_true_post+Sp-1)/(Se+Sp-1)
  
  # calculate the BC interval
  # Ntot=1, we calculate the BC for prevalence estimator
  # Ntot=N, we calculate the BC for case count estimator
  Ntot = 1
  N_iter = Ntot*p_true_post
  
  N_iter_lower = max(quantile(N_iter,c(0.025)),0)
  N_iter_upper = max(quantile(N_iter,c(0.975)),0)
  
  N_iter_interval_width = N_iter_upper - N_iter_lower
  
  N_iter_median = max(quantile(N_iter,c(0.5)),0)
  
  return(list(BC_lower = N_iter_lower,BC_upper = N_iter_upper,BC_width = N_iter_interval_width,
              BC_median = N_iter_median))
}


## A function to summarize simulation results

summary_stats_wald = function(est.N,est.sd,Npos){
  
  ## Parameters
  ## est.N: estimated N 
  ## est.sd: estimated standard error
  ## Npos: number of true diseased people
  
  # collect all estimates and sd, and calc wald type CI, width
  df = data.frame(est.N,est.sd)
  df2 = apply(df,1,coverage_wald,N_truth = Npos)
  col_names = names(unlist(df2)[1:4])
  df2 = t(matrix(unlist(df2),nrow=4))
  df2 = as.data.frame(df2)
  colnames(df2) = col_names
  
  results2 = unlist(apply(df2,2,mean))
  
  N.mean = mean(est.N,na.rm = T)
  N.sd = sd(est.N,na.rm = T)
  N.avgse = mean(est.sd,na.rm = T)
  N.width = as.numeric(results2[1])
  CI.pct = as.numeric(results2[2])*100
  
  rst.ls = list(N.mean = N.mean, N.sd = N.sd, N.avgse = N.avgse,
                N.width = N.width, CI.pct = CI.pct)
  
  return(unlist(rst.ls))
}


## A function to evaluate the coverage for wald-type confidence interval

coverage_wald = function(est.vec,N_truth){
  
  ## Parameters
  ## est.vec: estimation vector
  ## N_truth: number of true diseased people
  
  N = est.vec[1]
  N_sd = est.vec[2]
  
  upper = N + 1.96*N_sd
  lower = N - 1.96*N_sd
  
  if(is.na(upper)){   # remove NA case
    N_coverage = NA
    N_width = NA
    
  }else{
    if(upper >= N_truth & lower <= N_truth){
      N_coverage = 1
    }else{
      N_coverage = 0
    }
    N_width = upper - lower
  }
  
  return(list(width = N_width,coverage = N_coverage, lower = lower, upper = upper))
}


## A function to calculate three types of standard errors

RS_SE_compare = function(N, p_case, p.RS, Se, Sp){
  
  ## Parameters
  ## N: total population size
  ## p_case: true disease prevalence
  ## p.RS: the sampling rate for random samples
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # first N*p_case individual: case; remaining: health
  N_true = c(rep(1,round(N*p_case)),rep(0,round(N*(1-p_case)))) 
  
  # generate data set
  simu2 = simu_sym_RS(N,N_true,p.RS,Se,Sp)
  test2 = simu2$test
  testpos2 = simu2$testpos
  
  # RS method
  n_stream2 = sum(test2) 
  n_pos_stream2 = max(sum(testpos2),0.001)
  
  # calculate the estimate of test positive frequency ("pi" in paper) and the 
  # bias-corrected disease prevalence ("pi_c" in paper) with threshold [0, 1]
  p_star_RS = n_pos_stream2/n_stream2
  p_RS = max((p_star_RS+Sp-1)/(Se+Sp-1),0)
  
  # three types of standard errors
  fpc2 = n_stream2*(N-n_stream2)/(N*(n_stream2-1))
  V_pi_star0 = p_star_RS*(1-p_star_RS)/n_stream2
  V_pi_star1 = V_pi_star0*fpc2
  V_pi_star2 = (p_RS*Se*(1-Se) + (1-p_RS)*Sp*(1-Sp))/N
  
  p_RS_se1 = 1/(Se+Sp-1)*sqrt(V_pi_star0)
  p_RS_Se2 = 1/(Se+Sp-1)*sqrt(V_pi_star1)
  p_RS_se3 = 1/(Se+Sp-1)*sqrt(V_pi_star1+V_pi_star2)
  
  return(c(p_RS,p_RS_se1,p_RS_Se2,p_RS_se3))
}


## A function to evaluate the three types of averaged standard errors and "true variance"

RS_SE_compare_main = function(N, p_case, p.RS, Se, Sp){
  
  ## Parameters
  ## N: total population size
  ## p_case: true disease prevalence
  ## p.RS: the sampling rate for random samples
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  n_sim = 50000
  result = replicate(n=n_sim, RS_SE_compare(N, p_case, p.RS, Se, Sp))
  p_RS_sd=sd(result[1,])
  p_RS_se1=mean(result[2,])
  p_RS_Se=mean(result[3,])
  p_RS_se3=mean(result[4,])
  return(c(p_RS_sd,p_RS_se1,p_RS_Se,p_RS_se3))
}