##############################################################################
#                                                                            #
#                                                                            #
#  Simulation function used for "Enhanced Inference for Finite Population    #
#  Sampling-Based Prevalence Estimation with Misclassification Errors        #
#                                                                            #
#                                                                            #
##############################################################################

## A function to generate simulation results for the manuscript

find_p_p2 = function(N, p_case, p.RS,Se,Sp){

  ## Parameters
  ## N: total population size
  ## p_case: true disease prevalence
  ## p.RS: sampling rate for random sample
  ## Se: Sensitivity of the testing tool
  ## Sp: Specificity of the testing tool
  
  # set the number of the simulation runs
  n_sim = 5000
  
  # set the individual level true disease status
  N_true = c(rep(1,round(N*p_case)),rep(0,round(N*(1-p_case)))) 
  
  Npos = sum(N_true)
  Nneg = N-Npos
  True_val = Npos/N
  
  # store simulation results
  N_RS = rep(0,n_sim)
  N_RS_sd = rep(0,n_sim)
  N_RS_sd2 = rep(0,n_sim)
  
  BC_interval_coverage = rep(0,n_sim)
  BC_interval_coverage_025 = rep(0,n_sim)
  BC_interval_coverage_975 = rep(0,n_sim)
  BC_interval_length = rep(0,n_sim)
  BC_median = rep(0,n_sim)
  
  for(i in 1:n_sim){
    
    # simu RS steam
    simu.RS = simu_sym_RS(N,N_true,p.RS,Se,Sp)
    test.RS = simu.RS$test
    testpos.RS = simu.RS$testpos
    
    # summary for RS method
    n.RS = sum(test.RS) 
    n_pos.RS = max(sum(testpos.RS),0.001)
    
    # Ntot=1, we calculate for prevalence estimator
    # Ntot=N, we calculate for case count estimator
    Ntot = 1
    
    # calculate the estimate of test positive frequency ("pi" in paper) and the 
    # bias-corrected disease prevalence ("pi_c" in paper) with threshold [0, 1]
    p_star_RS = n_pos.RS/n.RS
    p_RS = max((p_star_RS+Sp-1)/(Se+Sp-1),0)
    N_RS[i] = Ntot*p_RS
    
    # V2(pi) in eqn.(2)
    fpc2 = n.RS*(N-n.RS)/(N*(n.RS-1)) 
    V_pi_star1 = p_star_RS*(1-p_star_RS)/n.RS*fpc2
    V_pi_star2 = (p_RS*Se*(1-Se) + (1-p_RS)*Sp*(1-Sp))/N
    
    # V3(pi) in eqn.(12) and insert into eqn.(6)
    se_pi_RS = 1/(Se+Sp-1)*sqrt(V_pi_star1+V_pi_star2)
    
    N_RS_sd[i] = Ntot*se_pi_RS
    
    # BC interval calculation
    BC_interval = RS_BC(N, n_pos.RS,n.RS,Se,Sp)
    if(BC_interval$BC_upper >= True_val && BC_interval$BC_lower <= True_val){
      BC_interval_coverage[i] = 1
    }else{
      BC_interval_coverage[i] = 0
    }  
    
    if(BC_interval$BC_upper < True_val){  # upper 2.5% prob
      BC_interval_coverage_975[i] = 1
    }else{
      BC_interval_coverage_975[i] = 0
    } 
    
    if(BC_interval$BC_lower > True_val){  # lower 2.5% prob
      BC_interval_coverage_025[i] = 1
    }else{
      BC_interval_coverage_025[i] = 0
    } 
    BC_interval_length[i] = BC_interval$BC_width
    
    BC_median[i] = BC_interval$BC_median
    
  }
  
  # summarize the simulation results
  var_list = c('RS')
  
  table_all = c(True_val,NA,NA,NA,NA)
  for(var_name in var_list){
    v1_est = eval(as.symbol(paste('N',var_name,sep = '_')))
    v1_sd = eval(as.symbol(paste('N',var_name,'sd',sep = '_')))
    rst_tmp = summary_stats_wald(v1_est,v1_sd,True_val)
    
    table_all = cbind(table_all,rst_tmp)
  }
  colnames(table_all) = c('N_true',paste('N',var_list,sep='_'))
  
  res.true = table_all[1,1]
  res.RS = table_all[,2]
  res.BC_median = mean(BC_median)
  res.BC_width = mean(BC_interval_length) 
  res.BC_pct = mean(BC_interval_coverage)*100
  
  return(list(res.N = N,res.true = res.true, res.Se = Se, res.Sp = Sp,
              res.RS = res.RS,res.BC_width = res.BC_width, res.BC_pct = res.BC_pct,
              res.BC_median = res.BC_median))
  
}