# -----------------------------------------------------------------------------------------------------------
# Research: A comparison of meta-analysis methods: Understanding the influence of experiments' 
#           statistical parameters
# File name: svc-functions.R
# File type: secondary R Script
# Date: March 2015
# R Script contributors: Oscar Dieste, Omar S. Gomez
# Purpose: functions to perform meta-analysis using SVC
# ----------------------------------------------------------------------------------------------------------
#


#
# The vote from a single study 
tell_vote <-
  function (control, treatment, level) {
    p <- t.test(x=control, y=treatment, conf.level=level)["p.value"]
    
    if (p > 0.05) return (0)
    else return(1)   
  }


#
# Auxiliary functions
calculate_nis <-
  function (ncs, nts)
    return ((ncs*nts)/(ncs+nts))

likelihood_function <-
  function (votes, nis, delta) (sum((votes*log(1-pnorm(-1*sqrt(nis)*delta))) + ((1-votes)*log(pnorm(-1*sqrt(nis)*delta)))))
    

tell_Di <-
  function (nis, delta) (sqrt(nis/(2*pi))*exp(-0.5*nis*(delta^2)))


tell_pi <-
  function (nis, delta) (1-pnorm(-1*sqrt(nis)*delta))


#
# Calculate the most probable effect size
svc_calculate_delta <-
  function (votes, ncs, nts){
    nis <- calculate_nis(ncs, nts)
    int <- 1
    
    #
    # coarse interaction
    max_l <- -1 * .Machine$double.xmax
    value_delta <- -10
    num_values <- (-1*2*value_delta*(1/int))+1
    deltas <- vector(length=num_values)
    deltas[1] = value_delta
    for (i in 2:num_values) deltas[i] <- deltas[i-1]+int
    for (i in 1:num_values) {
      l <- likelihood_function(votes, nis, deltas[i])
      if (!is.nan(l)) {
        if (max_l < l) {
          max_l <- l
          value_delta <- deltas[i]
        }
      }
    }
    
    #
    # finer iteration
    max_l <- -1 * .Machine$double.xmax
    value_delta <- value_delta - int
    num_values <- (100 * int * 2) + 1
    deltas <- vector(length=num_values)
    deltas[1] = value_delta
    for (i in 2:num_values) deltas[i] <- deltas[i-1]+0.01
    for (i in 1:num_values) {
      l <- likelihood_function(votes, nis, deltas[i])
      if (!is.nan(l)) {
        if (max_l < l) {
          max_l <- l
          value_delta <- deltas[i]
        }
      }
    }
    
    return(value_delta)
  }


#
# Accuracy calculation
tell_svc_accuracy <- 
  function (ncs, nts, delta, cutoff,d,list_effect_sizes) {
    nis <- calculate_nis(ncs, nts)
    svc_var <- 1/sum((tell_Di(nis, delta)^2)/(tell_pi(nis, delta)*(1-tell_pi(nis, delta))))
    left_side <- delta - sqrt(svc_var)*cutoff
    right_side <- delta + sqrt(svc_var)*cutoff 
    return(if((list_effect_sizes[d] >= left_side) && (list_effect_sizes[d] <= right_side)) 1 else 0)
  }

#
# Empirical power calculation
tell_svc_emp_power <- 
  function (ncs, nts, delta, cutoff) {
    nis <- calculate_nis(ncs, nts)
  
    svc_var <- 1/sum((tell_Di(nis, delta)^2)/(tell_pi(nis, delta)*(1-tell_pi(nis, delta))))
    
    left_side <- delta - sqrt(svc_var)*cutoff 

    return(if(0 < left_side) 1 else 0)
  }

