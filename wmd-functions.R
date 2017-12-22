# -----------------------------------------------------------------------------------------------------------
# Research: A comparison of meta-analysis methods: Understanding the influence of experiments' 
#           statistical parameters
# File name: wmd-functions.R
# File type: secondary R Script
# Date: March 2015
# R Script contributors: Oscar Dieste, Omar S. Gomez
# Purpose: functions to perform meta-analysis using WMD
# ----------------------------------------------------------------------------------------------------------
#


#
# tell_J performs the hedges' correction on d
tell_J <- function (N) (1-(3/((4*N)-9)))

#
# Functions that calculate the poopled variance and effect size
tell_pooled_variance <- 
  function (var_control, var_treatment, n_control, n_treatment) 
    ((((n_control-1)*var_control)+((n_treatment-1)*var_treatment))/(n_control+n_treatment-2))

tell_effect_size <- 
  function (mean_control, mean_treatment, var_control, var_treatment, n_control, n_treatment) {
    pooled_variance <- tell_pooled_variance(var_control, var_treatment, n_control, n_treatment)
    return(tell_J(n_control+n_treatment)*((mean_treatment-mean_control)/sqrt(pooled_variance)))
  }

#
# Variance of the estimador d
tell_variance_of_d <- 
  function (ds, ncs, nts) {
    results <- vector(length=length(ds))
    for (d in 1:length(ds)) results[d] <- ((nts[d]+ncs[d])/(nts[d]*ncs[d]))+((ds[d]^2)/(2*(nts[d]+ncs[d])))
    return(results)
  }

#
# Accuracy calculation
tell_wmd_accuracy <- 
  function (effect_sizes, ncs, nts, cutoff, effect) {
    weights <- 1/tell_variance_of_d(effect_sizes, ncs, nts)
    global_effect_size <- sum(weights*effect_sizes)/sum(weights)
    var_of_global_effect_size <- 1/sum(weights)
    left_side <- global_effect_size - sqrt(var_of_global_effect_size)*cutoff
    right_side <- global_effect_size + sqrt(var_of_global_effect_size)*cutoff   
    return(if((effect >= left_side) && (effect <= right_side)) 1 else 0)
}

#
# Empirical power calculation
tell_wmd_emp_power <- 
  function (effect_sizes, ncs, nts, cutoff) {
    weights <- 1/tell_variance_of_d(effect_sizes, ncs, nts)
    global_effect_size <- sum(weights*effect_sizes)/sum(weights)
    var_of_global_effect_size <- 1/sum(weights)
    left_side <- global_effect_size - sqrt(var_of_global_effect_size)*cutoff   
    return(if(0 < left_side) 1 else 0)
  }

#
# Between studies variance
tell_tau_squared_wmd <- 
  function (effect_sizes, weights, df) {
    q <- sum(weights*(effect_sizes^2))-((sum(weights*effect_sizes)^2)/sum(weights))
    c <- sum(weights) - (sum(weights^2)/sum(weights))
    if (q > df) return ((q-df)/c)
    else return (0)
  } 


#
# Random accuracy calculation
tell_random_wmd_accuracy <- 
  function (effect_sizes, ncs, nts, cutoff, effect) {
    weights <- 1/tell_variance_of_d(effect_sizes, ncs, nts)
    tau_squared <- tell_tau_squared_wmd (effect_sizes, weights, length(effect_sizes)-1)
    weights <- 1/(tell_variance_of_d(effect_sizes, ncs, nts) + tau_squared)
    global_effect_size <- sum(weights*effect_sizes)/sum(weights)
    var_of_global_effect_size <- 1/sum(weights)
    left_side <- global_effect_size - sqrt(var_of_global_effect_size)*cutoff
    right_side <- global_effect_size + sqrt(var_of_global_effect_size)*cutoff   
    return(if((effect >= left_side) && (effect <= right_side)) 1 else 0)
  }

#
# Random empirical power calculation
tell_random_wmd_emp_power <- 
  function (effect_sizes, ncs, nts, cutoff) {
    weights <- 1/tell_variance_of_d(effect_sizes, ncs, nts)
    tau_squared <- tell_tau_squared_wmd (effect_sizes, weights, length(effect_sizes)-1)
    weights <- 1/(tell_variance_of_d(effect_sizes, ncs, nts) + rep (tau_squared, length(weights)))
    global_effect_size <- sum(weights*effect_sizes)/sum(weights)
    var_of_global_effect_size <- 1/sum(weights)
    left_side <- global_effect_size - sqrt(var_of_global_effect_size)*cutoff   
    return(if(0 < left_side) 1 else 0)
  }
