# ----------------------------------------------------------------------------------------------------------------
# Research             : A procedure for selecting an appropriate meta-analysis method: A comparison of approaches
# File name            : rr-functions.R
# File type            : Auxiliary functions
# Purpose              : Functions to perform meta-analysis using response ratios
# Creation date        : March 2015 (later history available at GitHub)
# R Script contributors: Oscar Dieste, Omar S. Gomez
# ----------------------------------------------------------------------------------------------------------------
#


# Parametric & non-parametric response ratio
#
tell_rr <- function (mean_control, mean_treatment, uselogs) {
  if(uselogs) {
    return(log(mean_treatment/mean_control))
  } else {
    return(mean_treatment/mean_control)
  }
}


#
# Variance of the estimador prr
tell_variance_of_prr <- 
  function (mean_control, mean_treatment, var_control, var_treatment, n_control, n_treatment) 
    ((var_treatment/(n_treatment*mean_treatment^2))+(var_control/(n_control*mean_control^2)))

#
# Variance of the estimador nprr
tell_variance_of_nprr <- 
  function (n_control, n_treatment) ((n_control+n_treatment)/(n_control*n_treatment))

#
# Parametric & non-parametric response ratio accuracy calculation
tell_prr_accuracy <- 
  function (prrs, prr_vars, cutoff, control_mean, treatment_mean, uselogs) {
    rr <- tell_rr(control_mean, treatment_mean, uselogs)
    weights <- 1/prr_vars
    global_prr <- sum(weights*prrs)/sum(weights)
    var_of_global_prr <- 1/sum(weights)
    left_side <- global_prr - sqrt(var_of_global_prr)*cutoff
    right_side <- global_prr + sqrt(var_of_global_prr)*cutoff   
    return(if((rr >= left_side) && (rr <= right_side)) 1 else 0)
  }

tell_nprr_accuracy <- function (nprrs, nprr_vars, cutoff, control_mean, treatment_mean, uselogs) (tell_prr_accuracy (nprrs, nprr_vars, cutoff, control_mean, treatment_mean, uselogs))
  

#
# Parametric & non-parametric response ratio accuracy calculation empirical power calculation
tell_prr_emp_power <- 
  function (prrs, prr_vars, cutoff) {
    weights <- 1/prr_vars
    global_prr <- sum(weights*prrs)/sum(weights)
    var_of_global_prr <- 1/sum(weights)
    left_side <- global_prr - sqrt(var_of_global_prr)*cutoff
    return(if(0 < left_side) 1 else 0)
  }

tell_nprr_emp_power <- function (nprrs, nprr_vars, cutoff) (tell_prr_emp_power (nprrs, nprr_vars, cutoff))

#
# Between studies variance
tell_tau_squared_prr <- 
  function (prrs, weights, df) {
    q <- sum(weights*(prrs^2))-((sum(weights*prrs)^2)/sum(weights))
    c <- sum(weights) - (sum(weights^2)/sum(weights))
    if (q > df) return ((q-df)/c)
    else return (0)
  } 


#
# Random accuracy calculation
tell_random_prr_accuracy <- 
  function (prrs, prr_vars, cutoff, control_mean, treatment_mean, uselogs) {
    rr <- tell_rr(control_mean, treatment_mean, uselogs)
    weights <- 1/prr_vars
    tau_squared <- tell_tau_squared_prr (prrs, weights, length(prrs)-1)
    weights <- 1/(prr_vars + rep (tau_squared, length(weights)))
    global_prr <- sum(weights*prrs)/sum(weights)
    var_of_global_prr <- 1/sum(weights)
    left_side <- global_prr - sqrt(var_of_global_prr)*cutoff
    right_side <- global_prr + sqrt(var_of_global_prr)*cutoff   
    return(if((rr >= left_side) && (rr <= right_side)) 1 else 0)
  }

#
# Random empirical power calculation
tell_random_prr_emp_power <- 
  function (prrs, prr_vars, cutoff) {
    weights <- 1/prr_vars
    tau_squared <- tell_tau_squared_prr (prrs, weights, length(prrs)-1)
    weights <- 1/(prr_vars + rep (tau_squared, length(weights)))
    global_prr <- sum(weights*prrs)/sum(weights)
    var_of_global_prr <- 1/sum(weights)
    left_side <- global_prr - sqrt(var_of_global_prr)*cutoff
    return(if(0 < left_side) 1 else 0)
  }
