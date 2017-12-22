# -----------------------------------------------------------------------------------------------------------
# Research: A comparison of meta-analysis methods: Understanding the influence of experiments' 
#           statistical parameters
# File name: auxiliary_unequal.R
# File type: auxiliary R Script
# Date: March 2015
# R Script contributors: Oscar Dieste, Omar S. Gomez
# Purpose: auxiliary functions
# ----------------------------------------------------------------------------------------------------------
#

#
# Output the header
output_header <-
  function (f) {
    cat ("measure", "method", "d", "sdv", "tau", "nc", "nt", "num_exp", "value", file = f, sep = ",") 
    cat("\n", file=f)
  }

#
# Output a row to a file
output_line <-
  function (f, measure, method, d, sdv, t, nc, nt, list_experiments, list_values) {
    for (i in 1:length(list_experiments)){
      cat (measure,",",method,",",d,",",sdv,",",t,"'",nc[i,],"'",",","'",nt[i,],"'",list_experiments[i],",",list_values[i], file = f)
      cat("\n", file=f)
    }   
  }