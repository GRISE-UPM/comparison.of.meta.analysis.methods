# -----------------------------------------------------------------------------------------------------------
# Research: A comparison of meta-analysis methods: Understanding the influence of experiments' 
#           statistical parameters
# File name: ctrees_script.R
# File type: secondary R Script
# Date: March 2016
# R Script main contributor: Oscar Dieste
# R Script secondary contributor: Omar S. Gómez
# Purporse: Find out optimal regions of the parameter space for the meta-analysis methods
# ----------------------------------------------------------------------------------------------------------
#


# Library for classification trees
library(party)
library(Hmisc)


## Clasification Trees, set the mysql user and password 
print_ct <- function(measure,method,filename){
  con <- dbConnect(MySQL(), user="db_usr", password="password", dbname="dbsimulation", host="localhost")
  rs <- dbSendQuery(con, paste0("select * from methods where measure = '",measure,"' and method = '",method,"';"))
  ds <- fetch(rs, n=27000)
  
  ct <- ctree(value ~ d + sdv + tau + nc + num_exp, data = ds, controls = ctree_control(maxdepth = 4))  
  
  postscript(filename, width = 13.3, height = 6, paper = "special", family = "Helvetica")
  plot(ct, type = "simple")
  dev.off()
  
  huh <- dbHasCompleted(rs)
  dbClearResult(rs)
  dbDisconnect(con)
  
  return(TRUE)
}

m<-"rejection_rate"
print_ct( m, "wmd", "wmd_rejection_rate.png")
print_ct( m, "random-wmd", "random-wmd_rejection_rate.png")
print_ct( m, "noln_prr","noln_prr_rejection_rate.png")
print_ct( m, "ln_prr","ln_prr_rejection_rate.png")
print_ct( m, "noln_random-prr","noln_random-prr_rejection_rate.png")
print_ct( m, "ln_random-prr","ln_random-prr_rejection_rate.png")
print_ct( m, "noln_nprr","noln_nprr_rejection_rate.png")
print_ct( m, "ln_nprr","ln_nprr_rejection_rate.png")
print_ct( m, "svc","svc_rejection_rate.png")

m<-"accuracy"
print_ct( m, "wmd","wmd_accuracy.ps")
print_ct( m, "random-wmd","random-wmd_accuracy.png")
print_ct( m, "noln_prr","noln_prr_accuracy.png")
print_ct( m, "ln_prr","ln_prr_accuracy.png")
print_ct( m, "noln_random-prr","noln_random-prr_accuracy.png")
print_ct( m, "ln_random-prr","ln_random-prr_accuracy.png")
print_ct( m, "noln_nprr","noln_nprr_accuracy.png")
print_ct( m, "ln_nprr","ln_nprr_accuracy.png")
print_ct( m, "svc","svc_accuracy.png")

m<-"empirical_power"
print_ct( m, "wmd","wmd_empirical_power.png")
print_ct( m, "random-wmd","random-wmd_empirical_power.png")
print_ct( m, "noln_prr","noln_prr_empirical_power.png")
print_ct( m, "ln_prr","ln_prr_empirical_power.png")
print_ct( m, "noln_random-prr","noln_random-prr_empirical_power.png")
print_ct( m, "ln_random-prr","ln_random-prr_empirical_power.png")
print_ct( m, "noln_nprr","noln_nprr_empirical_power.png")
print_ct( m, "ln_nprr","ln_nprr_empirical_power.png")
print_ct( m, "svc","svc_empirical_power.png")
