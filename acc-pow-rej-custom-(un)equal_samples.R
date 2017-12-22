# -----------------------------------------------------------------------------------------------------------
# Research: A comparison of meta-analysis methods: Understanding the influence of experiments' 
#           statistical parameters
# File name: acc-pow-rej-custom-(un)equal_samples.R
# File type: main R Script
# Date: March 2016
# R Script main contributor: Oscar Dieste
# R Script secondary contributor: Omar S. Gómez
# Purporse: simulation of meta-analysis techniques for rejection rate, accuracy and emperical power
#           tested under different statistical parameters combinations.
# ----------------------------------------------------------------------------------------------------------
#


#
# Auxiliary functions

# Set simulation R-scripts working directory 


#setwd("/media/omar/Windows7_OS/SysPart/Default/Users/Omar/Dropbox/espoch/papers/tse_simulation/r_script/journal_script_simulation")


source('wmd-functions.R')
source('rr-functions.R')
source('svc-functions.R')
source('auxiliary_unequal.R')


# Generation of randonm samples sizes for each study (balanced studies)
random_sample <- 
  function (e, list_ss) {
    rs_v <- vector(length=length(e))
    for (x in 1:e){
      if(length(list_ss)==1)
        rs_v[x] <- list_ss[1]
      else
        rs_v[x] <- sample(list_ss,1)
    }  
    return(rs_v)
  }

# Simulation function 
# acc, pow, rej are boolean variables intended to be included or excluded in the simulation, 
# it similarly applies with wmd,random wmd,noln_prr,ln_prr,noln_random prr,ln_random prr,noln_nprr,ln_nprr,svc
# variables
montecarlo<-function(acc,pow,rej,m1,m2,m3,m4,m5,m6,m7,m8,m9,repetitions,alpha,cutoff,
                     list_experiments,list_unequal_sample_sizes,rnd_unbalanced,fixed_unbalanced_factor,fixed_unbalanced,
                     list_effect_sizes,real_control_mean,list_standard_deviations,list_taus,fname,mainDir,subDir,ps){
  #
  # set output directory
  
dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))
  
  #
  # Variables with varying sizes. They are predefined to the max length
  tmp_sample_effect_sizes <- vector(length=length(list_experiments))
  noln_tmp_sample_rrs     <- vector(length=length(list_experiments))
  ln_tmp_sample_rrs       <- vector(length=length(list_experiments))
  tmp_sample_votes        <- vector(length=length(list_experiments))
  
  noln_tmp_prr_vars       <- vector(length=length(list_experiments))
  noln_tmp_nprr_vars      <- vector(length=length(list_experiments))
  ln_tmp_prr_vars         <- vector(length=length(list_experiments))
  ln_tmp_nprr_vars        <- vector(length=length(list_experiments))
  
  
  wmd_accuracy            <- matrix (nrow=length(list_experiments), ncol=repetitions)
  random_wmd_accuracy     <- matrix (nrow=length(list_experiments), ncol=repetitions)
  noln_prr_accuracy       <- matrix (nrow=length(list_experiments), ncol=repetitions)
  noln_random_prr_accuracy<- matrix (nrow=length(list_experiments), ncol=repetitions)
  noln_nprr_accuracy      <- matrix (nrow=length(list_experiments), ncol=repetitions)
  ln_prr_accuracy         <- matrix (nrow=length(list_experiments), ncol=repetitions)
  ln_random_prr_accuracy  <- matrix (nrow=length(list_experiments), ncol=repetitions)
  ln_nprr_accuracy        <- matrix (nrow=length(list_experiments), ncol=repetitions)
  svc_accuracy            <- matrix (nrow=length(list_experiments), ncol=repetitions)
  tmp_accuracy            <- vector(length=length(list_experiments))
  
  
  #
  # Depending on wether d == 0 or not, *_emp_power can be used to calculate rejection rate or empirical power 
  wmd_emp_power           <- matrix (nrow=length(list_experiments), ncol=repetitions)
  random_wmd_emp_power    <- matrix (nrow=length(list_experiments), ncol=repetitions)
  noln_prr_emp_power      <- matrix (nrow=length(list_experiments), ncol=repetitions)
  noln_random_prr_emp_power<- matrix (nrow=length(list_experiments), ncol=repetitions)
  noln_nprr_emp_power     <- matrix (nrow=length(list_experiments), ncol=repetitions)
  ln_prr_emp_power        <- matrix (nrow=length(list_experiments), ncol=repetitions)
  ln_random_prr_emp_power <- matrix (nrow=length(list_experiments), ncol=repetitions)
  ln_nprr_emp_power       <- matrix (nrow=length(list_experiments), ncol=repetitions)
  svc_emp_power           <- matrix (nrow=length(list_experiments), ncol=repetitions)
  tmp_emp_power           <- vector(length=length(list_experiments))
  
  #TODO
  #
  # Open the output file
  if(!file.exists(fname))
       file.create(fname)

  f<-file(fname, "w")

  output_header(f)
  
  unequal_samples_c<-matrix(nrow=length(list_experiments),ncol=max(list_experiments))
  
 # Generate a matrix with random samples per each element of the experiments list
  for (i in 1:length(list_experiments)) {    
    aux<-random_sample(list_experiments[i],list_unequal_sample_sizes)
    for (j in 1:length(aux))
      unequal_samples_c[i,j]<-aux[j]
  }  
  
 unequal_samples_t<-unequal_samples_c

 #in case of true add a random unbalanced factor
  if(rnd_unbalanced){
    for (i in 1:length(list_experiments))     
      for (j in 1:list_experiments[i]){
        unequal_samples_c[i,j]<-unequal_samples_c[i,j]+round(unequal_samples_c[i,j]*runif(1,-fixed_unbalanced_factor,fixed_unbalanced_factor))
        unequal_samples_t[i,j]<-unequal_samples_t[i,j]+round(unequal_samples_t[i,j]*runif(1,-fixed_unbalanced_factor,fixed_unbalanced_factor))
      }
  }
 
  if(fixed_unbalanced){
    for (i in 1:length(list_experiments))     
      for (j in 1:list_experiments[i])
        unequal_samples_t[i,j]<-unequal_samples_t[i,j]+round(unequal_samples_t[i,j]*fixed_unbalanced_factor)    
  }
 
 
  #
  # Simulation loops
  for (d in 1:length(list_effect_sizes)) {
    
    for (sdv in 1:length(list_standard_deviations)) {
      
      for (t in 1:length(list_taus)) {
            
            real_treatment_mean <- real_control_mean+(list_effect_sizes[d]*list_standard_deviations[sdv])
            
            start_time <- Sys.time()
            
            
            for (e in 1:length(list_experiments)) {
              cat("with ",list_experiments[e]," experiments: \n")
              for (y in 1:list_experiments[e])  
                cat("     exp. [",y,"] with ", unequal_samples_c[e,y],",",unequal_samples_t[e,y],"(c,t) subjects \n")
              
              for (x in 1:repetitions) {
                
                #
                # This loop is for wmd, random_wmd and svc, because they don't cause problems with the signs
                for (y in 1:list_experiments[e]) {
                  
                  #
                  # Generation of the samples
                  heterogeneity <- rnorm(1, mean=0, sd=list_taus[t])
                  control <- rnorm(unequal_samples_c[e,y], mean=real_control_mean, sd=list_standard_deviations[sdv]) + heterogeneity
                  treatment <- rnorm(unequal_samples_t[e,y], mean=real_treatment_mean, sd=list_standard_deviations[sdv]) + heterogeneity
                  
                  #
                  # Calculation of the groups' means and variances
                  mean_control <- mean(control)
                  mean_treatment <- mean(treatment)
                  var_control <- var(control)
                  var_treatment <- var(treatment)
                  
                  #
                  # Effect size calculation
                  tmp_sample_effect_sizes [y] <- tell_effect_size (mean_control, mean_treatment, var_control, var_treatment, unequal_samples_c[e,y], unequal_samples_t[e,y]) 
                  
                  #
                  # Vote counting calculation
                  tmp_sample_votes [y] <- tell_vote (control, treatment, 1-alpha)     
                  
                  
                  uselogs <- FALSE
                  noln_tmp_sample_rrs [y] <- tell_rr (mean_control, mean_treatment, uselogs)
                  noln_tmp_prr_vars [y] <- tell_variance_of_prr (mean_control, mean_treatment, var_control, var_treatment, unequal_samples_c[e,y], unequal_samples_t[e,y])
                  noln_tmp_nprr_vars [y] <- tell_variance_of_nprr (unequal_samples_c[e,y], unequal_samples_t[e,y])
                }
                
                #
                # This loop is for prr, random_prr and nprr, because there are not ln's of negative numbers
                uselogs <- TRUE
                y<- 1
                while (y <= list_experiments[e]) {             
                  #
                  # Generation of the samples
                  heterogeneity <- rnorm(1, mean=0, sd=list_taus[t])
                  control <- rnorm(unequal_samples_c[e,y], mean=real_control_mean, sd=list_standard_deviations[sdv]) + heterogeneity
                  treatment <- rnorm(unequal_samples_t[e,y], mean=real_treatment_mean, sd=list_standard_deviations[sdv]) + heterogeneity
                  
                  #
                  # Calculation of the groups' means and variances
                  mean_control <- mean(control)
                  mean_treatment <- mean(treatment)
                  var_control <- var(control)
                  var_treatment <- var(treatment)
                  
                  #
                  # In case that a zero/negative value in either mean_control or mean_treatment shows up
                  # we skip that interaction for prr and nprr
                  #
                  # We count the number of times we skip to adjust the accuracy and emp_power calculations later
                  if ((mean_control > 0) && (mean_treatment > 0)) {
                    #
                    # Response ratio calculation
                    ln_tmp_sample_rrs [y] <- tell_rr (mean_control, mean_treatment, uselogs)
                    ln_tmp_prr_vars [y] <- tell_variance_of_prr (mean_control, mean_treatment, var_control, var_treatment, unequal_samples_c[e,y], unequal_samples_t[e,y])
                    ln_tmp_nprr_vars [y] <- tell_variance_of_nprr (unequal_samples_c[e,y], unequal_samples_t[e,y])
                    y <- y + 1
                  }
                }
     
                
                if (acc | pow | rej){
                  #
                  # WMD's accuracy & empirical power calculation
                  if(m1){
                    wmd_accuracy[e, x] <- tell_wmd_accuracy(tmp_sample_effect_sizes [1:list_experiments[e]], unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]], cutoff, list_effect_sizes[d])
                    wmd_emp_power[e, x] <- tell_wmd_emp_power(tmp_sample_effect_sizes [1:list_experiments[e]], unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]], cutoff)
                  }
                  
                  if(m2){
                    random_wmd_accuracy[e, x] <- tell_random_wmd_accuracy(tmp_sample_effect_sizes [1:list_experiments[e]], unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]], cutoff, list_effect_sizes[d])
                    random_wmd_emp_power[e, x] <- tell_random_wmd_emp_power(tmp_sample_effect_sizes [1:list_experiments[e]], unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]], cutoff)
                  }
                  
                  #
                  # PRR's accuracy & empirical power calculation (without logs)
                  if(m3){
                    uselogs <- FALSE
                    noln_prr_accuracy[e, x] <- tell_prr_accuracy (noln_tmp_sample_rrs[1:list_experiments[e]], noln_tmp_prr_vars[1:list_experiments[e]], cutoff, real_control_mean, real_treatment_mean, uselogs)
                    noln_prr_emp_power[e, x] <- tell_prr_emp_power (noln_tmp_sample_rrs[1:list_experiments[e]], noln_tmp_prr_vars[1:list_experiments[e]], cutoff)
                  }
                  
                  if(m5){
                    uselogs <- FALSE
                    noln_random_prr_accuracy[e, x] <- tell_random_prr_accuracy (noln_tmp_sample_rrs[1:list_experiments[e]], noln_tmp_prr_vars[1:list_experiments[e]], cutoff, real_control_mean, real_treatment_mean, uselogs)
                    noln_random_prr_emp_power[e, x] <- tell_random_prr_emp_power (noln_tmp_sample_rrs[1:list_experiments[e]], noln_tmp_prr_vars[1:list_experiments[e]], cutoff)   
                  }
                  
                  #
                  # PRR's accuracy & empirical power calculation (with logs)
                  if(m4){
                    uselogs <- TRUE
                    ln_prr_accuracy[e, x] <- tell_prr_accuracy (ln_tmp_sample_rrs[1:list_experiments[e]], ln_tmp_prr_vars[1:list_experiments[e]], cutoff, real_control_mean, real_treatment_mean, uselogs)
                    ln_prr_emp_power[e, x] <- tell_prr_emp_power (ln_tmp_sample_rrs[1:list_experiments[e]], ln_tmp_prr_vars[1:list_experiments[e]], cutoff)
                  }
                  
                  if(m6){
                    uselogs <- TRUE
                    ln_random_prr_accuracy[e, x] <- tell_random_prr_accuracy (ln_tmp_sample_rrs[1:list_experiments[e]], ln_tmp_prr_vars[1:list_experiments[e]], cutoff, real_control_mean, real_treatment_mean, uselogs)
                    ln_random_prr_emp_power[e, x] <- tell_random_prr_emp_power (ln_tmp_sample_rrs[1:list_experiments[e]], ln_tmp_prr_vars[1:list_experiments[e]], cutoff)
                  }
                  
                  #
                  # NPRR's accuracy & empirical power calculation (without logs)
                  if(m7){
                    uselogs <- FALSE
                    noln_nprr_accuracy[e, x] <- tell_nprr_accuracy (noln_tmp_sample_rrs[1:list_experiments[e]], noln_tmp_nprr_vars[1:list_experiments[e]], cutoff, real_control_mean, real_treatment_mean, uselogs)
                    noln_nprr_emp_power[e, x] <- tell_nprr_emp_power (noln_tmp_sample_rrs[1:list_experiments[e]], noln_tmp_nprr_vars[1:list_experiments[e]], cutoff)
                  }
                  
                  #
                  # NPRR's accuracy & empirical power calculation (with logs)
                  if(m8){
                    uselogs <- TRUE
                    ln_nprr_accuracy[e, x] <- tell_nprr_accuracy (ln_tmp_sample_rrs[1:list_experiments[e]], ln_tmp_nprr_vars[1:list_experiments[e]], cutoff, real_control_mean, real_treatment_mean, uselogs)
                    ln_nprr_emp_power[e, x] <- tell_nprr_emp_power (ln_tmp_sample_rrs[1:list_experiments[e]], ln_tmp_nprr_vars[1:list_experiments[e]], cutoff)
                  }            
                  
                  #
                  # SVC's accuracy & empirical power calculation
                  if(m9){
                    delta = svc_calculate_delta(tmp_sample_votes[1:list_experiments[e]], unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]])                 
                    svc_accuracy[e, x] <- tell_svc_accuracy(unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]], delta, cutoff,d,list_effect_sizes)
                    svc_emp_power[e, x] <- tell_svc_emp_power(unequal_samples_c[e,1:list_experiments[e]], unequal_samples_t[e,1:list_experiments[e]], delta, cutoff)
                  }      
                }
                
              }
            }
            
              aux_nc<-"rnd"
              aux_nt<-"rnd"
            
            if(fixed_unbalanced || length(list_unequal_sample_sizes)==1 ){
              aux_nc<-unequal_samples_c[1,1]
              aux_nt<-unequal_samples_t[1,1]
            }
            
            if (acc){
              txt_leg<-c()
              num_leg<-c()
              #
              # Set up the plot for accuracy
              if(ps)
                postscript(paste("Accuracy for d =", as.character(list_effect_sizes[d]), ";sdv =", as.character(list_standard_deviations[sdv]/real_control_mean*100), " percent; t =", as.character(list_taus[t]/real_control_mean*100), " percent; nc=", as.character(aux_nc), ";nt=", as.character(aux_nt), ".ps"),
                         width = 6, height = 6, paper = "special", family = "Helvetica")
              else
                png(paste("Accuracy for d =", as.character(list_effect_sizes[d]), ";sdv =", as.character(list_standard_deviations[sdv]/real_control_mean*100), " percent; t =", as.character(list_taus[t]/real_control_mean*100), " percent; nc=", as.character(aux_nc), "; nt=", as.character(aux_nt), ".png"))
              
              plot(0, 0, 
                   main=paste("Accuracy for d =", as.character(list_effect_sizes[d]), ";sdv=", as.character(list_standard_deviations[sdv]/real_control_mean*100), "%;t=", as.character(list_taus[t]/real_control_mean*100), "%;nc=", as.character(aux_nc), ";nt=", as.character(aux_nt)), 
                   type = "n", 
                   xlab = "Number of experiments", ylab = "Accuracy",
                   xlim = c(min(list_experiments), max(list_experiments)+4), ylim = c (0, 1))
              
              lines (c (min(list_experiments), max(list_experiments)), c (1-alpha, 1-alpha), lwd=3)
              text (max(list_experiments), (1-alpha), pos=4, paste("1-alpha=", as.character(1-alpha)))
              
              #
              # Draw the WMD's accuracy in the plot
              if(m1){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(wmd_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="1")
                output_line(f, "accuracy", "wmd", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"wmd")
                num_leg<-append(num_leg,"1")       
              }
              
              #
              # Draw the random WMD's accuracy in the plot
              if(m2){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(random_wmd_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="2")
                output_line(f, "accuracy", "random-wmd", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"random wmd")
                num_leg<-append(num_leg,"2") 
              }
              
              #
              # Draw the PRR's accuracy in the plot (without logs)
              if(m3){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(noln_prr_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="3")  
                output_line(f, "accuracy", "noln_prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"noln_prr")
                num_leg<-append(num_leg,"3") 
              }
              
              #
              # Draw the PRR's accuracy in the plot (with logs)
              if(m4){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(ln_prr_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="4")  
                output_line(f, "accuracy", "ln_prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"ln_prr")
                num_leg<-append(num_leg,"4") 
              }
              
              #
              # Draw the random PRR's accuracy in the plot (without logs)
              if(m5){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(noln_random_prr_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="5")
                output_line(f, "accuracy", "noln_random-prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"noln_random prr")
                num_leg<-append(num_leg,"5")
              }
              
              #
              # Draw the random PRR's accuracy in the plot (with logs)
              if(m6){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(ln_random_prr_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="6")
                output_line(f, "accuracy", "ln_random-prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"ln_random prr")
                num_leg<-append(num_leg,"6")
              }
              
              #
              # Draw the NPRR's accuracy in the plot (without logs)
              if(m7){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(noln_nprr_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="7")
                output_line(f, "accuracy", "noln_nprr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"noln_nprr")
                num_leg<-append(num_leg,"7")
              }
              
              #
              # Draw the NPRR's accuracy in the plot (with logs)
              if(m8){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(ln_nprr_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="8")
                output_line(f, "accuracy", "ln_nprr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"ln_nprr")
                num_leg<-append(num_leg,"8")
              }
              
              #
              # Draw the SVC's accuracy in the plot
              if(m9){
                for (i in 1:length(list_experiments)) tmp_accuracy[i] <- sum(svc_accuracy[i,])/repetitions
                lines(list_experiments, tmp_accuracy, type="b", pch="9")
                output_line(f, "accuracy", "svc", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_accuracy)
                txt_leg<-append(txt_leg,"svc")
                num_leg<-append(num_leg,"9")
              }
              
              #
              # Legend  
              legend ("bottomright", txt_leg, pch=num_leg) 
              dev.off()        
            }
            
            #
            # Set up the plot for rejection rate or empirical power      
            if (list_effect_sizes[d] == 0) {
              if(rej){
                #
                # Rejection rate
                if(ps)
                  postscript(paste("Rej.. rate for d =", as.character(list_effect_sizes[d]), "; sdv =", as.character(list_standard_deviations[sdv]/real_control_mean*100), " percent; t =", as.character(list_taus[t]/real_control_mean*100), " percent; nc=", as.character(aux_nc), ";nt=", as.character(aux_nt), ".ps"), 
                           width = 6, height = 6, paper = "special", family = "Helvetica")
                else
                  png(paste("Rej.. rate for d =", as.character(list_effect_sizes[d]), "; sdv =", as.character(list_standard_deviations[sdv]/real_control_mean*100), " percent; t =", as.character(list_taus[t]/real_control_mean*100), " percent; nc=", as.character(aux_nc), "; nt=", as.character(aux_nt), ".png"))
                
                plot(0, 0, 
                     main=paste("Rej. rate for d=", as.character(list_effect_sizes[d]), ";sdv=", as.character(list_standard_deviations[sdv]/real_control_mean*100), "%;t=", as.character(list_taus[t]/real_control_mean*100), "%;nc=", as.character(aux_nc), ";nt=", as.character(aux_nt)), 
                     type = "n", 
                     xlab = "Number of experiments", ylab = "Rejection rate",
                     xlim = c(min(list_experiments), max(list_experiments)+4), ylim = c (0, 1))
                lines (c (min(list_experiments), max(list_experiments)), c (alpha, alpha), lwd=3)
                text (max(list_experiments), (alpha), pos=4, paste("alpha=", as.character(alpha)))
                measure <- "rejection_rate"
                #  print("rej")
              }
            } else {
              if(pow){
                #
                # Empirical power
                if(ps)
                 postscript(paste("Emp. power for d =", as.character(list_effect_sizes[d]), "; sdv =", as.character(list_standard_deviations[sdv]/real_control_mean*100), " percent; t =", as.character(list_taus[t]/real_control_mean*100), " percent; nc=", as.character(aux_nc), ";nt=", as.character(aux_nt), ".ps"),
                           width = 6, height = 6, paper = "special", family = "Helvetica")
                else
                  png(paste("Emp. power for d =", as.character(list_effect_sizes[d]), "; sdv =", as.character(list_standard_deviations[sdv]/real_control_mean*100), " percent; t =", as.character(list_taus[t]/real_control_mean*100), " percent; nc=", as.character(aux_nc), "; nt=", as.character(aux_nt), ".png"))
                
                plot(0, 0, 
                     main=paste("Emp. power for d=", as.character(list_effect_sizes[d]), ";sdv=", as.character(list_standard_deviations[sdv]/real_control_mean*100), "%;t=", as.character(list_taus[t]/real_control_mean*100), "%;nc=", as.character(aux_nc), ";nt=", as.character(aux_nt)), 
                     type = "n", 
                     xlab = "Number of experiments", ylab = "Empirical power",
                     #xlim = c(min(list_experiments), max(list_experiments)+4), ylim = c (0.6, 1))
                     xlim = c(min(list_experiments), max(list_experiments)+4), ylim = c (0, 1))
                lines (c (min(list_experiments), max(list_experiments)), c (0.8, 0.8), lwd=3)
                text (max(list_experiments), (0.8), pos=4, paste("1-beta=", as.character(0.8))) 
                measure <- "empirical_power"
                #  print("pow")
              }
            }
            
            if((rej & list_effect_sizes[d] == 0) | (pow & list_effect_sizes[d] > 0) ){
              txt_leg<-c()
              num_leg<-c()
              #
              # Draw the WMD's accuracy/rej. rate in the plot
              if(m1){  
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(wmd_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="1") 
                output_line(f, measure, "wmd", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"wmd")
                num_leg<-append(num_leg,"1")
              }
              
              #
              # Draw the random WMD's accuracy/rej. rate in the plot
              if(m2){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(random_wmd_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="2")
                output_line(f, measure, "random-wmd", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"random wmd")
                num_leg<-append(num_leg,"2")
              }
              
              #
              # Draw the PRR's accuracy/rej. rate in the plot (without logs)
              if(m3){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(noln_prr_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="3")  
                output_line(f, measure, "noln_prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"noln_prr")
                num_leg<-append(num_leg,"3")
              }
              
              #
              # Draw the PRR's accuracy/rej. rate in the plot (with logs)
              if(m4){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(ln_prr_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="4")  
                output_line(f, measure, "ln_prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"ln_prr")
                num_leg<-append(num_leg,"4")
              }
              
              #
              # Draw the random PRR's accuracy/rej. rate in the plot (without logs)
              if(m5){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(noln_random_prr_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="5")
                output_line(f, measure, "noln_random-prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"noln_random prr")
                num_leg<-append(num_leg,"5")
              }
              
              #
              # Draw the random PRR's accuracy/rej. rate in the plot (with logs)
              if(m6){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(ln_random_prr_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="6")
                output_line(f, measure, "ln_random-prr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"ln_random prr")
                num_leg<-append(num_leg,"6")
              }
              
              #
              # Draw the NPRR's accuracy/rej. rate in the plot (without logs)
              if(m7){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(noln_nprr_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="7")
                output_line(f, measure, "noln_nprr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"noln_nprr")
                num_leg<-append(num_leg,"7")
              }
              
              #
              # Draw the NPRR's accuracy/rej. rate in the plot (with logs)
              if(m8){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(ln_nprr_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="8")
                output_line(f, measure, "ln_nprr", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"ln_nprr")
                num_leg<-append(num_leg,"8")
              }
              
              #
              # Draw the SVC's accuracy/rej. rate in the plot
              if(m9){
                for (i in 1:length(list_experiments)) tmp_emp_power[i] <- sum(svc_emp_power[i,])/repetitions
                lines(list_experiments, tmp_emp_power, type="b", pch="9")
                output_line(f, measure, "svc", list_effect_sizes[d], list_standard_deviations[sdv], list_taus[t], unequal_samples_c, unequal_samples_t, list_experiments, tmp_emp_power)
                txt_leg<-append(txt_leg,"svc")
                num_leg<-append(num_leg,"9")
              }
              
            }
            
            #
            # Legend
            if (list_effect_sizes[d] == 0) { #0=rej
              if (rej){
                legend ("topleft", txt_leg, pch=num_leg)
              }        
            } else {
              if (pow){
                legend ("bottomright", txt_leg, pch=num_leg)
              }
            }
            
            if((rej & list_effect_sizes[d] == 0) | (pow & list_effect_sizes[d] > 0) )
              dev.off()   
            print(Sys.time()-start_time)

      }#tau
    }#sdv     
  }#d
  close (f)
  return(TRUE)
}


#----  Starting variables initialization
#
# set output graphs directory
mainDir<-'~/Desktop/Dropbox/RStudio/Simulacion TOSEM (web)(refactorizada)'
subDir<-'plots'


#
# set output text file name
fname <- "output.txt"


#
# Basic simulation parameters
repetitions <- 1000
alpha       <- 0.05
cutoff      <- qnorm(1-(alpha/2))

#
# Specific simulation parameters
# plus some others
list_experiments            <- c (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

#if it contains only one value, sample size is the same in all experiments, otherwise, sample sizes of all experiments are randomized
list_sample_sizes           <- c(22) # reference values c(2, 4, 6, 8, 10, 15, 20, 30)

list_effect_sizes           <- c(0,0.5) # reference values c(0, 0.2, 0.5, 0.8, 1.2, 2)
real_control_mean           <- 10

#factor of debalancing
fixed_unbalanced_factor     <- 0.5

#when it is true all experiments are debalanced by the same fixed unbalanced factor
fixed_unbalanced            <- FALSE

#random factor to debalancing samples, when it is true, random factor takes values less or equal than fixed_unbalanced_factor
random_unbalanced           <- FALSE

#
# The standard deviation, we assume reference values of 10% - 50% - 100% - 200% - 400% of the control treatment
list_standard_deviations <- c(0.5*real_control_mean)  # reference values c(0.1*real_control_mean, 0.5*real_control_mean, real_control_mean, 2*real_control_mean, 4*real_control_mean)
#
# tau-squared's, to explore the effect of heterogeneity,
list_taus <- c(0.1*real_control_mean) #reference values c(0, 0.01*real_control_mean, 0.05*real_control_mean, 0.1*real_control_mean, 0.5*real_control_mean, real_control_mean, 2*real_control_mean, 4*real_control_mean)

#Output postscript or png graphics
ps<-FALSE #TRUE




#run the simulation
simulation <- montecarlo(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,repetitions,alpha,cutoff,list_experiments,
                         list_sample_sizes,random_unbalanced,fixed_unbalanced_factor,fixed_unbalanced,list_effect_sizes,real_control_mean,
                         list_standard_deviations,list_taus,fname,mainDir,subDir,ps)

rm(list=ls())

