
library(tidyverse)
library(reshape2)

##Used to calculate AUC directly from prevalence data. 


#Compile simulator
system("g++ mnrm_sir.cpp -o ./sir_sim -lgsl -lgslcblas")
sir_command_root <- "./sir_sim 10 1000000 20 3"

#Calculate EWS for times series data using spaero
get_ews <- function(df){
  
  ## Wrapper for spaero::get_stats, sets function variables and returns only stats as a data frame. 
  #spaero parameters are hardcoded in here
  ews <- function(x){
    return(
      as.data.frame(spaero::get_stats(x, center_trend ="local_constant",
                                      center_kernel = "uniform",
                                      center_bandwidth =100, 
                                      stat_trend = "local_constant", 
                                      stat_kernel = "uniform",
                                      stat_bandwidth=100, 
                                      lag = 1)$stats)
    )
  }
  
  emerge <- df
  names(emerge) <- c("Time", "Infected", "R0", "N", "T", "Run")
  
  emerge %>%
    group_by(Run) %>%
    do(ews = ews(.$Infected)) -> a
  bind_rows(a$ews) -> b
  emerge_ews <- cbind(emerge,b)
}

#calculate taus
get_taus <- function(df){
  df <- melt(df, id.vars = names(df)[1:6]) 
  df %>%
    group_by(variable, Run) %>%
    summarise(tau = stats::cor(value,Time, use="pairwise.complete.obs", method="kendall")) -> taus
  return(taus) 
}

#calculate AUC
get_auc <- function(df){
  calc_auc <- function(predictions, is_null){
    test <- !is.na(predictions)
    predictions <- predictions[test]
    is_null <- is_null[test]
    r <- rank(predictions)
    r1 <- sum(r[!is_null])
    n1 <- sum(!is_null)
    n2 <- sum(is_null)
    (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
  }
  df %>%
    group_by(variable) %>%
    summarise(AUC = calc_auc(tau, isnull)) -> out
  
  new_levels <-  c("Mean", "Variance", "Var. 1st Diff", 
                   "Index of Dis.", "Autocovar.", "Autocorr.",
                   "Decay time", "Coeff. Var", "Skewness",
                   "Kurtosis")
  out$variable <- factor(out$variable,levels(out$variable)[c(6,1,2,7,3,4,5,8,9,10)])
  #names(new_levels) <- levels(out$variable)
  levels(out$variable) <- new_levels
  
  
  
  return(out)
}




#Simulates 1000 replicates, N=1e6, 20yrs of data, rng seed = 1

calculate_auc <- function(infectious_period){
  
  sir_command <- paste(sir_command_root,c("-e","-n"),infectious_period)
  emerge_data <- read.table(text = system(sir_command[1], intern = TRUE))
  null_data <- read.table(text = system(sir_command[2], intern = TRUE))
  
  calculate_auc_agg <- function(aggregation_period){
    edf <- emerge_data
    ndf <- null_data
    if(aggregation_period == "-monthly"){
      edf <-  edf[edf$V1 %% 4 == 0, ]
      ndf <-  ndf[ndf$V1 %% 4 == 0, ]
    }
    
    emerge_ews <- get_ews(edf)
    not_ews <- get_ews(ndf)
    taus_test <- get_taus(emerge_ews)
    taus_null <- get_taus(not_ews)
    taus_test$isnull <- FALSE
    taus_null$isnull <- TRUE
    taus <- rbind(taus_test, taus_null)
    
    df <- get_auc(taus)
    df$`Infectious period` <- as.factor(infectious_period)
    df$`Aggregation period` <- as.factor(aggregation_period)
    return(df)
  }

  df <- do.call(rbind,lapply(c("-weekly","-monthly"),calculate_auc_agg))
  return(df)
}

auc_data <- do.call(rbind,lapply(c("-weekly", "-monthly"),calculate_auc))

write_csv(auc_data,"prevalence_auc.csv")
