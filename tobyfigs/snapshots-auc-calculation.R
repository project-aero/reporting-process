
library(tidyverse)
library(reshape2)

##Used to calculate AUC directly from snapshots data.


#load simulated cases
load("cases.RData")

#Calculate EWS for times series data using spaero
get_ews <- function(emerge){

  ## Wrapper for spaero::get_stats, sets function variables and returns only stats as a data frame.
  #spaero parameters are hardcoded in here
  ews <- function(x){
    return(
      as.data.frame(spaero::get_stats(x, center_trend ="local_constant",
                                      center_kernel = "uniform",
                                      center_bandwidth =35,
                                      stat_trend = "local_constant",
                                      stat_kernel = "uniform",
                                      stat_bandwidth=35,
                                      lag = 1)$stats)
    )
  }

  emerge %>%
    group_by(sim) %>%
    do(ews = ews(.$I)) -> a
  bind_rows(a$ews) -> b
  emerge_ews <- cbind(emerge,b)
}

#calculate taus
get_taus <- function(df){
  df <- melt(df, id.vars = c("time", "S", "I", "R", "N", "cases", "beta_t", "sim"))
  df %>%
    group_by(variable, sim) %>%
    summarise(tau = stats::cor(value,time, use="pairwise.complete.obs", method="kendall")) -> taus
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

  if (infectious_period == "-weekly"){
    i <- which(dplyr::near(process_des_mat$infectious_days, 7))
  } else {
    i <- which(dplyr::near(process_des_mat$infectious_days, 30))
  }
  emerge_data <- simulated_procs[[i]]$test
  null_data <- simulated_procs[[i]]$null

  calculate_auc_agg <- function(aggregation_period){
    edf <- emerge_data
    ndf <- null_data
    is_snap <- function(x, period){
      mod <- x %% period
      dplyr::near(mod, 0)
    }
    if(aggregation_period == "-monthly"){
      per <- 30
    } else {
      per <- 7
    }
    edf <- edf[is_snap(edf$time, per), ]
    ndf <- ndf[is_snap(ndf$time, per), ]

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

write_csv(auc_data,"snapshots-auc.csv")
