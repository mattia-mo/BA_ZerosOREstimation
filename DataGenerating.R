library(epitools)
library(tidyverse)
library(utils)

# reprod
set.seed(2409)

# function to create a fourfold table, based on the sample size, the true Odds Ratio value and 
# the proportion between exposed and non-exposed subjects
# arguments: 
# n:    sample size of the observed study (n_obs)
# OR:   true, known Odds Ratio (theta)
# prop: proportion between exposed and non-exposed

simuldata <-  function(n, OR, prop){
  
  times0 <- as.integer(n*prop)
  times1 <- n - times0
  
  x <- c(rep(0, times0), rep(1, times1))
  if(n%%2 == 1) x <- c(x,0)
  
  # drawing Y based on the binomial distribution
  y <- rbinom(n, 1, plogis(x*log(OR)))
  
  # 
  df <- data.frame(x = x, y = y)
  tab <- with(df, table(x,y))
  #creating a 2x2 table of X and Y
  if(all(y == 0)){
    tab <- cbind(tab, "1" = c(0,0))
  } else if(all(y == 1)){
    tab <- cbind("0" = c(0,0), tab)
  }
  new_order <- c("1","0")
  tab <- tab[new_order, new_order]
  return(tab)
}

# N as simulation sample size, here denoted as n_rep
n_rep <-  10000

# sample sizes of studies
n_obs <- c(50,100,500)
#n_names <- paste("n_", n_obs, sep = "")

# true OR values
or_sim <- c(1, 2, 10, 30)
#or_names <- paste("or_", or_sim, sep = "")

# proportion exposed vs non-exposed
prop_sim <- c(0.35, 0.5)
#prop_names <- paste("prop_", prop_sim, sep = "")

# estimation methods
methods <- c("midp", "fisher", "wald", "small")
names(methods) <- c("Mid-p", "Fisher", "Wald", "Small")
n_methods <- length(methods)

# error handling of OR methods in cases of zeros
# as oddsratio.midp is throwing an error if there is a zero included and not only returning NA
myTryCatch <- function(expr) {
  err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }))
  list(value=value, error=err)
}

# create process of estimations for each table 
onerep <- function(rep, n_obs, prop = 0.5, oddsratio) {
  
  # simulate the table
  table_of_rep <- simuldata(n_obs, oddsratio, prop)
  table_as_list <- list(table_of_rep)
  
  #check if zero in table
  zero <- any(table_of_rep == 0)
  
  # applying estimation methods on table
  
  or_wald <- oddsratio.wald(table_of_rep)$measure[2,1]
  lower_wald <- oddsratio.wald(table_of_rep)$measure[2,2]
  upper_wald <- oddsratio.wald(table_of_rep)$measure[2,3]
  
  or_fisher <- oddsratio.fisher(table_of_rep)$measure[2,1]
  lower_fisher <-  oddsratio.fisher(table_of_rep)$measure[2,2]
  upper_fisher <-  oddsratio.fisher(table_of_rep)$measure[2,3]
  
  or_small <- oddsratio.small(table_of_rep)$measure[2,1]
  lower_small <- oddsratio.small(table_of_rep)$measure[2,2]
  upper_small <- oddsratio.small(table_of_rep)$measure[2,3]
  
  err <- myTryCatch(oddsratio.midp(table_of_rep))
  if(is.null(err$error)) {
    or_midp <- oddsratio.midp(table_of_rep)$measure[2,1]
    lower_midp <- oddsratio.midp(table_of_rep)$measure[2,2]
    upper_midp <- oddsratio.midp(table_of_rep)$measure[2,3]
  } else {
    or_midp <- NA
    lower_midp <- NA
    upper_midp <- NA
  }
  
  # prepare output of the parameters, OR and CI
  out <- data.frame(
    rep = rep,
    or = oddsratio,
    n_obs = n_obs,
    prop = prop,
    table = NA,
    zero = zero,
    or_wald = or_wald, or_fisher = or_fisher, OR_small = or_small, or_midp = or_midp,
    lower_wald = lower_wald, lower_fisher = lower_fisher, lower_small = lower_small, lower_midp = lower_midp,
    upper_wald = upper_wald, upper_fisher = upper_fisher, upper_small = upper_small, upper_midp = upper_midp
  )
  out$table <- list(table_of_rep)
  return(out)
}

# create all possible scenarios based on the predefined parameters
scenarios <- expand_grid(or_sim, n_obs, prop_sim)
n_scenarios <- nrow(scenarios)

# function to conduct multiple repetitions of one scenario inclusive storing the states
multreps <- function(n_obs, or, prop, n_rep){
  
  # predefine the output
  estimates <- data.frame(matrix(ncol = 6+(n_methods*3), nrow = n_rep))
  x <- c("rep", "or", "n_obs", "prop", "table", "zero", 
         "or_wald", "or_fisher", "or_small", "or_midp", 
         "lower_wald", "lower_fisher", "lower_small", "lower_midp", 
         "upper_wald", "upper_fisher", "upper_small", "upper_midp")
  colnames(estimates) <- x
  states <- matrix(ncol = 627, nrow = n_rep)
  
  # Run all nsim reps
  for (r in 1:n_rep) {
    states[r, ] <- c(rep = r, .Random.seed)
    estimates[r, ] <- onerep(rep = r, n_obs = n_obs, oddsratio = or, prop = prop)
  }
  return(list(estimates, states))
}

# create multiple repetitions on all scenarios and store the results (states inclusive)
k <- estimates_all <- scenarios %>% 
  apply(.,1, function(scen){
    multreps(n_obs = scen[2], or = scen[1], prop = scen[3], n_rep = n_rep)
  })

estimates_all <- lapply(k, `[[`, 1) %>% 
  bind_rows() %>% 
  group_by(or, n_obs, prop) %>% 
  mutate(scenario = cur_group_id(), .before = "rep")
states_all <- lapply(k, `[[`, 2) %>%
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rename(rep = V1) %>% 
  mutate(scenario = rep(1:n_scenarios, each = n_rep), .before = "rep")
