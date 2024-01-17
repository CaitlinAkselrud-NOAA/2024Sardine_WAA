# Purpose: method for filling missing weight-at-age data in stock assessment, 
#          incl cohort/age/year effects
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 1/16/23
# 
# Modified: Caitlin Allen Akselrud 12.12.23
# contact: caitlin.allen_akselrud@noaa.gov
# Notes: set up as main code to call functions in other files
# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(TMB)
library(cowplot)
library(janitor)
library(magrittr)
library(gt)

# Compile and load in model -----------------------------------------------

wd <- here("src")
compile(here("src","GMRF_WAA.cpp"))
dyn.load(dynlib(here("src","GMRF_WAA")))

# Data --------------------------------------------------------------------
# FISHERY DATA
# Load in WAA matrix
waa_df_raw <- read_csv(here("data", "sardine_fishery_waa_kg.csv")) %>%
  janitor::clean_names()

# Load in std for WAA matrix
waa_std_df_raw <- read_csv(here("data", "sardine_fishery_waa_kg_sd.csv")) %>%
  janitor::clean_names()

# clean data

# fleets
# 1 = MexCal S1
# 2 = MexCal S2
# 3 = PNW
# 4 = AT survey
waa_df <- waa_df_raw %>% 
  # mutate(year = case_when(seas == 1 ~ number_yr + 0.0,
  #                         seas == 2 ~ number_yr + 0.5)) %>% 
  rename(year = model_year) %>% 
  mutate(fleet = case_when(unq == "MexCal 1" ~ 1,
                           unq == "MexCal 2" ~ 2,
                           unq == "PNW 1" ~ 3)) %>%
  dplyr::select(-x9, -x10) %>% 
  dplyr::select_if(is.numeric) #remove columns with notes


waa_std_df <- waa_std_df_raw %>% 
  # mutate(year = case_when(seas == 1 ~ number_yr + 0.0,
  #                         seas == 2 ~ number_yr + 0.5)) %>% 
  rename(year = model_year) %>%
  mutate(fleet = case_when(unq == "MexCal 1" ~ 1,
                           unq == "MexCal 2" ~ 2,
                           unq == "PNW 1" ~ 3)) %>%
  dplyr::select(-x9, -x10) %>% 
  dplyr::select_if(is.numeric) #remove cols  with notes


# data (survey) -----------------------------------------------------------

# SURVEY DATA
# Load in WAA matrix
waa_df_raw <- read_csv(here("data", "2024_survey_waa.csv")) %>% 
  janitor::clean_names()

# Load in std for WAA matrix
waa_std_df_raw <- read_csv(here("data", "2024_survey_waa_sd.csv")) %>% 
  janitor::clean_names() 
  
waa_df <- waa_df_raw %>% 
  rename(year = yr) %>% 
  dplyr::select(-seas, -sex, -bio_pattern, -birth_seas) %>% 
  mutate(fleet = 1)

waa_std_df <- waa_std_df_raw %>% 
  rename(year = yr) %>% 
  dplyr::select(-seas, -sex, -bio_pattern, -birth_seas) %>% 
  mutate(fleet = 1)

# Fill NA values ----------------------------------------------------------
# CA: updated 20.12.2023 to not fill NAs
# rpl_sd <- rep(1.111, times = length(names(waa_df)))
# names(rpl_sd) <- names(waa_df)
# 
# waa_std_df %<>%
#   replace_na(as.list(rpl_sd))
# # waa filled in by fleet inside get_fleet_dat fxn
# overall_mean_WAA <- waa_df %>%
#   colMeans(na.rm = T) %>%
#   as.list()

# Fixed inputs ------------------------------------------------------------

projection_time <- 2 #fishery: 9
proj_yr_start <- 2021 #fishery: 2015

model_name <- "2024_sardine_"
model_year <- 2024

# * fleet setup ---------------------------------------------------------------
n_fleets <- 1 # fishery: 3

# * other model settings --------------------------------------------------

newton_steps = 3

# Set up model functions ----------------------------------------
get_fleet_dat <- function(waa_dat, waa_std_dat, fleet_num)
{
  waa_f <- waa_dat %>% dplyr::filter(fleet == fleet_num)
  # mean_WAA <- waa_f %>%
  #   colMeans(na.rm = T) %>%
  #   as.list()
  # waa_f %<>%
  #   replace_na(mean_WAA)
  # waa_f %<>%
  #   replace_na(overall_mean_WAA) #if there is no fleet-specific mean, fill with overall mean
  
  waa_std_f = waa_std_dat %>% dplyr::filter(fleet == fleet_num)
  return(list(dat = waa_f, sd = waa_std_f))
}

borrow_data <- function(waa_dat, waa_std_dat, fleet_num_add, fleet_num_borrow, replace_by = "m", sd_repl_value = 1.111)
{
  waa_f <- waa_dat %>% dplyr::filter(fleet == fleet_num_add)
  if(replace_by == "y") #by year
  {
    waa_borrow <- waa_dat %>% 
      dplyr::filter(fleet %in% fleet_num_borrow) %>% 
      group_by(year) %>%
      summarize_all(.funs = mean, na.rm=T)
    all_na <- waa_f %>% keep(~all(is.na(.x))) %>% names
    replace_c <- waa_borrow %>% dplyr::select(year, all_of(all_na)) %>% 
      dplyr::filter(year %in% unique(waa_f$year))
    waa_f %<>% dplyr::select(-all_of(all_na)) %>% 
      full_join(replace_c)
  }else if(replace_by == "m") #general mean
  {
    waa_borrow <- waa_dat %>% 
      dplyr::filter(fleet %in% fleet_num_borrow) %>% 
      summarize_all(.funs = mean, na.rm=T) %>% 
      dplyr::select(-year)
    all_na <- waa_f %>% keep(~all(is.na(.x))) %>% names
    replace_c <- waa_borrow %>% dplyr::select(all_of(all_na))
    waa_f %<>% dplyr::select(-all_of(all_na)) %>%
      mutate(replace_c)
  }
  
  rpl_sd <- rep(sd_repl_value, times = length(names(waa_f)))
  names(rpl_sd) <- names(waa_df)
  
  rpl_sd_f <- waa_std_dat %>% 
    dplyr::filter(fleet == fleet_num) %>% 
    dplyr::select(all_of(all_na)) %>% 
    replace_na(as.list(rpl_sd))
  
  waa_std_f = waa_std_dat %>% 
    dplyr::filter(fleet == fleet_num) %>% 
    dplyr::select(-all_of(all_na)) %>%
    mutate(rpl_sd_f)
  
  return(list(dat = waa_f, sd = waa_std_f))
}

# TMB_setup <- function(proj_yrs = 2, waa_df, waa_std_df)
# {
#   # Number of projection years
#   n_proj_years <- proj_yrs
#   
#   # Years
#   years <- waa_df$year
#   
#   # Ages (goes from age 3 - 15+)
#   ages <- str_extract_all(names(waa_df), '[0-9]+') %>% 
#     unlist() %>% as.numeric()
#   
#   # Read in data weight at age matrix
#   X_at <- waa_df %>% 
#     select(starts_with("x")) %>% #keep only age cols
#     t() 
#   
#   # Create projection columns (append to X_at matrix)
#   proj_cols <- matrix(NA, nrow = length(ages), ncol = n_proj_years) 
#   
#   # Append NA for projection year
#   X_at <- cbind(X_at, proj_cols) 
#   
#   # Read in standard deviations for weight at age matrix
#   Xse_at <- waa_std_df %>% 
#     select(starts_with("x")) %>% #keep only age cols
#     t() 
#   
#   # Convert to CV
#   Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )
#   
#   # Now convert back to sd in lognormal space
#   Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))
#   
#   # Create an index for ages and years to feed into TMB, which helps construct the precision matrix
#   ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)), #CIA: does age need to be age - 1?? CHECK
#                                     "year" = seq_len(length(years) + n_proj_years) ))
#   
#   return(list(years = years,
#               ages = ages,
#               X_at = X_at,
#               Xsd_at = Xsd_at,
#               ay_Index = ay_Index,
#               n_proj_years = n_proj_years ))
# }
TMB_setup <- function(proj_yrs = 2, waa_df, waa_std_df)
{
  # Number of projection years
  n_proj_years <- proj_yrs
  
  # Years
  years <- waa_df$year
  
  # Ages (goes from age 3 - 15+)
  ages <- str_extract_all(names(waa_df), '[0-9]+') %>%
    unlist() %>% as.numeric()
  
  # Create projection columns (append to X_at matrix)
  proj_cols <- matrix(NA, nrow = length(ages), ncol = n_proj_years)
  
  # Read in standard deviations for weight at age matrix
  Xse_at <- waa_std_df %>%
    select(starts_with("x")) %>% #keep only age cols
    t()
  
  # Convert to CV
  Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )
  
  # Now convert back to sd in lognormal space
  Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))
  Xsd_at <- cbind(Xsd_at, proj_cols) # append NAs for projection
  
  # Read in data weight at age matrix
  X_at <- waa_df %>%
    select(starts_with("x")) %>% #keep only age cols
    t()
  
  # Append NA for projection year
  X_at <- cbind(X_at, proj_cols)
  X_at[is.na(Xsd_at)] <- NA
  
  # Create an index for ages and years to feed into TMB, which helps construct the precision matrix
  ay_Index <- as.matrix(expand.grid("age" = seq_len(length(ages)), #CIA: does age need to be age - 1?? CHECK
                                    "year" = seq_len(length(years) + n_proj_years) ))
  
  return(list(years = years,
              ages = ages,
              X_at = X_at,
              Xsd_at = Xsd_at,
              ay_Index = ay_Index,
              n_proj_years = n_proj_years ))
}

get_params_Y_at <- function(dat_in)
{
  yat <- array(0,dim=dim(dat_in$X_at))
  return(yat)
}

fact_design <- function(y = 0:1, c = 0:1, a = 0:1)
{
  # Generate full factorial design
  map_factorial <- tidyr::crossing(rho_y = y, rho_c = c, rho_a = a) %>% 
    data.frame()
  return(map_factorial)
}

run_model <- function(map_factorial, n.newton = 3, data, parameters)
{
  # Define number of extra newton steps we want to take
  
  # Empty list to store model objects
  models <- list()
  
  for(n_fact in 1:nrow(map_factorial)) {
    
    # Create empty map list object
    map <- list()
    
    # Extract combinations of parameters estimated here
    map_fact <- map_factorial[n_fact,]
    
    # Create our mapping list to feed into MakeADFun
    for(m in 1:length(names(map_fact))) {
      
      if(map_fact[1,names(map_fact)[m]] == 0) { # if factorial = 0, turn estimation off
        map[[m]] <- factor(NA)
        names(map)[m] <- names(map_fact)[m] # Name list object
      } else{ # factorial == 1
        map[[m]] <- factor(1)
        names(map)[m] <- names(map_fact)[m] # Name list object
      } # ifelse statement for constructing map list object
    } # end m loop
    
    map <- rlist::list.append(map, "ln_Linf" = factor(NA), 
                              "ln_beta" = factor(NA))
    
    # Now, make AD model function
    waa_model <- MakeADFun(data = data, parameters = parameters,
                           map = map, random = "ln_Y_at",
                           DLL = "GMRF_WAA", silent = FALSE)
    
    # Now, optimize the function
    waa_optim <- stats::nlminb(waa_model$par, waa_model$fn, waa_model$gr,  
                               control = list(iter.max = 1e5, eval.max = 1e5))
    
    # Take some additional newton steps to make sure we reach a minimum
    tryCatch(expr = for(i in 1:n.newton) {
      g = as.numeric(waa_model$gr(waa_optim$par))
      h = optimHess(waa_optim$par, fn = waa_model$fn, gr = waa_model$gr)
      waa_optim$par = waa_optim$par - solve(h,g)
      waa_optim$objective = waa_model$fn(waa_optim$par)
    }, error = function(e){e})
    
    # Save optimized model results
    waa_model$optim <- waa_optim
    
    # Get report
    waa_model$rep <- waa_model$report(waa_model$env$last.par.best)
    
    # Get sd report
    waa_model$sd_rep <- sdreport(waa_model)
    
    models[[n_fact]] <- waa_model
    
    print(n_fact)
    
  } # loop through to run multiple models
  return(models)
}

margAIC <- function(optim_model) {
  
  # Get number of parameters 
  k <- length(optim_model[["par"]])
  
  # Extract objective function
  nLL <- optim_model[["objective"]]
  
  # Calculate AIC
  margAIC <- 2*k + 2*nLL
  
  return(margAIC)
}

get_diagnostics <- function(fleet, fleet_name, model_name, models)
{
  # Generate full factorial design for naming models
  map_factorial <- tidyr::crossing(rho_y = 0:1, rho_c = 0:1, rho_a = 0:1) %>% 
    data.frame()
  
  # Model Checking ----------------------------------------------------------
  
  AIC_models <- vector() # store aic
  nlminb_conv <- vector() # store nlminb convergence diagnostic
  grad_conv <- vector() # store maximum gradient
  grad_name_conv <- vector() # store parameter with maximum gradient
  hess_conv <- vector() # whether or not Hessian is positive definite 
  par_values <- data.frame(rho_y = NA, rho_a = NA, 
                           rho_c = NA, log_sigma2 = NA) # parameter values
  par_sd_values <- data.frame(rho_y_sd = NA, rho_a_sd = NA, 
                              rho_c_sd = NA,  log_sigma2_sd = NA) # se values
  
  for(i in 1:length(models)) {
    
    # Convergence diagnostics
    AIC_models[i] <- margAIC(models[[i]]$optim) # get aic
    nlminb_conv[i] <- models[[i]]$optim$convergence # get nlminb convergence
    grad_conv[i] <- max(abs(models[[i]]$sd_rep$gradient.fixed)) # max gradient
    # Get max gradient parameter name
    grad_name_conv[i] <- names(models[[i]]$sd_rep$par.fixed)[which.max(abs(models[[i]]$sd_rep$gradient.fixed))]
    hess_conv[i] <- models[[i]]$sd_rep$pdHess # whether or not this is pd HESS
    
    # Grab parameter values here and store
    par_values[i,] <- models[[i]]$sd_rep$par.fixed[match(colnames(par_values),
                                                         names(models[[i]]$sd_rep$par.fixed))]
    
    # Get parameter std values here and store
    par_sd_values[i,] <- sqrt(
      diag(models[[i]]$sd_rep$cov.fixed)[match(colnames(par_values),
                                               names(diag(models[[i]]$sd_rep$cov.fixed)))]
    )
    
    if(i == length(models)) {
      model_diag <- data.frame(AIC = AIC_models, nlminb_conv = nlminb_conv,
                               max_grad_name = grad_name_conv,
                               max_grad = grad_conv, pd_Hess = hess_conv)
      
      # Bind parameter values and sd together
      model_diag <- cbind(model_diag, par_values, par_sd_values)
    } # when done w/ evaluating all models
    
  } # end i loop
  
  
  # Munge into plot format --------------------------------------------------
  
  # Create model names to differentiate models
  model_names <- map_factorial %>% 
    mutate(rho_y_lab = case_when(rho_y == 1 ~ "y"),
           rho_a_lab = case_when(rho_a == 1 ~ "a"),
           rho_c_lab = case_when(rho_c == 1 ~ "c")) %>% 
    dplyr::select(rho_y_lab, rho_a_lab, rho_c_lab) %>% 
    dplyr::rowwise() %>% 
    tidyr::unite('model', na.rm = TRUE)
  
  # Input model names above into model_diag df
  model_diag$model <- model_names$model
  # No correlation parameters estimated
  model_diag$model[model_diag$model == ""] <- "None"
  
  model_diag %<>% dplyr::filter(!is.na(AIC))
  
  # Pivot this dataframe longer 
  # Parmeter MLEs here only (doing ses and binding in 2 steps)
  model_pars_long <- model_diag %>% 
    # mutate(
    #   rho_a_trans = 2 / (1 + exp(-2 * rho_a)) - 1,
    #   rho_y_trans = 2 / (1 + exp(-2 * rho_y)) - 1,
    #   rho_c_trans = 2 / (1 + exp(-2 * rho_c)) - 1
    # ) %>% # transform rhos
    dplyr::select(-rho_a_sd, -rho_y_sd, -rho_c_sd, -log_sigma2_sd) %>% 
    pivot_longer(cols = c(rho_a, rho_c, rho_y, log_sigma2), 
                 names_to = "parameters",  values_to = "mle_val")
  
  # Get SE values now
  model_se_long <- model_diag %>% 
    dplyr::select(rho_a_sd, rho_c_sd, rho_y_sd, log_sigma2_sd) %>% 
    pivot_longer(cols = everything(),names_to = "sd",  values_to = "sd_val")
  
  # Now bind, these two together
  model_diag_long <- cbind(model_pars_long, model_se_long)
  
  # Create lwr and upr confidence intervals here
  model_diag_long <- model_diag_long %>% 
    mutate(lwr_95 = ifelse(parameters == "log_sigma2", exp(mle_val - (1.96 * sd_val)),
                           mle_val - (1.96 * sd_val)),
           upr_95 = ifelse(parameters == "log_sigma2", exp(mle_val + (1.96 * sd_val)),
                           mle_val + (1.96 * sd_val)),
           mle_val = ifelse(parameters == "log_sigma2", exp(mle_val), mle_val),
           dAIC = min(AIC) - AIC)
  
  # Calcualte wAIC
  wAIC <- model_diag_long %>% 
    group_by(model) %>% 
    summarize(dAIC = abs(mean(dAIC))) %>% 
    ungroup() %>% 
    mutate(wAIC = exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)))
  
  # Left join wAIC
  model_diag_long <- model_diag_long %>% 
    left_join(wAIC %>% dplyr::select(wAIC, model), by = c("model"))
  
  # Output to csv
  write.csv(model_diag_long, here::here("output", paste0(model_name, fleet_name, "_model_diag_vals.csv")))
  return(model_diag_long)
}

# conditional model -------------------------------------------------------

# * choose factorial design -----------------------------------------------
model_fact <- fact_design(y = 0:1, 
                          c = 0:1, 
                          a = 0:1)


for(i in 1:n_fleets)
{
  fleet_name <- paste0("fleet",i,"_")
  fleet_name <- paste0("survey_")
  fleet_num <- i
  
  waa_fleet <- get_fleet_dat(waa_dat = waa_df,
                             waa_std_dat = waa_std_df,
                             fleet_num = fleet_num)
  if(i == 3)
  {
    borrow <- c(1) # only mexcal s1, not both
    waa_fleet <- borrow_data(waa_dat = waa_df,
                             waa_std_dat = waa_std_df,
                             fleet_num_add = fleet_num,
                             fleet_num_borrow = borrow,
                             replace_by = "m", #m = overall mean; y = by year
                             sd_repl_value = 1.111)
    waa_fleet$dat %<>%
      relocate(x0, .after = year)
    waa_fleet$sd %<>%
      relocate(x0, .after = year)

  }
  
  dat_setup <- TMB_setup(proj_yrs = projection_time,      # this is proj time steps (2 in this code = 2 seasons, aka 1 year)
                         waa_df = waa_fleet$dat, 
                         waa_std_df = waa_fleet$sd)
  
  # dim check
  check <- dim(dat_setup$X_at) == dim(dat_setup$Xsd_at)
  if(check[1] == "FALSE" || check[2] == 'FALSE')
  {
    print("WARNING: matrix dimensions don't match, reconfig data setup")
  }
  
  # set conditional variance
  dat_setup$Var_Param <- 0 # Var_Param == 0 Conditional, == 1 Marginal
  
  # Now, input these components into a data list
  data_in <- dat_setup
  
  # Input parameters into a list
  parameters_in <- list( rho_y = 0, 
                         rho_a = 0,
                         rho_c = 0,
                         log_sigma2 = log(0.03), #
                         ln_L0 = log(9), #first length bin sardine
                         ln_Linf = log(28),  # last length bin for sardine
                         ln_k = log(0.15),
                         ln_alpha = log(3.5e-7), # Start alpha at a reasonable space 
                         # Starting value for alpha derived from a run where none of the rhos were estimated.
                         ln_beta = log(3), # Fix at isometric
                         ln_Y_at =  get_params_Y_at(data_in))
  
  
  # * run conditional model -------------------------------------------------
  
  models_cond <- run_model(map_factorial = model_fact, 
                           n.newton = newton_steps, 
                           data = data_in, 
                           parameters = parameters_in)
  
  save(models_cond, file = here("output", paste0(model_name, fleet_name, "cond_var_waa_models.RData")))
}

# marginal model ----------------------------------------------------------

# * run marginal model ----------------------------------------------------

# CA: doesn't converge for fleet 2 or 3

# for(i in 1:n_fleets)
# {
#   fleet_name <- paste0("fleet",i,"_")
#   fleet_num <- i
# 
#   waa_fleet <- get_fleet_dat(waa_dat = waa_df,
#                              waa_std_dat = waa_std_df,
#                              fleet_num = fleet_num)
# 
#   dat_setup <- TMB_setup(proj_yrs = projection_time,      # this is proj time steps (2 in this code = 2 seasons, aka 1 year)
#                          waa_df = waa_fleet$dat,
#                          waa_std_df = waa_fleet$sd)
# 
#   dat_setup$Var_Param <- 1
# 
#   data_in <- dat_setup
# 
#   parameters_in <- list( rho_y = 0,
#                          rho_a = 0,
#                          rho_c = 0,
#                          log_sigma2 = log(0.1), #
#                          ln_L0 = log(9), #first length bin sardine
#                          ln_Linf = log(28),  # last length bin for sardine
#                          ln_k = log(0.15),
#                          ln_alpha = log(3.5e-7), # Start alpha at a reasonable space
#                          # Starting value for alpha derived from a run where none of the rhos were estimated.
#                          ln_beta = log(3), # Fix at isometric
#                          ln_Y_at =  get_params_Y_at(data_in))
# 
#   models_marg <- run_model(map_factorial = model_fact,
#                            n.newton = newton_steps,
#                            data = data_in,
#                            parameters = parameters_in)
# 
#   save(models_marg, file = here("output", paste0(model_name, fleet_name, "marg_var_waa_models.RData")))
# }


# model selection ---------------------------------------------------------
get_model <- function(fleet_diag, model_fac)
{
  best_mod <- read_csv(here::here("output", fleet_diag)) %>% 
    dplyr::filter(AIC == min(AIC)) %>% 
    distinct(model) %>% 
    mutate(rho_y = if_else(str_detect(model, "y"), 1, 0),
           rho_c = if_else(str_detect(model, "c"), 1, 0),
           rho_a = if_else(str_detect(model, "a"), 1, 0))
  model_fac %<>% mutate(model_num = row_number())
  mod_index <- inner_join(model_fac, best_mod)
  return(mod_index$model_num)
}  

get_cond_waa <- function(fleet_index, cond_model, fleet_name, fleet_num, min_age, max_age, proj_t, proj_yr)
{
  proj_names = paste0("Proj_", seq(from = proj_yr, to = proj_yr+proj_t-1, by = 1))
  
  fleet_waa <- cond_model[[fleet_index]]$sd_rep$par.random %>% 
    matrix(nrow= (max_age - min_age +1)) %>% 
    exp()
  rownames(fleet_waa) = paste0("age_", seq(from = min_age, to = max_age, by = 1))
  colnames(fleet_waa) = c(cond_model[[fleet_index]]$env$data$years, proj_names)
  
  fleet_waa_ss <- fleet_waa %>%
    t() %>%
    as_tibble() %>% 
    mutate(Yr = c(models_cond[[fleet_index]]$env$data$years, proj_names),
           Seas = 1,
           Sex = 1,
           Bio_Pattern = 1,
           BirthSeas = 1,
           Fleet = fleet_num,
           comment = paste0("#wt_flt_", fleet_num)) %>%
    relocate(Yr, Seas, Sex, Bio_Pattern, BirthSeas, Fleet)
  
  heatmap(fleet_waa, 
          Rowv = NA, Colv = NA,
          main = paste("Fleet", fleet_num, "conditional weight-at-age"))
  
  write_csv(as.data.frame(fleet_waa), 
            file = here::here("output", paste0(model_name, fleet_name, "_waa_output.csv")))
  write_csv(fleet_waa_ss, 
            file = here::here("output", paste0(model_name, fleet_name, "_waa_ss.csv")))
  return(fleet_waa)
}

get_cond_waa_sd <- function(fleet_index, cond_model, fleet_name, fleet_num, min_age, max_age, proj_t, proj_yr)
{
  proj_names = paste0("Proj_", seq(from = proj_yr, to = proj_yr+proj_t-1, by = 1))
  fleet_waa_sd <- cond_model[[fleet_index]]$sd_rep$diag.cov.random %>% 
    matrix(nrow = (max_age - min_age +1))
  rownames(fleet_waa_sd) = paste0("age_", seq(from = min_age, to = max_age, by = 1))
  colnames(fleet_waa_sd) = c(cond_model[[fleet_index]]$env$data$years, proj_names)
  
  heatmap(fleet_waa_sd, 
          Rowv = NA, Colv = NA,
          main = paste("Fleet", fleet_num, "standard dev"))
  write_csv(as.data.frame(fleet_waa_sd), 
            file = here::here("output", paste0(model_name, fleet_name, "_waa_sd_output.csv")))
  return(fleet_waa_sd)
}
# * * fleet 1 -------------------------------------------------------------

load(here("output", paste0(model_year, "_sardine_fleet1_cond_var_waa_models.RData")))
f1_diag <- get_diagnostics(fleet = 1, 
                           fleet_name = "fleet1",
                           model_name = model_name,
                           models = models_cond) 

fleet1_index <- get_model(fleet_diag = paste0(model_year, "_sardine_fleet1_model_diag_vals.csv"),
                          model_fac = model_fact)
fleet1_waa <- get_cond_waa(fleet_index = fleet1_index,
                           cond_model = models_cond,
                           fleet_name = "fleet1",
                           fleet_num = 1,
                           min_age = 0,
                           max_age = 8,
                           proj_t = projection_time,
                           proj_yr = proj_yr_start)

# heatmap(fleet1_waa, Rowv = NA, Colv = NA)

# * std dev
fleet1_waa_sd <- get_cond_waa_sd(fleet_index = fleet1_index,
                                 cond_model = models_cond,
                                 fleet_name = "fleet1",
                                 fleet_num = 1,
                                 min_age = 0,
                                 max_age = 8,
                                 proj_t = projection_time,
                                 proj_yr = proj_yr_start)


# CA: you are here: plots; input vs ouput comparisons for validation
# CA: create working plot fxn for all fleets 


# cols = ages, rows = years
# exp() log space

# * * fleet 2 -------------------------------------------------------------

load(here("output", paste0(model_year, "_sardine_fleet2_cond_var_waa_models.RData")))
f2_diag <- get_diagnostics(fleet = 2, 
                           fleet_name = "fleet2",
                           model_name = model_name,
                           models = models_cond)

fleet2_index <- get_model(fleet_diag = paste0(model_year, "_sardine_fleet2_model_diag_vals.csv"),
                          model_fac = model_fact)
fleet2_waa <- get_cond_waa(fleet_index = fleet2_index,
                           cond_model = models_cond,
                           fleet_name = "fleet2",
                           fleet_num = 2,
                           min_age = 0,
                           max_age = 8,
                           proj_t = projection_time,
                           proj_yr = proj_yr_start)

# * sd
fleet2_waa_sd <- get_cond_waa_sd(fleet_index = fleet2_index,
                                 cond_model = models_cond,
                                 fleet_name = "fleet2",
                                 fleet_num = 2,
                                 min_age = 0,
                                 max_age = 8,
                                 proj_t = projection_time,
                                 proj_yr = proj_yr_start)

# fleet2_waa %>% image()

# * * fleet 3 -------------------------------------------------------------

load(here("output", paste0(model_year, "_sardine_fleet3_cond_var_waa_models.RData")))
f3_diag <- get_diagnostics(fleet = 3, 
                           fleet_name = "fleet3",
                           model_name = model_name,
                           models = models_cond)
fleet3_index <- get_model(fleet_diag = paste0(model_year, "_sardine_fleet3_model_diag_vals.csv"),
                          model_fac = model_fact)
fleet3_waa <- get_cond_waa(fleet_index = fleet3_index,
                           cond_model = models_cond,
                           fleet_name = "fleet3",
                           fleet_num = 3,
                           min_age = 0,
                           max_age = 8,
                           proj_t = projection_time,
                           proj_yr = proj_yr_start)
# * sd
fleet3_waa_sd <- get_cond_waa_sd(fleet_index = fleet3_index,
                                 cond_model = models_cond,
                                 fleet_name = "fleet3",
                                 fleet_num = 3,
                                 min_age = 0,
                                 max_age = 8,
                                 proj_t = projection_time,
                                 proj_yr = proj_yr_start)

# fleet3_waa %>% image()

# fleet3_waa_sd <- models_cond[[fleet3_index]]$sd_rep$diag.cov.random %>% 
#   matrix(nrow = 9)

# * * survey -------------------------------------------------------------

load(here("output", paste0(model_year, "_sardine_survey_cond_var_waa_models.RData")))
surv_diag <- get_diagnostics(fleet = -1, 
                           fleet_name = "survey",
                           model_name = model_name,
                           models = models_cond)
surv_index <- get_model(fleet_diag = paste0(model_year, "_sardine_survey_model_diag_vals.csv"),
                          model_fac = model_fact)
surv_waa <- get_cond_waa(fleet_index = surv_index,
                           cond_model = models_cond,
                           fleet_name = "survey",
                           fleet_num = -1,
                           min_age = 0,
                           max_age = 10,
                           proj_t = projection_time,
                           proj_yr = proj_yr_start)
# * sd
surv_waa_sd <- get_cond_waa_sd(fleet_index = surv_index,
                                 cond_model = models_cond,
                                 fleet_name = "survey",
                                 fleet_num = -1,
                                 min_age = 0,
                                 max_age = 10,
                                 proj_t = projection_time,
                                 proj_yr = proj_yr_start)
# model checks ------------------------------------------------------------

# CA: compare output to original input
# fleet1_waa --> model waa
# waa_df --> orig waa (all fleets)

compare_waa <- function(orig_dat, model_dat, n_proj, fleet_name, model_name)
{
  orig_dat %<>% dplyr::select(-fleet, -year) %>% t()
  rownames(orig_dat) = rownames(model_dat)
  colnames(orig_dat) = colnames(model_dat)[1:(length(colnames(model_dat))-n_proj)]
  
  comp <- orig_dat - model_dat[,1:(length(colnames(model_dat))-n_proj)]
  # heatmap(comp, Rowv = NA, Colv = NA)
  dat_comp <- comp %>% 
    as_tibble %>% 
    rownames_to_column("Age") %>% 
    mutate(Age = as.numeric(Age) -1) %>% 
    pivot_longer(-Age, names_to = "Year", values_to = "waa")
  
  p_comp <- ggplot(dat_comp, aes(Age, Year)) +
    geom_tile(aes(fill = waa)) +
    geom_text(aes( #label=formatC(waa, format = "e"))) +
      label = round(waa, 3))) +
    scale_fill_gradient(low = "white", high = "cadetblue") +
    theme_bw() +
    labs(title = paste("Difference in waa:", fleet_name), 
         subtitle = "Original data vs selected model")
  ggsave(filename = here::here("output", paste0(model_name = model_name, fleet_name, "diff_obsv_pred_waa.png")), 
         width = 10, height = 10)
  return(comp)
}

f1_compare <- compare_waa(orig_dat = (waa_df %>% dplyr::filter(fleet == 1)),
                          model_dat = fleet1_waa,
                          n_proj = projection_time,
                          fleet_name = "fleet1", 
                          model_name = model_name)

f2_compare <- compare_waa(orig_dat = (waa_df %>% dplyr::filter(fleet == 2)),
                          model_dat = fleet2_waa,
                          n_proj = projection_time,
                          fleet_name = "fleet2", 
                          model_name = model_name)

f3_compare <- compare_waa(orig_dat = (waa_df %>% dplyr::filter(fleet == 3)),
                          model_dat = fleet3_waa,
                          n_proj = projection_time,
                          fleet_name = "fleet3", 
                          model_name = model_name)

surv_compare <- compare_waa(orig_dat = (waa_df),
                          model_dat = surv_waa,
                          n_proj = projection_time,
                          fleet_name = "survey", 
                          model_name = model_name)
# additional output -------------------------------------------------------

# tables:
# 1 model diag values as table (arrange by min AIC) for each fleet
f1_diag %<>% mutate(fleet = "Fleet 1")
f2_diag %<>% mutate(fleet = "Fleet 2")
f3_diag %<>% mutate(fleet = "Fleet 3")
mod_results <- bind_rows(f1_diag, f2_diag, f3_diag) %>%
  dplyr::select(-nlminb_conv, -max_grad_name, -max_grad) %>% 
  group_by(fleet) %>%
  arrange(min(wAIC)) %>% 
  relocate(AIC, .before = dAIC) %>% 
  relocate(pd_Hess, .after = wAIC) %>% 
  gt()
# gtsave(data = mod_results, 
#        filename = "model_compare.png", 
#        path = here::here("output"))

mod1_results <- f1_diag %>%
  dplyr::select(-nlminb_conv, -max_grad_name, -max_grad, - fleet, -sd) %>% 
  arrange(min(wAIC)) %>% 
  relocate(AIC, .before = dAIC) %>% 
  relocate(pd_Hess, .after = wAIC) %>% 
  gt() %>% 
  tab_header(
    title = "Fleet 1",
    subtitle = "Factorial model results"
  ) %>% 
  cols_label(
    model = "Model",
    parameters = "Parameter",
    mle_val = "Parameter estimate",
    sd_val = "Standard deviation",
    lwr_95 = "95 CI lower",
    upr_95 = "95 CI upper",
    pd_Hess = "pos-def Hessian"
  )

mod2_results <- f2_diag %>%
  dplyr::select(-nlminb_conv, -max_grad_name, -max_grad, - fleet, -sd) %>% 
  arrange(min(wAIC)) %>% 
  relocate(AIC, .before = dAIC) %>% 
  relocate(pd_Hess, .after = wAIC) %>% 
  gt() %>% 
  tab_header(
    title = "Fleet 2",
    subtitle = "Factorial model results"
  ) %>% 
  cols_label(
    model = "Model",
    parameters = "Parameter",
    mle_val = "Parameter estimate",
    sd_val = "Standard deviation",
    lwr_95 = "95 CI lower",
    upr_95 = "95 CI upper",
    pd_Hess = "pos-def Hessian"
  )

mod3_results <- f3_diag %>%
  dplyr::select(-nlminb_conv, -max_grad_name, -max_grad, - fleet, -sd) %>% 
  arrange(min(wAIC)) %>% 
  relocate(AIC, .before = dAIC) %>% 
  relocate(pd_Hess, .after = wAIC) %>% 
  gt() %>% 
  tab_header(
    title = "Fleet 3",
    subtitle = "Factorial model results"
  ) %>% 
  cols_label(
    model = "Model",
    parameters = "Parameter",
    mle_val = "Parameter estimate",
    sd_val = "Standard deviation",
    lwr_95 = "95 CI lower",
    upr_95 = "95 CI upper",
    pd_Hess = "pos-def Hessian"
  )

gtsave(data = mod1_results, 
       filename = "model_compare_f1.png", 
       path = here::here("output"))

gtsave(data = mod2_results, 
       filename = "model_compare_f2.png", 
       path = here::here("output"))

gtsave(data = mod3_results, 
       filename = "model_compare_f3.png", 
       path = here::here("output"))


# 2 waa and waa sd by fleet for selected model
# 3 diff btw 2020 benchmark wt at age and new

# plots
# WAA heatmaps for each fleet 
# SD heatmaps for each fleet
