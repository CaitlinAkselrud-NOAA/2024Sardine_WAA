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


# Compile and load in model -----------------------------------------------

wd <- here("src")
compile(here("src","GMRF_WAA.cpp"))
dyn.load(dynlib(here("src","GMRF_WAA")))

# Data --------------------------------------------------------------------

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


# Fill NA values ----------------------------------------------------------

rpl_sd <- rep(1.111, times = length(names(waa_df)))
names(rpl_sd) <- names(waa_df)

waa_std_df %<>% 
  replace_na(as.list(rpl_sd))

# waa filled in by fleet inside get_fleet_dat fxn
overall_mean_WAA <- waa_df %>%
  colMeans(na.rm = T) %>% 
  as.list()

# Fixed inputs ------------------------------------------------------------

projection_time <- 2

model_name <- "2020_sardine_"

# * fleet setup ---------------------------------------------------------------
n_fleets <- 3

# * other model settings --------------------------------------------------

newton_steps = 3

# Set up model functions ----------------------------------------
get_fleet_dat <- function(waa_dat, waa_std_dat, fleet_num)
{
  waa_f <- waa_dat %>% dplyr::filter(fleet == fleet_num)
  mean_WAA <- waa_f %>%
    colMeans(na.rm = T) %>% 
    as.list()
  waa_f %<>%
    replace_na(mean_WAA)
  waa_f %<>%
    replace_na(overall_mean_WAA) #if there is no fleet-specific mean, fill with overall mean
  
  waa_std_f = waa_std_dat %>% dplyr::filter(fleet == fleet_num)
  return(list(dat = waa_f, sd = waa_std_f))
}

TMB_setup <- function(proj_yrs = 2, waa_df, waa_std_df)
{
  # Number of projection years
  n_proj_years <- proj_yrs
  
  # Years
  years <- waa_df$year
  
  # Ages (goes from age 3 - 15+)
  ages <- str_extract_all(names(waa_df), '[0-9]+') %>% 
    unlist() %>% as.numeric()
  
  # Read in data weight at age matrix
  X_at <- waa_df %>% 
    select(starts_with("x")) %>% #keep only age cols
    t() 
  
  # Create projection columns (append to X_at matrix)
  proj_cols <- matrix(NA, nrow = length(ages), ncol = n_proj_years) 
  
  # Append NA for projection year
  X_at <- cbind(X_at, proj_cols) 
  
  # Read in standard deviations for weight at age matrix
  Xse_at <- waa_std_df %>% 
    select(starts_with("x")) %>% #keep only age cols
    t() 
  
  # Convert to CV
  Xcv_at <- sqrt( (exp(Xse_at^2) - 1) )
  
  # Now convert back to sd in lognormal space
  Xsd_at <- sqrt((log((Xcv_at)^2 + 1))/(log(10)^2))
  
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

# conditional model -------------------------------------------------------

# * choose factorial design -----------------------------------------------
model_fact <- fact_design(y = 0:1, 
                          c = 0:1, 
                          a = 0:1)


for(i in 1:n_fleets)
{
  fleet_name <- paste0("fleet",i,"_")
  fleet_num <- i
  
  waa_fleet <- get_fleet_dat(waa_dat = waa_df,
                             waa_std_dat = waa_std_df,
                             fleet_num = fleet_num)
  
  dat_setup <- TMB_setup(proj_yrs = projection_time,      # this is proj time steps (2 in this code = 2 seasons, aka 1 year)
                         waa_df = waa_fleet$dat, 
                         waa_std_df = waa_fleet$sd)
  
  # set conditional variance
  dat_setup$Var_Param <- 0 # Var_Param == 0 Conditional, == 1 Marginal
  
  # Now, input these components into a data list
  data_in <- dat_setup
  
  # Input parameters into a list
  parameters_in <- list( rho_y = 0, 
                         rho_a = 0,
                         rho_c = 0,
                         log_sigma2 = log(0.1), #
                         ln_L0 = log(9), #first length bin sardine
                         ln_Linf = log(28),  # last length bin for sardine
                         ln_k = log(0.15),
                         ln_alpha = log(3.5e-7), # Start alpha at a reasonable space 
                         # Starting value for alpha derived from a run where none of the rhos were estimated.
                         ln_beta = log(3), # Fix at isometric
                         ln_Y_at = array(0,dim=dim(data_in$X_at))) 
  
  
  # * run conditional model -------------------------------------------------
  
  models_cond <- run_model(map_factorial = model_fact, 
                           n.newton = newton_steps, 
                           data = data_in, 
                           parameters = parameters_in)
  
  save(models_cond, file = here("output", paste0(model_name, fleet_name, "cond_var_waa_models.RData")))
}

# marginal model ----------------------------------------------------------

data_in$Var_Param <- 1

# * run marginal model ----------------------------------------------------

models_marg <- run_model(map_factorial = model_fact, 
                    n.newton = newton_steps, 
                    data = data_in, 
                    parameters = parameters_in)

save(models_marg, file = here("output", paste0(model_name, fleet, "marg_var_waa_models.RData")))


# model selection ---------------------------------------------------------

models = models_marg

