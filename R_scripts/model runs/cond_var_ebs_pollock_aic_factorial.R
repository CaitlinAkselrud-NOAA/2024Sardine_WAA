# Purpose: To run a series of triple separable models testing for evidence of
# cohort/age/year effects for EBS pollock w/ conditional variance
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 1/16/23
# 
# Modified: Caitlin Allen Akselrud 12.12.23
# contact: caitlin.allen_akselrud@noaa.gov
# Notes: set up for 2024 sardine benchmark assessment; run this code first
# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(TMB)
library(cowplot)
library(janitor)

# Compile and load in model
wd <- here("src")
compile(here("src","GMRF_WAA.cpp"))
dyn.load(dynlib(here("src","GMRF_WAA")))

# raw data

# Load in WAA matrix (only use fishery data)
# waa_df <- read_csv(here("data", "ebs_waa.csv")) %>% 
#   filter(source == "fishery") %>% 
#   dplyr::select(-source)

waa_df_raw <- read_csv(here("data", "2020_sardine_waa.csv")) %>% 
  janitor::clean_names()

# Load in std for WAA matrix
# waa_std_df <- read_csv(here("data", "ebs_waa_std.csv")) %>% 
#   filter(source == "fishery") %>% 
#   dplyr::select(-source)

waa_std_df_raw <- read_csv(here("data", "2020_sardine_waa_std.csv")) %>% 
  janitor::clean_names()
# CIA: **DUMMY FILE, put in actual std devs**

# clean data

# fleets
# 1 = MexCal S1
# 2 = MexCal S2
# 3 = PNW
# 4 = AT survey
waa_df <- waa_df_raw %>% 
  mutate(year = case_when(seas == 1 ~ number_yr + 0.0,
                          seas == 2 ~ number_yr + 0.5)) %>% 
  dplyr::filter(fleet > 0) %>% 
  dplyr::select_if(is.numeric) #remove columns with notes
  

waa_std_df <- waa_std_df_raw %>% 
  mutate(year = case_when(seas == 1 ~ number_yr + 0.0,
                          seas == 2 ~ number_yr + 0.5)) %>% 
  dplyr::filter(fleet > 0) %>% 
  dplyr::select_if(is.numeric) #remove cols  with notes

waa_f1 <- waa_df %>% dplyr::filter(fleet == 1)
waa_f2 <- waa_df %>% dplyr::filter(fleet == 2)
waa_f3 <- waa_df %>% dplyr::filter(fleet == 3)
waa_f4 <- waa_df %>% dplyr::filter(fleet == 4)

waa_sd_f1 <- waa_std_df %>% dplyr::filter(fleet == 1)
waa_sd_f2 <- waa_std_df %>% dplyr::filter(fleet == 2)
waa_sd_f3 <- waa_std_df %>% dplyr::filter(fleet == 3)
waa_sd_f4 <- waa_std_df %>% dplyr::filter(fleet == 4)


# Set up TMB data ----------------------------------------

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



# Set up TMB Model --------------------------------------------------------

dat_setup_f4 <- TMB_setup(proj_yrs = 2,      # this is proj time steps (2 in this code = 2 seasons, aka 1 year)
                          waa_df = waa_f4, 
                          waa_std_df = waa_sd_f4)

dat_setup_f4$Var_Param <- 0

# Now, input these components into a data list
data <- dat_setup_f4
  # years
  # ages
  # X_at
  # Xsd_at
  # ay_Index
  # n_proj_years # proj time: 2 seasons = 1 year in sardine mode
  # Var_Param = 0) # Var_Param == 0 Conditional, == 1 Marginal

# Input parameters into a list
# CIA: need to adjust for sardine
parameters <- list( rho_y = 0,
                    rho_a = 0,
                    rho_c = 0,
                    log_sigma2 = log(0.1),
                    ln_L0 = log(45),
                    ln_Linf = log(80),  # Fixed at arbitrary value
                    ln_k = log(0.15),
                    ln_alpha = log(3.5e-7), # Start alpha at a reasonable space 
                    # Starting value for alpha derived from a run where none of the rhos were estimated.
                    ln_beta = log(3), # Fix at isometric
                    ln_Y_at = array(0,dim=dim(X_at))) 


# Run factorial models ----------------------------------------------------

# Generate full factorial design
map_factorial <- tidyr::crossing(rho_y = 0:1, rho_c = 0:1, rho_a = 0:1) %>% 
  data.frame()

# Define number of extra newton steps we want to take
n.newton <- 3

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

save(models, file = here("output", "2020_sardine_cond_var_waa_models.RData"))

