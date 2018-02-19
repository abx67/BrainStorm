# Andrey Simonov
## Correction for initial conditions
## Run file
## October 2017

## load the libraries
library(data.table)
library(bayesm)
#library(xtable)
library(coda)
library(MASS)

rm(list = ls())

## set the run parameters 
## t = (3, 4, 5, 8, 10, 15, 25, 50, 100) # number of choices in a balanced panel

# simulate the data or use actual data? 
exercise_type <- 1 # 0 = simulate, 1 = estimate using old data, 2 = estimation using new data
  
# run on the grid or not?
flag_grid = 0

# balanced or unbalanced choice panel? 
flag_unb = TRUE # unbalanced panel 

# min number of purchases observed per household
threshold_number <- 3 

# should we simulate the markovian prices? 
flag_markovian_simulate <- TRUE

## some additional parameters
n_clust <- 100 # number of price clusters

## for simulation
t <- 900 # number of time periods to simulate the data
tburn <- 100 # number of time periods to drop
# we use the price process from the new data for estimation 
# transition probability matrix is computed using the median shopping frequency
# demand parameters are taken from the old (or new?) data estimates 

## Specifications for the forward simulation case with Markovian prices
draws_length <- 10 # number of steps for backward price simulation
number_price_simulations <- 3 # number of price simulations
flag_markovian <- TRUE # flag to run an actual markovian case

if (flag_grid == 0) {
  setwd("~/Dropbox/State Dependence and Initial Conditions/")
} 

###### LOAD/SIMULATE DATA #######

source('scripts/20.1.simulate_data.R')

###### PREPARE DATA FOR ESTIMATION AND ESTIMATE #######
source('scripts/20.2.reshape_data_for_estimation.R')

####### ANALYZE THE RESULTS ########
source('scripts/20.3.examine_results.R')

