# Andrey Simonov
## Correction for initial conditions
## Reshape data for estimation 
## October 2017

## parameters are defined in the 20.core_code.R file
## not to run separately -- sourced in the 20.core_code.R file

## find how many initial choices were the choices of the outside option
count_out <- array(1,n) # number of observation to change in case of unknown initial state

all_purchased <- TRUE
j <- 1
for (i in 1:n) {
  flag <- 0
  while (flag == 0) {
    if (y[j] == m) {
      count_out[i] <- count_out[i] + 1 # count one if person chose outside option
      j <- j + 1 # next observation if person chose outside option
    }
    if (y[j] != m) flag <- 1 # exit if not an outside option
    if (j == tcum[i]) {flag <- 1; all_purchased <- FALSE} # exit if last observation for this person
  }
  j <- tcum[i] + 1 # next person's first observation
}


store_data <- NULL

j1 <- 1
j2 <- 1
for (i in 1:n) {
  if (i == 1) { 
    # works only for 4 states!!!! grid_X0 specific assignment
    store_data[[i]] <- list(y=y[j1:tcum[i]], X=X[j2:(tcum[i]*m),], 
                            count_out = count_out[i], initial_price_id = initial_price_id[i], 
                            transition_matrix_median_steps = transition_matrix_median_steps,price_transition_states = as.matrix(price_transition_states),vec_price_states_probs = price_states_probs[,probs],draws_length=draws_length,number_price_simulations=number_price_simulations,flag_markovian=flag_markovian, flag_know_state=FALSE)
  }
  else {
    store_data[[i]] <- list(y=y[j1:tcum[i]], X=X[j2:(tcum[i]*m),], count_out = count_out[i], initial_price_id = initial_price_id[i]) # store the data
  }
  j1 <- tcum[i] + 1
  j2 <- tcum[i]*m + 1
}



## data with dropped observations -- to test the upper bound
## assignment works in all count_out > 0 
print(paste0("Number of units with known initial state: ", sum(count_out==0)))

store_data_dropped <- NULL
# we lose 12.3 percent of out data in this case 

## need extra inputs: s1 (initial choice state) and price when s1 was created

j1 <- 1 + count_out[1]
j2 <- 1 + count_out[1] * m
for (i in 1:n) {
  initial_price_state <- array(0,c(m*(m-1), m+1)) # matrix of all possible initial states with observed prices in count[i]-1
  temp_price_i <- X[(1:m) + m*(count_out[i]-1),m] # store price in the observation before the initial one
  for (j in 1:(m-1)) { # create a matrix with each possible state
    initial_price_state[(j-1)*m + 1:m,1:m] <- diag(m)
    initial_price_state[(j-1)*m + 1:m,m]   <- temp_price_i
    initial_price_state[(j-1)*m + j,m+1]   <- 1    
  }
  if (i == 1) { 
    store_data_dropped[[i]] <- list(y=y[j1:tcum[i]], X=X[j2:(tcum[i]*m),],
                                    count_out = count_out[i], initial_price_id = initial_price_id[i], initial_price_state_dropped = initial_price_state, s1 = y[j1-1],
                                    transition_matrix_median_steps = transition_matrix_median_steps,price_transition_states = as.matrix(price_transition_states),vec_price_states_probs = price_states_probs[,probs],draws_length=draws_length,number_price_simulations=number_price_simulations,flag_markovian=flag_markovian, flag_know_state=TRUE)
  }
  else {
    store_data_dropped[[i]] <- list(y=y[j1:tcum[i]], X=X[j2:(tcum[i]*m),], count_out = count_out[i], initial_price_id = initial_price_id[i], initial_price_state_dropped = initial_price_state, s1 = y[j1-1]) 
  }
  j1 <- tcum[i] + 1 + count_out[i+1]
  j2 <- tcum[i]*m + 1 + count_out[i+1]*m
}


### MCMC estimation

# load necessary functions

source('scripts/9_1_corrLike.R')
source('scripts/9_1_llmnl_initial_mcprice.R')

### 

nummix <- 1 # number of mixtures of normals
keep <- 5 # thinning parameter
R <- 1000 # number of MCMC draws
# a = rep(3/nummix, nummix) # prior on dirichlet distribution

#a = rep(.5/ncomp,ncomp)
#Amu = matrix(1/16,ncol=1)

Prior1=list(ncomp=nummix) # prior specification (number of mixtures)
Mcmc1=list(R=R,keep=keep) # MCMC specification
Data1=list(p=m,lgtdata=store_data) # data specification, p is number of alts
Data1_dropped=list(p=m,lgtdata=store_data_dropped) # data specification, p is number of alts

# specification (1)
out_wout <- rhierMnlRwMixture(Data=Data1,Prior=Prior1,Mcmc=Mcmc1) # results without correction
# specification (2)
out_wout_d <- rhierMnlRwMixture(Data=Data1_dropped,Prior=Prior1,Mcmc=Mcmc1) # results without correction
plot(apply(out_wout$betadraw[,6,],2,mean))
plot(apply(out_wout_d$betadraw[,6,],2,mean))
# specification (3)
random_seed <- sample.int(100000,1)
set.seed(random_seed)
Data1_dropped$lgtdata[[1]]$flag_markovian = FALSE
out_with_d <- rhierMnlRwMixture_inLike_X0_mcprice(Data=Data1_dropped,Prior=Prior1,Mcmc=Mcmc1) # results with correction
plot(apply(out_with_d$betadraw[,6,],2,mean))
# specification (5)
set.seed(random_seed)
Data1_dropped$lgtdata[[1]]$flag_markovian = TRUE
out_with_d_markovian <- rhierMnlRwMixture_inLike_X0_mcprice(Data=Data1_dropped,Prior=Prior1,Mcmc=Mcmc1) # results with correction
plot(apply(out_with_d_markovian$betadraw[,6,],2,mean))
# specification (4)
set.seed(random_seed)
Data1$lgtdata[[1]]$flag_markovian = FALSE
out_with <- rhierMnlRwMixture_inLike_X0_mcprice(Data=Data1,Prior=Prior1,Mcmc=Mcmc1) # results with correction
pdf(file = 'false_markov_desktop.pdf')
plot(apply(out_with$betadraw[,6,],2,mean))
dev.off()
# specification (6)
set.seed(random_seed)
Data1$lgtdata[[1]]$flag_markovian = TRUE
out_with_markovian <- rhierMnlRwMixture_inLike_X0_mcprice(Data=Data1,Prior=Prior1,Mcmc=Mcmc1) # results with correction
pdf(file = 'true_markov_desktop.pdf')
plot(apply(out_with_markovian$betadraw[,6,],2,mean))
dev.off()

if (exercise_type == 0) save(out_wout, out_wout_d, out_with, out_with_d, out_with_d_markovian, out_with_markovian, n, alpha, eta, gamma, Data1, Data1_dropped, nummix, random_seed, file = paste0("data/outputs/new_simulation_desktop_exercise_type_",exercise_type,".RData"))
if (exercise_type != 0) save(out_wout, out_wout_d, out_with, out_with_d, out_with_d_markovian, out_with_markovian, n, Data1, Data1_dropped, nummix, random_seed, file = paste0("data/outputs/new_simulation_desktop_exercise_type_",exercise_type,".RData"))

