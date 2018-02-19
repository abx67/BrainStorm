library(Rcpp)
library(bayesm)
library(devtools)#,lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library(RcppArmadillo)#,lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
betac=rnorm(10)
y=sample(1:8,10,replace = TRUE)
X=matrix(runif(900),90)
count_out=1
initial_price_id=5
initial_price_state_dropped=matrix(runif(720),72,10)*10
s1=5
transition_matrix_median_steps=matrix(runif(100),10)
price_transition_states=matrix(rnorm(100),10)
vec_price_states_probs=runif(10)
draws_length=5
number_price_simulations=10
flag_markovian=FALSE    #FALSE
flag_know_state=TRUE    #TRUE
#TRUE TRUE exact   #TRUE FALSE or FALSE FALSE decimal right 3 loc  #

sourceCpp("llmnl_initial_mcprice_bayesm.cpp")
source("9_1_llmnl_initial_mcprice2.r")
llmnl_initial_mcprice(betac,y,X,count_out,initial_price_id,initial_price_state_dropped,s1,
                       transition_matrix_median_steps,price_transition_states,
                       vec_price_states_probs,draws_length,number_price_simulations,
                       flag_markovian,flag_know_state)
rm(llmnl_initial_mcprice)

sourceCpp("test2.cpp")
bb(beta,initial_price_state_dropped)
