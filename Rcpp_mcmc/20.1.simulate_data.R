# Andrey Simonov
## Correction for initial conditions
## Simulate data 
## October 2017

## parameters are defined in the 20.core_code.R file
## not to run separately -- sourced in the 20.core_code.R file

## INPUTS to change:
# -- prices_actual: observed prices, currently based on prices in the Nielsen data. used to estimate the price process. 
#                   Note: code is for 5 choice alternatives, needs to be adjusted if it is different (creating ids section)
# -- median_purch_freq: median purchase frequency 
# -- n: number of observations. currently 1829, based on Nielsen data
# -- purchase_per_hh: number of purchases per household, currently empirical distribution from Nielsen data
# -- parameters eta, alpha, gamma: currently based on the results in Dube, Hitsch, Rossi (2010)

## load actual data to get the price realizations
if (exercise_type != 1) {
# new data
load("data/inputs/newdata.RData")
prices_actual <- unique(data_groups_bal[, list(groups, week_end, price, store_code_uc)])
# median purchase frequency
x <- unique(data_groups_bal[, list(id = household_code, period = as.Date(week_end))]); x_id <- unique(x[, id]); setkey(x, period)
purch_freq <- NULL
for (i in 1:length(x_id)) {x_tmp <- x[id == x_id[i], period]; purch_freq <- c(purch_freq,(x_tmp[2:length(x_tmp)] - x_tmp[2:length(x_tmp)-1])/7)}
purch_freq <- purch_freq[!is.na(purch_freq)]
median_purch_freq <- median(purch_freq)
}

if (exercise_type == 1) {
## old data
data <- fread("data/inputs/data.csv")
prices_actual <- unique(data[, list(groups = brand, week_end = period, price, store_code_uc= nhsout)])
# median purchase frequency
x <- unique(data[, list(id, period)]); x_id <- unique(x[, id])
purch_freq <- NULL
for (i in 1:length(x_id)) {x_tmp <- x[id == x_id[i], period]; purch_freq <- c(purch_freq,(x_tmp[2:length(x_tmp)] - x_tmp[2:length(x_tmp)-1]))}
median_purch_freq <- median(purch_freq)
}

## continue 
setkey(prices_actual, store_code_uc, week_end, groups)
prices_actual[, store_code_uc := NULL]
prices_actual[, week_end := NULL]

## data parameters
m <- prices_actual[, max(groups)] # number of choice alternatives (with the outside option)

## reduce the number of price states to 100
## cluster prices 
prices_all <- data.table(data.frame(t(array(prices_actual[, price], c(m,dim(prices_actual)[1]/m)))))
prices_clustering <- kmeans(prices_all, n_clust, iter.max = 1000)
prices_all[, id := prices_clustering$cluster]
# cluster means and price state ids
price_transition_states <- data.table(data.frame(prices_clustering$centers, id = 1:n_clust))
setkey(price_transition_states, id)

## estimate transition probability matrix for the markovian prices 

## create a transition probability matrix
transition_matrix <- array(0, c(n_clust,n_clust))
for (i in 2:dim(prices_all)[1]) transition_matrix[prices_all[i-1,id],prices_all[i,id]] <- transition_matrix[prices_all[i-1,id],prices_all[i,id]] + 1
transition_matrix <- transition_matrix/as.matrix(apply(transition_matrix, 1,sum))[,rep(1,dim(price_transition_states)[1])]

## define transition probability matri for a media frequency of consumer visits 
# function to change transition probability matrix to the median frequency of shopping day difference 
steps_t_matrix_fn <- function(t_matrix, n_steps){ y <- transition_matrix; for (i in 1:n_steps) {y <- y %*% transition_matrix}; return(y)}
# create price states transition probability matrix for median_purch_freq steps 
transition_matrix_median_steps <- steps_t_matrix_fn(transition_matrix, median_purch_freq)

## case of iid prices -- price state probabilities
prices_all[, ind := 1]
price_states_probs <- prices_all[, sum(ind),by='id']
setkey(price_states_probs, id)
price_states_probs[, probs := V1/sum(V1)]
price_states_probs[, V1 := NULL]

##### choice data -- simulate

if (exercise_type == 0) {
# data dimensions
n <- length(unique(data_groups_bal[, household_code])) # number of people
k <- m - 1 # number of product options

## number of purchases per person
purchase_per_hh <- data_groups_bal[groups != 5, sum(bought), by='household_code']
setkey(purchase_per_hh, household_code)
# sample the number of purchase observations from the empirical distribution
if (flag_unb) obs_number <- sample(purchase_per_hh[V1 >= threshold_number, V1], size = n, replace = T)
if (flag_unb==F) obs_number <- rep(threshold_number,n)

## parameters
eta <- rnorm(n, -1.18, 1.17) # price coef
alpha <- cbind(rnorm(n,-1.62, 2.9), 
               rnorm(n,0.39, 2.84), 
               rnorm(n,-1.32, 2.59), 
               rnorm(n,-0.58, 3.31)) # brand FEs
gamma = rnorm(n, 1, 1.103) # state dependence 

## simulate prices for t periods
price_index <- array(0, t)
price_index[1] <- which(rmultinom(1,1,price_states_probs[, probs]) == 1)
if (flag_markovian_simulate) { # Markovian case
  for (i in 2:t) price_index[i] <- which(rmultinom(1,1,transition_matrix_median_steps[price_index[i-1],]) == 1)
} else { # iid case
  for (i in 2:t) price_index[i] <- which(rmultinom(1,1,price_states_probs[, probs]) == 1)
}

price <- price_transition_states[price_index]
price <- as.matrix(price[, 1:4, with = F])

## choice simulation - copied from old code

state <- array(0,dim = c(n, t))
choice2 <- array(0, dim = c(n,t)) # matrix of choices
xbeta <- array(0, dim = c(n,m)) 
# first choice
xbeta[,1:k] <-  alpha + matrix(rep(eta,k), ncol = k)*matrix(rep(price[1,],n), nrow = n, byrow = T)
expxbeta <- exp(xbeta)
prob <- expxbeta / matrix(rep(apply(expxbeta,1,sum), m), ncol = m, byrow = F)
choice2[, 1] <- apply(prob, 1, function(x) which(rmultinom(1,1,x) == 1))
for (i in 2:t) {
  state[,i][choice2[,i-1] != m] <- choice2[,i-1][choice2[,i-1] != m] # if choice i-1 != outside option, update the state
  state[,i][choice2[,i-1] == m] <- state[,i-1][choice2[,i-1] == m] # if choice i-1 == outside option, keep the state
  xbeta[,1:k] <-  alpha + matrix(rep(eta,k), ncol = k)*matrix(rep(price[i,],n), nrow = n, byrow = T)
  for (i2 in 1:n) if (state[i2,i] != 0) xbeta[i2,state[i2,i]] <- xbeta[i2,state[i2,i]] + gamma[i2] # add state dependence
  expxbeta <- exp(xbeta)
  prob <- expxbeta / matrix(rep(apply(expxbeta,1,sum), m), ncol = m, byrow = F)
  choice2[, i] <- apply(prob, 1, function(x) which(rmultinom(1,1,x) == 1))
}

## check if outside option is chosen with a reasonable frequency
table(choice2)/n/t
table(choice2[1,])/t
table(choice2[2,])/t
table(choice2[3,])/t
table(choice2[4,])/t
table(choice2[5,])/t
table(choice2[6,])/t
table(choice2[7,])/t

## delete burn-in period
price <- price[(tburn+1):t,]
choice2 <- choice2[,(tburn+1):t]
t <- t - tburn
price_index <- price_index[(tburn+1):t]

## first time we observe the household
start_obs <- array(0,n) # where do we start observing the household
for (i in 1:n) start_obs[i] <- round(runif(1,1,t/10)) 

## keep exactly obs_number of purchases for each household
a <- array(0,n)
ppl_drop <- NULL # people with less than obs_number choices
for (i in 1:n) {
  flst <- 0
  while (flst < obs_number[i] & a[i] < (t - start_obs[i]-1)) {
    flst <- flst + (choice2[i, start_obs[i] + a[i]] != 5)
    a[i] <- a[i] + 1
  }
  if (flst < obs_number[i]) ppl_drop <- c(ppl_drop,i)
}

# remove people with less than obs_number purchase observations 
if (!is.null(ppl_drop)) {
  n <- n - length(ppl_drop)
  choice2 <- choice2[-ppl_drop, ]
  a <- a[-ppl_drop]
  start_obs <- start_obs[-ppl_drop]
  # update eta, alpha, gamma
  eta <- eta[-ppl_drop]
  gamma <- gamma[-ppl_drop]
  alpha <- alpha[-ppl_drop,]
}


## construct the data for estimation
ppl <- rep(c(1:n), a) # individuals-observations
N <- length(ppl) # number of observations total 

# choices
y <- NULL 
for (i in 1:n) y <- c(y,choice2[i,0:(a[i]-1) + start_obs[i]])

# X matrix
tmp <- cbind(diag(k),0)
X <- matrix(rep(tmp,N), ncol = k, byrow = T) # dummy columns
p <- NULL # prices for observations
for (i in 1:n)   p <- rbind(p, price[0:(a[i]-1) + start_obs[i],])
X <- cbind(X,as.vector(t(cbind(p,0)))) # prices
# state
acum <- cumsum(a)
state <- c(0,y[2:N-1]) # states are choices
state[acum[2:n-1]+1] <- 0 # for new people state are zero (first observations)
for (i in 1:N) {if (state[i] == m) state[i] <- state[i-1]} # if outside option - state is state before

tmp2 <- array(0,N*(k+1)) # map it to X matrix
tmp <- state + m*c(1:N-1) # locations of states
tmp2[tmp] <- 1
tmp2[seq(k+1, N*(k+1), k+1)] <- 0 # exclude outside options
X <- cbind(X, tmp2)

## initial price state
initial_price_id <- array(0, n)
for (i in 1:n) initial_price_id[i] <- price_index[start_obs[i]] # index of price in the first observation for each household

## mapping notation back to 13_run_newdata code - what to take the analysis from there
t <- a
tcum <- acum

#### output structure: 
# X, y -- observations and choices
# t, tcum  -- number of choices per person
# estimates of the prices process
}

if (exercise_type == 1) {
  m <- max(data[, choice]) # number of brands plus one (outside option)
  N <- dim(data)[1] / m # number of observations
  y <- data[seq(1,dim(data)[1],by=m),choice] # choices
  X <- data[,c(6:9,5,11), with = F]
  X[,price := price * 16] # change prices to price per 16 oz tube
  X <- as.matrix(X)
  
  n <- max(data[,id]) # number of people
  t <- as.vector(table(data[,id])) / m # number of time observations per person
  tcum <- cumsum(t)
  
  ## initial price state
  a <- data.table(data.frame(matrix(data[, price], ncol = m, byrow = T)))
  b <- unique(prices_all)
  setkey(b,id)
  b[, ind := NULL]
  prices_with_states <- merge(a, b, by.x = c("X1", "X2", "X3", "X4", "X5"), by.y = c("X1", "X2", "X3", "X4", "X5")) # match cluster ids to actual prices
  
  initial_price_id <- prices_with_states[c(1,(tcum + 1)[2:n-1]), id]
}

if (exercise_type == 2) {
  ppl <- match(data_groups_bal[, household_code], unique(data_groups_bal[, household_code]))
  m <- max(data_groups_bal[, groups])
  N <- length(ppl)/m
  n <- n_ppl <- max(ppl)
  
  y <- data_groups_bal[, bought*groups]
  y <- y[y != 0]
  
  ##### create a state variable
  tmp <- c(0, y[1:(N-1)])
  tmp[ppl[seq(1,N*m,by=m)] != c(0, ppl[seq(1,N*m,by=m)][1:(N-1)])] <- 0
  while (max(tmp) == m) {
    tmp2 <- tmp == m
    tmp[tmp2] <- tmp[c(tmp2[2:(N)], FALSE)]
  }
  tmp2 <- tmp + c(1:N-1)*m
  tmp2 <- tmp2[tmp != 0]
  tmp3 <- array(0, N*m)
  tmp3[tmp2] <- 1
  data_groups_bal[, state := tmp3]
  
  X <- data_groups_bal[, list(d1, d2, d3, d4, price, state)]
  X <- as.matrix(X)
  
  # number of observations ber household
  t <- as.vector(table(data_groups_bal[,household_code]))/m
  tcum <- cumsum(t)
  
  ## initial price state
  a <- data.table(data.frame(matrix(data_groups_bal[, price], ncol = m, byrow = T)))
  b <- unique(prices_all)
  setkey(b,id)
  b[, ind := NULL]
  prices_with_states <- merge(a, b, by.x = c("X1", "X2", "X3", "X4", "X5"), by.y = c("X1", "X2", "X3", "X4", "X5")) # match cluster ids to actual prices
  
  initial_price_id <- prices_with_states[c(1,(tcum + 1)[2:n-1]), id]
  
}
