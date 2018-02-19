rmultinomF=
  function(p) {
    return(sum(runif(1) > cumsum(p))+1)
  }

llmnl_initial_mcprice= 
  function(beta,y,X,count_out,initial_price_id,initial_price_state_dropped,s1,transition_matrix_median_steps,price_transition_states,vec_price_states_probs,draws_length,number_price_simulations,flag_markovian, flag_know_state) 
  { 
    
    ## data input types: 
    # beta = vector
    # y = vector
    # X = matrix
    # count_out = integer
    # initial_price_id = integer
    # initial_price_state_dropped = matrix
    # s1 = integer 
    # transition_matrix_median_steps = matrix
    # price_transition_states = matrix
    # vec_price_states_probs = vector
    # draws_length = integer 
    # number_price_simulations= integer
    # flag_markovian = boolean 
    # flag_know_state = boolean
    
      
    n=length(y) # number of choices
    j=nrow(X)/n # number of brands
    nvar <- ncol(X)
    
    ########## added to llmnl function ##############
    if (count_out != 0) { ### if the initial state is not known
      
    #### New August 2017 -- Markovian prices ### 
    if (flag_markovian) {
      X_forward <- cbind(diag(j)[,1:(j-1)],0,0)
      dif_probs <- array(0,c(number_price_simulations,j))
      for (j_markov in 1:number_price_simulations){
        # draw a price vector
        price_draw_this_iter <- array(0, draws_length)
        price_draw_this_iter[1] <- rmultinomF(transition_matrix_median_steps[initial_price_id,])
        #price_draw_this_iter[1] <- which(rmultinom(1,1,transition_matrix_median_steps[initial_price_id,])==1)
        for (i in 2:draws_length) {
          price_draw_this_iter[i] <- rmultinomF(transition_matrix_median_steps[price_draw_this_iter[i-1],])
          #price_draw_this_iter[i] <- which(rmultinom(1,1,transition_matrix_median_steps[price_draw_this_iter[i-1],])==1)
        }
        # reverse the simulated price vector
        price_draw_this_iter <- price_draw_this_iter[draws_length:1]
        # initial random state
        #X_forward[,nvar] <- c(rmultinom(1,1,rep(1/(j-1),j-1)),0)
        X_forward[rmultinomF(rep(1/(j-1),j-1)),nvar] <- 1
        X_forward[,j] <- price_transition_states[price_draw_this_iter[1],1:j]
        # forward simulate the choices
        for (i in 2:draws_length) {
          exp_xbeta_forward <- exp(X_forward%*%beta)
          prob_forward <- exp_xbeta_forward / sum(exp_xbeta_forward)
          #X_forward[,nvar] <- rmultinom(1,1,prob_forward)
          X_forward[rmultinomF(prob_forward),nvar] <- 1
          X_forward[,j] <- price_transition_states[price_draw_this_iter[i],1:j]
        }
        dif_probs[j_markov,] <- prob_forward
      }
      tmp_x <- apply(dif_probs,2,mean)
      # marginal initial state probability vector
      vec <- tmp_x[1:(j-1)]/sum(tmp_x[1:(j-1)])
    } else {
      ####### Older iid price case #########
      xbeta <- t(c(beta[1:(j-1)],0) + t(price_transition_states[,1:j] * beta[j]))
      Prob0 = array(0,c((j-1),(j-1)))
      for (j_state in 1:(j-1)) {
      xbeta_state <- xbeta
      xbeta_state[,j_state] <- xbeta[,j_state] + beta[nvar]
      xb0 <- exp(xbeta_state)
      prob0 <- xb0/apply(xb0,1,sum)
      tmp_prob <- vec_price_states_probs %*% prob0
      Prob0[,j_state] <- tmp_prob[1:(j-1)] # transposed transition matrix
      Prob0[j_state,j_state] <- Prob0[j_state,j_state] + tmp_prob[j] # if choose outside option, state stays the same
      }
      aaa = eigen(Prob0) 
      vec = abs(aaa$vectors[,1]) 
      if (is.complex(vec)) vec <- Re(vec) # for odd rare cases of i numbers
      vec = vec / sum(vec) # marginal probability
      # these are marginal probabilites of chosings states for average prices
    }
      if (!flag_know_state) { # if we are in the case when the data was not dropped and initial state is not known
      s0 <- rmultinomF(vec) # draw for the initial state
      X[1:(count_out*j),nvar] <- array(0,count_out*j)
      X[c(s0 + 0:(count_out-1)*j),nvar] <- array(1,count_out) # assign drawn initial state to all observations until the first brand is chosen
      }
    }
    ######################################################################################################
    
    Xbeta=X%*%beta
    Xbeta=matrix(Xbeta,byrow=T,ncol=j)
    ind=cbind(c(1:n),y)
    xby=Xbeta[ind]
    Xbeta=exp(Xbeta)
    iota=c(rep(1,j))
    denom=log(Xbeta%*%iota)
    
    if (flag_know_state) {
      xb_init = matrix(exp(initial_price_state_dropped%*%beta), ncol = j, byrow = T) # matrix of expxbeta by each possible (m-1) state
      prob_init = xb_init/apply(xb_init,1,sum) # matrix of choice probabilities
      vec_with_init = vec%*%prob_init # marginal probabilities
      return(sum(xby-denom) + log(vec_with_init[s1]))
    } else {
      return(sum(xby-denom))
    }
  }
