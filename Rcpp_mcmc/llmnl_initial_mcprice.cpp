#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


NumericMatrix my_cbind(NumericMatrix a, NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
};

NumericMatrix mmult(NumericMatrix m1, NumericMatrix m2){
  if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
    }
  }
  return out;
};

NumericMatrix transpose(NumericMatrix X){
  NumericMatrix ans (X.ncol(), X.nrow()); 
  for (int i = 0; i < X.nrow(); i++){ 
    for (int j = 0; j < X.ncol(); j++){ 
      ans(j,i) = X(i, j); 
    } 
  } 
  return(ans);
};

int my_rmultinomF(NumericVector a){
  Rcpp::Environment global = Rcpp::Environment::global_env();
  Function rmultinomF = global["rmultinomF"];
  NumericVector x = rmultinomF(a);
  return x[0];
};

NumericVector first_abs_eigenvec(NumericMatrix a){
  Rcpp::Environment global = Rcpp::Environment::global_env();
  Function getEigenVectors = global["getEigenVectors"];
  ComplexMatrix EigenVec = getEigenVectors(a);
  return Mod(EigenVec(_,0));
};
  
// [[Rcpp::export]]
NumericVector llmnl_initial_mcprice(NumericVector beta,NumericVector y,NumericMatrix X,
                                    int count_out,int initial_price_id,NumericMatrix initial_price_state_dropped,
                                    int s1,NumericMatrix transition_matrix_median_steps,NumericMatrix price_transition_states,
                                    NumericVector vec_price_states_probs,int draws_length,int number_price_simulations,
                                    bool flag_markovian,bool flag_know_state) {
  //envoronment
  //Rcpp::Environment global = Rcpp::Environment::global_env();
  //Function rmultinomF = global["rmultinomF"];

  //main
  int n=y.size(); // number of choices
  int j=X.nrow()/n; // number of brands
  int nvar = X.ncol();
  
  NumericVector vec;
  NumericMatrix mat_beta( beta.size() , 1 , beta.begin() ); //convert beta to matrix
  
  // added to llmnl function 
  if (count_out != 0) { // if the initial state is not known
    
    // New August 2017 -- Markovian prices ### 
    if (flag_markovian) {
      NumericMatrix X_forward(j,j+1);
      std::fill(X_forward.begin(),X_forward.end(),0);
      for(int i=0;i<j-1;i++){X_forward(i,i)=1;}
      NumericMatrix dif_probs(number_price_simulations,j);
      std::fill(dif_probs.begin(),dif_probs.end(),0);
      for(int j_markov=0;j_markov<number_price_simulations;j_markov++){
        // draw a price vector
        NumericVector price_draw_this_iter(draws_length,0.0);
        price_draw_this_iter[draws_length-1] = my_rmultinomF(transition_matrix_median_steps(initial_price_id,_));
        for (int i=1;i<draws_length;i++){
          price_draw_this_iter[draws_length-i-1] = my_rmultinomF(transition_matrix_median_steps(price_draw_this_iter[i],_));   //check later
        }
        // price_draw_this_iter <- price_draw_this_iter[draws_length:1] //no need anymore
        
        //# initial random state
        NumericVector temp(j-1,1/(j-1));
        X_forward(my_rmultinomF(temp),nvar) = 1;
        Range ran(0,j-1);
        X_forward(_,j-1) = price_transition_states(price_draw_this_iter[0],_);//require column = j
        //# forward simulate the choices
        NumericMatrix prob_forward;
        for (int i=1;i<draws_length;i++){
          NumericMatrix exp_xbeta_forward(j,1);
          //NumericMatrix mat_beta( beta.size() , 1 , beta.begin() ); //convert beta to matrix
          exp_xbeta_forward = mmult(X_forward,mat_beta);
          
          prob_forward = exp_xbeta_forward / sum(exp_xbeta_forward); //differ from original(mat/vec)
          X_forward(my_rmultinomF(prob_forward),nvar) = 1;    //attention to mat type although work
          X_forward(_,j) = price_transition_states(price_draw_this_iter[i],_);
        }
        dif_probs(j_markov,_) = prob_forward(_,0);
      }
      NumericVector tmp_x;
      for(int i=0;i<dif_probs.ncol();i++){    //correspond to apply(,mean)
        tmp_x[i]=mean(dif_probs(i,_));
      }
      //# marginal initial state probability vector
      NumericVector vec;
      Range ran2(0, j-1-1);
      vec = tmp_x[ran2]/sum(tmp_x[ran2]);
    } else {
      //####### Older iid price case #########
      NumericMatrix xbeta(price_transition_states.ncol(),price_transition_states.nrow());
      NumericVector tmp_beta = beta;  tmp_beta[j-1] = 0;
      NumericMatrix tmp_price_beta = transpose(price_transition_states*beta[j-1]);
      for (int i=0;i<tmp_price_beta.ncol();i++){
        xbeta(_,i)=tmp_price_beta(_,i)+tmp_beta;
        }
      NumericMatrix Prob0(j-1,j-1);
      std::fill(Prob0.begin(),Prob0.end(),0);
      for(int j_state=0;j_state<j-1;j_state++) {
        NumericMatrix xbeta_state=xbeta;
        xbeta_state(_,j_state) = xbeta(_,j_state) + beta[nvar];
        NumericMatrix xb0(xbeta_state.nrow(),xbeta_state.ncol());
        for(int i=0;i<xb0.ncol();i++){xb0(_,i) = exp(xbeta_state(_,i));}
        NumericMatrix prob0(xb0.nrow(),xb0.ncol());  //case matter??
        for(int i=0;i<xb0.nrow();i++){prob0(i,_) = xb0(i,_)/sum(xb0(i,_));}
        NumericMatrix tmp_prob;
        NumericMatrix mat_price_states_probs(1,vec_price_states_probs.size(),vec_price_states_probs.begin());
        tmp_prob = mmult(mat_price_states_probs,prob0);
        Range ran(0,j-1-1);
        NumericVector tmp_prob_vec = wrap(tmp_prob(_,ran));
        Prob0(_,j_state) =  tmp_prob_vec; //# transposed transition matrix
        Prob0(j_state,j_state) = Prob0(j_state,j_state) + tmp_prob[j]; //# if choose outside option, state stays the same
      }
      //aaa = eigen(Prob0) 
      //vec = abs(aaa$vectors[,1]);

      //if (is.complex(vec)) vec <- Re(vec) # for odd rare cases of i numbers
      //    vec = vec / sum(vec) # marginal probability
      vec = first_abs_eigenvec(Prob0);
      vec = vec / sum(vec);
      //# these are marginal probabilites of chosings states for average prices
    }
    if (!flag_know_state) {// # if we are in the case when the data was not dropped and initial state is not known
      int s0 = my_rmultinomF(vec); //# draw for the initial state
      Range ran(1,count_out*j);
      for(int i=0;i<count_out*j;i++) X(i,nvar) = 0;
      for(int i=0;i<(count_out-1)*j;i++) X(s0+i-1,nvar) = 1;  //# assign drawn initial state to all observations until the first brand is chosen
    }
  }
  //######################################################################################################

  NumericMatrix Xbeta = mmult(X,mat_beta);
  //Xbeta=matrix(Xbeta,byrow=T,ncol=j)
  NumericMatrix ind(2,n);
  ind(_,1)=y;
  for(int i=0;i<n;i++) ind(i,0)=i;
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


