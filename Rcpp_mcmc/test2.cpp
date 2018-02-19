#include <bayesm.h>
using namespace std;

int my_rmultinomF(vec p) {
  vec r(1,fill::randu);
  vec cp = cumsum(p);
  int counter=0;
  for(int i=0;i<cp.size();i++){
    if(r[0]>cp[i]){
      counter += 1;
    }
  }
  if(counter<100&&counter>10){counter=20;}
  else if(counter <=10){counter=5;}
  return(counter+1);
};

// [[Rcpp::export]]
double llmnl_initial_mcprice(vec const& beta,vec const& y,mat X,
                             int count_out,int initial_price_id,mat initial_price_state_dropped,
                             int s1,mat transition_matrix_median_steps,mat price_transition_states,
                             vec const& vec_price_states_probs,int draws_length,int number_price_simulations,
                             bool flag_markovian,bool flag_know_state) {
  
  

  
  
  if (flag_know_state) {
    //xb_init = matrix(exp(initial_price_state_dropped%*%beta), ncol = j, byrow = T) //later reshape
    mat xy(100,100);
    xy = initial_price_state_dropped*beta;
    //cout<<xy<<endl;
    return xy(0,0);
  }
  else{return 0;}
}
