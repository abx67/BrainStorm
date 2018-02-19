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
  // if(counter<100&&counter>10){counter=20;}
  // else if(counter <=10){counter=5;}
  return(counter+1);
};

// [[Rcpp::export]]
double llmnl_initial_mcprice(vec const& beta,vec const& y,mat X,
                                    int count_out,int initial_price_id,mat initial_price_state_dropped,
                                    int s1,mat transition_matrix_median_steps,mat price_transition_states,
                                    vec const& vec_price_states_probs,int draws_length,int number_price_simulations,
                                    bool flag_markovian,bool flag_know_state) {
  
  
  //main
  int n=y.size(); // number of choices
  int j=X.n_rows/n; // number of brands
  int nvar = X.n_cols;
  
  vec Vec(j-1);  //in .r file it was called vec
  
  // added to llmnl function 
  if (count_out != 0) { // if the initial state is not known
    
    // New August 2017 -- Markovian prices ### 
    if (flag_markovian) {
      mat X_forward(j,j+1,fill::zeros);
      for(int i=0;i<j-1;i++){X_forward(i,i)=1;}
      mat dif_probs(number_price_simulations,j,fill::zeros);
      for(int j_markov=0;j_markov<number_price_simulations;j_markov++){
        // draw a price vector
        vec price_draw_this_iter(draws_length,fill::zeros);
        mat transition_matrix_median_steps_t = transition_matrix_median_steps.t();
        price_draw_this_iter[draws_length-1] = my_rmultinomF(transition_matrix_median_steps_t.col(initial_price_id-1));
        for (int i=1;i<draws_length;i++){
          price_draw_this_iter[draws_length-i-1] = my_rmultinomF(transition_matrix_median_steps_t.col(price_draw_this_iter[i]));   //check later
        }
        // price_draw_this_iter <- price_draw_this_iter[draws_length:1] //no need anymore
        //# initial random state
        ///////return price_draw_this_iter[0];
        vec temp(j-1);
        temp.fill(1.0/(j-1));
        X_forward(my_rmultinomF(temp)-1,nvar-1) = 1;    //nvar and j+1
        X_forward.col(j-1) = price_transition_states( span(0, j-1), price_draw_this_iter[0] );
        //# forward simulate the choices
        vec prob_forward(j);

        for (int i=1;i<draws_length;i++){
          vec exp_xbeta_forward(j);
          exp_xbeta_forward = exp(X_forward*beta);
          prob_forward = exp_xbeta_forward / sum(exp_xbeta_forward);
          X_forward(my_rmultinomF(prob_forward)-1,nvar-1) = 1;
          //error notes:incompatible matrix dimensions: 10x1 and 1x10
          mat X_forward_tmp = price_transition_states(price_draw_this_iter[i]-1,span(0, j-1));
          X_forward.col(j-1) = X_forward_tmp.t();
        }
        dif_probs.row(j_markov) = prob_forward.t();
        //return dif_probs;
      }
      vec tmp_x(j);
      for(int i=0;i<j;i++){    //correspond to apply(,mean)
        tmp_x[i]=mean(dif_probs.col(i));
      }
      //# marginal initial state probability vector
      Vec = tmp_x(span(0,j-1-1))/sum(tmp_x(span(0,j-1-1)));

    } else {
      //####### Older iid price case #########
      mat xbeta(price_transition_states.n_cols,price_transition_states.n_rows);
      vec tmp_beta = beta(span(0,j-1));  tmp_beta[j-1] = 0;
      mat tmpmat = price_transition_states.cols(0,j-1) * beta[j];
      tmpmat = tmpmat.t();
      for(int i=0;i<j;i++){xbeta.row(i) = tmp_beta[i]+tmpmat.row(i);}
      xbeta = xbeta.t();
      mat Prob0(j-1,j-1,fill::zeros);
      for(int j_state=0;j_state<j-1;j_state++) {
        mat xbeta_state=xbeta;
        xbeta_state.col(j_state) = xbeta.col(j_state) + beta[nvar];
        mat xb0(xbeta_state.n_rows,xbeta_state.n_cols);
        xb0 = exp(xbeta_state);
        mat prob0(xb0.n_rows,xb0.n_cols);  //case matter??
        for(int i=0;i<xb0.n_rows;i++){prob0.row(i) = xb0.row(i)/sum(xb0.row(i));}
        //prob0: same size xb0 or xbeta_state or xbeta or t(price_transition_states); 
        mat tmp_prob = vec_price_states_probs.t() * prob0;    //careful!!
        
        Prob0.col(j_state) =  tmp_prob(0,span(0,j-1-1)).t(); //# transposed transition matrix

        Prob0(j_state,j_state) = Prob0(j_state,j_state) + tmp_prob(0,j-1); //# if choose outside option, state stays the same
      }

      cx_vec eigval;
      cx_mat eigvec;

      eig_gen(eigval, eigvec, Prob0);
      
      Vec = real(eigvec.col(0));
      Vec = Vec / sum(Vec);
      //# these are marginal probabilites of chosings states for average prices
    }
    if (!flag_know_state) {// # if we are in the case when the data was not dropped and initial state is not known
      int s0 = my_rmultinomF(Vec); //# draw for the initial state

      for(int i=0;i<count_out*j;i++) {

        X(i,nvar-1) = 0;     //due to rewrite X, cannot input type const& X

      }

      for(int i=0;i<(count_out-1)*j;i++) {
        X(s0+i-1,nvar-1) = 1;  //# assign drawn initial state to all observations until the first brand is chosen
      }
    }
  }

  //######################################################################################################
  mat Xbeta_tmp = X*beta;
  Xbeta_tmp.reshape(j,Xbeta_tmp.n_cols*Xbeta_tmp.n_rows/j);
  mat Xbeta = Xbeta_tmp.t();    //size: n*j;  row=n??
  //Xbeta=matrix(Xbeta,byrow=T,ncol=j)

  //ind=cbind(c(1:n),y)
  vec xby(n);
  for(int i=0;i<n;i++){xby[i] = Xbeta(i,y[i]-1);}
  Xbeta=exp(Xbeta);
  vec iota(j,fill::ones);
  vec denom=log(Xbeta*iota); //n*1

    if (flag_know_state) {
      mat xb_init_tmp(initial_price_state_dropped.n_rows,1);
      xb_init_tmp = exp(initial_price_state_dropped*beta); //matrix of expxbeta by each possible (m-1) state
      //xb_init = matrix(exp(initial_price_state_dropped%*%beta), ncol = j, byrow = T) //later reshape
      xb_init_tmp.reshape(j,xb_init_tmp.n_cols*xb_init_tmp.n_rows/j);
      mat xb_init = xb_init_tmp.t();
      mat prob_init(xb_init.n_rows,xb_init.n_cols); // matrix of choice probabilities
      for(int i=0;i<xb_init.n_rows;i++){prob_init.row(i) = xb_init.row(i)/sum(xb_init.row(i));}
      mat vec_with_init_tmp(1,prob_init.n_cols);
      vec_with_init_tmp = Vec.t()*prob_init; //marginal probabilities
      vec vec_with_init = vec_with_init_tmp.t();
      return(sum(xby-denom) + log(vec_with_init[s1-1]));
    } else {
      return(sum(xby-denom));
    }
}
