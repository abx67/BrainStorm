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
  return(counter+1);
};

// [[Rcpp::export]]

double bb(vec const& beta,mat initial_price_state_dropped){
  mat xy = initial_price_state_dropped*beta;

  return xy[0];
}
