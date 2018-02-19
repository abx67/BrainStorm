Hi Fan, 

Nice meeting you today!

Here are the files: 

(1) 9_1_llmnl_initial_mcprice.R ！ new likelihood function file, to be written in Rcpp 
(2) 9_1_corrLike.R ！ old overall function, just for your reference
(3) rhierMnlRwMixture_rcpp_loop.cpp ！ new c++ functions that need to be adjusted (3 new functions created, one based on  9_1_llmnl_initial_mcprice.R, and two that just call these functions with all the inputs)
(4) rhierMnlRwMixture_rcpp.R ！ new overall function that calls functions defined by (3).  

Also, here is a link to the bayesm package: https://cran.r-project.org/web/packages/bayesm/index.html
GitHub repo for this package: https://github.com/cran/bayesm

Just to recap, the task is to create three new Rcpp files that are in (3), using function defined in (1) instead of llmnl_con function there, and make sure that these new functions properly use objects defined in the bayesm package. 

Feel free to shoot me an email if you have some questions! 

################################################################
sourceCpp("optimal2.cpp")
Rcpp::Environment global = Rcpp::Environment::global_env();
Rcpp::Function mutualinfo = global["mutualinfo"];

R CMD BATCH --args xxx.R
cd ~/'iCloud Drive (Archive)'/Documents/GitHub/Rcpp_mcmc
################################################################

#9_1_corrLike.R
rhierMnlRwMixture_inLike_X0_mcprice=
  function(Data,Prior,Mcmc)
  {

#9_1_llmnl_initial_mcprice.R

llmnl_initial_mcprice= 
	function(beta,y,X,count_out,initial_price_id,initial_price_state_dropped,s1,transition_matrix_median_steps,price_transition_states,vec_price_states_probs,draws_length,number_price_simulations,flag_markovian, flag_know_state) 

#rhierMnlRwMixture_rcpp_loop.cpp
double llmnl_con(vec const& betastar, vec const& y, mat const& X, vec const& SignRes = NumericVector::create(0)){

List rhierMnlRwMixture_rcpp_loop(List const& lgtdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  double nu, mat const& V, double s,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta,  vec const& a, vec oldprob, mat oldbetas, vec ind, vec const& SignRes){