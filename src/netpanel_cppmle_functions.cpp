#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <numbers>
using namespace arma;

// [[Rcpp::export]]
arma::mat omega(double theta,
          double n,
          double t){
  
  arma::vec l = arma::ones(t); // the vector of ones
  arma::mat In = eye(n,n); // the standard N x N identity matrix
  arma::mat It = eye(t,t); // the standard T x T identity matrix
  arma::mat llt = l * l.t(); // per the formula
  // the omega matrix per the formula
  arma::mat omega = theta*(arma::kron(llt,In)) + arma::kron(It,In); 
  return omega;
}


// [[Rcpp::export]]
arma::mat netfilter(arma::mat W, double rho){
  arma::mat A = eye(W.n_rows,W.n_rows) - rho*W; 
  return A;
}

// [[Rcpp::export]]
arma::vec compute_v(arma::mat A,
              arma::vec beta,
              arma::mat X,
              arma::vec Y){
  
  arma::vec resid = A*Y - X*beta; // the classical residual formula
  return resid;
}


// [[Rcpp::export]]
arma::vec compute_beta( arma::mat A,
                  arma::mat omegainv,
                  arma::mat X,
                  arma::vec Y){
  
  arma::vec beta = solve(X.t()*omegainv*X, X.t()*omegainv*A*Y);
  return beta;
}

// [[Rcpp::export]]
double compute_sigma2(arma::vec v,
                      arma::mat omegainv){
  double sigma = arma::as_scalar(v.t() * omegainv * v);
  sigma = sigma/v.size();
  return sigma;
}


// [[Rcpp::export]]
double pnam_ll_rule_search(arma::vec Y,
                      arma::mat X,
                      arma::mat W,
                      double n,
                      double t,
                      double rho,
                      double theta){
 double pi = 2 * acos(0.0); 
 double c = -n*t/2*log(2*pi);
 arma::mat omegaest = omega(theta,n,t);
 
 double logdet_val;
 double sign;
 arma::log_det(logdet_val, sign, omegaest);
 
 arma::mat netA = netfilter(W,rho);
 
 double jacobian;
 double signjacob;
 arma::log_det(jacobian, signjacob, netA);

 arma::mat invomega = omegaest.i(); 
 
 arma::vec beta = compute_beta(netA,invomega,X,Y);
 arma::vec v = compute_v(netA, beta,X,Y); 
 double sigma2 = compute_sigma2(v,invomega); 
 
 double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + (-0.5)*logdet_val + 
                             jacobian - (1/(2*sigma2))*(v.t()*invomega*v));
   
 return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_rule(arma::mat netA,
                    double sigma2, 
                    arma::vec v,
                    arma::mat omega,
                    double n,
                    double t){

  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi) + -n*t/2*log(2*sigma2);
  
  arma::mat invomega = omega.i();
  
  double logdet_val;
  double sign;
  arma::log_det(logdet_val, sign, omega);

  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  
  double ll = arma::as_scalar(c + (-0.5)*logdet_val + jacobian - (1/(2*sigma2))*(v.t()*invomega*v));
  return ll;
}



// the following functions are to be used with random unit level effects and correlated errors


// [[Rcpp::export]]
arma::mat omega_rulewae(double theta,
                  double n,
                  double t,
                  arma::mat W1,
                  double lambda){
  
  arma::vec l = arma::ones(t); // the vector of ones
  arma::mat In = eye(n,n); // the standard N x N identity matrix
  arma::mat It = eye(t,t); // the standard T x T identity matrix
  arma::mat llt = l * l.t(); // per the formula
  arma::mat B = eye(W1.n_rows,W1.n_rows) - lambda*W1; 
  arma::mat Binv = B.i(); 
  arma::mat BBt= Binv * Binv.t(); 
  // the omega matrix per the formula
  arma::mat omega = theta*(arma::kron(llt,In)) + BBt; 
  return omega;
}



// [[Rcpp::export]]
double pnam_ll_rulewae_search(arma::vec Y,
                              arma::mat X,
                              arma::mat W,
                              arma::mat W1,
                              double n,
                              double t,
                              double rho,
                              double theta,
                              double lambda){
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi);
  arma::mat omegaest = omega_rulewae(theta,n,t,W1,lambda);
  
  double logdet_val;
  double sign;
  arma::log_det(logdet_val, sign, omegaest);
  
  arma::mat netA = netfilter(W,rho);
  
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);

  arma::mat invomega = omegaest.i(); 
  
  arma::vec beta = compute_beta(netA,invomega,X,Y);
  arma::vec v = compute_v(netA, beta,X,Y); 
  double sigma2 = compute_sigma2(v,invomega); 
  
  double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + (-0.5)*logdet_val + 
                            jacobian- (1/(2*sigma2))*(v.t()*invomega*v));
  

  return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_rulewae(arma::mat netA,
                       double sigma2, 
                       arma::vec v,
                       arma::mat omega,
                       double n,
                       double t){
  
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi) + -n*t/2*log(2*sigma2);
  
  arma::mat invomega = omega.i();
  
  double logdet_val;
  double sign;
  arma::log_det(logdet_val, sign, omega);
  
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  
  double ll = arma::as_scalar(c + (-0.5)*logdet_val + jacobian  - (1/(2*sigma2))*(v.t()*invomega*v));
  return ll;
}


// the following functions are to be used with dynamic network autocorrelation models

// [[Rcpp::export]]
arma::mat netfilter_dpnam(arma::mat W, 
                    arma::mat L,
                    double rho,
                    double delta,
                    double phi){
  arma::mat A = eye(W.n_rows,W.n_rows) - rho*W - delta*L - phi*L*W; 
  return A;
}

// [[Rcpp::export]]
arma::vec compute_beta_norandom(arma::mat X, 
                          arma::mat Y,
                          arma::mat A){
  
  arma::vec beta = solve(X.t()*X, X.t()*A*Y);
  return beta;
}

// [[Rcpp::export]]
double compute_sigma2_norandom(arma::vec v){
  double sigma = arma::as_scalar(v.t() * v);
  sigma = sigma/v.size();
  return sigma;
}


// [[Rcpp::export]]
double pnam_ll_dynam_search(arma::vec Y,
                            arma::mat X,
                            arma::mat W,
                            arma::mat L,
                            double n,
                            double t,
                            double rho,
                            double delta,
                            double phi){
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi);
  
  arma::mat netA = netfilter_dpnam(W,L,rho,delta,phi);
  
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  
  
  arma::vec beta = compute_beta_norandom(X,Y,netA);
  arma::vec v = netA*Y - X*beta; 
  double sigma2 = compute_sigma2_norandom(v); 
  double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + jacobian - (1/(2*sigma2))*(v.t()*v));
  return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_dynam(arma::mat netA,
                     double sigma2, 
                     arma::vec v,
                     double n,
                     double t){
  
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi) + -n*t/2*log(sigma2);
  
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  
  double ll = arma::as_scalar(c + jacobian - (1/(2*sigma2))*(v.t()*v));
  return ll;
}


// the following functions are to be used with dynamic network autocorrelation models with autocorrelated errors



















// [[Rcpp::export]]
arma::mat B_dpnam(arma::mat B, 
                 double lambda){
  arma::mat Ba = eye(B.n_rows,B.n_rows) - lambda*B; 
  arma::mat Bi = Ba.i(); 
  return Bi;
}



// [[Rcpp::export]]
arma::vec compute_beta_norandom_wae(arma::mat X, 
                              arma::mat Y,
                              arma::mat A,
                              arma::mat Lambda){
  
  arma::vec beta = solve(X.t()*Lambda*X, X.t()*Lambda*A*Y);
  return beta;
}

// [[Rcpp::export]]
double compute_sigma2_norandom_wae(arma::vec v,
                               arma::mat Lambda){
  double sigma = arma::as_scalar(v.t() *Lambda* v);
  sigma = sigma/v.size();
  return sigma;
}


// [[Rcpp::export]]
double pnam_ll_dynam_wae_search(arma::vec Y,
                            arma::mat X,
                            arma::mat W,
                            arma::mat L,
                            arma::mat B,
                            double n,
                            double t,
                            double rho,
                            double delta,
                            double phi,
                            double lambda){
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi);
  
  arma::mat netA = netfilter_dpnam(W,L,rho,delta,phi);
  arma::mat netB = B_dpnam(B, lambda);
  
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  
  double jacobianB;
  double signjacobB;
  arma::log_det(jacobianB, signjacobB, netB.i());
  arma::mat Lambda = netB.t()*netB; 
  
  arma::vec beta = compute_beta_norandom_wae(X,Y,netA, Lambda);
  arma::vec e = netB.i()*(netA*Y - X*beta); 
  double sigma2 = compute_sigma2_norandom_wae(e,Lambda); 

  double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + jacobian + jacobianB+ - (1/(2*sigma2))*(e.t()*e));
  return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_dynam_wae(arma::mat netA,
                     arma::mat netB,
                     double sigma2, 
                     arma::vec e,
                     double n,
                     double t){
  
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi) + -n*t/2*log(sigma2);
  
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  
  double jacobianB;
  double signjacobB;
  arma::log_det(jacobianB, signjacobB, netB.i());
  arma::mat Lambda = netB.t()*netB; 
  
  double ll = arma::as_scalar(c + jacobian + jacobianB - (1/(2*sigma2))*(e.t()*e));
  return ll;
}

