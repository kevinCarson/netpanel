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
arma::mat netfilter(Rcpp::List W, arma::vec rho){
  // extracting the size of the first block-diagonal matrix
  arma::mat NT = Rcpp::as<arma::mat>(W[0]); 
  arma::mat A = eye(NT.n_rows,NT.n_rows); // the NT x NT identity matrix
  int nNets = W.length();
  for(int i = 0; i < nNets; i ++){
    A = A - rho[i]*Rcpp::as<arma::mat>(W[i]); 
  }
  return A;
}

// [[Rcpp::export]]
arma::vec gls_compute_v(arma::mat A,
              arma::vec beta,
              arma::mat X,
              arma::vec Y){
  
  arma::vec resid = A*Y - X*beta; // the classical residual formula
  return resid;
}


// [[Rcpp::export]]
arma::vec gls_compute_beta(arma::mat A,
                       arma::mat omegainv,
                       arma::mat X,
                       arma::vec Y,
                       bool random){
  arma::vec beta; 
  if(random){
    beta = solve(X.t()*omegainv*X, X.t()*omegainv*A*Y);
  }else{
    beta = solve(X.t()*X, X.t()*A*Y);
  }
  return beta;
  
}

// [[Rcpp::export]]
double gls_compute_sigma2(arma::vec v,
                          arma::mat omegainv,
                          bool random){
  double sigma;
  if(random){
   sigma = arma::as_scalar(v.t() * omegainv * v);
    sigma = sigma/v.size();
  }else{
     sigma = arma::as_scalar(v.t() * v);
    sigma = sigma/v.size();
  }
  return sigma;
}












// the following functions are to be used with random unit level effects


// [[Rcpp::export]]
double pnam_ll_random_find(arma::vec Y,
                      arma::mat X,
                      Rcpp::List W,
                      double n,
                      double t,
                      arma::vec rho,
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
 
 arma::vec beta = gls_compute_beta(netA,invomega,X,Y,true);
 arma::vec v = gls_compute_v(netA, beta,X,Y); 
 double sigma2 = gls_compute_sigma2(v,invomega,true); 
 
 double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + (-0.5)*logdet_val + 
                             jacobian - (1/(2*sigma2))*(v.t()*invomega*v));
   
 return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_random(arma::mat netA,
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







// the following functions are to be used with fixed effects panel network autocorrelation models

// [[Rcpp::export]]
double pnam_ll_fixed_find(arma::vec Y,
                           arma::mat X,
                           Rcpp::List W,
                           double n,
                           double t,
                           arma::vec rho){
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi);
  arma::mat netA = netfilter(W,rho);
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  arma::vec beta = gls_compute_beta(netA,netA,X,Y,false);
  arma::vec v = gls_compute_v(netA,beta,X,Y); 
  double sigma2 = gls_compute_sigma2(v,netA,false); 
  double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) +  jacobian - (1/(2*sigma2))*(v.t()*v));
  return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_fixed(arma::mat netA,
                      double sigma2, 
                      arma::vec v,
                      double n,
                      double t){
  
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi) + -n*t/2*log(2*sigma2);
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  double ll = arma::as_scalar(c + jacobian - (1/(2*sigma2))*(v.t()*v));
  return ll;
}



// [[Rcpp::export]]
arma::mat create_B(arma::mat W2, 
                   double lambda){
  arma::mat B = eye(W2.n_rows,W2.n_rows) - lambda*W2; 
  arma::mat Bi = B.i(); 
  return Bi;
}

// [[Rcpp::export]]
arma::mat create_LAMBDA(arma::mat B){
  // per the formula: LAMBDA: (Bt)^(-1)^B(-1)
  arma::mat LAMBDA = B.t().i()*B.i(); 
  return LAMBDA;
}


// [[Rcpp::export]]
arma::vec gls_compute_e(arma::mat A,
                        arma::vec beta,
                        arma::mat B,
                        arma::mat X,
                        arma::vec Y){
  arma::vec resid = B.i()*(A*Y - X*beta); // the classical residual formula
  return resid;
}


// [[Rcpp::export]]
arma::vec gls_compute_beta_error(arma::mat A,
                                 arma::mat LAMBDA,
                                 arma::mat X,
                                 arma::vec Y){
  // for fixed/dynamic models, the lambda will be different 
  arma::vec beta =  solve(X.t()*LAMBDA*X, X.t()*LAMBDA*A*Y); 
  return beta;
  
}

// [[Rcpp::export]]
double gls_compute_sigma2_error(arma::vec v,
                                arma::mat omegainv){
  // for dynamic and fixed models, we will need to change omega to the different value
  double sigma;
  sigma = arma::as_scalar(v.t() * omegainv * v);
  sigma = sigma/v.size();
  return sigma;
}


// [[Rcpp::export]]
double pnam_ll_fixed_error_find(arma::vec Y,
                                arma::mat X,
                                Rcpp::List W,
                                double n,
                                double t,
                                arma::vec rho,
                                arma::mat W2,
                                double lambda){
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi);
  arma::mat netA = netfilter(W,rho);
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  arma::mat B = create_B(W2,lambda);
  double jacobianerror;
  double signjacoberror;
  arma::log_det(jacobianerror, signjacoberror, B.i());
  arma::mat LAMBDA = create_LAMBDA(B);
  arma::vec beta = gls_compute_beta_error(netA,LAMBDA,X,Y);
  arma::vec e = gls_compute_e(netA,beta,B,X,Y); 
  double sigma2 = gls_compute_sigma2_error(e,LAMBDA); 
  double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + jacobianerror + jacobian - (1/(2*sigma2))*(e.t()*e));
  return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_fixed_error(arma::mat netA,
                           arma::mat B,
                           double sigma2, 
                           arma::vec e,
                           double n,
                           double t){
  
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi) + -n*t/2*log(2*sigma2);
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  double jacobianerror;
  double signjacoberror;
  arma::log_det(jacobianerror, signjacoberror, B.i());
  double ll = arma::as_scalar(c + jacobian + jacobianerror - (1/(2*sigma2))*(e.t()*e));
  return ll;
}



// [[Rcpp::export]]
arma::mat omega_error(double theta,
                        double n,
                        double t,
                        arma::mat W2,
                        double lambda){
  
  arma::vec l = arma::ones(t); // the vector of ones
  arma::mat In = eye(n,n); // the standard N x N identity matrix
  arma::mat It = eye(t,t); // the standard T x T identity matrix
  arma::mat llt = l * l.t(); // per the formula
  arma::mat B = eye(W2.n_rows,W2.n_rows) - lambda*W2; 
  arma::mat Binv = B.i(); 
  arma::mat BBt= Binv * Binv.t(); 
  // the omega matrix per the formula
  arma::mat omega = theta*(arma::kron(llt,In)) + BBt; 
  return omega;
}


// [[Rcpp::export]]
double pnam_ll_random_error_find(arma::vec Y,
                           arma::mat X,
                           Rcpp::List W,
                           double n,
                           double t,
                           arma::vec rho,
                           double theta,
                           arma::mat W2,
                           double lambda){
  double pi = 2 * acos(0.0); 
  double c = -n*t/2*log(2*pi);
  arma::mat omegaest = omega_error(theta,n,t,W2,lambda);
  double logdet_val;
  double sign;
  arma::log_det(logdet_val, sign, omegaest);
  arma::mat netA = netfilter(W,rho);
  double jacobian;
  double signjacob;
  arma::log_det(jacobian, signjacob, netA);
  arma::mat invomega = omegaest.i(); 
  arma::vec beta = gls_compute_beta(netA,invomega,X,Y,true);
  arma::vec v = gls_compute_v(netA, beta,X,Y); 
  double sigma2 = gls_compute_sigma2(v,invomega,true); 
  double ll = arma::as_scalar(c -(n*t/2)*log(sigma2) + (-0.5)*logdet_val + 
                              jacobian - (1/(2*sigma2))*(v.t()*invomega*v));
  return -ll; 
}


// [[Rcpp::export]]
double pnam_ll_random_error(arma::mat netA,
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

