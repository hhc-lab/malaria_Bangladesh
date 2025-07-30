// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <math.h>
#include <RcppEigen.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
List malaria_ode_cpp(double t, NumericVector x, List params) {
  
  // define parameters
  double numpatch = params["numpatch"];
  double c = params["c"];
  double b = params["b"];
  double a = params["a"];
  double mu = params["mu"];
  double tau = params["tau"];
  double r = params["r"];
  std::vector<double> m = params["m"];
  std::vector<double> H = params["H"];
  NumericMatrix pij = params["pij"];
  
  // define model compartments
  double *X = &x[0];
  //double *C = &x[numpatch];
  std::vector<double> k(numpatch);
  
  for(int i = 0; i < numpatch; i++){
    double temp1 = 0;
    double temp2 = 0;
    
    for(int j = 0; j < numpatch; j++){
      temp1 += X[j] * H[j] * pij(j,i);
      temp2 += H[j] * pij(j,i);
    }
    
    k[i] = temp1 / temp2;
  }
  
  //std::vector<double> dX(numpatch, 0);
  //std::vector<double> dC(numpatch, 0);
  std::vector<double> rtnVec(numpatch*2, 0);
  
  for(int i = 0; i < numpatch; i++) {
    
    for(int j = 0; j < numpatch; j++){
      //dC[i] += pij(i,j) * m[j] * pow(a,2) * b * c * exp(-mu*tau) * k[j] / (a*c*k[j] +mu);
      rtnVec[i+numpatch-1] += pij(i,j) * m[j] * pow(a,2) * b * c * exp(-mu*tau) * k[j] / (a*c*k[j] +mu);
    }
    rtnVec[i+numpatch-1] = rtnVec[i+numpatch-1] * (1 - X[i]);
    rtnVec[i] = rtnVec[i+numpatch-1] - r*X[i];
    //dC[i] = dC[i] * (1 - X[i]);
    //dX[i] = dC[i] - r*X[i];
  }
  
  return List::create(rtnVec);
}
