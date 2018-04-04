#include <RcppArmadillo.h>
#include <cmath>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
// [[Rcpp::depends(RcppArmadillo)]]

using std::log;
using std::exp;
using std::max;
using std::min;
using std::abs;
using std::sqrt;
using std::pow;
using namespace Rcpp; 

void R_init_markovchain(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);	
}

extern void dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda,
                   double *w, double *work, int *lwork, int *info );

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov){
    int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - center;
    }
    return sum((x_cen * cov.i()) % x_cen, 1);
}

// [[Rcpp::export]]
arma::vec ei(arma::mat M) {
    return arma::eig_sym(M);
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma (arma::mat x,  arma::rowvec mean,  arma::mat sigma){
    arma::vec distval = Mahalanobis(x, mean, sigma);
    arma::vec ei_sigma = ei(sigma);
    double logdet = sum(arma::log(ei_sigma));
    arma::vec logretval = -((x.n_cols * log2pi + logdet + distval)/2) ;
    return(logretval);
}

// [[Rcpp::export]]
arma::mat rmvnorm_arma(int n,
                  const arma::vec& mu,
                  const arma::mat& Sigma) {
    unsigned int p = Sigma.n_cols;
    Rcpp::NumericVector draw = Rcpp::rnorm(n*p);
    arma::mat Z = arma::mat(draw.begin(), n, p, false, true);
    arma::mat Y = arma::repmat(mu, 1, n).t()+Z * arma::chol(Sigma);
    return Y;
}

// **********************************************************//
//                  Transpose                                //
// **********************************************************//
// [[Rcpp::export]]
arma::mat transpose (arma::mat x) {
	return x.t();
}

// **********************************************************//
//              Likelihood evaluation of Timepart            //
// **********************************************************//
// [[Rcpp::export]]
double Timepartsum (NumericMatrix mumat, double sigma_tau,
                   IntegerVector senders, NumericVector timestamps){
    int D = senders.size();
    double timesum = 0;
    for (unsigned int d = 0; d < D; d++) {
        int a_d = senders[d] - 1;
        timesum += R::dlnorm(timestamps[d], mumat(d, a_d), sigma_tau, TRUE);
        for (unsigned int i = 0; i < mumat.ncol(); i++) {
        if (i != a_d) {
            timesum += R::plnorm(timestamps[d], mumat(d, i), sigma_tau, FALSE, TRUE);
    	  }
        }
    }
    return timesum;
}

// **********************************************************//
//                    Call rmultinom from R                  //
// **********************************************************//
// [[Rcpp::export]]
IntegerVector callRMultinom (NumericVector x) {
    NumericVector x1 = x / sum(x);
    int n = x1.size();
    IntegerVector d(n);
    R::rmultinom(1, x1.begin(), n, d.begin());
    return d;
}

// **********************************************************//
//                    Multinomial Sampler                    //
// **********************************************************//
// [[Rcpp::export]]
int multinom_vec (NumericVector x) {
	NumericVector x1 = x / sum(x);
	int n = x1.size();
    IntegerVector d(n);
    R::rmultinom(1, x1.begin(), n, d.begin());
	for (unsigned int j = 0; j < n; j++) {
		if (d[j] == 1) {
			return j+1;
		}
	}
	return 0;
}

// **********************************************************//
//               Normalizing constant of Gibbs               //
// **********************************************************//
// [[Rcpp::export]]
double normalizer (arma::vec lambda_da) {
	return log(arma::prod(lambda_da+1)-1);
}

// **********************************************************//
//                   Posterior for Edgepart                  //
// **********************************************************//
// [[Rcpp::export]]
double Edgepartsum (List lambda, List u) {
	double out = 0;
	int D = lambda.size();
	IntegerMatrix samp = u[0];
	int A = samp.nrow();
	for (unsigned int d = 0; d < D; d++) {
		NumericMatrix lambda_d = lambda[d];
		NumericMatrix u_d = u[d];
		for (unsigned int a = 0; a < A; a++) {
			double num = sum(lambda_d(a, _) * u_d(a, _));
			NumericVector lambda_da = exp(lambda_d(a, _));
			lambda_da[a] = 0;
			double denom = normalizer(lambda_da);
			out += num - denom;
		}
	}
	return out;
}

// **********************************************************//
//                     Calculate mu matrix                   //
// **********************************************************//
// [[Rcpp::export]]
arma::mat mu_cpp (arma::cube Y, arma::rowvec eta) {
	int D = Y.n_rows;
	int A = Y.n_cols;
	int Q = Y.n_slices;
	arma::mat mu = arma::zeros(D, A);
	for (unsigned int d = 0; d < D; d++) {
		for (unsigned int a = 0; a < A; a++) {
			arma::vec Y_da = Y.subcube(d, a, 0, d, a, Q-1);
			mu(d, a) = sum(eta * Y_da);
		}
	}
	return mu;
}

// **********************************************************//
//                     Calculate lambda list                 //
// **********************************************************//
// [[Rcpp::export]]
arma::mat lambda_cpp (arma::cube X, arma::rowvec beta, double delta) {
	int A = X.n_rows;
	int P = X.n_slices;
	arma::mat lambda_d = arma::zeros(A, A);
	for (unsigned int a = 0; a < A; a++) {
		for (unsigned int r = 0; r < A; r++) {
			if (r != a) {
				arma::vec X_dar = X.subcube(a, r, 0, a, r, P-1);
				lambda_d(a,r) = delta + sum(beta * X_dar);
			}
		}
	}
	return lambda_d;
}
