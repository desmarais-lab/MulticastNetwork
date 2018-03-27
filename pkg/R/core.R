#' @useDynLib MulticastNetwork
#' @import stats
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom reshape melt
#' @importFrom coda mcmc geweke.diag
#' @importFrom combinat permn
#' @importFrom mgcv uniquecombs
#' @importFrom lubridate wday hour
#' @importFrom LaplacesDemon dhalfcauchy rhalfcauchy
#' @importFrom truncnorm rtruncnorm dtruncnorm
#' @importFrom MCMCpack rinvgamma dinvgamma

tvapply = function(...) transpose(vapply(...))

#' @title gibbs.measure.support
#' @description List out the support of Gibbs measure
#'
#' @param n length of the vector to be sampled 
#'
#' @return  a 2^n x n binary matrix representing the support of the binary Gibbs measure in n elements
#'
#' @export
gibbs.measure.support = function(n) {
	gibbs.support = rbind(rep(1, n))
	for(i in 1:(n-1)){
		gibbs.mat.i = do.call('rbind',permn(c(rep(1, i), rep(0, n-i))))
		gibbs.support = rbind(gibbs.support, gibbs.mat.i)
	}
	out = as.matrix(unique(gibbs.support))
	return(out)
}

#' @title r.gibbs.measure
#' @description List out the support of Gibbs measure
#'
#' @param lambda.i a vector of coefficients according to which each element
#' @param delta a positive real valued parameter that controls the density penalty 
#' @param support support of Gibbs measure
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
r.gibbs.measure <- function(lambda.i, delta, support) {
	#gibbsNormalizer = prod(exp(delta+lambda.i)+1)-1
	logitNumerator = vapply(1:nrow(support), function(s) {
		sum((delta+lambda.i)*support[s,])
		}, c(1))		
	samp = lmultinom(logitNumerator)
	return(support[samp,])	
}

#' @title adaptive.MH 
#' @description adaptive Metropolis Hastings to maintain target acceptance rate
#'
#' @param sigma.Q proposal distribution variance parameter
#' @param accept.rates acceptance rate from previous iteration
#' @param target target acceptance rate
#' @param update.size size of update to be adjusted 
#' @param tol tolerance level to determine if acceptance rate is too high
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
adaptive.MH = function(sigma.Q, accept.rates, target = 0.1, update.size, tol = 0.8) {
	for (i in 1:length(sigma.Q)) {
	  if (accept.rates[i] < target) {
			sigma.Q[i] = sigma.Q[i]-update.size[i]
		}
		if (accept.rates[i] > (target+tol)) {
			sigma.Q[i] = sigma.Q[i]+update.size[i]
		}		
	}
	return(sigma.Q)
}

#' @title Inference
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param sigma.Q proposal distribution variance parameter
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param initial list of initial values user wants to assign including (delta, b, eta, cd, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

Inference = function(edge, node, sigma.Q, prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner,
                          netstat, timestat, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
  # trim the edge so that we only model edges after 384 hours
  A = length(node)
  D = length(edge)
  timestamps = vapply(edge, function(d) { d[[3]] }, c(1))
  senders = vapply(edge, function(d) { d[[1]] }, c(1))
  edge.trim = which_num(384*timeunit, timestamps-timestamps[1]):D
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])/timeunit
  timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
  # initialization
  convergence = c()
  netstat = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = D, ncol = sum(timestat))
  if (sum(timestat) > 0) {
    Sys.setenv(TZ = tz)
    time_ymd = as.POSIXct(timestamps, tz = getOption("tz"), origin = "1970-01-01")
    if (timestat[1] > 0) {
      days = vapply(time_ymd, function(d) {wday(d)}, c(1))
      days[days==1] = 8
      timemat[,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      it = 1
    }
    if (timestat[2] > 0) {
      hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
      timemat[,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
    }     
  }
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  Q = length(prior.eta[[1]])
  proposal.var1 = diag(P)
  proposal.var2 = diag(Q)
  if (length(initial) == 0) {
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rinvgamma(1, prior.tau[1], prior.tau[2])
    b.old = rmvnorm_arma(1, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(1, prior.eta[[1]], prior.eta[[2]])
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow = A, ncol = A)
      diag(u[[d]]) = 0
      u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
  } else {
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = matrix(NA, nrow = P, ncol = Inner[1])
  etamat = matrix(NA, nrow = Q, ncol = Inner[2])
  deltamat = rep(NA, Inner[1])
  sigma_taumat = rep(NA, Inner[3])
  mu = matrix(0, nrow = D, ncol = A)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in 1:D) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)
      }
  }
  timeinterval = timefinder(timestamps, edge.trim, timeunit)
  X = list()
  for (d in edge.trim) {
    X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, A, timeunit, netstat)
  }
  #start outer iteration
  for (o in 1:Outer) {
    # Data augmentation
    for (d in edge.trim) {
        lambda = MultiplyXB(X[[d]], b.old)
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = lmultinom(probij)-1
            }
        }
    }
  # adaptive M-H   
    #if (o > 1) {
    #	accept.rates[1] = accept.rates[1]/Inner[1]
    #	accept.rates[2] = accept.rates[2]/Inner[2]
    #  accept.rates[3] = accept.rates[3]/Inner[1]
    #  accept.rates[4] = accept.rates[1]/Inner[3]
    #	sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    #}
    #accept.rates = rep(0, 4)
   
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    			 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], b.old, u[[max.edge]], delta)
    for (inner in 1:Inner[1]) {
      b.new = rmvnorm_arma(1, b.old, sigma.Q[1]*proposal.var1)
	  delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], b.new, u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
      #  accept.rates[1] = accept.rates[1]+1
      }
        bmat[,inner] = b.old
        deltamat[inner] = delta
    }
    
    if (sum(timestat) > 0) {
      	mu = eta.old[1] + matrix(timemat[edge.trim,] %*% eta.old[2:3], nrow = length(edge.trim), ncol = A)
    } else {
		mu = matrix(eta.old[1], length(edge.trim), A)
    }
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum(mu, sqrt(sigma_tau), senders[edge.trim], timeinc[edge.trim])
      for (inner in 1:Inner[2]) {
         eta.new = rmvnorm_arma(1, eta.old, sigma.Q[2]*proposal.var2)
         if (sum(timestat) > 0) {
      		mu = eta.new[1] + matrix(timemat[edge.trim,] %*% eta.new[2:3], nrow = length(edge.trim), ncol = A) 	
   		 } else {
			mu = matrix(eta.new[1], length(edge.trim), A)
    		}
         prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
         post.new2 = Timepartsum(mu, sqrt(sigma_tau), senders[edge.trim], timeinc[edge.trim])
         loglike.diff = prior.new2+post.new2-prior.old2-post.old2
         if (log(runif(1, 0, 1)) < loglike.diff) {
             eta.old = eta.new
             prior.old2 = prior.new2
             post.old2 = post.new2
             #accept.rates[2] = accept.rates[2]+1
     	}
          etamat[,inner] = eta.old
    }
    if (sum(timestat) > 0) {
      	mu = eta.old[1] + matrix(timemat[edge.trim,] %*% eta.old[2:3], nrow = length(edge.trim), ncol = A)
    } else {
    	mu = matrix(eta.old[1], length(edge.trim), A)
    }
     prior.old3 = log(dinvgamma(sigma_tau, prior.tau[1], prior.tau[2]))
     post.old3 = post.old2
     for (inner in 1:Inner[3]) {
      #sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      #while (sigma_tau.new > 10) {
      #	sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      #}	
      sigma_tau.new = exp(rnorm(1, log(sigma_tau), sqrt(sigma.Q[3])))
      prior.new3 = log(dinvgamma(sigma_tau.new, prior.tau[1], prior.tau[2]))
      post.new3 =  Timepartsum(mu, sqrt(sigma_tau.new), senders[edge.trim], timeinc[edge.trim])
      # loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   # log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   # prior.new3+post.new3-prior.old3-post.old3
      loglike.diff = prior.new3+post.new3-prior.old3-post.old3
     if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
      #   accept.rates[3] = accept.rates[3]+1
     }
        sigma_taumat[inner] = sigma_tau
     }
   	convergence[o] = post.old1 + post.old3
  }
  chain.final = list(b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, edge.trim = edge.trim, convergence = convergence)
  return(chain.final)
}

#' @title GenerateDocs
#' @description Generate a collection of documents according to the generative process of IPTM
#'
#' @param nDocs number of documents to be generated
#' @param node vector of node id's (ID starting from 1)
#' @param b coefficients for recipients
#' @param eta coefficients for timestamps
#' @param delta tuning parameter for the number of recipients
#' @param sigma_tau variance parameter for the timestamps
#' @param support support of latent recipients
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.data edges before 384 hours that is used to calculate initial history of interactions
#' @param backward Logigal indicating whether we are generating backward samples (if FALSE -> forward)
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return generated data including (author, recipients, timestamp, words)
#'
#' @export
GenerateDocs = function(nDocs, node, b, eta, delta, sigma_tau, support, netstat, timestat,
                        base.data = NULL, backward = FALSE, base = FALSE, timeunit = 3600,
                        tz = "America/New_York") {
  A = length(node)
  netstat = as.numeric(c("degree", "dyadic", "triadic" ) %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  t.d = ifelse(base, 0, base.data$edge[[1]][[3]]+384*timeunit)
  edge = base.data$edge
  base.length = length(edge)
  timemat = matrix(0, nrow = nDocs+base.length, ncol = sum(timestat))
  timestamps = rep(NA, nDocs+base.length)
  senders = rep(NA, nDocs+base.length)
  timeinterval = list()
  if (base.length > 0) {
  	timestamps[1:base.length] = vapply(edge, function(d) { d[[3]] }, c(1))
  	senders[1:base.length] = vapply(edge, function(d) { d[[1]] }, c(1))
  	timeinterval[1:base.length] = lapply(1:base.length, function(d) timefinder_vec(timestamps[1:d], d, timeunit))
  }

  if (!base) {
    if (sum(timestat) > 0) {
      Sys.setenv(TZ = tz)
      time_ymd = as.POSIXct(vapply(base.data$edge, function(d) {d[[3]]}, c(1)), tz = getOption("tz"), origin = "1970-01-01")
      if (timestat[1] > 0) {
        days = vapply(time_ymd, function(d) {wday(d)}, c(1))
        days[days==1] = 8
        timemat[1:base.length,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
        it = 1
      }
      if (timestat[2] > 0) {
        hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
        timemat[1:base.length,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
      }
    }
  }
  u = list()
  for (d in 1:nDocs) {
    u[[base.length+d]] = matrix(0, A, A)
    if (t.d >= 384*timeunit) {
       timeinterval = timefinder_vec(timestamps[1:(base.length+d-1)], base.length+d-1, timeunit)
       X = Netstats_cpp(edge, timestamps[1:(base.length+d-1)], timeinterval, senders[1:(base.length+d-1)], A, timeunit, netstat)
       lambda = MultiplyXB(X, b)
    } else {
       lambda = matrix(0, A, A)
    }
    for (i in node) {
      u[[base.length+d]][i,-i] = r.gibbs.measure(lambda[i,-i], delta, support)
    }
    if (sum(timestat) > 0) {
      	 mu = eta[1] + sum(timemat[base.length+d, ] * eta[2:3])
    } else {
      	 mu = eta[1]
    }
    timevec = rlnorm(A, mu, sqrt(sigma_tau))*timeunit
    i.d = which(timevec == min(timevec))
    j.d = which(u[[base.length+d]][i.d,] == 1)
    t.d = t.d+timevec[i.d]
    senders[base.length+d] = i.d
    timestamps[base.length+d] = t.d
    edge[[base.length+d]] = list(author = i.d, recipients = j.d, timestamp = t.d)
    if (t.d <= exp(38.7) & sum(timestat) > 0) {
        Sys.setenv(TZ = tz)
     	it = 0
      time_ymd = as.POSIXct(edge[[base.length+d]][[3]], tz = getOption("tz"), origin = "1970-01-01")
      if (timestat[1] > 0) {
      	it = it + 1
        days = vapply(time_ymd, function(d) {wday(d)}, c(1))
        days[days==1] = 8
        timemat[base.length+d,it] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      }
      if (timestat[2] > 0) {
      	it = it + 1
        hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
        timemat[base.length+d,it] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
      }
    }
  }
  if (base == TRUE & t.d > 384*timeunit) {
    cutoff = which_num(384*timeunit, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1)))-1
    edge = edge[1:cutoff]
    u = u[1:cutoff]
  }
  return(list(edge = edge, base = base.length, u = u, b = b, eta = eta, delta = delta, sigma_tau = sigma_tau))
} 


#' @title GiR_stats
#' @description Calculate several statistics from samples generated from forward or backward sampling
#'
#' @param GiR_sample one sample from generative process
#' @param V number of unique words
#' @param K number of topics
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#'
#' @return A vector of statistics calculated from one GiR sample
#'
#' @export

GiR_stats = function(GiR_sample, V, K, timeunit = 3600) {
  edge = GiR_sample$edge
  if (GiR_sample$base > 0)  {
  	edge = edge[-(1:GiR_sample$base)]
  }
  GiR_stats = c()
  D = length(edge)
  P = length(GiR_sample$b)
  Q = length(GiR_sample$eta)
  GiR_stats[1:P] = GiR_sample$b
  GiR_stats[(P+1):(P+Q)] = GiR_sample$eta
  GiR_stats[P+Q+1] = GiR_sample$delta
  GiR_stats[P+Q+2] = GiR_sample$sigma_tau
  GiR_stats[P+Q+3] = mean(vapply(1:D, function(d) {length(edge[[d]][[2]])}, c(1)))
  GiR_stats[P+Q+4] = mean(vapply(2:D, function(d) {edge[[d]][[3]]-edge[[d-1]][[3]]}, c(1))/timeunit) 			
  return(GiR_stats)
}


#' @title GiR_PP_Plots
#' @description Generate PP-plots (probability-probability plot) for Gettig it Right test
#'
#' @param Forward_stats statistics obtained from forward sampling
#' @param Backward_stats statistics obtained from backward sampling
#'
#' @return PP-plots for different GiR statistics of interest
#'
#' @export
GiR_PP_Plots = function(Forward_stats, Backward_stats) {
  nms = colnames(Forward_stats)
  
  for (i in 1:ncol(Forward_stats)) {
    all = c(Backward_stats[, i], Forward_stats[, i])
    
    quantiles = 100
    if (grepl("b_", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("eta_", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("delta", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("sigma_tau", nms[i]) ) {
      quantiles = 1000
    }
    
    uniqueValues = quantile(all,seq(0, 1, length = quantiles))
    qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
  		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(Forward_stats[, i] <= uniqueValues[j])
  		qx2[j] = mean(Backward_stats[, i] <= uniqueValues[j])
  	}
    
    qqplot(x = qx1,
           y = qx2,
           ylim = c(0, 1),
           xlim = c(0, 1),
           ylab = "Backward",
           xlab = "Forward",
           col = "blue",
           pch = 19,
           cex = 0.25,
           main = nms[i],
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    
      if (nrow(Forward_stats) > 10000) {
       thinning2 = seq(from = floor(nrow(Forward_stats) / 10), to = nrow(Forward_stats), length.out = 10000)
       Forward_test2 = Forward_stats[thinning2, i]
       Backward_test2 = Backward_stats[thinning2, i]
       } else {
        Forward_test2 = Forward_stats[, i]
        Backward_test2 = Backward_stats[, i]    	
      }
    text(paste("Backward Mean:", round(mean(Backward_stats[, i]), 4),
                "\nForward Mean:", round(mean(Forward_stats[, i]), 4),
                "\nt-test p-value:", round(t.test(Backward_test2, Forward_test2)$p.value, 4),
                "\nMann-Whitney p-value:", round(wilcox.test(Backward_test2, Forward_test2)$p.value, 4)),
                x = 0.65, y = 0.15, cex = 0.4)
  }
}      
          
     
#' @title GiR
#' @description Getting it Right test for the IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param nDocs number of documents to be generated per each sample
#' @param node vector of node id's (ID starting from 1)
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.data artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR = function(Nsamp, nDocs, node, prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.data, generate_PP_plots = TRUE) {
  
  A = length(node)
  netstat2 = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat2 = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat2[1]+2*netstat2[2]+4*netstat2[3])
  Q = length(prior.eta[[1]])
  support = gibbs.measure.support(A-1)
  
  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P+Q+5)
  colnames(Forward_stats) = c(paste0("b_",1:P), paste0("eta_",1:Q), "delta", "sigma_tau", 
                              "Mean_recipients", "Mean_timediff")
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = rmvnorm_arma(1, prior.b[[1]], prior.b[[2]])
    eta = rmvnorm_arma(1, prior.eta[[1]], prior.eta[[2]])
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rinvgamma(1, prior.tau[1], prior.tau[2])
    Forward_sample = GenerateDocs(nDocs, node, b, eta, delta, sigma_tau,
                     support, netstat, timestat, base.data = base.data, backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample)
  }
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs(nDocs, node, b, eta, delta, sigma_tau,
                    support, netstat, timestat, base.data = base.data, backward = FALSE, base = FALSE)
  for (i in 1:Nsamp) { 
    if (i %% 1 == 0) {cat("Backward Sampling", i, "\n")}
    Inference_samp = Inference(Backward_sample$edge, node, sigma.Q,
                     prior.b, prior.delta,prior.eta, prior.tau, Outer, Inner,
                     netstat, timestat, initial = NULL)
    b = Inference_samp$b[,ncol(Inference_samp$b)]
    eta = Inference_samp$eta[,ncol(Inference_samp$eta)]
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma_tau = Inference_samp$sigma_tau[length(Inference_samp$sigma_tau)]
    Backward_sample = GenerateDocs(nDocs, node, b, eta, delta, sigma_tau,
                      support, netstat, timestat, base.data = base.data, backward = TRUE, base = FALSE)
    Backward_stats[i, ] = GiR_stats(Backward_sample)
  }
  				
  if (generate_PP_plots) {
    par(mfrow=c(5,6), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          
      

     
#' @title Schein
#' @description Schein (preliminary GiR) test for the IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param nDocs number of documents to be generated per each sample
#' @param node vector of node id's (ID starting from 1)
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.data artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
Schein = function(Nsamp, nDocs, node, prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.data, generate_PP_plots = TRUE) {
  
  A = length(node)
  netstat2 = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat2 = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat2[1]+2*netstat2[2]+4*netstat2[3])
  Q = length(prior.eta[[1]])
  support = gibbs.measure.support(A-1)
  
  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P+Q+4)
  colnames(Forward_stats) = c(paste0("b_",1:P), paste0("eta_",1:Q), "delta", "sigma_tau",
                            "Mean_recipients", "Mean_timediff")
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
					 
  for (i in 1:Nsamp) { 
  	if (i %% 100 == 0) {cat("Sampling", i, "\n")}
    b = rmvnorm_arma(1, prior.b[[1]], prior.b[[2]])
    eta = rmvnorm_arma(1, prior.eta[[1]], prior.eta[[2]])
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rinvgamma(1, prior.tau[1], prior.tau[2])
    Forward_sample = GenerateDocs(nDocs, node, b, eta, delta, sigma_tau,
                     support, netstat, timestat, base.data = base.data, backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample)
  	initial = list(delta = delta, sigma_tau = sigma_tau, b = b, eta = eta, sigma.Q = sigma.Q, u = Forward_sample$u)
    Inference_samp = Inference(Forward_sample$edge, node, sigma.Q, prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner, netstat, timestat, initial = initial)
    b = Inference_samp$b[,ncol(Inference_samp$b)]
    eta = Inference_samp$eta[,ncol(Inference_samp$eta)]
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma_tau = Inference_samp$sigma_tau[length(Inference_samp$sigma_tau)]
    Backward_sample = GenerateDocs(nDocs, node, b, eta, delta, sigma_tau,
                      support, netstat, timestat, base.data = base.data, backward = TRUE, base = FALSE)
    Backward_stats[i, ] = GiR_stats(Backward_sample)
 }
   browser()
  if (generate_PP_plots) {
    par(mfrow=c(4,4), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  } 
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          
      
   
