usethis::use_package("dqrng")
usethis::use_package("pracma")
usethis::use_package("nloptr")

#' General use function
#' @description zeroremove is used to remove zeros and NAs in a vector of reals or of integers. This is a simple utility function.
#'
#' @param x vector of real
#'
#' @return real
#' @export
#'
#' @examples
#' zeroremove(c(2.1,0.0,3.0))
zeroremove=function(x){ y=x[x!=0]; z=y[!is.na(y)]; z}


#' Maximum likelihood estimation of \eqn{\theta} based on Ewens sampling formula
#'
#' @description optim.ewens computes the maximum likelihood estimate of Ewens \eqn{\theta} parameter
#' for a given species abundance distribution \eqn{(n_1, ..., n_k)}
#' where \eqn{n_i} is the number of organisms of species i in the sample
#' and k is the total number of species in the sample (for all \eqn{i, n_i > 0}).
#' The input should be a vector
#'
#' @param input_abundances a vector of integers
#'
#' @return theta the \eqn{\theta} parameter of Ewens sampling formula, a real value
#' @return logl the maximal value of the log-likelihood, a real value
#' @export
#'
#' @examples
#' input_abundances=c(234,87,34,5,4,3,3,2,2,1,1,1,1)
#' optim.ewens(input_abundances)
optim.ewens = function(input_abundances){

  J=sum(input_abundances)
  k=length(input_abundances)

  #Alternative calculation of theta: theta such that  k=kest(theta)
  #kest=function(theta){theta*digamma(theta+J)-theta*digamma(theta)-k}
  #out=uniroot(kest,lower=1,upper=J, tol = 1e-6)
  #theta=out$root

  # -loglikelihood function (log-transform of the Ewens sampling formula)
  logL = function(theta){
    value=k*log(theta)-lgamma(theta+J)+lgamma(theta)+lgamma(J+1)-lgamma(k+1)-sum(log(input_abundances))
    -value
  }

  # optimize the logL function. Initial value for theta is theta=1.0
  # theta is bounded by 1.0 and J
  optp = nloptr::nloptr(
    x0 = 1.0,
    eval_f = logL,
    lb=c(0.0),
    ub=c(as.numeric(J)),
    opts = list("algorithm" = "NLOPT_LN_SBPLX")
  )
  # return MLE of theta, and value of the logL function at MLE_theta
  return(list(theta=optp$solution[1],logl=--optp$objective))

}

#' Maximum likelihood estimation of local immigration rates \eqn{(m_1,...,m_K)} based on the K-deme sampling formula
#'
#' @description optim.multideme computes the maximum likelihood estimate for the multi-deme parameters
#' for a given species abundance matrix \eqn{(n_{1j}, ..., n_{kj})} for deme j (\eqn{j\in\{1,...,K\}})
#' where \eqn{n_{ij}} is the number of organisms of species i in deme j
#' and k is the total number of species (for all \eqn{i, n_i \geq 0}).
#' The formulation in this description is that of species abundance distributions
#' but is valid for any partition in a subdivided setting (e.g., word frequencies in multiple books).
#'
#' @param input_abundance_matrix a matrix \eqn{n_{ij}} of abundances for species i in deme (local site) j
#'
#' @return I rescaled immigration rate, \eqn{I=\frac{m}{1-m}(J-1)}, a vector of real numbers
#' @return m immigration rate, a vector of real numbers
#' @return J local community size, a vector of integers
#' @return k number of classes/species, a vector of integers
#' @export
#'
#' @examples
#' input_abundances1=c(44,37,34,5,4,3,3,2,2,1,1,1,1)
#' input_abundances2=c(240,20,48,2,21,1,3,2,5,2,0,1,1)
#' input_abundance_matrix=rbind(input_abundances1,input_abundances2)
#' optim.multideme(input_abundance_matrix)
optim.multideme = function(input_abundance_matrix){

  input_abundance=colSums(input_abundance_matrix) ## sum over species
  x=input_abundance/sum(input_abundance) ## regional species abundance distribution
  demes=dim(input_abundance_matrix)[1]

  I=rep(0,demes)
  m=rep(0,demes)
  J=rep(0,demes)
  k=rep(0,demes)
  for(j in 1:demes){
    local_abundance=as.numeric(input_abundance_matrix[j,])
    k[j]=length(zeroremove(local_abundance))
    J0=sum(zeroremove(local_abundance));J[j]=J0
    list1=which(local_abundance != 0)
    x1=as.numeric(x[list1]);local_abundance1=local_abundance[list1]

    # -loglikelihood function
    logL = function(I){
      value=sum(lgamma(I*x1+local_abundance1)-lgamma(I*x1))-lgamma(I+J0)+lgamma(I)
      -value
    }
    if(k[j]>1 && J[j]>1){
      optp = nloptr::nloptr(
        x0 = 1.0,
        eval_f = logL,
        lb=c(0.0),
        ub=c(as.numeric(J[j])),
        opts = list("algorithm" = "NLOPT_LN_SBPLX", maxeval=2000)
      )
      I[j]=optp$solution[1]
    }
    print(paste("Deme: ", j, "..."))

  }
  m=I/(J-1+I)

  # optimize the logL function. Initial value for theta is theta=1.0
  # theta is bounded by 1.0 and J
  # return MLE of theta, and value of the logL function at MLE_theta
  return(list(I=I,m=m,J=J,k=k))

}

#' Maximum likelihood estimation of \eqn{\theta} and \eqn{\sigma} based on Pitman sampling formula
#'
#' @description optim.pitman computes the maximum likelihood estimate of Pitman \eqn{(\theta,\sigma)} parameters
#' for a given species abundance distribution \eqn{(n_1, ..., n_k)}
#' where \eqn{n_i} is the number of organisms of species i in the sample
#' and k is the total number of species (for all \eqn{i, n_i > 0})
#' The input should be a vector and initial values of \eqn{(\theta,\sigma)}
#' with the argument init_vals=c(theta_init,sigma_init)
#' with the condition theta_init>0.0, and 1.0 > sigma_init > 0.0.
#'
#' @param input_abundances a vector of integers
#' @param init_vals initial values of (theta,sigma), a vector of two real numbers
#'
#' @return theta parameter \eqn{\theta}, a real value
#' @return sigma parameter \eqn{\sigma}, a real value
#' @return logl, maximal log-likelihood, a real value
#' @export
#'
#' @examples
#' optim.pitman(c(234,87,34,5,4,3,3,2,2,1,1,1,1),c(10.0,0.1))
optim.pitman = function(input_abundances,init_vals=c(10.0,0.1)){

  J=sum(input_abundances)
  k=length(input_abundances)

  # -loglikelihood function (log-transform of the Pitman sampling formula)
  logL = function(x){
    theta=x[1]
    sigma=x[2]
    value=k*log(sigma)-k*lgamma(1-sigma)+lgamma(theta)-lgamma(theta+J)+lgamma(theta/sigma+k)-lgamma(theta/sigma)+sum(lgamma(input_abundances-sigma))
    -value
  }

  # optimize the logL function. Initial values for theta and sigma are user-supplied
  # theta is bounded by 1.0 and J
  # sigma is bounded by 0.0 and 1.0
  optp = nloptr::nloptr(
    x0 = init_vals,
    eval_f = logL,
    lb=c(1.0,0.0),
    ub=c(as.numeric(J),1.0),
    opts = list("algorithm" = "NLOPT_LN_SBPLX", maxeval=2000)
  )

  # return MLE of (theta,sigma) and value of the logL function at MLE_(theta,sigma)
  return(list(theta=optp$solution[1],sigma=optp$solution[2],logl=-optp$objective))

}

###### GENERATING URN MODELS   ######################

#' Generation of a neutral partition given parameter \eqn{\theta} using Hoppe urn model
#'
#' @description generate.hoppe.urn0 creates a rank-abundance distribution based on Hoppe's urn scheme
#' for a given parameter \eqn{\theta} and sampling size J.
#' The process starts at zero species, corresponding to the single black ball with weight \eqn{\theta}.
#' If the black ball is picked, a new class is created, else the picked colored ball is duplicated
#' (its abundance increases by one unity).
#' Added a much faster random number generator (dqrng package).
#'
#' @param theta parameter \eqn{\theta}, a real value
#' @param J  number of individuals, an integer
#'
#' @return output_abundances is a vector of integers
#' @return k is the expected number of classes, a real number
#' @export
#'
#' @examples
#' generate.hoppe.urn0(11.3,234373)
generate.hoppe.urn0 = function(theta,J) {

  if (J < 1) {
    stop("Error generate.hoppe.urn0:: ",
         "J must be > 0")
  }
  if (theta < 0.0 || theta > J) {
    stop("Error generate.hoppe.urn0:: ",
         "theta>J or theta<0.0")
  }

  spec_ct = 0
  species  = rep(NaN, J)
  ## Draw all the random numbers at once (not in the loop)
  ## replace runif rng by dqrunif rng
  listrunif=dqrng::dqrunif(J, 0, 1)
  listrunif2=dqrng::dqrunif(J, 0, 1)

  for (j in 1:J) { # loop over each individual
    if (listrunif[j] <= theta / (theta + j - 1)) {
      spec_ct = j + 1
      species[j] = j
    } else {
      prior_lineage = pracma::ceil(listrunif2[j] * (j - 1))
      species[j] = species[prior_lineage]
    }
  }

  x = c(table(species))
  output_abundances = c()
  for (i in seq_along(x)) {
    output_abundances[i] = x[[i]]
  }
  output_abundances = sort(output_abundances, decreasing = TRUE)
  return(output_abundances)
}

#' Generation of a neutral partition given parameter \eqn{\theta} using the Griffith-Engen-McCloskey representation
#'
#' @description generate.hoppe.urn creates a rank-abundance distribution for a given \eqn{\theta} and J.
#' Beta-distributed random variables are picked from distribution \eqn{Beta(1,\theta)},
#' and the ressulting partition is one multinomial draw of J individuals drawn according to the GEM random variables.
#'
#' @param theta parameter \eqn{\theta}, a real value
#' @param J  number of individuals, an integer
#'
#' @return output_abundances is a vector of integers
#' @return k is the expected number of classes, a real number
#' @export
#'
#' @examples
#' generate.hoppe.urn(11.3,234373)
generate.hoppe.urn = function(theta,J){

  if (J < 1) {
    stop("Error generate.hoppe.urn0:: ",
         "J must be > 0")
  }
  if (theta < 0.0 || theta > J) {
    stop("Error generate.hoppe.urn0:: ",
         "theta>J or theta<0.0")
  }

  ## estimate k from theta and J
  kest=theta*(digamma(theta+J)-digamma(theta))
  ## Ensure that there are more classes than kest (set arbitrarily at 2*kest)
  k2=2.0*kest

  ## Draw k2 Beta(1,theta)-distributed random numbers Ws
  Ws=stats::rbeta(k2, 1, theta, ncp = 0)

  ## create Zs: product of (1-Ws) from 1 to i
  Zs=Ws;Zs[1]=1;for(i in 2:k2) Zs[i]=Zs[i-1]*(1-Ws[i])

  ## Create Ps: probability density of the GEM model
  Ps=Ws;for(i in 1:k2) Ps[i]=Zs[i]*Ws[i]

  ## Multinomial draw of J individuals from probability density Ps
  x=stats::rmultinom(1,J,prob=Ps)

  ## Filter and sort the species abundance distribution
  abundance=if(length(which(x==0)!=0)) x[-which(x==0)]
  output_abundance = sort(abundance, decreasing = TRUE)

  ## Output the species abundance distribution and the expected class number kest
  return(list(abundance=output_abundance,k=kest))
}


#' Generation of a Pitman partition given parameters \eqn{\theta} and \eqn{\sigma} from the Griffith-Engen-McCloskey representation
#'
#' @description The Pitman urn scheme is as follows: starting at zero species, corresponding to the single black ball with weight \eqn{\theta}
#' a species i, of k species total, is selected with probability:
#' \eqn{(n_i-\sigma)/(n+\theta)},
#' \eqn{\sigma} is strictly between 0.0 and 1.0, n is the current number of balls
#' the probability to pick the black ball i \eqn{(k\sigma+\theta)/(n+\theta)}.
#' This function generate.pitman.urn creates a rank-abundance distribution for a given theta,sigma and J objects.
#' Beta-distributed random variables are picked from distribution \eqn{Beta(1-\sigma,\theta+k\sigma)},
#' and the partition is the multinomial draw of J individuals drawn according to the GEM random variables
#'
#' @param theta parameter \eqn{\theta},  a real value
#' @param sigma parameter \eqn{\sigma},  a real value
#' @param J is the number of individuals, an integer
#'
#' @return output_abundances is a vector of integers
#' @return k is the expected number of classes, a real number
#' @export
#'
#' @examples
#' generate.pitman.urn(11.3,0.1,234373)
generate.pitman.urn = function(theta,sigma,J){

  if (theta < 0.0 || theta > J) {
    stop("Error generate.hoppe.urn0:: ",
         "theta>J or theta<0.0")
  }
  if (sigma <= 0.0 || sigma >= 1.0) {
    stop("Error generate.hoppe.urn0:: ",
         "sigma>=1.0 or sigma<=0.0")
  }

  ## Number of species (classes) k inferred from Piman model
  kest=(exp(lgamma(theta+1)-lgamma(theta+sigma)+lgamma(J+theta+sigma)-lgamma(J+theta))-theta)/sigma

  ## Ensure that there are more classes than kest (set arbitrarily at 2*kest)
  k2=2.0*kest

  ## create Ws: Beta random variable that depends on the sequence order
  Ws=rep(0,k2); notWs=rep(0,k2) ## see Eq (156) Pitman 2006
  for(i in 1:k2) {Ws[i]=stats::rbeta(1, 1-sigma, theta+i*sigma, ncp = 0);notWs[i]=1.0-Ws[i];}

  ## Define Zs[i] = notWs[1]*...notWs[i-1], with Zs[1]=1
  Zs=rep(1.0,k2); for(i in 2:k2) Zs[i]=Zs[i-1]*notWs[i]

  ## Create Ps: probability density of the GEM-Pitman model
  Ps=rep(0,k2); for(i in 1:k2) Ps[i]=Zs[i]*Ws[i]

  ## Multinomial draw of J individuals from probability density Ps
  x=stats::rmultinom(1,J,prob=Ps)

  abundance=if(length(which(x==0)!=0)) x[-which(x==0)]
  output_abundance = sort(abundance, decreasing = TRUE)

  ## Output the species abundance distribution and the expected class number kest
  return(list(abundance=output_abundance,k=kest))

}

