###library(dqrng)
###library(pracma)
###library(nloptr)
usethis::use_package("dqrng")
usethis::use_package("pracma")
usethis::use_package("nloptr")
## General use function
#' zeroremove
#'
#' @param x vector of real
#'
#' @return real
#' @export
#'
#' @examples
#' zeroremove(c(2.1,0.0,3.0))
zeroremove=function(x){ y=x[x!=0]; z=y[!is.na(y)]; z}

###### EWENS SAMPLING FORMULA        ######################
###### MAXIMUM LIKELIHOOD ESTIMATION ######################

### Compute the maximum likelihood estimate of Ewens theta parameter
### for a given species abundance distribution {n_1, ..., n_k}
### where n_i is the number of organisms of species i in the sample
### and k is the total number of species (for all i, n_i > 0)
### The input should be a vector
#' optim.ewens
#'
#' @param input_abundances a vector of integers
#'
#' @return theta a real value
#' @return log-likelihood a real value
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
    opts = list("algorithm" = "NLOPT_LN_SBPLX", maxeval=2000)
  )
  # return MLE of theta, and value of the logL function at MLE_theta
  return(list(theta=optp$solution[1],logl=--optp$objective))

}

### Compute the maximum likelihood estimate for the multi-deme parameters
### for a given species abundance matrix {n_j1, ..., n_jk} for deme j
### where n_ji is the number of organisms of species i in deme j
### and k is the total number of species (for all i, n_i >= 0)
#' optim.multideme
#'
#' @param input_abundance_matrix a matrix of abundances
#'
#' @return I a vector of real numbers
#' @return m a vector of real numbers
#' @return J a vector of integers
#' @return k a vector of real numbers
#' @export
#'
#' @examples
#' input_abundances1=c(234,87,34,5,4,3,3,2,2,1,1,1,1)
#' input_abundances2=c(240,20,48,2,21,1,3,2,5,2,1,1,1)
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
    local_abundance=input_abundance_matrix[,j]
    k[j]=length(zeroremove(local_abundance))
    J[j]=sum(zeroremove(local_abundance))
    # -loglikelihood function

    logL = function(I){
      value=sum(lgamma(I*x+local_abundance)-lgamma(I*x))-lgamma(I+J[j])+lgamma(I)
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

  }
  m=I/(J-1+I)

  # optimize the logL function. Initial value for theta is theta=1.0
  # theta is bounded by 1.0 and J
  # return MLE of theta, and value of the logL function at MLE_theta
  return(list(I=I,m=m,J=J,k=k))

}

### Compute the maximum likelihood estimate of Pitman (theta,sigma) parameters
### for a given species abundance distribution {n_1, ..., n_k}
### where n_i is the number of organisms of species i in the sample
### and k is the total number of species (for all i, n_i > 0)
### The input should be a vector and initial values of (theta,sigma)
### The argument should be init_vals=c(theta_init,sigma_init)
#' optim.pitman
#'
#' @param input_abundances a vector of integers
#' @param init_vals a vector of two real numbers
#'
#' @return theta a real value
#' @return sigma a real value
#' @return log-likelihood a real value
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

### Create a rank-abundance distribution based on Hoppe's urn scheme for a given theta
### and sampling size nb_balls = J
### Starting at zero species, corresponding to the single black ball with weight theta.
### If picked, a new class is created, else the picked ball is duplicated
### (abundance increases by one unity)

### Added a much faster random number generator
### https://cran.r-project.org/web/packages/dqrng/vignettes/dqrng.html
### and simplified only for the Ewens sampling formula
#' generate.hoppe.urn0
#'
#' @param theta real
#' @param J integer
#'
#' @return output_abundances a vector of integers
#' @return k a real number
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

#### Generate Hoppe urn scheme, ultrafast option
#### Based on the broken stick formulation
#' generate.hoppe.urn
#'
#' @param theta real
#' @param J integer
#'
#' @return output_abundances a vector of integers
#' @return k a real number
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

#### Generate Pitman urn scheme
### Create a rank-abundance distribution based on Pitman's urn scheme for a given (theta,sigma)
### and sampling size nb_balls = J
### Starting at zero species, corresponding to the single black ball with weight theta.
## a species i, of k species total, is selected with probability:
## (abundance[i]-sigma)/(n+theta),
## sigma is non-null, n is the current number of balls
## the probability to pick the black ball i (k*sigma+theta)/(n+theta)

#### Generate Pitman urn scheme, ultrafast option
#### Based on the broken stick formulation, exact
#' generate.pitman.urn
#'
#' @param theta real
#' @param sigma real
#' @param J integer
#'
#' @return output_abundances a vector of integers
#' @return k a real number
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

