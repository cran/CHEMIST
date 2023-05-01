#' @title Generation of Artificial Data
#'
#' @description This function shows the demonstration of data generation based on some
#' specific and commonly used settings, including exponential family distributed
#' potential outcomes, error-prone treatments, and covariates. In this function,
#' users can specify different magnitudes of measurement error and relationship
#' between outcome, treatment, and covariates.
#'
#' @param X The input of n x p dimensional matrix of true covariates, where n is
#' sample size and p is number of covariates. Users can customize the data structure
#' and distribution.
#'
#' @param alpha A vector of the parameters that reflects the relationship between
#' treatment model and covariates. The dimension of \code{alpha} should be equal to the
#' dimension of \code{beta}. If \code{alpha} and \code{beta} have the same nonzero
#' components, then we call them Xc (covariates associated with both outcome and treatment).
#' If components in \code{alpha} are zero but the same components in \code{beta} are
#' nonzero, we call them Xp (covariates associated with outcome only), If components
#' in \code{alpha} are nonzero but the same components in \code{beta} are zero, we
#' call them Xi (covariates associated with treatment only).
#' For example, if \code{alpha = c(2,2,0,0,1,1)} and \code{beta = c(3,3,1,1,0,0)}, then
#' the first two components are Xc, the middle two components are Xp, and the last two
#' components are Xi.
#'
#' @param beta A vector of the parameters that reflects the relationship between
#' outcome and covariates. The dimension of \code{alpha} should be equal to the
#' dimension of \code{beta}. If \code{alpha} and \code{beta} have the same nonzero
#' components, then we call them Xc (covariates associated with both outcome and treatment).
#' If components in \code{alpha} are zero but the same components in \code{beta} are
#' nonzero, we call them Xp (covariates associated with outcome only), If components
#' in \code{alpha} are nonzero but the same components in \code{beta} are zero, we
#' call them Xi (covariates associated with treatment only).
#' For example, if \code{alpha = c(2,2,0,0,1,1)} and \code{beta = c(3,3,1,1,0,0)}, then
#' the first two components are Xc, the middle two components are Xp, and the last two
#' components are Xi.
#'
#' @param theta The scalar of the parameter used to link outcome and
#' treatment.
#'
#' @param a A weight of \code{cov_e} in the measurement error model
#' W = cov_e*a + X + e, where W is observed covariates with measurement error,
#' X is actual covariates, and e is noise term with covaraince matrix \code{cov_e}.
#'
#' @param sigma_e \code{sigma_e} is the common diagonal entries of covariance
#' matrix in the measurement error model.
#'
#' @param e_distr Distribution of the noise term in the classical measurement
#' error model. The input "normal" refers to the normal distribution with mean
#' zero and covariance matrix with diagonal entries \code{sigma_e}. The scalar input "v"
#' represents t-distribution with degree of freedom v.
#'
#' @param num_pi Settings of misclassification probability with option 1 or 2.
#' \code{num_pi = 1} gives that pi_01 equals  pi_10, and \code{num_pi = 2} refers to that
#' pi_01 is not equal to pi_10.
#'
#' @param delta The parameter that determines number of treatment with measurement
#' error. \code{delta = 1} has equal number of treatment with and without measurement
#' error. We set \code{default = 0.5} since it has smaller number of treatment who has
#' measurement error.
#'
#' @param linearY The boolean option that determines the relationship between
#' outcome and covariates. \code{linearY = TRUE} gives linear relationship with a
#' vector of parameters \code{alpha}, \code{linearY = FALSE} refers to non linear
#' relationship between outcome and covariates, where the sin function is specified on
#' Xc and the exponential function is specified on Xp.
#'
#' @param typeY The outcome variable with exponential family distribution
#' "binary", "pois" and "cont". \code{typeY = "binary"} refers to binary random
#' variables, \code{typeY = "pois"} refers to Poisson random variables, and
#' \code{typeY = "cont"} refers to normally distributed random variables.
#'
#' @return \item{Data}{A n x (p+2) matrix of the original data without measurement error,
#' where n is sample size and the first p columns are covariates with the order being
#' Xc (the covariates associated with both treatment and outcome),
#' Xp (the covariates associated with outcome only),
#' Xi (the covariates associated with treatment only),
#' Xs (the covariates independent of outcome and treatment),
#' the last second column is treatment, and the last column is outcome.}
#'
#' @return \item{Error_Data}{A n x (p+2) matrix of the data with measurement error in covariates
#' and treatment, where n is sample size and the first p columns are covariates
#' with the order being
#' Xc (the covariates associated with both treatment and outcome),
#' Xp (the covariates associated with outcome only),
#' Xi (the covariates associated with treatment only),
#' Xs (the covariates independent of outcome and treatment),
#' the last second column is treatment, and the last column is outcome.}
#'
#' @return \item{Pi}{A n x 2 matrix containing two misclassification probabilities pi_10 =
#' P(Observed Treatment = 1 | Actual Treatment = 0) and pi_01 =
#' P(Observed Treatment = 0 | Actual Treatment = 1) in columns.}
#'
#' @return \item{cov_e}{A covariance matrix of the measurement error model.}
#'
#' @examples
#' ##### Example 1: A multivariate normal continuous X with linear normal Y #####
#'
#' ## Generate a multivariate normal X matrix
#' mean_x = 0; sig_x = 1; rho = 0
#' Sigma_x = matrix( rho*sig_x^2,nrow=120 ,ncol=120 )
#' diag(Sigma_x) = sig_x^2
#' Mean_x = rep( mean_x, 120 )
#' X = as.matrix( mvrnorm(n = 60,mu = Mean_x,Sigma = Sigma_x,empirical = FALSE) )
#'
#' ## Data generation setting
#' ## alpha: Xc's scale is 0.2 0.2 and Xi's scale is 0.3 0.3
#' ## so this refers that there is 2 Xc and Xi
#' ## beta: Xc's scale is 2 2 and Xp's scale is 2 2
#' ## so this refers that there is 2 Xc and Xp
#' ## rest with following setup
#' Data_fun <- Data_Gen(X, alpha = c(0.2,0.2,0,0,0.3,0.3), beta = c(2,2,2,2,0,0)
#' , theta = 2, a = 2, sigma_e = 0.75, e_distr = 10, num_pi = 1, delta = 0.8,
#' linearY = TRUE, typeY = "cont")
#'
#' ##### Example 2: A uniform X with non linear binary Y #####
#'
#' ## Generate a uniform X matrix
#' n = 50; p = 120
#' X = matrix(NA,n,p)
#' for( i in 1:p ){ X[,i] = sample(runif(n,-1,1),n,replace=TRUE ) }
#' X = scale(X)
#'
#' ## Data generation setting
#' ## alpha: Xc's scale is 0.1 and Xi's scale is 0.3
#' ## so this refers that there is 1 Xc and Xi
#' ## beta: Xc's scale is 2 and Xp's scale is 3
#' ## so this refers that there is 1 Xc and Xp
#' ## rest with following setup
#' Data_fun <- Data_Gen(X, alpha = c(0.1,0,0.3), beta = c(2,3,0)
#' , theta = 1, a = 2, sigma_e = 0.5, e_distr = "normal", num_pi = 2, delta = 0.5,
#' linearY = FALSE, typeY = "binary")
#'
#' @export
#' @importFrom stats "rgamma" "rt" "rbinom" "rpois" "rnorm"
#' @importFrom MASS "mvrnorm"
#' @importFrom LaplacesDemon "rmvt"

Data_Gen <- function(X,alpha,beta,theta,a,sigma_e,e_distr="normal",num_pi,delta,linearY,typeY){
  Data = NULL
  X=X
  n=nrow(X)
  p=ncol(X)
  delta=delta
  theta=theta
  a=a
  sigma_e=sigma_e
  var.list2=c()

  scale_list = as.data.frame(rbind(alpha,beta))

  colnames(scale_list)[ which(scale_list[1,]*scale_list[2,]!=0 ) ] = "Xc"
  colnames(scale_list)[ which(scale_list[1,]==0 ) ] = "Xp"
  colnames(scale_list)[ which(scale_list[2,]==0 ) ] = "Xi"

  pC = sum(colnames(scale_list)=="Xc" )  # pC: associate with both
  pP = sum(colnames(scale_list)=="Xp" )  # pP: associate with outcome
  pI = sum(colnames(scale_list)=="Xi" )  # pI: associate with treatment

  scale_1 = c( unlist( scale_list[1,colnames(scale_list)=="Xc"] ) )
  scale_2 = c( unlist( scale_list[1,colnames(scale_list)=="Xi"] ) )
  scale_3 = c( unlist( scale_list[2,colnames(scale_list)=="Xc"] ) )
  scale_4 = c( unlist( scale_list[2,colnames(scale_list)=="Xp"] ) )

  pS = p - (pC+pI+pP) # pS: associate with neither

  var.list2 = c( paste( "Xc", 1:pC, sep ="") , paste( "Xp", 1:pP, sep ="") ,
                 paste( "Xi", 1:pI, sep =""), paste( "XS", 1:pS, sep =""))

  Beta0 = rgamma(n,5,5)-delta
  Gamma0 = rgamma(n,5,5)-delta #delta=1 has equal num 1 0  ,default 0.5 since has small pi10
  pi_10 = 1-exp(Beta0)/ (1+exp(Beta0))
  pi_01 = 1-exp(Gamma0)/ (1+exp(Gamma0))

  if(num_pi==1){pi_10=pi_01}

  #generate data
  colnames(X) = var.list2

  #generate measurement error for X\
  W = matrix(NA,n,p)
  cov_e = matrix(0,p,p)
  diag(cov_e) = sigma_e
  mean_e = rep(0,p)
  if(e_distr == "normal"){
    e = mvrnorm(n = n, mu = mean_e, Sigma = cov_e)
    W = sigma_e*a + X + e

  }
  else{
    e =  rmvt(n = n, mu = mean_e, S=cov_e, df=e_distr)
    W = sigma_e*a + X + e
  }


  # set associate with treatment
  prob =  rowSums( cbind( X[,grepl("Xc", colnames(X))]*scale_1 , X[,grepl("Xi", colnames(X))]*scale_2 ) )
  prob = exp(prob)/(1+exp(prob))
  D = (prob >0.5)*1
  #D = rbinom(n, size = 1,prob)
  D_star = rep(0,n)
  for( i in 1:n){
    pi_matrix = matrix( c(1-pi_10[i],pi_10[i],pi_01[i],1-pi_01[i]) ,nrow=2,ncol=2)
    prob_star = pi_matrix %*% rbind(prob[i],1-prob[i])
    D_star[i] = (prob_star[1,]>prob_star[2,])*1
  }

  if(linearY == TRUE){
    CP_vec = rowSums( cbind( X[,grepl("Xc", colnames(X))]*scale_3 , X[,grepl("Xp", colnames(X))]*scale_4) )
  }
  else{
    CP_vec =  rowSums( cbind(sin(X[,grepl("Xc", colnames(X))])*scale_3 , exp(X[,grepl("Xp", colnames(X))])*scale_4))
  }

  #set associate with outcome
  if(typeY == "binary") {
    prob = exp(D*theta + CP_vec  ) /(1+exp(D*theta + CP_vec ))
    #y = rbinom(n, size = 1,prob )
    y = (prob >0.5)*1

    #TRUE ATE
    prob1 = exp(theta + CP_vec  ) /(1+exp(theta + CP_vec ))
    prob0 = exp(CP_vec) /(1+exp(CP_vec ))
    #ATE =mean(rbinom(n,1,prob1)) -mean(rbinom(n,1,prob0))
    ATE =mean((prob1 >0.5)*1) -mean((prob0 >0.5)*1)
  }
  else if(typeY == "pois"){
    lam = abs(D*theta + CP_vec )
    y = rpois(n,lam)

    #TRUE ATE
    lam1 = abs(theta + CP_vec )
    lam0 = abs(CP_vec )
    ATE = mean(rpois(n,lam1)) -mean(rpois(n,lam0))
  }
  else {
    y = D*theta + CP_vec +rnorm(n,0,1)

    #TRUE ATE
    ATE = mean(theta+CP_vec +rnorm(n,0,1))-mean(CP_vec +rnorm(n,0,1))
  }

  Data = as.data.frame(cbind(X,D,y) )
  colnames(Data) <- c(var.list2,"D","Y")
  Error_Data = as.data.frame(cbind(W,D_star,y) )
  colnames(Error_Data) <-  c(var.list2,"D","Y")
  return( list(Data=Data, Error_Data = Error_Data, Pi = cbind(pi_10,pi_01), cov_e = cov_e) )
}

