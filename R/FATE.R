#' @title Estimation of ATE under high-dimensional error-prone data
#'
#' @description This function aims to estimate ATE by selecting informative
#' covariates and correcting for measurement error in covariates and
#' misclassification in treatments. The function FATE reflects the strategy of
#' estimation method: Feature screening, Adaptive lasso, Treatment  adjustment,
#' and Error correction for covariates.
#'
#' @param Data A n x (p+2) matrix of the data, where n is sample size and the first
#' p columns are covariates with the order being
#' Xc (the covariates associated with both treatment and outcome),
#' Xp (the covariates associated with outcome only),
#' Xi (the covariates associated with treatment only),
#' Xs (the covariates independent of outcome and treatment),
#' the last second column is treatment, and the last column is outcome.
#'
#' @param cov_e Covariance matrix in the measurement error model.
#'
#' @param Consider_D Feature screening with treatment effects accommodated.
#' \code{Conidser_D = TRUE} refers to feature screening with A and (1-A) incorporated.
#' \code{Consider_D = FALSE} will not multiply with A and (1-A).
#'
#' @param pi_10 Misclassifcation probability is
#' P(Observed Treatment = 1 | Actual Treatment = 0).
#'
#' @param pi_01 Misclassifcation probability is
#' P(Observed Treatment = 0 | Actual Treatment = 1).
#'
#' @return \item{ATE}{A value of the average treatment effect.}
#'
#' @return \item{wAMD}{A weighted absolute mean difference.}
#'
#' @return \item{Coef_prop_score}{A table containing coefficients of propensity score.}
#'
#' @return \item{Kersye_table}{The selected covariates by feature screening.}
#'
#' @return \item{Corr_trt_table}{A summarized table containing corrected treatment.}
#'
#' @examples
#' ##### Example 1: Input the data without measurement correction #####
#'
#' ## Generate a multivariate normal X matrix
#' mean_x = 0; sig_x = 1; rho = 0; n = 50; p = 120
#' Sigma_x = matrix( rho*sig_x^2 ,nrow=p ,ncol=p )
#' diag(Sigma_x) = sig_x^2
#' Mean_x = rep( mean_x, p )
#' X = as.matrix( mvrnorm(n ,mu = Mean_x,Sigma = Sigma_x,empirical = FALSE) )
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
#' ## Extract Ori_Data, Error_Data, Pi matrix, and cov_e matrix
#' Ori_Data=Data_fun$Data
#' Pi=Data_fun$Pi
#' cov_e=Data_fun$cov_e
#' Data=Data_fun$Error_Data
#' pi_01 = pi_10 = Pi[,1]
#'
#' ## Input data into model without error correction
#' Model_fix = FATE(Data, matrix(0,p,p), Consider_D = FALSE, 0, 0)
#'
#' ##### Example 2: Input the data with measurement correction #####
#'
#' ## Input data into model with error correction
#' Model_fix = FATE(Data, cov_e, Consider_D = FALSE, Pi[,1],Pi[,2])
#'
#' @export
#' @importFrom stats "sd" "binomial" "coef"


FATE <- function(Data, cov_e, Consider_D, pi_10, pi_01){
  n = nrow(Data)
  pi_10 = pi_10
  pi_01 = pi_01


  if( Consider_D == TRUE){
    kersye <- choose_p_D(Data, cov_e)
  }
  if( Consider_D == FALSE){
    kersye <- choose_p(Data, cov_e)
  }
  kersye <- kersye[sort(names(kersye))]
  var.list = names( kersye )

  # find the position of choose p for cov_e
  p = length(var.list)
  pos = c()
  for(i in 1:length(var.list)){ pos = c(pos, which( colnames(Data)==var.list[i]) )   }
  Data <- Data[ , c(var.list, "D", "Y") ]
  cov_e <- cov_e[pos,pos]

  colnames(Data)[length(Data)-1] <- "A"
  n.p=n+p

  ## set vector of possible lambda's to try
  # set lambda values for grid search.

  lambda_vec =  c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
  names(lambda_vec) = as.character(lambda_vec)

  ## lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2

  ## get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*(gamma_convergence_factor - lambda_vec + 1)
  names(gamma_vals) = names(lambda_vec)

  ## normalize covariates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] = Data[,var.list] - Temp.mean
  temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] = Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))

  # weight for each variable for refined selection
  weight = kersye
  #betaXY <- weight

  mins = min(weight)
  maxs = max(weight)
  betaXY = (weight - mins) / (maxs-mins)
  betaXY[which(betaXY==0)] <- 10^-7
  #betaXY <- weight/ max(weight)

  ## want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec =coeff_XA=NULL
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list

  ## set the possible lambda2 value (taken from Zou and Hastie (2005))
  S_lam=c(0,10^c(-2,-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,1))

  ## want to save ATE, wAMD and propensity score coefficients for each lambda2 value
  WM_N=M_N=S_wamd=rep(NA,length(S_lam))
  M_mat=matrix(NA,length(S_lam),p)
  colnames(M_mat)= var.list

  ## fix A's measurements error
  matrix_A = as.data.frame( matrix( Data$A,n,3))
  colnames(matrix_A) = c("no fix A","fix A","get change")

  Data$A <- as.vector( ( Data$A - pi_10 ) / ( 1 - pi_10 - pi_01) )
  Data$A = (Data$A > 0.5)*1

  matrix_A[,2] = Data$A
  matrix_A[,3] = ifelse(matrix_A[,1]==matrix_A[,2] , " ", "fix" )

  ## fix X's measurement error
  W =  as.matrix( Data[ , c(var.list) ] )
  mu_W = colMeans(W);length(mu_W)
  sigma_W = (1/(n-1)) * t(W-mu_W)%*%(W-mu_W)  ; dim(sigma_W)
  #sigma_W = cov(W)

  for(i in 1: n){
    Data[i, c(var.list) ] =  mu_W + t( sigma_W-cov_e  ) %*% solve(sigma_W) %*% (W[i,]-mu_W)
  }


  for (m in 1:length(S_lam)) {
    ## create augmented A and X
    lambda2=S_lam[m]
    I=diag(1,p,p)
    Ip=sqrt(lambda2)*I
    Anp=c(Data$A,rep(0,p))
    Xnp=matrix(0,n+p,p)
    X=Data[,var.list]
    for (j in 1:p){
      Xnp[,j]=c(X[,j],Ip[,j])
    }
    newData=as.data.frame(Xnp)
    names(newData)=var.list
    newData$A=Anp

    ## want to save ATE, wAMD and propensity score coefficients for each lambda value
    ATE_try=ATE = wAMD_vec =coeff_XA=NULL;
    ATE_try=ATE = wAMD_vec = rep(NA, length(lambda_vec))
    names(ATE) = names(wAMD_vec) = names(lambda_vec)
    coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list

    ## run GOAL with lqa using augmented data
    # weight model with all possible covariates included, this is passed into lasso function

    for( lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]

      # create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n.p^(il),al.weights = abs(betaXY)^(-ig) )
      # run outcome-adaptive lasso model with appropriate penalty
      X=as.matrix(newData[var.list]);y=as.vector(newData$A);
      logit_oal = lqa.default ( X,y, penalty=oal_pen, family=binomial(logit) )
      # save propensity score coefficients    !!!!adaptive elastic net
      coeff_XA[var.list,lil] = (1+lambda2)*coef(logit_oal)[var.list]
      # generate propensity score
      Data[,paste("f.pA",lil,sep="")]=
        expit(as.matrix(cbind(rep(1,n),Data[var.list]))%*%as.matrix((1+lambda2)*coef(logit_oal)))

      # create inverse probability of treatment weights for ATE
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
      # save ATE estimate for this lambda value
      ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)


    } # close loop through lambda values

    # print out wAMD for all the lambda values evaluated
    wAMD_vec
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)
    # print out ATE corresponding to smallest wAMD value
    ATE[tt][[1]]
    # save the coefficients for the propensity score that corresponds to smallest wAMD value
    GOAL.lqa.c=coeff_XA[,tt]
    names(GOAL.lqa.c) = var.list

    # check which covariates are selected
    #M_mat[m,]=ifelse(abs(coeff_XA[,tt])> 10^(-8),1,0)
    M_mat[m,]=coeff_XA[,tt]
    # save the ATE corresponding to smallest wAMD value

    M_N[m]=ATE[tt][[1]]

    # save the smallest wAMD value
    WM_N[m]=wAMD_vec[tt][[1]]
  }

  # find the optimal lambda2 value that creates the smallest wAMD

  ntt= which.min(WM_N)

  # save the ATE corresponding the optimal lambda2 value

  GOAL_ATE =M_N[ntt]

  # save the propensity score coefficients corresponding the optimal lambda2 value

  GOAL.lqa.c=M_mat[ntt,]

  Coefficients_propensity_score_table = as.data.frame( matrix(GOAL.lqa.c,length(GOAL.lqa.c),3) )
  colnames(Coefficients_propensity_score_table) = c("Covariate","Coefficients Propensity Score","Significantness")
  Coefficients_propensity_score_table[,1] = var.list
  Coefficients_propensity_score_table[,3] = ifelse(abs(Coefficients_propensity_score_table[,3])> 10^(-8),"***"," ")

  # save wAMD value corresponding the optimal lambda2 value

  wAMD_vec=WM_N[ntt]

  return( list( ATE = GOAL_ATE , wAMD = wAMD_vec,
                Coef_prop_score = Coefficients_propensity_score_table,
                Kersye_table = kersye , Corr_trt_table = matrix_A  ) )
}

