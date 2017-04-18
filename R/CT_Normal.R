#' SSR for adaptive multi-arm design with normal outcomes
#'
#' Function to calculate conditional power and perform SSR for adaptive multi-arm designs with normal outcomes. The adaptive design is either based on inverse normal combination function or Fisher's combination function.
#' @param combF Combination function, either "Inverse Normal" or "Fisher".
#' @param a Number of experimental treatments (a = 2 or 3)
#' @param sel_rule Treatment selection rule at interim analysis. There are four treatment selection rule, 1) sel_tr = 1, select the best treatment; 2) sel_tr = 2, keep all promising treatments that do not cross either efficacy boundaries or futility boundaries; 3) sel_tr = 3, select the best treatment only if the mean difference (treatment - control if side = "upper", control - treatment if side = "lower") is above a threshold delta; 4) sel_tr = 4, keep all promising treatments where the mean difference is above delta.
#' @param n1 First stage sample size per arm (default is NULL).
#' @param N2 Total sample size per arm at the end of the trial (default is NULL).
#' @param initialSS Whether initial sample size needs to be calculated. If initial sample size does not need to be calculated, then n1 and N2 have to be provided and initialSS is set to 0. If initial sample size need to be calculated, then initialSS can take numeric value 1 (initial sample size will be calculated using the method for traditional fixed two-arm design) or 2 (initial sample size will be calculated using the method for group sequential two-arm design).
#' @param r21 Allocation ratio of N2 to n1 (r21 = N2/n1). If initialSS = NULL, then it should also be NULL. If initialSS = 1 or 2, then it should be a numeric value that is larger than 1.
#' @param rmax The ratio of the maximum sample size allowed per arm (N2_max) to N2 (rmax = N2_max/N2). It should be a numeric value that is larger than 1.
#' @param alpha One sided familywise error rate (default is 0.025).
#' @param beta Type II error. 1-beta is used as the desired power and conditional power.
#' @param ushape The shape of the efficacy boundaries (default is "obf"). Either a function specifying the shape or one of "pocock", "obf" (the default), "triangular" and "fixed".
#' @param lfix Fixed lower boundary on z-scale (default is 0). lfix must be smaller than Phi(1-alpha)/2 to ensure that it is smaller than the upper boundary.
#' @param theta Planned absolute mean difference of experimental treatment to control. It is used for initial sample size caculation and SSR based on the planned treatment effect.
#' @param sigma Known standard deviation of normal distribution.
#' @param x0_1 First stage sample mean of control arm. It has to be a one column matrix with the same number of rows as x1.
#' @param x1 First stage sample mean of experimental treatment arms. It has to be a matrix with each row represents one trial, and the number of colum should be equal to a.
#' @param ConditionalPower Methods for conditional power. If ConditionalPower = "planned", then conditional power will be calculated using the planned treatment effect (theta). If ConditionalPower = "observed", then conditional power will be calcuated using the maximum observed treatment effect at interim analysis.
#' @param ssr Whether sample size re-estimation should be performed. If ssr = 0, sample size re-estimation will not be performed. If ssr = 1, sample size will be performed if conditional power is lower than 1-beta.
#' @param delta The minimum treatment effect that is of interest. The default is NULL, but when sel_rule = 3 or 4, a numeric value has to be specified.
#' @param side Side for conditional power (default is "upper"). Either "upper" or "lower".
#' @details This function calculates conditional power and perform sample size re-estimation (SSR) based on the conditional power approach for two-stage adaptive multi-arm design with normal outcomes. It is assumed that sample size is equally allocated among arms and the observations from each arm follow a normal distribution with a common known variance. At interim analysis, a trial will be stopped if any treatments cross the efficacy boundaries, and any treatments cross the futility boundary will be dropped. Conditional power is calculated under the global alternative hypothesis where it is assumed all the treatment effects are equal to planned treatment effect (theta) or maximum observed treatment effect. If a trial stops at interim analysis, conditional power(CP), new second stage sample size per arm (n2new) and new conditional power(cp_new) will not be calculated and will all be returned as NA. If a trial can continue to the second stage and conditional power is larger than 1-beta, then n2new and cp_new will not be calculated and will be returned as NA. Multiple trials can be assessed at the same time.
#' @return A list containing the following components:
#' \item{EfficacyBound}{The efficacy boundaries for both interim analysis and final analysis.}
#' \item{FutilityBound}{The futility boundaries for both interim analysis and final analysis.}
#' \item{stop_efficacy}{A matrix indicating which arms are stopped for efficacy. 1 indicates stopping for efficacy and 0 otherwise.}
#' \item{stop_futility}{A matrix indicating which arms are stopped for futility. 1 indicates stopping for efficacy and 0 otherwise.}
#' \item{sel_tr}{A matrix indicating which experimental treatments can continue to the second stage.}
#' \item{cp}{Conditional power at interim analysis. If a trial stops at interim analysis, conditional power will not be calculated and returned as NA.}
#' \item{n1}{First stage sample size per arm.}
#' \item{n2new}{New second stage sample size per arm.}
#' \item{cp_new}{Conditional power based on the new second stage sample size.}
#' @references
#' Bauer P, Kohne K: Evaluation of experiments with adaptive interim analyses. Biometrics 1994, 50(4):1029-1041.
#' @references
#' Bretz F, Koenig F, Brannath W, Glimm E, Posch M: Adaptive designs for confirmatory clinical trials. Statistics in Medicine 2009, 28(8):1181-1217.
#' @examples
#'  # two-stage three experiments versus control design, initial sample size was given and all the promising treatments are selected at the interim analysis
#' CT_Normal(combF = "Inverse Normal", a = 3, sel_rule = 2, n1=72, N2=144, initialSS = 1, r21 = NULL, rmax= 1.5, alpha = 0.025, ushape = "obf" ,lfix=0, theta=2, sigma = 6, x0_1=matrix(2.2,nrow=2,ncol=1), x1=matrix(c(2,2.2,3.1,3.7,3.5,2.5),nrow=2,ncol=3), beta=0.2,ConditionalPower = "observed",ssr = 1,delta = NULL)
#' @export

CT_Normal = function(combF,a , sel_rule, n1=NULL, N2=NULL, initialSS = 2, r21, rmax, alpha=0.025, ushape = "obf" ,lfix=0, theta, sigma, x0_1, x1, beta, ConditionalPower, ssr, delta = NULL,side = "upper"){

  if(!require("MAMS")){
    print("trying to install MAMS")
    install.packages("MAMS")
    if(require(MAMS)){
      print("MAMS installed and loaded")
    } else {
      stop("could not install MAMS")
    }
  }
  if (side != "upper" & side != "lower"){
    stop("side should be either upper or lower")
  }
  Bf = function(p1,combF){
    if (combF == "Inverse Normal"){
      (qnorm(1-alpha2) - w1*qnorm_matrix(1-p1))/w2
    } else if (combF == "Fisher"){
      ce1 = alpha2/p1
      ce1[ce1>1]=1
      qnorm_matrix(1 - ce1)
    } else {
      stop("combination function should be Inverse Normal or Fisher")
    }
  }

  ce_ct = function(Calpha, alpha1, alpha0, p1,combF){
    pp1 =  p1* (p1 > 10^(-6)) * (p1 <= 0.9999999) + (p1 <= 10^(-6))*(10^(-6)) + (p1 >= 0.9999999)*0.9999999
    if (combF == "Inverse Normal"){
      ce1 = 1 - pnorm_matrix( (1/w2) * (qnorm(1-Calpha) - w1*qnorm_matrix(1 - pp1)))
    }
    if (combF == "Fisher"){
      ce1 = Calpha/p1
      ce1[ce1>1]=1
    }
    ce = ce1 * (p1 > alpha1) * (p1 < alpha0) + (p1 <= alpha1)
    return(ce)
  }


  z_stat = function (x, y, n, sigma){
    if (class(n)=="numeric"){
      n = matrix(n,nrow = nrow(x), ncol = 1)
    }
    matrix1 = matrix(1,nrow=1,ncol=ncol(x))

    z = ((x-(y%*%matrix1))/(sqrt(2)*sigma)) * (sqrt(n)%*%matrix1)
    if (side == "lower"){
     z = -z
    }
    return(z)
  }


  pnorm_matrix = function (z_matrix){
    if (class(z_matrix)=="numeric"){
      z_matrix = matrix(z_matrix,nrow = wh, ncol = 1)
    }
    probz = NULL
    for (i in 1:nrow(z_matrix)){
      prob1 = pnorm(z_matrix[i,],lower.tail = T)
      probz = rbind(probz,prob1)
    }
    rownames(probz)=NULL
    return(probz)
  }
  qnorm_matrix = function(p_matrix){
    if (class(p_matrix)=="numeric"){
      p_matrix = matrix(p_matrix,nrow = wh, ncol = 1)
    }
    quantile_matrix = NULL
    for (i in 1:nrow(p_matrix)){
      quantile1 = qnorm(p_matrix[i,])
      quantile_matrix = rbind(quantile_matrix,quantile1)
    }
    rownames(quantile_matrix)=NULL
    return(quantile_matrix)
  }

  ### matrix_summary function
  matrix_summary = function(m1,margin=1, method){
    if (match(method,c("sum","max","min","mean"),nomatch = 0) == 0)
      stop("wrong method")

    if (!(margin %in% c(1,2)))
      stop("margin should be 1 or 2")

    if (margin == 1){
      if (method == "sum"){
        ms = matrix(apply(m1,1,sum), ncol=1)
      } else if (method == "max"){
        ms = matrix(apply(m1,1,max), ncol=1)
      } else if (method == "min"){
        ms = matrix(apply(m1,1,min), ncol=1)
      } else {
        ms = matrix(apply(m1,1,mean), ncol=1)
      }
    } else {
      if (method == "sum"){
        ms = matrix(apply(m1,2,sum), nrow=1)
      } else if (method == "max"){
        ms = matrix(apply(m1,2,max), nrow=1)
      } else if (method == "min"){
        ms = matrix(apply(m1,2,min), nrow=1)
      } else {
        ms = matrix(apply(m1,2,mean), nrow=1)
      }
    }
    return(ms)
  }


  DunnettT = function (Z, select = matrix(1, ncol=ncol(Z),nrow = nrow(Z))) {
    if (sum(is.na(Z)) > 0) {
      stop("No missing values allowed in Z")
    }

    treats <- ncol(Z)
    hyp.comb <- list(NULL)
    nhyp.comb <- vector(length = treats)
    for (i in 1:treats) {
      comb.dat <- combn(1:treats, i)
      hyp.comb[[i]] <- comb.dat
      nhyp.comb[i] <- ncol(hyp.comb[[i]])
      rownames(hyp.comb[[i]]) <- 1:i
      hypcol <- NULL
      for (j in 1:nhyp.comb[i]) {
        ihypcol <- paste("H", paste(hyp.comb[[i]][, j], collapse = ""),
                         sep = "")
        hypcol <- append(hypcol, ihypcol)
      }
      colnames(hyp.comb[[i]]) <- hypcol
    }

    pvalue.all = NULL
    for (i1 in 1:nrow(Z)){

      pdunnett.test <- list(NULL)
      zscores <- list(NULL)
      int_dunnett <- function(x, z, k) {
        ((pnorm(sqrt(2) * z + x))^k) * dnorm(x)
      }
      for (i in 1:treats) {
        ptest <- NULL
        if (select[i1,i] == 0) {
          Z[i1,i] <- -Inf
        }
        for (j in 1:nhyp.comb[i]) {
          kselect <- sum(select[i1,c(hyp.comb[[i]][, j])])
          Zmax <- max(Z[i1,c(hyp.comb[[i]][, j])])
          if (kselect == 0) {
            dunnet_integral <- 0
            F_Zmax <- 1
          }
          else {
            dunnett_integral <- integrate(int_dunnett, -Inf,
                                          Inf, z = Zmax, k = kselect)
            F_Zmax <- 1 - dunnett_integral$value
          }
          ptest <- append(ptest, F_Zmax)
        }
        pdunnett.test[[i]] <- matrix(ptest, nrow = 1, ncol = length(ptest))
        colnames(pdunnett.test[[i]]) <- colnames(hyp.comb[[i]])
        rownames(pdunnett.test[[i]]) <- 1
        zscores[[i]] <- qnorm(1 - pdunnett.test[[i]])
      }
      pvalues = matrix(unlist(pdunnett.test),nrow=1)
      pvalue.all = rbind(pvalue.all,pvalues)
    }

    pnames = NULL
    for (i in 1:treats){
      pnames = append(pnames,colnames(hyp.comb[[i]]))
    }
    colnames(pvalue.all)=pnames
    return(pvalue.all)
  }


  p_dunmc = function(z_in,sel_tr) {
    p_dunnett_ih = DunnettT(Z =z_in,select = sel_tr)
    return(p_dunnett_ih)
  }

  goon_rule = function(rule=sel_rule, z1=z1,delta = delta ){
    rowm1 = matrix(1, nrow = 1, ncol = a)
    if (rule == 1){
      sel_tr = ( z1 >= (matrix_summary(z1,1,"max")%*%rowm1) )
    } else if (rule == 2){
      sel_tr = matrix(1, nrow = wh, ncol = a)
    } else if (rule == 3){
      sel_tr = ( z1 >= (matrix_summary(z1,1,"max")%*%rowm1) )
      if (side == "lower"){
        sel_tr = sel_tr * ( (-x1 + (x0_1 %*% rowm1))>=delta )
      } else if (side == "upper"){
        sel_tr = sel_tr * ( (x1 - (x0_1 %*% rowm1))>=delta )
      }
    } else if (rule == 4){
      if (side == "lower"){
        sel_tr = 1 * ( (-x1 + (x0_1 %*% rowm1))>=delta )
      } else if (side == "upper"){
        sel_tr = 1 * ( (x1 - (x0_1 %*% rowm1))>=delta )
      }
    }
    return(sel_tr)
  }

  fisher.alpha1 = function(alpha, alpha0){

    Calpha = exp(-0.5*qchisq(p = 1-alpha, df = 4))
    x_L = Calpha
    x_U = alpha
    f1 = function(Ca,a0,a1,a) {
      a1 + Ca*(log(a0) - log(a1)) - a
    }
    x_m = uniroot(f = f1, interval = c(x_L, x_U), Ca = Calpha,a0 = alpha0, a = alpha,tol = 1e-07)
    alpha1 = x_m$root
    return(c(alpha1, Calpha))
  }

  ce_ct = function(Calpha, alpha1, alpha0, p1,combF = "Inverse Normal"){
    pp1 =  p1* (p1 > 10^(-6)) * (p1 <= 0.9999999) + (p1 <= 10^(-6))*(10^(-6)) + (p1 >= 0.9999999)*0.9999999
    if (combF == "Inverse Normal"){
      ce1 = 1 - pnorm_matrix( (1/w2) * (qnorm(1-Calpha) - w1*qnorm_matrix(1 - pp1)))
    }
    if (combF == "Fisher"){
      ce1 = Calpha/p1
      ce1[ce1>1]=1
    }
    ce = ce1 * (p1 > alpha1) * (p1 < alpha0) + (p1 <= alpha1)
    return(ce)
  }
  macro_ct_ce_mcp = function(combF="Inverse Normal", Calpha, alpha1, alpha0 = alpha0,z1 = z1){
    nc2 = ncol(p1_dunnett_ih)
    p1_dunnett_ih1 = p1_dunnett_ih[,1:a]


    if (a == 2){
      p1_dunnett_ih2 = matrix(p1_dunnett_ih[,nc2],ncol=1)
      colnames(p1_dunnett_ih2) = "H12"
    } else if (a == 3) {
      p1_dunnett_ih2 = matrix(p1_dunnett_ih[,(a+1):(nc2-1)],ncol=3)
      p1_dunnett_ih3 = matrix(p1_dunnett_ih[,nc2],ncol=1)
      colnames(p1_dunnett_ih3) = "H123"
    }

    ce.ct.ih1 = ce_ct(Calpha, alpha1,alpha0,p1 = p1_dunnett_ih1,combF)
    ce.ct.ih2 = ce_ct(Calpha, alpha1,alpha0,p1 = p1_dunnett_ih2,combF)
    ce.ct.ih = list(ce.ct.ih1 = ce.ct.ih1, ce.ct.ih2 = ce.ct.ih2)

    if (a > 2){
      ce.ct.ih3 = ce_ct(Calpha, alpha1,alpha0,p1 = p1_dunnett_ih3,combF)
      ce.ct.ih = list(ce.ct.ih1 = ce.ct.ih1, ce.ct.ih2 = ce.ct.ih2, ce.ct.ih3 = ce.ct.ih3)
    }


    stop_rej = (ce.ct.ih1 >= 1)
    stop_f = (ce.ct.ih1 <=0)

    i3 = 0
    stop_rej_ih2 = matrix(nrow = wh, ncol = a)
    stop_f_ih2 = matrix(nrow = wh, ncol = a)
    if (a ==2){
      stop_rej_ih2 = (ce.ct.ih2 >= 1) %*% matrix(1,nrow=1,ncol=2)
      stop_f_ih2 = (ce.ct.ih2 <= 0) %*% matrix(1,nrow=1,ncol=2)
      stop_rej = stop_rej * stop_rej_ih2
      stop_f = stop_f + stop_f_ih2
    } else if (a == 3){
      ce.ct.ih2.rej = (ce.ct.ih2 >= 1)
      ce.ct.ih2.frej = (ce.ct.ih2 <= 0)
      for (i1 in 1:(a - 1)){
        beg_i2 = i1 + 1
        for (i2 in (i1+1):a){
          i3 = i3 + 1
          stop_rej_ih2[,i3] = ce.ct.ih2.rej[,i1] * ce.ct.ih2.rej[,i2]
          stop_f_ih2[,i3] = ce.ct.ih2.frej[,i1] + ce.ct.ih2.frej[,i2]
        }
      }
      stop_rej = stop_rej * stop_rej_ih2
      stop_f = stop_f + stop_f_ih2
      stop_rej = stop_rej * ( (ce.ct.ih3 >= 1) %*% matrix(1,nrow = 1, ncol = 3) )
      stop_f = stop_f + ( (ce.ct.ih3 <= 0) %*% matrix(1,nrow = 1, ncol = 3) )
    }
    stop_f = (stop_f >= 1) + 0

    nostop = ((stop_f + stop_rej) <= 0) + 0


    St1Decision = list(ct.ce.dunmc.stopf = stop_f, ct.ce.dunmc.stopr = stop_rej, ct.ce.dunmc.nostop = nostop )
    return(list(ce = ce.ct.ih, St1Decision = St1Decision ))
  }

  J = 2
  alpha0 = 1-pnorm(lfix)


  if (is.null(r21)){
    r21 = N2/n1
  }



  plan = mams(K = 1, J=J, alpha = alpha, power = 1-beta, r = c(1,r21), r0 = c(1,r21), p=pnorm(theta/sqrt(2)/sigma),p0=0.5, u.shape = ushape, l.shape = "fixed",lfix = lfix)
  if (initialSS == 1){  #sample size using fixed design
    n1 = ceiling( 2*( (qnorm(1-alpha)+qnorm(1-beta))*sigma/theta )^2 / r21)
    N2 = ceiling(n1*r21)
  } else if (initialSS == 2){
    n1 = plan$n
    N2 = ceiling(n1*r21)
  }

  r21 = N2/n1
  n2 = N2 - n1
  N2max = ceiling(rmax*N2)
  rN2max =  N2max/n1
  n2max = N2max - n1


  w1 = sqrt(1/r21)
  w2 = sqrt(1 - 1/r21)


  if (combF == "Inverse Normal"){
    alpha1 = 1-pnorm(plan$u[1])
    alpha2 = 1-pnorm(plan$u[2])
  } else if (combF == "Fisher"){
    Fi.alpha = fisher.alpha1(alpha, alpha0)
    alpha1 = Fi.alpha[1]
    alpha2 = Fi.alpha[2]
  } else{
    stop("combF can only be Inverse Normal or Fisher ")
  }



  z1 = z_stat(x=x1,y=x0_1,n=n1,sigma=sigma)
  p1 = pnorm_matrix(-z1)
  wh = nrow(z1)
  p1_dunnett_ih = p_dunmc(z_in = z1,sel_tr = matrix(1,nrow = wh,ncol = a))
  S1.decision = macro_ct_ce_mcp(combF=combF, Calpha =alpha2, alpha1 = alpha1, alpha0,z1)

  sel_tr = goon_rule(rule=sel_rule, z1 = z1, delta = delta )
  sel_tr = sel_tr * S1.decision$St1Decision$ct.ce.dunmc.nostop
  sel_tr =  (1 - (((matrix_summary(S1.decision$St1Decision$ct.ce.dunmc.stopr,1,"sum") > 0) %*% matrix(1,nrow=1,ncol=a)) >= 1)) * sel_tr
  sum_sel_tr = matrix_summary(m1 = sel_tr,margin =1,method ="sum")
  i_nostop = which(sum_sel_tr >=1)
  i_stop = which(sum_sel_tr == 0)

  cp_comb = matrix(nrow = wh, ncol = 1)
  cp_new = matrix(nrow=wh, ncol=1)

  if (ConditionalPower == "planned"){
    n2_new_combF = function(ssr=ssr,combF = combF){
      p1ih_allT = p1_dunnett_ih[,ncol(p1_dunnett_ih)]

      B = Bf(p1=matrix(p1ih_allT,nrow=wh,ncol=1),combF=combF)
      if (sum(sum_sel_tr) == 0) {
        n2_new = matrix(0,nrow = wh, ncol = 1)
      } else if (sum(sum_sel_tr) >= 1){
        n2_new = matrix(nrow = wh, ncol = 1)

        if (ssr==0){
          n2_new = n2 * (sum_sel_tr >=1)
          for (i in i_nostop){
            k1 = sum_sel_tr[i,]
            cov1 = matrix(1/2,nrow=k1,ncol=k1)
            diag(cov1)=1
            bi = B[i,]
            cp_comb[i,] = round(1 - pmvnorm(lower=rep(-Inf, k1),upper=rep(bi-theta*sqrt(n2/2)/sigma,k1),mean=rep(0,k1),sigma=cov1)[1],2)
          }
        }else if (ssr == 1){
          for (i in i_nostop){
            k1 = sum_sel_tr[i,]
            cov1 = matrix(1/2,nrow=k1,ncol=k1)
            diag(cov1)=1
            bi = B[i,]
            cp_comb[i,] = round(1 - pmvnorm(lower=rep(-Inf, k1),upper=rep(bi-theta*sqrt(n2/2)/sigma,k1),mean=rep(0,k1),sigma=cov1)[1],2)
            if (cp_comb[i,]>=1-beta){
              n2_new[i,]= n2
            } else {
              if (sum_sel_tr[i,] == 1){
                m2 = 2*(sigma/ theta)^2 * (bi-qnorm(beta))^2
                n2_new[i,] = min(ceiling(m2), n2max)
                cp_new[i,] = round(1-pnorm(bi - sqrt(n2_new[i,]/2)*theta/sigma),2)

              } else if (sum_sel_tr[i,] > 1) {
                f4 = function(x,k1,bi,theta,sigma,cov1){
                  beta - pmvnorm(lower=rep(-Inf, k1),upper=rep(bi-theta*sqrt(x/2)/sigma,k1),mean=rep(0,k1),sigma=cov1)[1]
                }
                cpdiff_n2max = round(f4(x=n2max,k1=k1,bi=bi,theta=theta,sigma=sigma,cov1=cov1),2)

                if (cpdiff_n2max <= 0){
                  n2_new[i,]=n2max
                  cp_new[i,] = cpdiff_n2max + 1 - beta
                }else{
                  root1 = NULL
                  try(root1 <- ceiling(uniroot(f=f4,interval=c(n2,n2max),k1=k1,bi=bi,theta=theta,sigma=sigma,cov1=cov1,tol=1)$root),silent=T)
                  if (is.null(root1)){
                    n2_newi = NULL
                    for (m2new in ((n2max-1):n2)){
                      if (m2new == n2){
                        n2_newi = n2+1
                        break
                      } else {
                        cpdiff = round(f4(x=m2new,k1=k1,bi=bi,theta=theta,sigma=sigma,cov1=cov1),2)
                        if (cpdiff<0){
                          n2_newi=m2new+1
                          break
                        }else if(cpdiff==0){
                          n2_newi=m2new
                          break
                        }else {
                          next
                        }
                      }
                    }
                    n2_new[i,]=n2_newi
                  } else {
                    n2_new[i,] = root1
                  }
                  cp_new[i,] = round(f4(x=n2_new[i,],k1=k1,bi=bi,theta=theta,sigma=sigma,cov1=cov1),2) + 1 - beta
                }
              }
            }
          }
          n2_new[i_stop,]=0
        }
      }
      return(list(n2_new = n2_new,cp = cp_comb, cp_new =cp_new))

    }
  } else if (ConditionalPower == "observed") {
    n2_new_combF = function(ssr=ssr,combF = combF){
      p1ih_allT = p1_dunnett_ih[,ncol(p1_dunnett_ih)]

      B = Bf(p1=matrix(p1ih_allT,nrow=wh,ncol=1),combF=combF)
      Otheta = matrix_summary(z1*sel_tr,1,"max")*sigma*sqrt(2/n1)

      if (sum(sum_sel_tr) == 0) {
        n2_new = matrix(0,nrow = wh, ncol = 1)
      } else if (sum(sum_sel_tr) >= 1){
        n2_new = matrix(nrow = wh, ncol = 1)

        if (ssr==0){
          n2_new = n2 * (sum_sel_tr >=1)
          for (i in i_nostop){
            k1 = sum_sel_tr[i,]
            cov1 = matrix(1/2,nrow=k1,ncol=k1)
            diag(cov1)=1
            cp_comb[i,] = round(1 - pmvnorm(lower=rep(-Inf, k1),upper=rep(B[i,]-Otheta[i,]*sqrt(n2/2)/sigma,k1),mean=rep(0,k1),sigma=cov1)[1],2)
          }
        }else if (ssr == 1){
          for (i in i_nostop){
            k1 = sum_sel_tr[i,]
            cov1 = matrix(1/2,nrow=k1,ncol=k1)
            diag(cov1)=1
            cp_comb[i,] = round(1 - pmvnorm(lower=rep(-Inf, k1),upper=rep(B[i,]-Otheta[i,]*sqrt(n2/2)/sigma,k1),mean=rep(0,k1),sigma=cov1)[1],2)
            if (cp_comb[i,]>=1-beta){
              n2_new[i,]= n2
            } else {
              if (sum_sel_tr[i,] == 1){
                m2 = 2*(sigma/ Otheta[i,])^2 * (B[i,]-qnorm(beta))^2
                n2_new[i,] = min(ceiling(m2), n2max)
                cp_new[i,] = round(1-pnorm(B[i,] - sqrt(n2_new[i,]/2)*Otheta[i,]/sigma),2)

              } else if (sum_sel_tr[i,] > 1) {
                bi = B[i,]

                f4 = function(x,k1,bi,theta,sigma,cov1){
                  beta - pmvnorm(lower=rep(-Inf, k1),upper=rep(bi-theta*sqrt(x/2)/sigma,k1),mean=rep(0,k1),sigma=cov1)[1]
                }

                cpdiff_n2max = round(f4(x=n2max,k1=k1,bi=bi,theta=Otheta[i,],sigma=sigma,cov1=cov1),2)
                if (cpdiff_n2max <= 0){
                  n2_new[i,]=n2max
                  cp_new[i,] = cpdiff_n2max + 1 - beta
                }else{
                  root1 = NULL
                  try(root1 <- ceiling(uniroot(f=f4,interval=c(n2,n2max),k1=k1,bi=bi,theta=Otheta[i,],sigma=sigma,cov1=cov1,tol=1)$root),silent=T)
                  if (is.null(root1)){
                    n2_newi = NULL
                    for (m2new in ((n2max-1):n2)){
                      if (m2new == n2){
                        n2_newi = n2+1
                        break
                      } else {
                        cpdiff = round(f4(x=m2new,k1=k1,bi=bi,theta=Otheta[i,],sigma=sigma,cov1=cov1),2)
                        if (cpdiff<0){
                          n2_newi=m2new+1
                          break
                        }else if(cpdiff==0){
                          n2_newi=m2new
                          break
                        }else {
                          next
                        }
                      }
                    }
                    n2_new[i,]=n2_newi

                  } else {
                    n2_new[i,] = root1
                  }
                  cp_new[i,] = round(f4(x=n2_new[i,],k1=k1,bi=bi,theta=Otheta[i,],sigma=sigma,cov1=cov1),2)
                }
              }
            }
          }
          n2_new[i_stop,]=0
        }
      }
      return(list(n2_new = n2_new,cp = cp_comb, cp_new = cp_new))
    }
  }
  CT_ssr = n2_new_combF(ssr = ssr, combF = combF)

  output = list(EfficacyBound = c(alpha1,alpha2),FutilityBound = alpha0, stop_efficacy= S1.decision$St1Decision$ct.ce.dunmc.stopr, stop_futility= S1.decision$St1Decision$ct.ce.dunmc.stopf, sel_tr = sel_tr,cp = CT_ssr$cp, n1 = n1,  n2new = CT_ssr$n2_new, cp_new = CT_ssr$cp_new)
  return(output)
}
