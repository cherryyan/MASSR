#' SSR for flexible group sequential multi-arm design with normal outcomes
#'
#' Calculate conditional power and perform SSR at the interim analysis for flexible group sequential multi-arm design with normal outcomes. If conditional power is less than desired, perform SSR and calculate the new conditional power based on the new sample size.
#' @param a Number of experimental treatments (a = 2 or 3)
#' @param sel_rule Treatment selection rule at interim analysis. There are four treatment selection rule, 1) sel_tr = 1, select the best treatment; 2) sel_tr=2, keep all promising treatments that do not cross either efficacy boundaries or futility boundaries; 3) sel_tr = 3, select the best treatment only if the mean difference to control is above a threshold delta; 4) sel_tr = 4, keep all promising treatments where the mean difference to control is above delta.
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
#' @param ConditionalPower Methods for conditional power. If ConditionalPower = "planned", then conditional power will be calculated using the planned treatment effect (theta). If ConditionalPower = "observed", then conditional power is calcuated using the maximum observed treatment effect at interim analysis.
#' @param ssr Whether sample size re-estimation should be performed. If ssr = 0, sample size re-estimation will not be performed. If ssr = 1, sample size will be performed if conditional power is lower than 1-beta.
#' @param delta The minimum treatment effect that is of interest. The default is NULL, but when sel_rule = 3 or 4, a numeric value has to be specified.
#' @param side Side for conditional power (default is "upper"). Either "upper" or "lower".
#' @details This function calculates conditional power and perform sample size re-estimation (SSR) based on the conditional power approach for two-stage flexible group sequential multi-arm study design (as described in Magirr et al (2014)) with normal outcomes. It is assumed that sample size is equally allocated among arms and the observations from each arm follow a normal distribution with a common known variance. At interim analysis, a trial will be stopped if any treatments cross the efficacy boundaries and any treatments cross the futility boundary will be dropped. Conditional power is calculated under the global alternative hypothesis where it is assumed all the treatment effects are equal to planned treatment effect (theta) or maximum observed treatment effect. If a trial stops at interim analysis, conditional power(CP), new second stage sample size per arm (n2new) and new conditional power(cp_new) will not be calculated and will all be returned as NA. If a trial can continue to the second stage and conditional power is larger than 1-beta, then n2new and cp_new will not be calculated and will be returned as NA. Multiple trials can be assessed at the same time.
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
#' \item{NoRoot}{If new sample size is successful found, then NoRoot = NULL. If new sample size cann't be found, then NoRoot will return the index of the trial for which the new sample size is not found.}
#' @references
#' Magirr D, Jaki T, Whitehead J: A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection. Biometrika 2012, 99(2):494-501.
#' @references
#' Magirr D, Stallard N, Jaki T: Flexible sequential designs for multi-arm clinical trials. Statistics in Medicine 2014, 33(19):3269-3279.
#' @examples
#' two-stage three experiments versus control design, initial sample size was calculated based on group sequential design and the best treatment is selected at the interim analysis
#' FGS_Normal(a = 3, sel_rule = 2, n1=NULL, N2=NULL, initialSS = 2, r21 = 2, rmax= 1.5, alpha = 0.025, ushape = "obf", lfix=0, theta=2, sigma = 6, x0_1=matrix(2.2,nrow=2,ncol=1), x1=matrix(c(2,2.2,3.1,3.0,3.5,2.5), nrow=2,ncol=3), beta=0.2, ConditionalPower = "observed",ssr = 1,delta = NULL)
#' @export
FGS_Normal = function(a, sel_rule, n1=NULL, N2=NULL, initialSS = 2, r21, rmax, alpha = 0.025, ushape = "obf" ,lfix=0, theta, sigma, x0_1, x1, beta, ConditionalPower, ssr, delta = NULL, side = "upper"){

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
    if (rule == 1){      #select the best
      sel_tr = ( z1 >= (matrix_summary(z1,1,"max")%*%rowm1) )
    } else if (rule == 2){  #select all promising
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


  J = 2

  early_test = 1
  alpha0 = 1-pnorm(lfix)

  rowm1 = matrix(1,nrow=1,ncol=a)
  if (is.null(r21)){
    r21 = N2/n1
  }


  plan = mams(K = 1, J=J, alpha = alpha, power = 1-beta, r = c(1,r21), r0 = c(1,r21), p=pnorm(theta/sqrt(2)/sigma),p0=0.5, u.shape = ushape, l.shape = "fixed",lfix = lfix)
  if (initialSS == 1){
    n1 = ceiling( 2*( (qnorm(1-alpha)+qnorm(1-beta))*sigma/theta )^2 / r21)
    N2 = ceiling(n1*r21)
  } else if (initialSS == 2){
    n1 = plan$n
    N2 = ceiling(n1*r21)
  }
  r21 = N2/n1
  n2 = N2 - n1
  N2max = ceiling(rmax*N2)
  rN2max = N2max/n1
  n2max = N2max - n1




  z1 = z_stat(x=x1,y=x0_1,n=n1,sigma=sigma)
  p1 = pnorm_matrix(-z1) #first stage p-value
  wh = nrow(z1)
  p1_dunnett_ih = p_dunmc(z_in = z1,sel_tr = matrix(1,nrow = wh,ncol = a)) #return a matrix with every column contains the p-value for every intersection hypothesis

  if (sel_rule == 1| sel_rule ==3){
    selection ="select_best"
  } else if (sel_rule == 2|sel_rule ==4){
    selection = "all_promising"
  }



  tSDM = step_down_mams(nMat = matrix(plan$n * plan$rMat[1,],nrow=J,ncol=a+1), alpha_star = plan$alpha_star, lb=lfix,selection = selection)
  u = matrix(unlist(tSDM$u),nrow=2^a-1,ncol=J,byrow = T)
  u1 = u[,1]
  u2 = u[,2]
  l1 = tSDM$l[[1]][1]
  if (a == 2){
    u1_1 = u1[c(1,2)]
    u1_2 = u1[3]
  } else if (a == 3){
    u1_1 = u1[c(1,2,4)]
    u1_2 = u1[c(3,5,6)]
    u1_3 = u1[7]
    z1_ih3 = matrix_summary(z1,1,"max")
  }

  stop_rej = (z1 >= u1_1[1] )
  stop_f = (z1 <= l1 )


  z_ih2f = function(z){
    z_ih2 = NULL
    for (i1 in 1:(a - 1)){
      for (i2 in (i1+1):a){
        zmax2 = matrix_summary(matrix(z[,c(i1,i2)],ncol=2),1,"max") #get test statistics for intersection hypothesis invovle 2 T arms
        z_ih2 = cbind(z_ih2, zmax2)
      }
    }
    return(z_ih2)
  }

  z1_ih2 = z_ih2f(z1)

  i3 = 0
  stop_rej_ih2 = matrix(nrow = wh, ncol = a)
  stop_f_ih2 = matrix(nrow = wh, ncol = a)
  if (a == 2){
    stop_rej_ih2 = (z1_ih2 >= u1_2[1]) %*%  rowm1
    stop_f_ih2 = (z1_ih2 <= l1) %*% rowm1
    stop_rej = stop_rej * stop_rej_ih2
    stop_f = stop_f + stop_f_ih2
  } else if (a == 3){
    ih2.rej = (z1_ih2 >= u1_2[1])
    ih2.frej = (z1_ih2 <= l1)
    for (i1 in 1:(a-1)){
      for (i2 in (i1+1):a){
        i3 = i3+1
        stop_rej_ih2[,i3] = ih2.rej[,i1] * ih2.rej[,i2]
        stop_f_ih2[,i3] = ih2.frej[,i1] + ih2.frej[,i2]
      }
    }
    stop_rej = stop_rej * stop_rej_ih2
    stop_f = stop_f + stop_f_ih2
    stop_rej = stop_rej * ( (z1_ih3 >= u1_3) %*% matrix(1,nrow = 1, ncol = 3) )
    stop_f = stop_f + ( (z1_ih3 <= l1) %*% matrix(1,nrow = 1, ncol = 3) )
  }
  stop_f = (stop_f >= 1) + 0

  nostop = ((stop_f + stop_rej) <= 0) + 0


  sel_tr = goon_rule(rule=sel_rule, z1 = z1, delta = delta )
  sel_tr = sel_tr*nostop
  sel_tr =  (1 - (((matrix_summary(stop_rej,1,"sum") > 0) %*% matrix(1,nrow=1,ncol=a)) >= 1)) * sel_tr
  sum_sel_tr = matrix_summary(m1 = sel_tr,margin =1,method ="sum")
  i_nostop = which(sum_sel_tr >=1)
  i_stop = which(sum_sel_tr == 0)

  if (ConditionalPower == "planned"){
    N2_newF = function(ssr=ssr,rN2max =rN2max){
      u2last = u2[2^a - 1]

      rejrH0 = function(u2last=u2last,z1k=z1k,r21=r21){
        (u2last - z1k*sqrt(1/r21))/sqrt(1-1/r21)
      }

      rejrHA = function (u2last=u2last,z1k=z1k,r21=r21){
        (u2last - z1k*sqrt(1/r21)-theta*(r21-1)*sqrt(n1/2/r21)/sigma)/sqrt(1-1/r21)
      }
      if (sum(sum_sel_tr) == 0) {
        N2_new = n1
        cp = NA
        NoRoot = NA
        cp_new = NA
      } else if (sum(sum_sel_tr) >= 1){
        N2_new = matrix(nrow=wh, ncol=1)
        cp = matrix(nrow = wh, ncol = 1)
        if (ssr==0){
          N2_new = N2 * (sum_sel_tr >=1)
          for (i in i_nostop){
            z1k = z1[i,which(sel_tr[i,]==1)]
            rejr = rejrHA(u2last,z1k,r21)
            b1 = length(rejr)
            cov1 = matrix(1/2,nrow = b1,ncol = b1)

            diag(cov1) = 1
            cp[i,] = round(1-pmvnorm(lower=rep(-Inf, b1),upper = rejr,mean=rep(0,b1), sigma=cov1)[1],2)
          }
        }else if (ssr == 1){
          cp_new = matrix(nrow = wh, ncol = 1)
          for (i in i_nostop){
            z1k = z1[i,which(sel_tr[i,]==1)]
            rejr = rejrHA(u2last,z1k,r21)
            b1 = length(rejr)
            cov1 = matrix(1/2,nrow = b1,ncol = b1)

            diag(cov1) = 1
            cp[i,] = round(1-pmvnorm(lower=rep(-Inf, b1),upper = rejr,mean=rep(0,b1), sigma=cov1)[1],2)
            if (cp[i,]>= 1-beta){
              N2_new[i,] = N2
            }else {
              zce =rejrH0(u2last, z1k, r21)
              if (sum_sel_tr[i,] == 1){
                m2 = 2*((zce - qnorm(beta))*sigma/theta)^2 + n1
                N2_new[i,] = min(ceiling(m2),N2max)
                cp_new[i,] = round(1-pnorm(zce - sqrt((N2_new[i,] - n1)/2)*theta/sigma),2)
              } else if (sum_sel_tr[i,] > 1){
                ce = 1- pmvnorm(lower=rep(-Inf, b1),upper = rejrH0(u2last,z1k,r21),mean=rep(0,b1), sigma=cov1)[1]

                ff1 = function(x,b1,z1k,n1,cov1,ce,r21new){
                  1 - ce - pmvnorm(lower=rep(-Inf, b1),upper = rejrH0(x,z1k,r21new),mean=rep(0,b1), sigma=cov1)[1]
                }

                ff2 = function(y){
                  u2last_new = uniroot(ff1,c(0,10),b1=b1,z1k=z1k,n1=n1,cov1=cov1,ce=ce,r21new=y/n1)$root
                  CPdiff = beta - pmvnorm(lower=rep(-Inf, b1),upper = rejrHA(u2last_new,z1k,r21 = y/n1), mean=rep(0,b1), sigma=cov1)[1] #beta - CP = 0
                  return(CPdiff)
                }
                cpdiff_N2max = round(ff2(N2max),2)
                if (cpdiff_N2max <= 0){
                  N2_new[i,] = N2max
                  cp_new[i,] = cpdiff_N2max + 1 - beta
                } else {
                  root1 = NULL
                  try(root1 <- ceiling(uniroot(ff2,c(N2,N2max),tol=1)$root),silent=T)
                  if (is.null(root1)){
                    N2_newi = NULL
                    for (M2new in ((N2max-1):N2)){
                      if(M2new == N2){
                        N2_newi = N2+1
                        break
                      }else{
                        cpdiff = round(ff2(M2new),2)
                        if (cpdiff<0){
                          N2_newi=M2new+1
                          break
                        }else if(cpdiff==0){
                          N2_newi=M2new
                          break
                        }else {
                          next
                        }
                      }
                    } #end of for
                    if(is.null(N2_newi)){
                      NoRoot = c(NoRoot,i)
                    }else{
                      N2_new[i,]=N2_newi
                    }

                  }else {
                    N2_new[i,] = root1
                  }
                  cp_new[i,] = round(ff2(N2_new[i,]),2) + 1 - beta
                }
              }
            }
          }
        }
        N2_new[i_stop,] = n1
      }

      return(list(N2_new = N2_new, cp = cp,cp_new = cp_new, NoRoot))
    }
  }else if (ConditionalPower == "observed") {
    N2_newF = function(ssr=ssr,rN2max =rN2max){
      u2last = u2[2^a - 1]
      Otheta = matrix_summary(z1*sel_tr,1,"max")*sigma*sqrt(2/n1)

      rejrH0 = function(u2last=u2last,z1k=z1k,r21=r21){
        (u2last - z1k*sqrt(1/r21))/sqrt(1-1/r21)
      }

      rejrHA = function (u2last=u2last,z1k=z1k,r21=r21,theta){
        (u2last - z1k*sqrt(1/r21)-theta*(r21-1)*sqrt(n1/2/r21)/sigma)/sqrt(1-1/r21)
      }

      if (sum(sum_sel_tr) == 0) {
        N2_new = n1
        cp = NA
        NoRoot = NA
        cp_new = NA
      } else if (sum(sum_sel_tr) >= 1){

        N2_new = matrix(nrow=wh, ncol=1)
        cp = matrix(nrow = wh, ncol = 1)
        if (ssr==0){
          N2_new = N2 * (sum_sel_tr >=1)
          for (i in i_nostop){
            z1k = z1[i,which(sel_tr[i,]==1)]
            rejr = rejrHA(u2last,z1k,r21,theta= Otheta[i,])
            b1 = length(rejr)
            cov1 = matrix(1/2,nrow = b1,ncol = b1)

            diag(cov1) = 1
            cp[i,] = round(1-pmvnorm(lower=rep(-Inf, b1),upper = rejr,mean=rep(0,b1), sigma=cov1)[1],2)
          }

        }else if (ssr == 1){
          cp_new = matrix(nrow=wh, ncol=1)
          for (i in i_nostop){
            z1k = z1[i,which(sel_tr[i,]==1)]
            rejr = rejrHA(u2last,z1k,r21,theta= Otheta[i,])
            b1 = length(rejr)
            cov1 = matrix(1/2,nrow = b1,ncol = b1)

            diag(cov1) = 1
            cp[i,] = round(1-pmvnorm(lower=rep(-Inf, b1),upper = rejr,mean=rep(0,b1), sigma=cov1)[1],2)
            if (cp[i,]>= 1-beta){
              N2_new[i,] = N2

            }else {
              zce =rejrH0(u2last, z1k, r21)
              if (sum_sel_tr[i,] == 1){
                m2 = 2*((zce - qnorm(beta))*sigma/Otheta[i,])^2 + n1
                N2_new[i,] = min(ceiling(m2),N2max)
                cp_new[i,] = round(1- pnorm(zce - sqrt((N2_new[i,] - n1)/2)*Otheta[i,]/sigma),2)
              } else if (sum_sel_tr[i,] > 1){
                ce = 1- pmvnorm(lower=rep(-Inf, b1),upper = rejrH0(u2last,z1k,r21),mean=rep(0,b1), sigma=cov1)[1]

                ff1 = function(x,b1,z1k,n1,cov1,ce,r21new){
                  1 - ce - pmvnorm(lower=rep(-Inf, b1),upper = rejrH0(x,z1k,r21new),mean=rep(0,b1), sigma=cov1)[1]
                }

                ff2 = function(y,theta){
                  u2last_new = uniroot(ff1,c(0,10),b1=b1,z1k=z1k,n1=n1,cov1=cov1,ce=ce,r21new=y/n1)$root
                  CPdiff = beta - pmvnorm(lower=rep(-Inf, b1),upper = rejrHA(u2last_new,z1k,r21 = y/n1,theta), mean=rep(0,b1), sigma=cov1)[1]
                  return(CPdiff)
                }
                cpdiff_N2max = round(ff2(N2max,theta=Otheta[i,]),2)
                if (cpdiff_N2max <= 0){
                  N2_new[i,] = N2max
                  cp_new[i,] = cpdiff_N2max + 1 - beta
                } else {
                  root1 = NULL
                  try(root1 <- ceiling(uniroot(ff2,c(N2,N2max),theta=Otheta[i,],tol=1)$root),silent=T)
                  if (is.null(root1)){
                    N2_newi = NULL
                    for (M2new in ((N2max-1):N2)){
                      if(M2new == N2){
                        N2_newi = N2+1
                        break
                      }else{
                        cpdiff = round(ff2(M2new,theta=Otheta[i,]),2)
                        if (cpdiff<0){
                          N2_newi=M2new+1
                          break
                        }else if(cpdiff==0){
                          N2_newi=M2new
                          break
                        }else {
                          next
                        }
                      }
                    } #end of for
                    if(is.null(N2_newi)){
                      NoRoot = c(NoRoot,i)
                    }else{
                      N2_new[i,]=N2_newi
                    }

                  }else {
                    N2_new[i,] = root1
                  }
                  cp_new[i,] = round(ff2(N2_new[i,], theta=Otheta[i,]),2) + 1 - beta

                }
              }
            }
          }
        }
        N2_new[i_stop,] = n1
      }

      return(list(N2_new = N2_new, cp = cp,cp_new = cp_new, NoRoot))
    }
  }

  NoRoot = NULL
  N2_ssr = N2_newF(ssr=ssr,rN2max =rN2max)

  output = list(EfficacyBound = tSDM$u, FutilityBound = tSDM$l, stop_efficacy= stop_rej, stop_futility= stop_f, sel_tr = sel_tr, cp = N2_ssr$cp,n1 = n1,n2new = N2_ssr$N2_new - n1, cp_new = N2_ssr$cp_new, NoRoot = N2_ssr$NoRoot)
  return(output)
}

