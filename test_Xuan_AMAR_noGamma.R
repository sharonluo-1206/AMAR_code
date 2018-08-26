

#C.LIB.NAME <- "Xuan_AMAR_bictool_backup_v2.so";
#dyn.load(C.LIB.NAME);

library("mnormt");
library("nleqslv");
library("MCMCpack");
library("mvtnorm");
library("numDeriv");
library("rootSolve");


###########################################################################################
###### the pdf of dirichlet multinomial distribution update on 08/23/2015
dirmulti <- function(eta, n, xi) {
     nttl <- sum(n);
     lengthn <- length(n);
     rst <- factorial(nttl)/prod(factorial(n))*gamma(1/eta)/gamma(nttl+1/eta)*prod(sapply(1:lengthn, function(x) gamma(n[x]+xi[x]/eta)/gamma(xi[x]/eta)));
     rst;
}

# find the empirical Bayes estimate of lambda
# optimize(dirmulti, c(0,100), n=emp.v.n, xi=ksi, maximum=TRUE)


#############################
### example for using this file

xuan.all <- function(test){
    SIMU.MARGIN  <- 3;
    SIMU.SIZE <- 100;

    ##values of V
    V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1),
                        c(0,1), c(0,1), c(0,1));
    ##cutoffs of V multinormal setup
    V.CUTS <- list(-0.2, -0.3, -0.4, -0.3, -0.1,
                   -0.4, -0.3, -0.3);
    ##multivariate normal variance matrix with equal correlation
    V.SIGMA       <- matrix(0.4, length(V.CATS), length(V.CATS));
    diag(V.SIGMA) <- 1;
    ##set V consts
    V.INX.ALL     <- get.v.all(V.CATS);
    V.PROB.ALL.MN <- get.prob.v.mnorm(V.SIGMA, V.CUTS);
    ##Parameters≈
    NTIME <- 3; #not including t=0 baseline    
    ### the transitonal model parameters
    BETA  <- c(0.5, 0.25, 0.3);
    ### the imputation model parameters
    ALPHA <- matrix(c(0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7,0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7),8,NTIME);
    GAMMA <- 0.3;
    ### the MDM parameters, use 4 v
    phi.v <- c(0.6, 0.7, 0.5, 0.4, 0, 0, 0, 0);
    PHI <- c(-3.5, phi.v, 0.5, 0);
    ##get delta for simulation
    DELTA.ALL<- get.Deltat(NTIME, ALPHA, BETA, GAMMA, V.INX.ALL, V.PROB.ALL.MN, SIMU.MARGIN);
    DELTA.ALL;

    ##-----------------------------------------------------------------------------------------------
    ##           SIMULATION
    ##-----------------------------------------------------------------------------------------------
    v     <- simu.V(SIMU.SIZE, V.INX.ALL, V.PROB.ALL.MN);
    ### simulate the Y based on the imputation model
    ### then simulate the R, and basing on R set corresponding Y as missing
    y.full     <- simu.Y(DELTA.ALL, ALPHA, GAMMA, v);
    ### in this AMAR case, the SIMU.R.RST is a matrix with dim(nsize, ntime), 1 as missing, 0 as observed
    r.ind <- simu.R(y.full, v, PHI);
    ### gives how many Y are observed for each subjects
    r <- NTIME-sapply(1:SIMU.SIZE, function(x) sum(r.ind[x,]));
    ### The observed Y outcome aftering setting NA based on R
    y.obs  <- set.missing.y(y.full, r.ind);
    ### empirical counts of every possible v
    emp.v.n <- get.v.inx(v, V.CATS, counts=TRUE);
    write.table(cbind(v, r, y.obs), file=paste("Size", SIMU.SIZE, "data1_vry.txt", sep=""));

    ##-----------------------------------------------------------------------------------------------
    ##           Theta Estimation
    ##-----------------------------------------------------------------------------------------------
    # for estimation of theta
    probv.NITER <- (10)^4;
    probv.NBURN <- (10)^3;
    eta.init <- 1;
    eta.lowbound <- 0;
    eta.upbound <- 10;
    
    prvdan.rst <- prv.dan(probv.NITER, probv.NBURN, v, v.log.fx, V.CATS, SIMU.SIZE, true.theta);
    ### objects(pvdan.rst) == all.theta, mean.theta, eta, mse.theta, pchi.theta

    prvslice.rst <- prv.slice(V.INX.ALL, emp.v.n, prob.NITER, probv.NBURN, eta.init, eta.lowbound, eta.upbound, true.theta);
    ### objects(prvslice.rst) ==  all.theta, mean.theta, mse.theta, pchi.theta, eta
    ### there are also prv.slice.v4, prv.slice.v5, prv.slice.v7 for different dimension of V

    ### This method is not good choice, better not use. It gains small improvement with hugh computation time
    prvslice2way.rst <- prv.slice.new.twoway(emp.v.n, V.INX.ALL, probv.NITER, probv.NBURN, inits=NULL, eta.upbound, sigma.lambda, V.CATS, true.theta);
    ### objects(prv.slice.new.twoway) == all.theta, mean.theta, mse.theta, pchi.theta, all.eta

    prvdir.rst <- xuan.pvdirichlet.fast.1127(k,data.i, probv.NITER, probv.NBURN, max.h, a.alpha, b.alpha, V.INX.ALL, V.CATS, true.theta);
    #### objects(xuan.pvdirichlet.fast.1127) == theta, n.hitmax
    ### for more information returned, check xuan.pvdirichlet.fast.1002, but note there are 2 version of theta which is going to take a huge physical memory. If you use it, make necessary modification
    
    ##-----------------------------------------------------------------------------------------------
    ##           Standard Uncongenial Analysis
    ##-----------------------------------------------------------------------------------------------
    impt.par <- uncong.impt(v, r, y.obs);
    est.Alpha <- impt.par$est.Alpha;
    est.Gamma <- impt.par$est.Gamma;
    est.Delta <- impt.par$est.Delta;
    uncong.std.rst <- uncong.trans.0630(nimpt, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v, r);
    #rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   
    ##-----------------------------------------------------------------------------------------------
    ##           Empirical Uncongenial Analysis
    ##-----------------------------------------------------------------------------------------------


}



### define expit function
expit <- function(x) {
   exp(x)/(1+exp(x));
 }

logit <- function(x) {
   log(x/(1-x));
 }


### based on BETA, calculate the transitional probability
get.prob.5 <- function(x) {
     rst <- numeric(5);
     rst[1] <- expit(x[1]);
     rst[2] <- expit(x[1]+x[2]);
     rst[3] <- expit(x[1]+x[2]+x[3]);
     rst[4] <- expit(x[1]+2*x[2]);
     rst[5] <- expit(x[1]+2*x[2]+x[3]);
     rst;
}


### get the standard error for mle p1-p5 based on mle beta[1:3] add on 6/19/2015
### AMAR case, matching the imputation model and marginal conditional model
get.se.p <- function(beta, hes) {
 
     var.beta <- (-1)*solve(hes);
     
     eta <- numeric(5);
     eta[1] <- beta[1];
     eta[2] <- sum(beta[1:2]);
     eta[3] <- sum(beta[1:3]);
     eta[4] <- beta[1]+2*beta[2];
     eta[5] <- beta[1]+2*beta[2]+beta[3];    

     var.eta <- numeric(5);
     var.eta[1] <- var.beta[1,1];
     var.eta[2] <- var.beta[1,1]+var.beta[2,2]+2*var.beta[1,2];
     var.eta[3] <- var.beta[1,1]+var.beta[2,2]+var.beta[3,3]+2*var.beta[1,2]+2*var.beta[1,3]+2*var.beta[2,3];
     var.eta[4] <- var.beta[1,1]+4*var.beta[2,2]+4*var.beta[1,2];
     var.eta[5] <- var.beta[1,1]+4*var.beta[2,2]+var.beta[3,3]+2*var.beta[1,3]+4*var.beta[1,2]+4*var.beta[2,3];

     var.p <- numeric(5);
     var.p <- sapply(1:5, function(x) var.eta[x]*(exp(eta[x]/2)+exp(-1*eta[x]/2))^(-4));
     se.p <- sqrt(var.p);

     se.p
}







#### Calculate the Pearson's Chi-Square, used to measure the accuracy of theta estimation
probv.pchi <- function(theta, true.theta) {
  ntheta <- length(theta);
  p.chi <- sum(sapply(1:ntheta, function(x)  (true.theta[x]-theta[x])^2/true.theta[x]));
  p.chi
}




### Based on BETA(0,1,2) get the Marginal probability of P(Yt=1) for t=1,2,3 in AMAR case
get.marg.prob <- function(x) {
   rst <- numeric(3);
   rst[1] <- expit(x[1]);
   rst[2] <- expit(x[1]+x[2])*(1-rst[1])+expit(x[1]+x[2]+x[3])*rst[1];
   rst[3] <- expit(x[1]+2*x[2])*(1-rst[2])+expit(x[1]+2*x[2]+x[3])*rst[2];
   rst;
}

### Based on BETA(0,1,2) get the Conditional probability of P(Yt=1|Y_(t-1)=y) for t=1,2,3 and y=0,1 in AMAR case
get.cond.prob <- function(x) {
     rst <- numeric(5);
     rst[1] <- expit(x[1]);
     rst[2] <- expit(x[1]+x[2]);
     rst[3] <- expit(x[1]+x[2]+x[3]);
     rst[4] <- expit(x[1]+2*x[2]);
     rst[5] <- expit(x[1]+2*x[2]+x[3]);
     rst;
}


get.marginv <- function(V.INX.ALL, true.theta) {
   n.v <- ncol(V.INX.ALL);
   p.v1 <- NULL;
   p.v11 <- matrix(0, n.v, n.v);
   p.v01 <- matrix(0, n.v, n.v);
   p.v10 <- matrix(0, n.v, n.v);
   p.v00 <- matrix(0, n.v, n.v);
  
   for(i in 1:n.v) {
      p.v1[i] <- sum(true.theta[which(V.INX.ALL[,i]==1)]); }

   for(i in 1:n.v) {
      for(j in 1:n.v) {
         p.v11[i,j] <- sum(true.theta[intersect(which(V.INX.ALL[,i]==1),which(V.INX.ALL[,j]==1))]);
         p.v01[i,j] <- sum(true.theta[intersect(which(V.INX.ALL[,i]==0),which(V.INX.ALL[,j]==1))]);
         p.v10[i,j] <- sum(true.theta[intersect(which(V.INX.ALL[,i]==1),which(V.INX.ALL[,j]==0))]);
         p.v00[i,j] <- sum(true.theta[intersect(which(V.INX.ALL[,i]==0),which(V.INX.ALL[,j]==0))]);
       }
    }    

}

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
###### Different ways to generate Theta

###  (1) Simulation of Theta Based ON 3-way interaction existance in loglinear model
gen.theta.3wayll <- function(V.CATS, V.INX.ALL, n.v, msd.vec){
   # main factor
   lambda1 <- runif(n.v,0,1);
   
   # second order factor 
   lengthv <- length(lambda1);
   lambda2 <- matrix(0, lengthv, lengthv);
   for(i in 1:(lengthv-1))
   {  for(j in (i+1):lengthv) 
      { lambda2[i,j] <- rnorm(1,msd.vec[1], msd.vec[2]);}
   }
   for(i in 2:lengthv)
   { 
      for(j in 1:(i-1))
      { lambda2[i,j] <- lambda2[j,i];}
   }
   
   # three order factor
   lambda3 <- array(0, dim=c(n.v, n.v, n.v));
   for(i in 1:(n.v-2)) {
    for(j in (i+1):(n.v-1)) {
      for(k in (j+1):n.v) {
         lambda3[i,j,k] <- rnorm(1,msd.vec[3], msd.vec[4]);}
    }
   }  
   
   ###########################
   lambdasum <- NULL;
   muvall <- NULL;
   nvall <- nrow(V.INX.ALL);
   
   for(i in 1:nvall) {
     # find the index of V having value=1
     index1 <- which(V.INX.ALL[i,]==1);
   
     # main factor sum
     oneinter <- sum(lambda1[index1]);
   
     # second order interaction sum
     lengthv1 <- length(index1);
     twointer <- 0;
     if(lengthv1 < 2) { twointer <- 0;}
     else
     {  for(j in 1:(lengthv1-1))
        {
           for(k in (j+1):lengthv1) { twointer=twointer+lambda2[index1[j], index1[k]];}
         }
      }
   
     # three order interaction sum
     #threeinter <- rnorm(1, 0, 1);  
     threeinter <- 0;
     if(lengthv1 < 3) { threeinter <- 0;}
     else
     {
        for(j in 1:(lengthv1-2))
          {  for(k in (j+1):(lengthv1-1))
               { for(l in (k+1):(lengthv1)) {
                   threeinter=threeinter+lambda3[j,k,l]; }
               }
           }   
      }
     lambdasum[i] <- oneinter+twointer+threeinter;  
     muvall[i] <- exp(oneinter+twointer+threeinter);
   }
   
   # return the value of pvall as theta, which sums to 1
   pvall <- muvall/sum(muvall);
   pvall;
 }


### (2) Simulation of Theta Based ON Sum of 2 Independent Multinomial Product
gen.theta.sumIndep <- function(V.INX.ALL, weight, phi) {
  ### weight: length L, gives the weight of each multinomial distribution
  ### phi: row L, col 8, gives the L different versions of 8 marginals for the Indep Multinomial Product 
   theta <- rep(0, nrow(V.INX.ALL));
   
   for(j in 1:nrow(V.INX.ALL)) {
      for(h in 1:nrow(phi)) {
         theta[j] <- theta[j]+weight[h]*prod(sapply(1:ncol(V.INX.ALL), function(x) phi[h,x]^(1-V.INX.ALL[j,x])*(1-phi[h,x])^V.INX.ALL[j,x]));
       }
   }
   #return the Theta value
   theta;
}

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

#### From Chenguang Wang for CTQ data

##get combinations of vall in its raw values
get.v.all <- function(lst.v) {
    ##inx
    inx.v  <- NULL;
    for (i in 1:length(lst.v)) {
        cur.v <- lst.v[[i]];
        if (1 == i) {
            inx.v <- cbind(cur.v);
        } else {
            left.v.inx <- rep(1:nrow(inx.v), each=length(cur.v));
            left.v     <- cbind(inx.v[left.v.inx,]);
            right.v    <- rep(cur.v, nrow(inx.v));
            inx.v      <- cbind(left.v, right.v);
        }
    }
    inx.v
}

##calculate the prob of Vall
##based on multinomial
get.prob.v <- function(lst.pv, eta, standardize=FALSE) {
    prob.v <- 1/eta;

    for (i in 1:length(lst.pv)) {
        cur.p  <- lst.pv[[i]];
        prob.v <- kronecker(prob.v, cur.p);
    }
    if (standardize) prob.v <- prob.v / sum(prob.v);
    prob.v
}


##calculate the prob of Vall based on multivariate normal
##i.e simulate multivariate normal then integrate to get probability with given cutoffs
get.prob.v.mnorm <- function(sigma, lst.cutoffs=NULL, lst.mar.p=NULL) {

    if (is.null(lst.cutoffs)) {
        lst.cutoffs <- mapply(qnorm, lst.mar.p);
    }

    ##boundaries for integration
    lb.v <- NULL;
    ub.v <- NULL;

    for (i in 1:length(lst.cutoffs)) {
        cur.v <- lst.cutoffs[[i]];

        if (1 == i) {
            lb.v <- cbind(c(-Inf, cur.v));
            ub.v <- cbind(c(cur.v, +Inf));
        } else {
            left.v.inx <- rep(1:nrow(lb.v), each=length(cur.v)+1);

            left.lb    <- cbind(lb.v[left.v.inx,]);
            left.ub    <- cbind(ub.v[left.v.inx,]);

            right.lb   <- rep(c(-Inf,cur.v), nrow(lb.v));
            right.ub   <- rep(c(cur.v,+Inf), nrow(ub.v));

            lb.v       <- cbind(left.lb, right.lb);
            ub.v       <- cbind(left.ub, right.ub);
        }
    }

    prob.v <- NULL;
    for (i in 1:nrow(lb.v)){
        cur.p  <- sadmvn(lower = lb.v[i,], upper=ub.v[i,], rep(0, length(lst.cutoffs)), sigma);
        prob.v <- c(prob.v, cur.p);
    }

    prob.v <- prob.v/sum(prob.v);
    prob.v
}

##translate raw v values to index values
##i.e. 10 30 20 to 1,2,3
get.trans.v <- function(raw.v, v.inx) {
    pv      <- length(v.inx);
    nlen    <- nrow(raw.v);

    trans.v <- NULL;
    for (i in 1:pv){
        cur.vinx <- v.inx[[i]];
        cur.v    <- rep(NA, nlen);
        for (j in 1:length(cur.vinx)) {
            cur.v[which(cur.vinx[j] == raw.v[,i])] <- j;
        }
        trans.v <- cbind(trans.v, cur.v);
    }

    trans.v
}


####By Xuan: transfer a realization of (V,R) into a simple value that indicated the position in all realization
### for ease of getting the corresponding prob(V,R) from V.PROB.ALL.MN
### it can only apply to ONE realization of (V,R) at each time
vtosimvalue <- function(v, listcat) {
   value <- 0;
   vlength <- length(v);
   for (i in 1:vlength) {
       if (i==vlength) {
          svalue <- v[vlength];
        } else {
         svalue <- v[i];
         for (j in (i+1):vlength)
         { svalue <- svalue*length(listcat[[j]]);}
       }
   value <- value + svalue;
 }
 value  
 }



###### for AMAR case  Pr(Yt=1|Y_(t-1)=y, v) for y=0, 1 and 1 v
AMAR.imputatemdl <- function(delta, tt, v, alpha.t, gamma) {
  rst <- c(0,0);
  nalpha <- length(alpha.t);

 # cat("alpha=", alpha.t, ", v=", v, ", delta=", delta, "\n");
  ttmp <- sum(sapply(c(1:nalpha), function(x) v[x]*alpha.t[x]));
 # cat("ttmp=", ttmp, "\n");
  if (tt==1) {
      rst[1]= expit(ttmp+delta[1]);
      ### Y0=0 with prob=1, hence rst[2]=0 corresponds to Y0=1
      rst[2] <- 0;}
  else {  
      rst[1] <- expit(ttmp+delta[1]); # delta[t,1] corresponds to prevy=0
      rst[2] <- expit(ttmp+gamma+delta[2]);}  # delta[t,2] corresponds to prevy=1
  rst
}




### AMAR case
### Calculate the E(Y_t|V) from the sum_y P(Y_t|V,Y_(t-1))*P(Y_(t-1)=y|v) for all realization of V 
### prev.pyt.gvn.v = Pr(Y_(t-1)=1|v) for all v
AMAR.eygvnvonimputate <- function(delta, t, vall, prev.pyt.gvn.v, alpha.t, gamma) {
   rst <- array(0);
   vrow <- dim(vall)[1];

 ##  cat("for T=", t, ", delta=", delta, "\n");
 ### for every possible v
   for (i in 1:vrow)
   {
      v <- vall[i,];
      if (t==1) { eygvnvr <- AMAR.imputatemdl(delta, t, v, alpha.t, gamma)[1]; 
      } else {
      ###=(Pr(Yt|Y_(t-1)=0, V),Pr(Yt|Y_(t-1)=1,V))
       cur.pyt <- AMAR.imputatemdl(delta, t, v, alpha.t, gamma);
       ### Pr(Y_(t-1)=1|V)
       cur.prev.pyt <- prev.pyt.gvn.v[i];
       eygvnvr=cur.pyt[1]*(1-cur.prev.pyt)+cur.pyt[2]*cur.prev.pyt; 
     } 
     rst[i] <- eygvnvr;
  }     
rst  
}


###### for AMAR case  Pr(Yt=1|Y_(t-1)=y, v) for y=0, 1 and all v
AMAR.imputatemdl.all <- function(delta, tt, vall, alpha.t, gamma) {
  rst <- matrix(0, nrow(vall), 2);
  nalpha <- length(alpha.t);

  ttmp <- vall%*%alpha.t;
  rst[,1] <- expit(ttmp+delta[1]);
  if(tt>1) rst[,2] <- expit(ttmp+gamma+delta[2]);
  rst
}



AMAR.eygvnvonimputate.fast <- function(delta, t, vall, prev.pyt.gvn.v, alpha.t, gamma) {
   rst <- array(0);
   vrow <- dim(vall)[1];

   eygvnvr.all <- AMAR.imputatemdl.all(delta, t, vall, alpha.t, gamma);
   if(t==1)  {
       rst <- eygvnvr.all[,1];
   }else{
       rst <- eygvnvr.all[,1]*(1-prev.pyt.gvn.v)+eygvnvr.all[,2]*prev.pyt.gvn.v
   }
   rst;
}


############################################################################################################
###   AMAR.eygvnvonimputate.parallel <- function(delta, t, vall, prev.pyt.gvn.v, alpha.t, gamma) {
###      rst <- array(0);
###      vrow <- dim(vall)[1];
###   
###    ##  cat("for T=", t, ", delta=", delta, "\n");
###    ### for every possible v
###      rst <- foreach (i = 1:vrow, .combine='c', .export=c('AMAR.imputatemdl', 'expit', 'logit')) %dopar% {
###         v <- vall[i,];
###         if (t==1) { eygvnvr <- AMAR.imputatemdl(delta, t, v, alpha.t, gamma)[1]; 
###         } else {
###         ###=(Pr(Yt|Y_(t-1)=0, V),Pr(Yt|Y_(t-1)=1,V))
###          cur.pyt <- AMAR.imputatemdl(delta, t, v, alpha.t, gamma);
###          ### Pr(Y_(t-1)=1|V)
###          cur.prev.pyt <- prev.pyt.gvn.v[i];
###          eygvnvr=cur.pyt[1]*(1-cur.prev.pyt)+cur.pyt[2]*cur.prev.pyt; 
###        } 
###        eygvnvr;
###     }     
###   rst  
###   }
###   


##### for AMAR case Pr(v|y_(t-1)=y) for all v and for y=0, and y=1
### based on Pr(V) and Pr(Y_(t-1)=1|V)
AMAR.vgvnprevy <- function(probv, eygvnv) {  
  vrow <- length(probv);
  rst <- matrix(0, vrow, 2);
    eymarginal <- sum(probv*eygvnv); 
    if (eymarginal==1) { rst[,1] <- 0; rst[,2] <- probv;
    } else if (eymarginal==0) { rst[,1] <- probv; rst[,2] <- 0;
    } else { rst[,1] <- probv*(1-eygvnv)/(1-eymarginal);
            rst[,2] <- probv*eygvnv/eymarginal;
    }
  rst    
}

##solve Delta_t for t=0, .., T
##params: alpha, beta, gamma, V and P(V)
get.Deltat <- function(NTime, alpha, beta, gamma, vall, probv, margin=MARGINAL.MODEL) {
    ##all delta t
    rst.delta <- matrix(0,NTime, 2);
    nvall<- nrow(vall);
 
    ##Pr(Y0=1|V)=0 as Y0=1 with prob 1
    prev.pyt.gvn.v <- rep(0, nvall);

    ##Pr(V|Y0=y) for (all v & y=0) as col 1, and (all v & y=1) as col 2
    ##hence the 1st column=probv, and the 2nd column=0 as Y0=1 with prob 0
    v.gvn.prevy0 <- probv;
    v.gvn.prevy1 <- rep(0, nvall);
    
    for (tt in 1:NTime) {
     alpha.t <- alpha[,tt];

     if(tt==1) { 
       rst1 <- nleqslv.fn.t1(tt, alpha.t, beta, gamma, vall, v.gvn.prevy0, v.gvn.prevy1, margin);
       rst.delta[1,1] <- rst1;
     }else{  rst2 <- nleqslv.fn.t2(tt, alpha.t, beta, gamma, vall, v.gvn.prevy0, v.gvn.prevy1, margin);
        rst.delta[tt,] <- rst2;
     }
          
      ### Pr(Yt=1|V) for all v
        prev.pyt.gvn.v <- AMAR.eygvnvonimputate.fast(rst.delta[tt, ], tt, vall, prev.pyt.gvn.v, alpha.t, gamma);
        ### Pr(V|Yt) for all v and Yt in (0,1)
        v.gvn.prevy <- AMAR.vgvnprevy(probv, prev.pyt.gvn.v);
        v.gvn.prevy0 <- v.gvn.prevy[,1];
        v.gvn.prevy1 <- v.gvn.prevy[,2];
        #cat("For T=",tt, ", Delta=", rst.delta[tt,], ".\n");
     } 
    rst.delta
}


# for t=1, delta dim=1 as Y0=0 with prob 1
get.diff.t1 <- function(delta, tt, alpha.t, beta, gamma, vall, v.gvn.prevy0, v.gvn.prevy1, margin) {
    transpose.vall <- t(vall);
    vec.vall <- as.vector(transpose.vall);   
    rst <- .C("XuanAMAR_transitional_getdiff_t1", as.double(delta), as.double(0), as.double(alpha.t), as.double(beta), as.double(gamma), as.integer(tt), as.double(vec.vall), as.double(v.gvn.prevy0), as.double(v.gvn.prevy1), as.integer(ncol(vall)), as.integer(nrow(vall)), as.integer(margin), as.double(0))[[13]];
    stopifnot(!is.nan(rst));
    rst
}


# for t>=2, delta dim=2 for prevy=0, 1
get.diff.t2 <- function(delta, tt, alpha.t, beta, gamma, vall, v.gvn.prevy0, v.gvn.prevy1, margin) {
    rst <- c(0,0);
    transpose.vall <- t(vall);
    vec.vall <- as.vector(transpose.vall);
    tmp <- c(0,0);
    rst <- .C("XuanAMAR_transitional_getdiff_t2", as.double(delta), as.double(alpha.t), as.double(beta), as.double(gamma), as.integer(tt), as.double(vec.vall), as.double(v.gvn.prevy0), as.double(v.gvn.prevy1), as.integer(ncol(vall)), as.integer(nrow(vall)), as.integer(margin), as.double(tmp))[[12]];
     rst
}


##use nonlinear equation solving method
nleqslv.fn.t1 <- function(...) {
    cur.sol <- nleqslv(0, get.diff.t1, jac=NULL,..., method="Newton");
  
    if (cur.sol$termcd > 3) {
        print("Delta t is not solved....");
        stopifnot(1);
    }
    cur.sol$x;
}


nleqslv.fn.t2 <- function(...) {
    delta.start <- c(0,0);
    cur.sol <- nleqslv(delta.start, get.diff.t2, jac=NULL, ..., method="Newton");
 
    if (cur.sol$termcd > 3) {
        print("Delta t is not solved....");
        stopifnot(1);
    }
    cur.sol$x;
}

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

############# Simulate the Data set #############


##simu auxilary covariates. The model is as follows:
## p(V1=v1, V2=v2, ..., Vp=vp) = eta * p(V1=v1)p(V2=v2)...p(Vp=vp)
## eta is the shrinkage parameter
simu.V <- function(sizen, inx.v, prob.v) {

    ##simulation
    smp.v <- rmultinom(sizen, 1, prob.v);
    id.v  <- apply(smp.v, 2, function(x) {which(1==x)});
    rst.v <- inx.v[id.v,];

    rst.v
}

### for AMAR case, simulate complete Y based on the Imputation model
simu.Y <- function(delta, alpha, gamma, simu.v) {

   ntime <- nrow(delta);
   sizen <- nrow(simu.v);
   lengthv <- ncol(simu.v);
   rst <- matrix(NA, sizen, ntime);

   ## Y0=0 with Prob=1
   prev.y <- rep(0,sizen);

   for (l in 1:ntime) {
      cur.py <- expit(delta[l,prev.y+1]+apply(simu.v,1, function(x) { sum(alpha[,l]*x) })+gamma*prev.y);
      cur.y <- sapply(cur.py, function(x) {rbinom(1,1,x)});
      rst[,l] <- cur.y;
      prev.y <- cur.y;
    }
   rst
}


### MDM model
simu.mdm <- function(phi, cur.v, cur.y, t) {
    rst <- sum(phi*c(1,cur.v, cur.y[t-1], cur.y[t]));
    rst
 } 



### for AMAR case, simulate the R based on MDMs
simu.R <- function(yall, vall, phi, f.mdm=simu.mdm) {
    nsize <- nrow(yall);
    ntime <- ncol(yall);

    rst.r <- matrix(NA, nsize, ntime);
    for (i in 1:nsize) {
      cur.v <- vall[i,];
      cur.y <- yall[i,];

      #time 1 y always observed
      rst.r[i,1] <- 0;

      for (t in 2:ntime) {
         logit.pr <- f.mdm(phi, cur.v, cur.y, t);
         pr <- expit(logit.pr);
         cur.r <- rbinom(1,1,pr);

         ### monotone missing, once missing all following missing
         if (1==cur.r) {
           rst.r[i,t:ntime] <-1;
           break;
         }else{
           rst.r[i,t] <- 0
         }
      }   
    }
   rst.r 
}


### for AMAR case, basing on R vector setting corresponding Y as missing
set.missing.y <- function(y.full, r.ind) {

   nsub <- nrow(y.full);
   y.obs <- y.full;
   for(i in 1:nsub) {
       y.obs[i, which(r.ind[i,]==1)] <- NA;  ### r.ind ==1 meaning missing!
   }
   y.obs;
}


##---------------------for covariates --------------------
##counts numbers for each category with covariates data v
##v.inx gives all values for each v
##need a translation from raw V to 111,112,121,122,..
##the final category index corresponds to the id
get.v.inx <- function(vall, v.cats, counts=FALSE) {

    ##remove rows with missing NA
    mis.inx <- which(apply(is.na(vall), 1, sum) > 0);
    if (length(mis.inx) > 0) {
        vall    <- vall[-mis.inx,];
    }

    ##number of columns
    pv <- length(v.cats);

    ##translate v
    trans.v <- get.trans.v(vall, v.cats);

    ##v categories
    ncats   <- mapply(length, v.cats);

    bs <- 1;
    if (pv > 1) {
        for (i in pv:2) {
            bs <- c(prod(ncats[i:pv]), bs);
        }
    }

    ##minus 1 for cols 2 and above
    alt.v <- trans.v;
    if (pv > 1) {
        alt.v[, 1:(pv-1)] <- alt.v[, 1:(pv-1)] - 1;
    }
    rst <- apply(alt.v, 1, function(x){ sum(x*bs) });

    if (counts) {
        rst.counts <- NULL;
        for (i in 1:prod(ncats)) {
            rst.counts <- c(rst.counts, length(which(i == rst)));
        }
        rst <- rst.counts;
    }

    rst
}



### This function is used to calculate the number of rows in matrix vr.all that are equal to vector vr
matchn <- function(vr, vr.all)
{ 
   ttl <- 0;  
   samplesize <- nrow(vr.all);
   for (i in 1:samplesize)   { ttl <- ttl+isTRUE(all.equal(vr,vr.all[i,], tolerance=0))*1;}
   ttl
}

############ End of Simulation of Data set #########################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

#### For Uncongenial Analysis


### Find the index value (which row) of a realization of vr in the matrix of vr.all, which is all category combination of vr
findindex <- function(vr, vr.all){
    ttlrow <- length(vr.all)/length(vr);
    indexvr <- NULL;
    lengthvr <- length(vr);
    for (i in 1:ttlrow) {
       tmp.vr <- vr.all[((i-1)*lengthvr+1):(i*lengthvr)];
       if (isTRUE(all.equal(c(vr), c(tmp.vr), tolerance=0)))  indexvr=i;
     }
    indexvr
}


#### based on the estimated imputation model, imputate the missing ys
SIMU.Y.AMAR.xuan <- function(delta, alpha, gamma, simu.v, simu.r, simu.y) {
    sizen <- nrow(simu.v);

    indexR1 <- which(simu.r==1);
    indexR2 <- which(simu.r==2);

    cur.py <- expit(delta[3]+simu.v[indexR2,]%*%alpha+gamma*simu.y[indexR2,2]);
    #wrong size cur.y <- rbinom(1,1,cur.py);
    cur.y <- sapply(cur.py,function(x) rbinom(1,1,x)); #update 3/19
    simu.y[indexR2, 3] <- cur.y;

    cur.py <- expit(delta[2]+simu.v[indexR1,]%*%alpha+gamma*simu.y[indexR1,1]);
    #cur.y <- rbinom(1,1,cur.py);
    cur.y <- sapply(cur.py,function(x) rbinom(1,1,x)); #update 3/19
    simu.y[indexR1, 2] <- cur.y;
    cur.py <- expit(delta[3]+simu.v[indexR1,]%*%alpha+gamma*simu.y[indexR1,2]);
    cur.y <- sapply(cur.py, function(x) rbinom(1,1,x));
    simu.y[indexR1, 3] <- cur.y;

    simu.y
  }


###########################################################################
###########################################################################
###########################################################################

### For Congenial Analysis AMAR -- HMC

### AMAR case for congenial analysis, each time has two deltas
### imputation model for observed outcome, only return probability Pr(Y=1|all) from imputation model
### works on single observation
###obs.imputatemdl <- function(t,v,r,y, delta, alpha, gamma) {
###  rst <- NULL;
###  
###  if (t==1) {
###    rst= expit(delta[1,1]+(v)%*%alpha[,1]);
###  }else if (t==2) {
###    if (r==1)  { rst=0; }
###    else { rst <- expit(delta[2,y[1]+1]+(v)%*%alpha[,2]+gamma*y[1]); }
###  } else {
###    if (r<=2) { rst=0; }
###    else { rst <- expit(delta[3,y[2]+1]+(v)%*%alpha[,3]+gamma*y[2]); }
###  }
###  rst
###}
###

###
###### calculate for the observed outcome, the logit(pr(Y=1)) part:= eta
###### if outcome is missing for this R, let eta=0
###old.obs.etainimputatemdl <- function(v,r,y,delta,alpha1,gamma,alpha3) {
###  rst <- 0;
###
###  n.obs <- nrow(v);
###
###  eta.matrix <- matrix(0,n.obs,ncol(y)+1);
###  
### for(i in 1:n.obs)
### {
###   if (r[i]==1)
###     { eta.matrix[i,1] <- delta[1]+(v[i,])%*%alpha1[,1];
###       eta.matrix[i,4] <- y[i,1]*eta.matrix[i,1]-log(1+exp(eta.matrix[i,1])); }
###   else if (r[i]==2)
###     { eta.matrix[i,1] <- delta[1]+(v[i,])%*%alpha1[,1]+alpha3[1];
###       eta.matrix[i,2] <- delta[2]+(v[i,])%*%alpha1[,2]+gamma[1]*y[i,1];
###       eta.matrix[i,4] <- y[i,1]*eta.matrix[i,1]-log(1+exp(eta.matrix[i,1]))+y[i,2]*eta.matrix[i,2]-log(1+exp(eta.matrix[i,2]));}
###   else
###     { eta.matrix[i,1] <- delta[1]+(v[i,])%*%alpha1[,1]+alpha3[2];
###       eta.matrix[i,2] <- delta[2]+(v[i,])%*%alpha1[,2]+gamma[1]*y[i,1];
###       eta.matrix[i,3] <- delta[3]+(v[i,])%*%alpha1[,3]+gamma[2]*y[i,2];
###       eta.matrix[i,4] <- y[i,1]*eta.matrix[i,1]-log(1+exp(eta.matrix[i,1]))+y[i,2]*eta.matrix[i,2]-log(1+exp(eta.matrix[i,2]))+y[i,3]*eta.matrix[i,3]-log(1+exp(eta.matrix[i,3]));}
### }
###  
###  
###  eta.matrix;
###}
###


### AMAR case
### the similar function with obs.imputatemdl, only that this one works on all observations at once
### it gives Pr(Y_t=1|Y_(t-1),X,V,R) for each observation, and for all t
### returning matrix with dim=(n.obs, n.T)
new.all.obs.imputatemdl <- function(v, r, y, delta, alpha, gamma) {
   n.obs <- nrow(y);
   n.T <- ncol(y);
   rst <- matrix(0, n.obs, n.T);
   
   rst[,1] <- expit(delta[1,1]+v%*%alpha[,1]);
   indexr2 <- which(r>=2);
   rst[indexr2,2] <- expit(delta[2,y[indexr2,1]+1]+v[indexr2,]%*%alpha[,2]+gamma*y[indexr2,1]);
   indexr3 <- which(r>=3);
   rst[indexr3,3] <- expit(delta[3,y[indexr3,2]+1]+v[indexr3,]%*%alpha[,3]+gamma*y[indexr3,2]);
   rst
}
 
### AMAR case
### gives for each subject, the contribution in the observed likelihood as in 1207 P4 above
### give log( sum_obs(y) p^y*(1-p)^(1-y)) for each subject
new.obs.likelihood <- function(v,r,y,delta,alpha,gamma) {
  rst <- 0;
  n.obs <- nrow(v);

  ### get the pr(Y_t=1|Y_(t-1),V,R), bounded in [0,1] with end pointes included (n*3) dimension
  prob.impt <- new.all.obs.imputatemdl(v,r,y,delta, alpha, gamma);  
  ind.sum.loglike <- NULL;

  ### bound loglikelihood from -Inf
  bound.loglike <- matrix(0, n.obs, ncol(y));

  for(j in 1:ncol(y)) {
  
  indx1 <- intersect(which(r>=j), which(prob.impt[,j]+y[,j]==1));
  bound.loglike[indx1,j] <- -100;

  indxmis <- which(r<j);
  indx2 <- intersect(which(r>=j),which(prob.impt[,j]==y[,j]));
  indx2 <- union(indx2, indxmis);
  bound.loglike[indx2,j] <- 0;

  union1 <- union(indx1, indx2);
  indxnormal <- setdiff(1:n.obs, union1);
  bound.loglike[indxnormal,j] <- y[indxnormal,j]*logit(prob.impt[indxnormal,j])+log(1-prob.impt[indxnormal,j]);
  }
  
  ind.sum.loglike <- sapply(1:n.obs, function(x) sum(bound.loglike[x,1:r[x]]));
  ind.sum.loglike;
}





### This is the correct U function for gradient calculation with grad()
new.U1 <- function(u.q, u.v, u.r, u.y, u.sigma) {
  
  # parameters in transition model
  u.Beta <- u.q[1:3];
  u.Alpha <- u.q[4:11];
  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Gamma <- u.q[12];
 
  logPrior <- sum(-u.q^2/(2*u.sigma^2));

  u.Delta <- get.Deltat(ncol(u.y), u.fullAlpha, u.Beta, u.Gamma, V.INX.ALL, est.theta, SIMU.MARGIN)
    
  Eta.matrix <- new.obs.likelihood(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
  logLike = sum(Eta.matrix);
  
  U= -(logLike + logPrior);
  cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
  U
}


#####################################################################################################################
#####################################################################################################################

### AMAR case for slice sampling of Eta
### copy from Daniels & Chenguang


##vall: the observed/simulated V data
##niter: number of iterations
##log.fx: loglikelihood function
##inits: initial values
##v.inx:categories for all V auxillary covariates
slice.v.par <- function(vall, niter, log.fx, v.inx,
                        prior.neta=NULL, inits=NULL) {

    bound.eta <- 10;

    ##v categories
    ncats  <- mapply(length, v.inx);
    ntheta <- prod(ncats);
    npv    <- length(v.inx);

    ##init values
    if (is.null(inits)) {
        for (i in 1:npv) {
            ##pi_is init at 1/L_i
            inits <- c(inits, rep(1/ncats[i], ncats[i]));
        }
        ##eta init at 1/N
        inits <- c(inits, 1/nrow(vall));
    }

    ##get count data
    vcount <- get.v.inx(vall, v.inx, counts=TRUE);

    ##initialize results
    rst.theta <- array(NA, dim=c(ntheta, niter));
    rst       <- array(NA, dim=c(niter+1, sum(ncats)+1));
    rst[1, ]  <- inits;
    for (s in 2:(niter+1)) {
        ##browser();
        rst[s,] <- rst[s-1,];
        ##gibbs theta L
        pv.eta      <- v.par.vec2lst(rst[s,], ncats);
        lst.prob.v  <- pv.eta[[1]];
        eta         <- pv.eta[[2]];
        prob.v      <- get.prob.v(lst.prob.v, eta);
        cur.dir.par <- vcount + prob.v;
        next.theta  <- rdirichlet(1, cur.dir.par);
        rst.theta[, s-1] <- next.theta;

        ##Gibbs sampling for pi's
        cur.start <- 1;
        for (i in 1:npv) {
            cur.end     <- (cur.start-1) + ncats[i];
            next.start  <- cur.end + 1;

            ##random choose one pi to be 1-sum of the rest pis
            cur.fix.inx <- sample(cur.start:cur.end, 1);
            for (j in cur.start:cur.end) {
                ##fixed pi will not be sampled
                if (cur.fix.inx == j) next;

                ##find bound
                up.bound <- rst[s,j] + rst[s,cur.fix.inx];

                ##slice sampling
                slc.pi <- slice.onestep(next.theta, rst[s,], j, cur.fix.inx, up.bound,
                                        log.fx, ncats, prior.neta);

                ##update
                rst[s, j] <- slc.pi;
                rst[s, cur.fix.inx] <- up.bound - slc.pi;
            }

            cur.start <- next.start;
        }

        ##Gibbs/Slice sampling for eta
        slc.eta <- slice.onestep(next.theta, rst[s,], ncol(rst), NULL,
                                 bound.eta, log.fx, ncats, prior.neta);
        rst[s, ncol(rst)] <- slc.eta;

    }

    list(post.pieta=rst, post.prob.v=rst.theta);
}


##one iteration of slice sampling
##parinx is the current parameter for slice sampling
##fix.inx is the one represents 1-the rest
##if fix.inx is NULL, that's for ETA
##nvcats is number of categories for all auxillary covariates
slice.onestep <- function(ndata, curpar, parinx, fix.inx, upbound, log.fx, ..., max.iter=80000) {
    ##browser();
    bound <- c(0, upbound);
    x0    <- curpar[parinx];

    ##get log f(x0)
    log.fx0 <- log.fx(ndata, curpar, ...);

    ##get auxiliary y~unif(0,fx0), that is
    log.y   <- log.fx0 - rexp(1);

    ##sampling from slice
    flag   <- FALSE;
    n.iter <- 0;
    while (!flag & (n.iter<max.iter)) {
        x1             <- runif(1, bound[1], bound[2]);
        curpar[parinx] <- x1;

        if (!is.null(fix.inx)) {
            curpar[fix.inx]<- upbound - x1;
        }
        log.fx1 <- log.fx(ndata, curpar, ...);

        if (log.y < log.fx1) {
            flag <- TRUE;
        } else if (x1 < x0) {
            bound[1] <- x1;
        } else if (x1 > x0) {
            bound[2] <- x1;
        }

        n.iter <- n.iter + 1;
    }

    if (n.iter >= max.iter) {
        print("slice failed");
    }
    x1;
}



##multinomial log likelihood for auxillary covariates
##last curpar is eta
##curpar is a vector, need to make it a list first and then
##call get.prob.v
v.log.fx <- function(vcount, curpar, nvcats, neta=NULL) {
    ##convert to lst of prob v
    pv.eta     <- v.par.vec2lst(curpar, nvcats);
    lst.prob.v <- pv.eta[[1]];
    eta        <- pv.eta[[2]];
    prob.v     <- get.prob.v(lst.prob.v, eta);

    logD   <- lgamma(sum(prob.v)) - sum(lgamma(prob.v));
    logV   <- as.vector(log(vcount));

    if (!is.null(neta)) {
        ##uniform shrinkage prior
        logeta <- -2 * log(neta * eta + 1);
    } else {
        logeta <- 0;
    }

    ##bound away from -Inf
    logV[which(is.infinite(logV))] <- -1000;
    rst  <- logD + sum((prob.v-1)*logV) + logeta;
      
    stopifnot(!is.infinite(rst))
    rst
}



##convert vector to list for parameters for
##auillary covariates V
v.par.vec2lst <- function(vpar, nvcats, addeta=TRUE) {
    lst.prob.v <- list();
    tmp.start  <- 1;
    for (i in 1:length(nvcats)) {
        cur.n <- nvcats[i];
        if (cur.n > 0) {
            lst.prob.v[[i]] <- vpar[tmp.start:sum(nvcats[1:i])];
        } else {
            lst.prob.v[[i]] <- NULL;
        }
        tmp.start <- 1 + sum(nvcats[1:i]);
    }

    if (addeta) {
        ##rst: prob.v and eta
        rst <- list(lst.prob.v, vpar[length(vpar)]);
    } else {
        rst <- lst.prob.v;
    }

    rst
}

## change from list to vector
## this may be done by R function unlist
v.par.lst2vec <- function(lst.prob.v, eta=NULL) {
    rst <- NULL;
    for (i in 1:length(lst.prob.v)) {
        rst <- c(rst, lst.prob.v[[i]]);
    }

    rst <- c(rst, eta);
    rst
}



#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
####     
####     ### PROBLEMATIC !!!! updated in 06/07/2014
####     #### Congenial Analysis for Posterior sampling of Eta
####     ### HMC for eta only
####     U.eta <- function(u.eta, u.ksi, u.theta, nsize) {
####     
####       nrowv <- length(u.theta);
####       logPrior <- -2*log(u.eta*nsize+1);
####       
####       logtheta <- log(u.theta);
####        logtheta[which(is.infinite(logtheta))] <- -1000;
####       p1 <- lgamma(1/u.eta);
####       p2 <- sum(sapply(1:nrowv, function(x) -lgamma(u.ksi[x]/u.eta)));
####     
####       indx.thetagt0 <- which(u.theta>0);
####       p3 <- sum(sapply(1:nrowv, function(x) logtheta[x]*(u.ksi[x]/u.eta-1)));
####       logLike <- p1+p2+p3;
####       
####       E= -(logLike+logPrior);
####       #cat("eta=", u.eta, "logPrior=", logPrior, "logLike p1=", p1, "p2=", p2, "p3=", p3, "loglike=", logLike, "hence U=", E, ".\n");
####       E
####     }
####     
####     ### PROBLEMATIC !!!! updated in 06/07/2014
####     grad.U.eta <- function(g.eta, g.ksi, g.theta, nsize) {
####     
####       nrowv <- length(g.theta);
####       dLogPrior <- -2*nsize/(g.eta*nsize+1);
####     
####       logtheta <- log(g.theta);
####          logtheta[which(is.infinite(logtheta))] <- -1000;
####       dLogLike <- -digamma(1/g.eta)/(g.eta^2) + sum(sapply(1:nrowv, function(x) -g.ksi[x]*logtheta[x]/(g.eta^2)+g.ksi[x]/(g.eta^2)*digamma(g.ksi[x]/g.eta)));
####       g <- -(dLogPrior+dLogLike);
####       g
####     }  
####     
####     ### PROBLEMATIC !!!! updated in 06/07/2014
####     hmc.eta.standard <- function(ksi, eta, n.iter, epsilon, L, diag.sigma.p, emp.v.n, lowbound.eta, upbound.eta, true.theta) {
####     
####       post.samp <- NULL;
####       post.samp.theta <- matrix(0, n.iter+1, length(emp.v.n));
####     
####       current.eta <- eta;
####       current.theta <- rdirichlet(1, emp.v.n+ksi/current.eta);
####       
####       nsize <- sum(emp.v.n);
####       n.accept <- 0;
####     
####       if(current.eta<0.01) { epsilon <- (10)^(-5); }
####        
####       for (i in 1:n.iter) {
####         if(i==1) {
####         current.U <- U.eta(current.eta, ksi, current.theta, nsize);
####         ### both grad.U.eta and grad(U.eta) give the same gradient
####         current.g <- grad.U.eta(current.eta, ksi, current.theta, nsize);
####         #current.g <- grad(U.eta, x=current.eta, u.ksi=ksi, u.theta=current.theta, nsize=nsize);
####         post.samp[1] <- current.eta;
####         post.samp.theta[1,] <- current.theta;
####         }
####         
####         samp <- eta.getSamp(current.eta, current.U, current.g, epsilon, L, diag.sigma.p, ksi, current.theta, nsize, upbound.eta, lowbound.eta, emp.v.n);
####     
####         current.eta <- samp$q;
####         current.update <- samp$Update;
####         
####         current.theta <- rdirichlet(1, emp.v.n+ksi/current.eta);
####         current.U <- U.eta(current.eta, ksi, current.theta, nsize);
####         current.g <- grad.U.eta(current.eta, ksi, current.theta, nsize);
####         #current.g <- grad(U.eta, x=current.eta, u.ksi=ksi, u.theta=current.theta, nsize=nsize);
####         
####         post.samp[i+1] <- current.eta;
####         post.samp.theta[i+1,] <- current.theta;
####     
####         n.accept <- n.accept+1*current.update;
####         acceptrate <- n.accept/i;
####     
####         avg.pchi <- probv.pchi(apply(post.samp.theta[1:(i+1),],2,mean), true.theta);
####         crnt.pchi <- probv.pchi(current.theta, true.theta);
####         
####         if(current.update==0) { cat("finish iteration", i, ", eta==",  current.eta, ", pchi===(", avg.pchi, ",", crnt.pchi, ").  Reject with acceptrate=", acceptrate, "\n");
####         }else {
####           cat("finish iteration", i, ",  eta==", current.eta, ", pchi===(", avg.pchi, ",", crnt.pchi, "). Accept with acceptrate=", acceptrate, "\n");
####         }
####      }
####     
####     list(samp=post.samp, samp.theta=post.samp.theta, acceptrate=acceptrate)
####     }
####     
####     ### PROBLEMATIC !!!! updated in 06/07/2014
####     new.hmc.eta.standard <- function(ksi, eta, n.iter, epsilon, L, diag.sigma.p, emp.v.n, lowbound.eta, upbound.eta, true.theta) {
####     
####       post.samp <- NULL;
####       post.samp.theta <- matrix(0, n.iter+1, length(emp.v.n));
####     
####       current.eta <- eta;
####       current.theta <- rdirichlet(1, emp.v.n+ksi/current.eta);
####       
####       nsize <- sum(emp.v.n);
####       n.accept <- 0;
####       for (i in 1:n.iter) {
####         if(i==1) {
####           logtheta <- log(current.theta);
####           logtheta[which(is.infinite(logtheta))] <- -1000;
####         cat("length theta==0)===", length(which(current.theta==0)), " , sum==", sum(logtheta*emp.v.n), "\n");
####         current.U <- U.eta(current.eta, ksi, current.theta, nsize)+sum(logtheta*emp.v.n);
####         ### both grad.U.eta and grad(U.eta) give the same gradient
####           current.g = grad.U.eta(current.eta, ksi, current.theta, nsize);
####         #current.g <- grad(U.eta, x=current.eta, u.ksi=ksi, u.theta=current.theta, nsize=nsize);
####         post.samp[1] <- current.eta;
####         post.samp.theta[1,] <- current.theta;
####         cat("for iter 1, c.U==", current.U, "\n");
####         }
####         
####         samp <- new.eta.getSamp(current.eta, current.U, current.g, epsilon, L, diag.sigma.p, ksi, current.theta, nsize, upbound.eta, lowbound.eta, emp.v.n);
####     
####         current.eta <- samp$q;
####         current.update <- samp$Update;
####         current.theta <- samp$theta;
####         current.U <- samp$U;
####         current.g <- samp$g;
####     
####         post.samp[i+1] <- current.eta;
####         post.samp.theta[i+1,] <- current.theta;
####     
####         n.accept <- n.accept+1*current.update;
####         acceptrate <- n.accept/i;
####     
####         if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
####             {
####               cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
####               list(samp=post.samp, acceptrate=acceptrate)
####               break
####             }
####         avg.pchi <- probv.pchi(apply(post.samp.theta[1:(i+1),],2,mean), true.theta);
####         crnt.pchi <- probv.pchi(current.theta, true.theta);
####         
####         if(current.update==0) { cat("finish iteration", i, ", eta== ", current.eta, ", pchi===(", avg.pchi, ",", crnt.pchi, ").  Reject with acceptrate=", acceptrate, "\n");
####         }else {
####           cat("finish iteration", i, ",  eta==", current.eta, ",  pchi===(", avg.pchi, ",", crnt.pchi, "). Accept with acceptrate=", acceptrate, "\n");
####         }
####      }
####     
####     list(eta=post.samp, theta=post.samp.theta, acceptrate=acceptrate)
####     }  
####     
####     ### PROBLEMATIC !!!! updated in 06/07/2014
####     new.eta.getSamp <- function(q, U, g, et, L, diag.sigma.p, ksi, theta, nsize, upbound.eta, lowbound.eta, emp.v.n){
####       current.q <- q;
####       current.g <- g;
####       current.U <- U;
####       p <- rnorm(1,0,diag.sigma.p^2);
####       current.p <- p;
####     
####       current.K <- current.p^2/(2*diag.sigma.p^2);
####       current.H <- current.U+current.K;
####     
####       et <- runif(1, et*0.8, et*1.2);
####       if(current.q < 0.01) { et <- (10)^(-5); }
####     
####       for(leap in 1:L) {
####          p = p-et*g/2;
####          q = q+et*p/(diag.sigma.p^2);
####          
####          while(q<lowbound.eta||q>upbound.eta)
####          { if(q<lowbound.eta){
####            q=2*lowbound.eta-q;
####            p=-p;
####            }
####            if(q>upbound.eta) {
####            q=2*(upbound.eta)-q;
####            p=-p;
####            }
####          }
####          #g <- grad(U.eta, x=q, u.ksi=ksi, u.theta=theta, nsize=nsize);
####          g = grad.U.eta(q, ksi, theta, nsize);
####          p = p-et*g/2;
####        }
####     
####       proposed.theta <- rdirichlet(1, emp.v.n+ksi/q);
####         logtheta <- log(proposed.theta);
####         logtheta[which(is.infinite(logtheta))] <- -1000;
####       proposed.U = U.eta(q, ksi, proposed.theta, nsize)+sum(logtheta*emp.v.n);
####       proposed.K = p^2/(2*diag.sigma.p^2);
####       proposed.H = proposed.U+proposed.K;
####     
####       acceptProb=min(1, exp(current.H-proposed.H));
####      
####       #cat("\n q=", q, ", sqr(diff.q)=", (current.q-q)^2, ", c.U-p.U=", current.U-proposed.U, ", c.K-p.K=", current.K-proposed.K, ", accept prob=", acceptProb, "\n");
####       
####       ### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
####       if(runif(1,0,1) < acceptProb){
####         #proposed.g <-  grad(U.eta, x=q, u.ksi=ksi, u.theta=proposed.theta, nsize=nsize);
####         proposed.g = grad.U.eta(q, ksi, proposed.theta, nsize);
####         list(q = q, Update=1, theta=proposed.theta, U=proposed.U, g=proposed.g);  # accept
####        }else{
####         update.theta <- rdirichlet(1, emp.v.n+ksi/current.q);
####          logtheta <- log(update.theta);
####         logtheta[which(is.infinite(logtheta))] <- -1000;
####         current.U <- U.eta(current.q, ksi, update.theta, nsize)+sum(logtheta*emp.v.n);
####         #current.g <-  grad(U.eta, x=current.q, u.ksi=ksi, u.theta=update.theta, nsize=nsize);
####         current.g = grad.U.eta(current.q, ksi, update.theta, nsize);
####         list(q = current.q, Update=0, theta=update.theta, U=current.U, g=current.g);  # reject
####        }
####     
####     }
####     
####     eta.getSamp <- function(q, U, g, et, L, diag.sigma.p, ksi, theta, nsize, upbound.eta, lowbound.eta, emp.v.n){
####       current.q <- q;
####       current.g <- g;
####       current.U <- U;
####       p <- rnorm(1,0,diag.sigma.p^2);
####       current.p <- p;
####     
####       current.K <- current.p^2/(2*diag.sigma.p^2);
####       current.H <- current.U+current.K;
####     
####       et <- runif(1, et*0.8, et*1.2);
####        if(current.q < 0.1) { et <- (10)^(-5); }
####     
####      
####       for(leap in 1:L) {
####          p = p-et*g/2;
####          q = q+et*p/diag.sigma.p^2;
####     
####          while(q<lowbound.eta||q>upbound.eta)
####          { if(q<lowbound.eta){
####             q=2*lowbound.eta-q;
####             p=-p;
####            }
####            if(q>upbound.eta) {
####             q=2*(upbound.eta)-q;
####             p=-p;
####            }
####          }  
####          g = grad.U.eta(q, ksi, theta, nsize);
####          #g <- grad(U.eta, x=q, u.ksi=ksi, u.theta=theta, nsize=nsize);
####          p = p-et*g/2;
####          #cat("After p= ", p, ", diff.q= ", q-current.q, ", g= ", g, "\n");   
####        }
####       proposed.U = U.eta(q, ksi, theta, nsize);
####       proposed.K = p^2/(2*diag.sigma.p^2);
####       proposed.H = proposed.U+proposed.K;
####     
####       acceptProb=min(1, exp(current.H-proposed.H));
####      
####       cat("\n sqr(diff.q)=", (current.q-q)^2, ", c.U-p.U=", current.U-proposed.U, ", c.K-p.K=", current.K-proposed.K, ", accept prob=", acceptProb, "\n");
####     
####       ### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
####       if(runif(1,0,1) < acceptProb){
####         list(q = q, Update=1); # accept
####        }else{
####         list(q = current.q, Update=0); # reject
####        }
####     
####     }
####     
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

### For Random Walk Metropolis Hasting for Congenial Analysis in AMAR, xuan

rwmh.abg.AMAR.xuan <- function(log.fx, init=NULL, iteration=2500, ys, vs, rs, prob.v.dist, VALL, SIMU.MARGIN, best.scale, prior.std=10, filename, tune=FALSE){
  
  npars <- length(init);
  cur.par <- init;
  nTime <- ncol(ys);
  
  nxt.delta <- NA;
  nxt.ll <- NULL; ### temporarily keep the likelihood

   
  plotindx <- c(100, 2500, 5000, 7500, 10000);
  true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
  
  rstall.mcmc <- array(NA, dim=c(iteration, length(init)));

  for (i in 1: iteration) {
     for (j in 1:npars) {
       # based on current q, get Delta, Alpha, Beta, Gamma
       alpha <- cbind(cur.par[4:11], cur.par[4:11], cur.par[4:11]);
       beta <- cur.par[1:3];
       gamma <- cur.par[12];

       if(i==1&j==1)
       { nxt.delta = get.Deltat(nTime, alpha, beta, gamma, VALL, prob.v.dist, SIMU.MARGIN); }

       rb0 <- rwmh.onestep.AMAR(log.fx, cur.par, j, best.scale[j], cur.pos=nxt.ll, ys, vs, rs, prior.std, prob.v.dist, cur.delta=nxt.delta, VALL, SIMU.MARGIN);
       nxt.par <- rb0$par;
       nxt.ll <- rb0$ll;
       cur.par[j] <- nxt.par;
       rstall.mcmc[i,j] <- nxt.par;
       nxt.delta <- rb0$delta;
     }
   ###  write.table(rstall.mcmc[1:i,], file=paste(filename,"_RWMH.txt", sep=""), append=FALSE);
   ###

   if (length(which(plotindx==i))>0)
    { jpeg(paste(filename,"_p1Niter_",i,".jpg", sep=""));
      par(mfrow=c(2,2));
      
      for(j in 1:3) {
        plot(1:i, rstall.mcmc[1:i,j], main=paste("BETA",j), type="l");
        abline(true.q[j], 0, col="blue");
      }
      dev.off();

      jpeg(paste(filename,"_p2Niter_",i,".jpg", sep=""));
      par(mfrow=c(3,3));
      for(j in 4:12) {
        if(j < 12) {
         plot(1:i, rstall.mcmc[1:i,j], main=paste("ALPHA",j-3), type="l");
         abline(true.q[j], 0, col="blue");
       } else {
         plot(1:i, rstall.mcmc[1:i,12], main="GAMMA", type="l");
         abline(true.q[12], 0, col="blue");
       }  
      }
      dev.off();
    }
   ###  
    if(i >=2) {
    nupdate <- npars-sum(sapply(1:npars, function(x) I(rstall.mcmc[i-1,x]==rstall.mcmc[i,x])));
    cat("\n Finish iteration", i, "update ", nupdate, " var with abg=", rstall.mcmc[i,], "\n");
    }
    else{
    cat("\n Finish iteration", i, "with abg=", rstall.mcmc[i,], "\n"); }
  }
  #thin
  #rstall.mcmc <- rstall.mcmc[seq(thin, by=thin, length.out=iteration),];
  ##burnin
  #rstall.mcmc <- rstall.mcmc[-(1:nburnin),];

  list(mcmc=rstall.mcmc, scale=best.scale)
}


### random walk metropolis hasting
## scale is the size of random walk

rwmh.onestep.AMAR <- function(log.fx, cur.par, parinx, scale, cur.pos=NULL, ys, vs, rs, prior.std, prob.v.dist, cur.delta, VALL, SIMU.MARGIN) {

   # if current log(prior*obs) is null, calculate it
   if (is.null(cur.pos)) {
      cur.pos <- log.fx(cur.par, ys, vs, rs, cur.delta, prior.std, prob.v.dist, VALL, SIMU.MARGIN);
   }

     rst.b0 <- NULL;
     current.b0 <- cur.par[parinx];
     next.b0 <- current.b0+scale*rnorm(1);
   
     next.par <- cur.par;
     next.par[parinx] <- next.b0;
   
     alpha <- cbind(next.par[4:11],next.par[4:11],next.par[4:11]);
     beta <- next.par[1:3];
     gamma <- next.par[12];

     next.delta <- get.Deltat(ncol(ys), alpha, beta, gamma, VALL, prob.v.dist, SIMU.MARGIN);
     next.pos <- log.fx(next.par, ys, vs, rs, next.delta, prior.std, prob.v.dist, VALL, SIMU.MARGIN);

     prob <- min(1,exp(next.pos-cur.pos));

    if(is.infinite(next.pos)) { print("Inf of posterior log likelihood ..."); }
  
  if(runif(1)<prob) { list(par=next.b0, rate=1, scale=scale, ll=next.pos, delta=next.delta);}
  else { list(par=current.b0, rate=0, scale=scale, ll=cur.pos, delta=cur.delta); }
 }


### get posterior of alpha, beta, gamma
get.log.post.abg <- function(vecpar, ys, vs, rs, delta=NULL, prior.std, prob.v.dist, VALL, SIMU.MARGIN) {
    
  if(is.null(delta))
    { alpha <- cbind(vecpar[4:11], vecpar[4:11], vecpar[4:11]);
      delta=get.Deltat(ncol(ys), alpha, vecpar[1:3], vecpar[12], VALL, prob.v.dist, SIMU.MARGIN); }
  prior.abg <- sum(dnorm(c(vecpar), 0, prior.std, log=TRUE));
  logl <- abg.obs.loglike(vecpar, ys, vs, rs, delta);
  #cat("prior.abg=", prior.abg, ", logl=", logl, "\n");
  rst <- prior.abg+logl;
  rst
}



abg.obs.loglike <- function(q, ys, vs, rs, u.Delta) {

  n.param = length(q);
  n.obs = nrow(vs);

  u.Beta <- q[1:3];
  u.Alpha <- q[4:11];
  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Gamma <- q[12];
  
  Eta.matrix <- obs.likelihood(vs, rs, ys, u.Delta, u.fullAlpha, u.Gamma);
  logLike = sum(Eta.matrix[1:n.obs]);
  
  logLike
}


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

### obtain the theta, and mse.theta basing on Daniels & Chenguang's original method
### assume independence among covariates

prv.dan <- function(niter, nburn, v, v.log.fx, V.CATS, ttln, true.theta)
{

  post.v.par <- slice.v.par(v, niter, v.log.fx, V.CATS, ttln);
  prob.v.dist <- post.v.par$post.prob.v;

  mean.theta.dan <- apply(prob.v.dist[,(nburn+1):niter], 1, mean);
  mse.dan <- t(true.theta-mean.theta.dan)%*%(true.theta-mean.theta.dan);
  pchi.probv.rst <- probv.pchi(mean.theta.dan, true.theta);

    eta.inx <- ncol(post.v.par$post.pieta);

  list(all.theta=prob.v.dist, all.eta=post.v.par$post.pieta[,eta.inx], mean.theta=mean.theta.dan, mse.theta=mse.dan, pchi.theta=pchi.probv.rst);
}

### DON'T use!!!! Update 04/01-- this method is problematic. DONT use!!!!
### similary approach is prv.slice --- difference is the sampling method for eta (1) hmc or (2) slice sampling
### obtain the theta, and mse.theta basing on Loglinear + HMC for theta
#prv.hmc <- function(niter, nburn, eptime, leapFrog, diag.sigma.p, Eta.lowbound, Eta.upbound, eta.init, V.INX.ALL, emp.v.n, true.theta) {

#  v.loglinear.dat <- as.data.frame(cbind(V.INX.ALL, emp.v.n));
 # colnames(v.loglinear.dat) <- c("X1","X2","X3","X4","X5","X6","X7","X8","Freq");
#
#  ### To make the sequence of fitted.values of this log linear model is the same with VR.INX.ALL
#  ### need to fit the model with factor sequence as: R, X8, X7, X6, ..., X1
###  attach(v.loglinear.dat);
#  v.loglinear <- glm(Freq ~ .+.*., family=poisson, data=v.loglinear.dat);
#  #v.loglinear <- glm(Freq ~ .+.*.+.*.*., family=poisson, data=v.loglinear.dat);
#  ### the fitted values gives the estimated frequence for each category, hence for the probability of each category,
#  ### we need to divide by the total number of frequence
#  ksi <- (v.loglinear$fitted.values)/sum(emp.v.n);
#
#  cat("finish estimation of ksi \n");
#  ### samp gives the sampling obtained by using HMC
#  samp <- hmc.eta.standard(ksi, eta.init, niter, eptime, leapFrog, diag.sigma.p, emp.v.n, Eta.lowbound, Eta.upbound);
#  
#  samp.theta <- samp$samp.theta;
#  samp.eta <- samp$samp;
#  mean.theta.hmc <- sapply(1:length(ksi), function(x) mean(samp.theta[nburn:niter,x]));
#  
#  ### mean square error of hmc method of theta
#  mse.hmc <- t(true.theta-mean.theta.hmc)%*%(true.theta-mean.theta.hmc);
#  pchi.probv.rst <- probv.pchi(mean.theta.hmc, true.theta);
#
#  list(all.theta=samp.theta, mean.theta=mean.theta.hmc, mse.theta=mse.hmc, pchi.theta=pchi.probv.rst, all.eta=samp, rate=samp$acceptrate);
#}
#


### update on 2015/02/12 -- change the colnames definition with loop, so fit for all dimension of V
### obtain the theta, and mse.theta basing on Loglinear as ksi from modeling + Slice Sampling for theta
prv.slice <- function(V.INX.ALL, emp.v.n, niter, nburn, eta.init, eta.lowbound, eta.upbound, true.theta) { 
  v.loglinear.dat <- as.data.frame(cbind(V.INX.ALL, emp.v.n));
  #colnames(v.loglinear.dat) <- c("X1","X2","X3","X4","X5","X6","X7","X8","Freq");
  colnames(v.loglinear.dat) <- c( paste("V", 1:ncol(V.INX.ALL),sep=""),"Freq");
  
  ### To make the sequence of fitted.values of this log linear model is the same with VR.INX.ALL
  ### need to fit the model with factor sequence as: R, X8, X7, X6, ..., X1
  attach(v.loglinear.dat);
  v.loglinear <- glm(Freq ~ .+.*. , family=poisson, data=v.loglinear.dat);
  ### fitted values gives the estimated frequence for each category, hence for the probability of each category,
  ### we need to divide by the total number of frequence
  ksi <- (v.loglinear$fitted.values)/sum(emp.v.n);

  #cat("finish estimation of ksi \n");

   pv <- slice.v.par.AMAR(ksi, emp.v.n, niter, eta.init, eta.lowbound, eta.upbound, true.theta);
   pv.theta <- pv$post.theta;

   mean.theta.slice <- apply(pv.theta[(nburn+1):niter,],2,mean);
   mse.slice <- t(true.theta-mean.theta.slice)%*%(true.theta-mean.theta.slice);
   pchi.probv.rst <- probv.pchi(mean.theta.slice,true.theta);

   list(all.theta=pv.theta, mean.theta=mean.theta.slice, mse.theta=mse.slice, pchi.theta=pchi.probv.rst, all.eta=pv$post.eta);
}


prv.slice.v5 <- function(V.INX.ALL, emp.v.n, niter, nburn, eta.init, eta.lowbound, eta.upbound, true.theta) { 
  v.loglinear.dat <- as.data.frame(cbind(V.INX.ALL, emp.v.n));
  colnames(v.loglinear.dat) <- c("X1","X2","X3","X4","X5","Freq");

  ### To make the sequence of fitted.values of this log linear model is the same with VR.INX.ALL
  ### need to fit the model with factor sequence as: R, X8, X7, X6, ..., X1
  attach(v.loglinear.dat);
  v.loglinear <- glm(Freq ~ .+.*. , family=poisson, data=v.loglinear.dat);
  ### fitted values gives the estimated frequence for each category, hence for the probability of each category,
  ### we need to divide by the total number of frequence
  ksi <- (v.loglinear$fitted.values)/sum(emp.v.n);

  #cat("finish estimation of ksi \n");

   pv <- slice.v.par.AMAR(ksi, emp.v.n, niter, eta.init, eta.lowbound, eta.upbound, true.theta);
   pv.theta <- pv$post.theta;

   mean.theta.slice <- apply(pv.theta[(nburn+1):niter,],2,mean);
   mse.slice <- t(true.theta-mean.theta.slice)%*%(true.theta-mean.theta.slice);
   pchi.probv.rst <- probv.pchi(mean.theta.slice,true.theta);

   list(all.theta=pv.theta, mean.theta=mean.theta.slice, mse.theta=mse.slice, pchi.theta=pchi.probv.rst, all.eta=pv$post.eta);
}

prv.slice.v4 <- function(V.INX.ALL, emp.v.n, niter, nburn, eta.init, eta.lowbound, eta.upbound, true.theta) { 
  v.loglinear.dat <- as.data.frame(cbind(V.INX.ALL, emp.v.n));
  colnames(v.loglinear.dat) <- c("X1","X2","X3","X4", "Freq");

  ### To make the sequence of fitted.values of this log linear model is the same with VR.INX.ALL
  ### need to fit the model with factor sequence as: R, X8, X7, X6, ..., X1
  attach(v.loglinear.dat);
  v.loglinear <- glm(Freq ~ .+.*. , family=poisson, data=v.loglinear.dat);
  ### fitted values gives the estimated frequence for each category, hence for the probability of each category,
  ### we need to divide by the total number of frequence
  ksi <- (v.loglinear$fitted.values)/sum(emp.v.n);

  #cat("finish estimation of ksi \n");

   pv <- slice.v.par.AMAR(ksi, emp.v.n, niter, eta.init, eta.lowbound, eta.upbound, true.theta);
   pv.theta <- pv$post.theta;

   mean.theta.slice <- apply(pv.theta[(nburn+1):niter,],2,mean);
   mse.slice <- t(true.theta-mean.theta.slice)%*%(true.theta-mean.theta.slice);
   pchi.probv.rst <- probv.pchi(mean.theta.slice,true.theta);

   list(all.theta=pv.theta, mean.theta=mean.theta.slice, mse.theta=mse.slice, pchi.theta=pchi.probv.rst, all.eta=pv$post.eta);
}


prv.slice.v7 <- function(V.INX.ALL, emp.v.n, niter, nburn, eta.init, eta.lowbound, eta.upbound, true.theta) { 
  v.loglinear.dat <- as.data.frame(cbind(V.INX.ALL, emp.v.n));
  colnames(v.loglinear.dat) <- c("X1","X2","X3","X4","X5","X6","X7", "Freq");

  ### To make the sequence of fitted.values of this log linear model is the same with VR.INX.ALL
  ### need to fit the model with factor sequence as: R, X8, X7, X6, ..., X1
  attach(v.loglinear.dat);
  v.loglinear <- glm(Freq ~ .+.*. , family=poisson, data=v.loglinear.dat);
  ### fitted values gives the estimated frequence for each category, hence for the probability of each category,
  ### we need to divide by the total number of frequence
  ksi <- (v.loglinear$fitted.values)/sum(emp.v.n);

  #cat("finish estimation of ksi \n");

   pv <- slice.v.par.AMAR(ksi, emp.v.n, niter, eta.init, eta.lowbound, eta.upbound, true.theta);
   pv.theta <- pv$post.theta;

   mean.theta.slice <- apply(pv.theta[(nburn+1):niter,],2,mean);
   mse.slice <- t(true.theta-mean.theta.slice)%*%(true.theta-mean.theta.slice);
   pchi.probv.rst <- probv.pchi(mean.theta.slice,true.theta);

   list(all.theta=pv.theta, mean.theta=mean.theta.slice, mse.theta=mse.slice, pchi.theta=pchi.probv.rst, all.eta=pv$post.eta);
}


############################################################################################################

### Used in prv.slice (updated on 07/02)
### For prv.slice method (use the ksi from fitting loglinear model + Dirichlet) --xi method
### AMAR case, Congenial Analysis for Eta using slice sampling
slice.v.par.AMAR <- function(ksi, vcounts, niter, init.eta, eta.lowbound, eta.upbound, true.theta) {

  #bound.eta <- sum(vcounts); 
  rst.eta <- NULL;
  rst.theta <- matrix(0, (niter+1), length(vcounts));

  ### current eta and theta
  c.eta <- init.eta;
  c.theta <- rdirichlet(1, vcounts+ksi/init.eta);

  rst.eta <- c.eta;
  rst.theta[1,] <- c.theta;
  
  for(s in 1:niter){

    ## Gibbs/Slice sampling for eta
    c.eta <- slice.onestep.AMAR(c.eta, ksi, c.theta, vcounts, eta.lowbound, eta.upbound);
    ## posterior sampling for theta
    c.theta <- rdirichlet(1, vcounts+ksi/c.eta);
    #cat("length of c.theta=", length(c.theta), "dim of rst.theta=", dim(rst.theta), "\n");
    rst.theta[s+1,] <- c.theta;
    rst.eta <- c(rst.eta, c.eta);
    m <- apply(rst.theta[1:(s+1),],2,mean);
    cat("finish iter =====", s, " with eta=", c.eta, ", and pchi===", probv.pchi(m,true.theta), ", and c.pchi=", probv.pchi(rst.theta[s+1,], true.theta), "\n");
  }

  list(post.eta=rst.eta, post.theta=rst.theta);
}
    

### Used in prv.slice (updated on 07/02)
### Gibbs/Slice sampling for eta
slice.onestep.AMAR <- function(curpar, ksi, theta, vcounts, lowbound, upbound, max.iter=10000){
  # bound for eta
  bound <- c(lowbound, upbound);
  x0 <- curpar;

  ## current f(x)
  log.fx0 <- eta.log.fx(ksi, theta, x0, vcounts);

 ### get auxiliary u ~ unif(0, f(x))
  log.u <- log.fx0-rexp(1);
  
  ### sampling from slice
  flag <- FALSE;
  n.iter <- 0;
  while(!flag & (n.iter<max.iter)) {
      x1 <- runif(1, bound[1], bound[2]);
      #cat("bound=", bound, "\n");
      log.fx1 <- eta.log.fx(ksi, theta, x1, vcounts);
      if (log.u < log.fx1) { flag <- TRUE; }
      else if (x1 < x0) { bound[1] <- x1; }
      else if (x1 > x0) { bound[2] <- x1; }
            n.iter <- n.iter+1;
    }

if(n.iter>=max.iter) { cat("slice failed with bounds as (", bound, "\n" ); }
x1;
}



### Used in prv.slice (updated on 07/02)
eta.log.fx <- function(ksi, theta, eta, vcounts){
  nttl <- sum(vcounts);

  logprior <- -2*log(eta*nttl+1);

  logtheta <- log(theta);
  ### for theta=0, bound from -Inf, assign a really large number -500
  logtheta[which(theta==0)]= -1000;

  loglike <- lgamma(1/eta)+sum(sapply(1:length(vcounts), function(x) -lgamma(ksi[x]/eta)))+sum(sapply(1:length(vcounts), function(x) logtheta[x]*(ksi[x]/eta-1))); 
  
  rst <- loglike+logprior;
  #cat("logprior=", logprior, "loglike=", loglike, "log(post.eta)=", rst, "\n");
  rst
}

############################################################################################################
############################################################################################################



###   ### obtain the theta, and mse.theta basing on Loglinear (only two way to save computation time) + Slice Sampling for theta
###   ### not use ksi from fitting loglinear model, but assume such a structure of two way interaction exists and the parameters unknown
###   ###   prv.slice.new.twoway <- function(vcounts, v.inx.all, niter=probv.NITER, nburn=probv.NBURN, inits, eta.upbound, sigma.lambda, v.cats, true.theta) {  
###     
###      pv <- slice.v.par.AMAR.twoway0702.v2(vcounts, niter, v.cats, inits=NULL, v.inx.all, sigma.lambda, eta.upbound, true.theta);
###      #pv.eta <- pv$post.eta;
###      pv.theta <- pv$post.theta;
###   
###      mean.theta.slice <- sapply(1:length(emp.v.n), function(x) mean(pv.theta[(nburn+1):niter,x]));
###      mse.slice <- t(true.theta-mean.theta.slice)%*%(true.theta-mean.theta.slice);
###      pchi.probv.rst <- probv.pchi(mean.theta.slice,true.theta);
###   
###      list(all.theta=pv.theta, mean.theta=mean.theta.slice, mse.theta=mse.slice, pchi.theta=pchi.probv.rst, all.eta=pv$post.eta);
###   }
###   

###   ### Used on prv.slice.new.twoway (updated on 07/02)
###   ##vall: the observed/simulated V data
###   ##niter: number of iterations
###   ##log.fx: loglikelihood function
###   ##inits: initial values
###   ##v.inx:categories for all possible V auxillary covariates (=V.CATS)
###   ##vall: is the observed V (=SIMU.V.RST)
###   slice.v.par.AMAR.twoway <- function(vcounts, niter, v.cats, inits=NULL, v.inx.all, sigma.lambda, eta.upbound, true.theta) {
###     
###       ##v categories
###       #ncats  <- mapply(length, v.cats); # v.inx = V.CATS, ncats gives the number of levels for each V variable
###       #ntheta <- prod(ncats); # length of theta
###       ntheta <- nrow(v.inx.all);
###       npv    <- length(v.cats); # number of covariates V=8
###       # number of (lambda's + eta)
###       nparlength <- choose(npv,1)+choose(npv,2)+1;
###   
###       ##init values
###       if (is.null(inits)) {
###           ## randomly choose lambda as initial value
###           inits <- c(inits, rnorm((nparlength-1), 0, sigma.lambda));        
###           ##eta init at 1/N=1/(sample size)
###           inits <- c(inits, 1/sum(vcounts));
###       }
###   
###       ##initialize results
###       rst.theta <- array(NA, dim=c(niter, ntheta));
###       rst       <- array(NA, dim=c(niter+1, nparlength));
###       rst[1, ]  <- inits;
###       for (s in 2:(niter+1)) {
###           ##browser();
###           rst[s,] <- rst[s-1,];
###           ##gibbs theta L
###                ### v.par.vec2lst.slice() transfer the parameter vector into 3 parts as (lambda1, labmda2, eta);
###           lambda.eta      <- v.par.vec2lst.slice(rst[s,], npv);
###           lambda1 <- lambda.eta$lambda1;
###           lambda2 <- lambda.eta$lambda2;
###           eta <- lambda.eta$eta;
###           #cat("current eta=", eta, "\n");
###              ### get.v.prob.slice() based on the lambda1, and lambda2 get the standardized ksi value (which sum to 1);
###           ksi <- get.v.prob.slice(lambda1, lambda2, v.inx.all);
###           cur.dir.par <- vcounts + ksi/eta;
###           next.theta  <- rdirichlet(1, cur.dir.par);
###           rst.theta[s-1, ] <- next.theta;
###   
###           ##Gibbs sampling for lambda's
###           for (i in 1:(nparlength-1)) {
###              slc.pi <- slice.onestep.new2(next.theta, rst[s,], npv, vcounts, sigma.lambda, v.inx.all, i, eta.upbound);
###              rst[s,i] <- slc.pi;
###              #cat("iter====",s, "finish slice sampling for par======", i, "\n");
###           }
###   
###           ##Gibbs/Slice sampling for eta
###           slc.eta <- slice.onestep.new2(next.theta, rst[s,], npv, vcounts, sigma.lambda, v.inx.all, nparlength, eta.upbound);
###           rst[s,nparlength] <- slc.eta;
###           if(s>5){
###           m <- apply(rst.theta[1:(s-1),],2,mean);
###           cat("finish iter====",s, "with eta=", slc.eta, ", with pchi===", probv.pchi(m,true.theta), ", *** c.pchi=", probv.pchi(rst.theta[s-1,],true.theta), "\n");
###         }
###       }
###   
###       list(post.eta=rst[,nparlength], post.theta=rst.theta);
###   }
###   



###   ### Updated on 07/02
###   slice.v.par.AMAR.twoway0702 <- function(vcounts, niter, v.cats, inits=NULL, v.inx.all, sigma.lambda, eta.upbound, true.theta) {
###    
###       ntheta <- nrow(v.inx.all);
###       npv    <- length(v.cats); # number of covariates V=8
###       # number of (lambda's + eta)
###       nparlength <- choose(npv,1)+choose(npv,2)+1;
###   
###       ##init values
###       if (is.null(inits)) {
###           ## randomly choose lambda as initial value
###           inits <- c(inits, rnorm((nparlength-1), 0, sigma.lambda));        
###           ##eta init at 1/N=1/(sample size)
###           inits <- c(inits, 1);
###       }
###   
###       ##initialize results
###       rst.theta <- array(NA, dim=c(niter, ntheta));
###       rst       <- array(NA, dim=c(niter+1, nparlength));
###       rst[1, ]  <- inits;
###           for (s in 2:(niter+1)) {
###           rst[s,] <- rst[s-1,];
###           cat("begin iter==", s, "finish update ");
###           #update theta
###           ### v.par.vec2lst.slice() transfer the parameter vector into 3 parts as (lambda1, labmda2, eta);
###           lambda.eta      <- v.par.vec2lst.slice(rst[s,], npv);
###           lambda1 <- lambda.eta$lambda1;
###           lambda2 <- lambda.eta$lambda2;
###           eta <- lambda.eta$eta;
###              ### get.v.prob.slice() based on the lambda1, and lambda2 get the standardized ksi value (which sum to 1);
###           ksi <- get.v.prob.slice(lambda1, lambda2, v.inx.all);
###           cur.dir.par <- vcounts + ksi/eta;
###           next.theta  <- rdirichlet(1, cur.dir.par);
###           rst.theta[s-1, ] <- next.theta;
###           
###           ##Gibbs sampling for lambda's
###           for(i in 1:npv) {
###              slc.pi <- slice.onestep.new0702(next.theta, rst[s,], lambda1, lambda2, eta, 2, i, NULL, npv, vcounts, sigma.lambda, v.inx.all, i, eta.upbound);
###              lambda1[i] <- slc.pi;
###              rst[s,i] <- slc.pi;
###              cat( i, ", ");
###           }
###   
###           prevn <- npv;
###           for(i in 1:(npv-1)) {
###               for(j in (i+1):npv) {
###                  slc.pi <- slice.onestep.new0702(next.theta, rst[s,], lambda1, lambda2, eta, 3, i, j, npv, vcounts, sigma.lambda, v.inx.all, i, eta.upbound);
###                  lambda2[i,j] <- slc.pi;
###                  rst[s, prevn+1] <- slc.pi;
###                  prevn <- prevn+1;
###               }
###               cat( "(", i, j, ")," );
###          }     
###   
###           ##Gibbs/Slice sampling for eta
###           slc.eta <- slice.onestep.new0702(next.theta, rst[s,], lambda1, lambda2, eta, 1, NULL, NULL, npv, vcounts, sigma.lambda, v.inx.all, nparlength, eta.upbound);
###           eta <- slc.eta;
###           rst[s, nparlength] <- slc.eta;
###           
###           if(s>5){
###           m <- apply(rst.theta[1:(s-1),],2,mean);
###           cat("finish iter====",s, "with eta=", slc.eta, ", with pchi===", probv.pchi(m,true.theta), ", *** c.pchi=", probv.pchi(rst.theta[s-1,],true.theta), "\n");
###          }      
###       }
###       
###       list(post.eta=rst[,nparlength], post.theta=rst.theta);
###   }
###   
###   

###   ### Updated on 07/02, with intercept added
###   slice.v.par.AMAR.twoway0702.v2 <- function(vcounts, niter, v.cats, inits=NULL, v.inx.all, sigma.lambda, eta.upbound, true.theta) {
###    
###       ntheta <- nrow(v.inx.all);
###       npv    <- length(v.cats); # number of covariates V=8
###       # number of (intercept+lambda's + eta)
###       nparlength <- 1+choose(npv,1)+choose(npv,2)+1;
###   
###       ##init values
###       if (is.null(inits)) {
###           ## randomly choose lambda as initial value
###           inits <- c(inits, rnorm((nparlength-1), 0, sigma.lambda));        
###           ##eta init at 1/N=1/(sample size)
###           inits <- c(inits, 1);
###       }
###   
###       ##initialize results
###       rst.theta <- array(NA, dim=c(niter, ntheta));
###       rst       <- array(NA, dim=c(niter+1, nparlength));
###       rst[1, ]  <- inits;
###           for (s in 2:(niter+1)) {
###           rst[s,] <- rst[s-1,];
###           cat("begin iter==", s, "finish update ");
###           #update theta
###           ### v.par.vec2lst.slice() transfer the parameter vector into 3 parts as (lambda1, labmda2, eta);
###           lambda.eta      <- v.par.vec2lst.slice.v2(rst[s,], npv);
###           
###           lambda0 <- lambda.eta$lambda0;
###           lambda1 <- lambda.eta$lambda1;
###           lambda2 <- lambda.eta$lambda2;
###           eta <- lambda.eta$eta;
###           
###              ### get.v.prob.slice() based on the lambda1, and lambda2 get the standardized ksi value (which sum to 1);
###           ksi <- get.v.prob.slice.v2(lambda0, lambda1, lambda2, v.inx.all);
###           cur.dir.par <- vcounts + ksi/eta;
###           next.theta  <- rdirichlet(1, cur.dir.par);
###           rst.theta[s-1, ] <- next.theta;
###           
###           ##Gibbs sampling for lambda's
###               slc.pi <- slice.onestep.new0702.v2(next.theta, rst[s,], lambda0, lambda1, lambda2, eta, 1, 0, NULL, npv, vcounts, sigma.lambda, v.inx.all, eta.upbound);
###              lambda0 <- slc.pi;
###              rst[s,1] <- slc.pi;
###   
###           for(i in 1:(npv)) {
###              slc.pi <- slice.onestep.new0702.v2(next.theta, rst[s,], lambda0, lambda1, lambda2, eta, 2, i, NULL, npv, vcounts, sigma.lambda, v.inx.all, eta.upbound);
###              lambda1[i] <- slc.pi;
###              rst[s,i+1] <- slc.pi;
###              #cat( i, ", ");
###           }
###   
###           prevn <- npv+1;
###           for(i in 1:(npv-1)) {
###               for(j in (i+1):npv) {
###                  slc.pi <- slice.onestep.new0702.v2(next.theta, rst[s,], lambda0, lambda1, lambda2, eta, 3, i, j, npv, vcounts, sigma.lambda, v.inx.all, eta.upbound);
###                  lambda2[i,j] <- slc.pi;
###                  rst[s, prevn+1] <- slc.pi;
###                  prevn <- prevn+1;
###               }
###               #cat( "(", i, j, ")," );
###          }     
###   
###           ##Gibbs/Slice sampling for eta
###           slc.eta <- slice.onestep.new0702.v2(next.theta, rst[s,], lambda0, lambda1, lambda2, eta, 0, NULL, NULL, npv, vcounts, sigma.lambda, v.inx.all, eta.upbound);
###           eta <- slc.eta;
###           rst[s, nparlength] <- slc.eta;
###   
###           if(s<=5) { cat("finish iter===", s, "\n"); }
###           if(s>5){
###           m <- apply(rst.theta[1:(s-1),],2,mean);
###           cat("finish iter====",s, "with eta=", slc.eta, ", with pchi===", probv.pchi(m,true.theta), ", *** c.pchi=", probv.pchi(rst.theta[s-1,],true.theta), "\n");
###          }      
###       }
###       
###       list(post.eta=rst[,nparlength], post.theta=rst.theta);
###   }
###   
###   



###   
###   ### the loglikelihood adjust to new slice
###   eta.log.fx.slice <- function(curpar, theta, nvcat, vcounts, sigma.lambda, V.INX.ALL){
###     nttl <- sum(vcounts);
###     npar <- length(curpar);
###     lambdavec <- curpar[1:(npar-1)];
###     
###     lambda.rst <- v.par.vec2lst.slice(curpar,nvcat);
###     lambda1 <- lambda.rst$lambda1;
###     lambda2 <- lambda.rst$lambda2;
###     eta <- lambda.rst$eta;
###   
###     ksi <- get.v.prob.slice(lambda1, lambda2, V.INX.ALL);
###     
###     logprioreta <- -2*log(eta*nttl+1);
###     logpriorlambda <- sum(-lambdavec^2/(2*sigma.lambda^2));
###   
###     logprior <- logprioreta+logpriorlambda;
###   
###     logtheta <- log(theta);
###     ### for theta=0, bound from -Inf, assign a really large number -500
###     logtheta[which(theta==0)]= -500;
###   
###     loglike <- lgamma(1/eta)+sum(sapply(1:length(vcounts), function(x) -lgamma(ksi[x]/eta)))+sum(sapply(1:length(vcounts), function(x) logtheta[x]*(ksi[x]/eta-1))); 
###   
###     rst <- loglike+logprior;
###     #cat("logprior=", logprior, "loglike=", loglike, "log(post.eta)=", rst, "\n");
###     rst
###   }
###   


##########
##one iteration of slice sampling
##parinx is the current parameter for slice sampling
##fix.inx is the one represents 1-the rest
##if fix.inx is NULL, that's for ETA
##nvcats is number of categories for all auxillary covariates
slice.onestep.new <- function(theta, curpar, nvcat, vcounts, sigma.lambda, v.inx.all, parinx, eta.upbound, max.iter=80000) {
    nlengthpar <- length(curpar);
    ### initially only bound for eta is non-negative, lambda's don't have bounds
    if(parinx<nlengthpar) { bound <- c(-100, 100); 
    }else { bound <- c(0, eta.upbound);}

    init.bound <- bound;
     
    # par gives which parameter are being slice sampled
    x0    <- curpar[parinx];

    ##get log f(x0)
    log.fx0 <- eta.log.fx.slice(curpar, theta, nvcat, vcounts, sigma.lambda, v.inx.all);
    ##get auxiliary y~unif(0,fx0), that is
    log.y   <- log.fx0 - rexp(1);

    ##sampling from slice
    flag   <- FALSE;
    n.iter <- 0;
    while (!flag & (n.iter<max.iter)) {
        x1 <- runif(1, bound[1], bound[2]);
        curpar[parinx] <- x1;
        log.fx1 <- eta.log.fx.slice(curpar, theta, nvcat, vcounts, sigma.lambda, v.inx.all);
         #cat("c.log=", log.y, "p.log=", log.fx1, "\n");

        if (log.y < log.fx1) {
            flag <- TRUE;
        } else if (x1 < x0) {
            bound[1] <- max(init.bound[1],x1);
        } else if (x1 > x0) {
            bound[2] <- min(x1,init.bound[2]);
        }
       # if(parinx==nlengthpar) { cat("bound for eta=", bound, "\n"); }
        n.iter <- n.iter + 1;
    }

    if (n.iter >= max.iter) {
        print("slice failed");
    }
    x1;
}


###   ### update on 07/02
###   ### change the vector of lambdas into lambda1, lambda2, eta
###   v.par.vec2lst.slice <- function(lambdapar, nv){
###     tmp.start <- 1;
###   
###     lengthinter <- choose(nv,1);
###     tmp.end <- tmp.start+lengthinter-1;
###     lambda1 <- NULL;
###     lambda1 <- lambdapar[tmp.start:tmp.end];
###   
###     tmp.start <- tmp.end+1;
###     # obtain number of two-way interactions
###     lengthinter <- choose(nv,2);
###     tmp.end <- tmp.start+lengthinter-1;
###     lambda2vec <- lambdapar[tmp.start:tmp.end];
###     indx <- 1;
###     lambda2 <- matrix(0, nv, nv);
###     for(i in 1:(nv-1)) {
###       for(j in (i+1):nv) {
###         lambda2[i,j] <- lambda2vec[indx];
###         indx <- indx+1;
###         # cat("for(", i, j, "), lambd2==", lambda2[i,j], ", indx==", indx, "*****");
###       }
###      }
###     for(i in 2:8) {
###        for(j in 1:(i-1)) {
###             lambda2[i,j] <- lambda2[j,i];
###         }
###     }   
###     eta <- lambdapar[tmp.end+1];
###     list(lambda1=lambda1, lambda2=lambda2, eta=eta);
###   }
###   

###   ### update on 07/02-- with intercept lambda added
###   ### change the vector of lambdas into lambda1, lambda2, eta
###   v.par.vec2lst.slice.v2 <- function(lambdapar, nv){
###     lambda0 <- lambdapar[1];
###     
###     tmp.start <- 2;
###     lengthinter <- choose(nv,1);
###     tmp.end <- tmp.start+lengthinter-1;
###     lambda1 <- NULL;
###     lambda1 <- lambdapar[tmp.start:tmp.end];
###   
###     tmp.start <- tmp.end+1;
###     # obtain number of two-way interactions
###     lengthinter <- choose(nv,2);
###     tmp.end <- tmp.start+lengthinter-1;
###     lambda2vec <- lambdapar[tmp.start:tmp.end];
###     indx <- 1;
###     lambda2 <- matrix(0, nv, nv);
###     for(i in 1:(nv-1)) {
###       for(j in (i+1):nv) {
###         lambda2[i,j] <- lambda2vec[indx];
###         indx <- indx+1;
###         # cat("for(", i, j, "), lambd2==", lambda2[i,j], ", indx==", indx, "*****");
###       }
###      }
###     for(i in 2:nv) {
###        for(j in 1:(i-1)) {
###             lambda2[i,j] <- lambda2[j,i];
###         }
###     }   
###     eta <- lambdapar[tmp.end+1];
###     list(lambda0=lambda0, lambda1=lambda1, lambda2=lambda2, eta=eta);
###   }
###   

###   ## updated on 07/02 
###   v.par.lst2vec.slice <- function(lambda1, lambda2, eta) {
###   
###        nv <- length(lambda1);
###        rst <- c(lambda1);
###        for(i in 1:(nv-1)) {
###           for(j in (i+1):nv) {
###              rst <- c(rst, lambda2[i,j]);
###           }
###        }
###        rst <- c(rst,eta);
###        rst;
###   }
###   
###   ### udpate on 07/02 -- with intercept of lambda added
###    v.par.lst2vec.slice.v2 <- function(lambda0, lambda1, lambda2, eta) {
###   
###        nv <- length(lambda1);
###        rst <- c(lambda0, lambda1);
###        for(i in 1:(nv-1)) {
###           for(j in (i+1):nv) {
###              rst <- c(rst, lambda2[i,j]);
###           }
###        }
###        rst <- c(rst,eta);
###        rst;
###   }
###   

###   #checked on 07/02
###   ### based on lambda1, lambda2, lambda3, get estimated \xi
###   get.v.prob.slice <- function(lambda1, lambda2, V.INX.ALL){
###      muvall <- NULL;
###      pvall <- NULL;
###      nvall <- nrow(V.INX.ALL);
###   
###     for(i in 1:nvall) {
###     # find the index of V having value=1
###     index1 <- which(V.INX.ALL[i,]==1);
###     
###     # main factor sum
###     oneinter <- sum(lambda1[index1]);
###   
###     # second order interaction sum
###     lengthv1 <- length(index1);
###     twointer <- 0;
###     if(lengthv1 < 2) { twointer <- 0;}
###     else
###     {  for(j in 1:(lengthv1-1))
###        {
###           for(k in (j+1):lengthv1) { twointer=twointer+lambda2[index1[j], index1[k]];}
###         }
###      }
###     
###     muvall[i] <- exp(oneinter+twointer);
###     } 
###   
###     # pvall gives the probability of P(v), by standardizing muvall
###     pvall <- muvall/sum(muvall);
###     pvall
###   }
###   
###   
###    ### upated on 07/02 -- for case with intercept lambda added
###    ### based on lambda1, lambda2, lambda3, get estimated \xi
###   get.v.prob.slice.v2 <- function(lambda0, lambda1, lambda2, V.INX.ALL){
###      muvall <- NULL;
###      pvall <- NULL;
###      nvall <- nrow(V.INX.ALL);
###   
###     for(i in 1:nvall) {
###     # find the index of V having value=1
###     index1 <- which(V.INX.ALL[i,]==1);
###     
###     # main factor sum
###     oneinter <- sum(lambda1[index1]);
###   
###     # second order interaction sum
###     lengthv1 <- length(index1);
###     twointer <- 0;
###     if(lengthv1 < 2) { twointer <- 0;}
###     else
###     {  for(j in 1:(lengthv1-1))
###        {
###           for(k in (j+1):lengthv1) { twointer=twointer+lambda2[index1[j], index1[k]];}
###         }
###      }
###     
###     muvall[i] <- exp(lambda0+oneinter+twointer);
###     } 
###   
###     # pvall gives the probability of P(v), by standardizing muvall
###     pvall <- muvall/sum(muvall);
###     pvall
###   }
###       
###   
###   slice.onestep.new2 <- function(theta, curpar, nvcat, vcounts, sigma.lambda, v.inx.all, parinx, eta.upbound, max.iter=100) {
###       nlengthpar <- length(curpar);
###       ### initially only bound for eta is non-negative, lambda's don't have bounds
###       if(parinx<nlengthpar) { bound <- c(-50, 50); 
###       }else { bound <- c(0, eta.upbound);}
###   
###       init.bound <- bound;
###   
###     lambda.rst <- v.par.vec2lst.slice(curpar,nvcat);
###     lambda1 <- lambda.rst$lambda1;
###     lambda2 <- lambda.rst$lambda2;
###     eta <- lambda.rst$eta;
###   
###        
###       # par gives which parameter are being slice sampled
###       x0    <- curpar[parinx];
###       ##get log f(x0)
###       log.fx0 <- eta.log.fx.slice2(curpar,lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, v.inx.all);
###       ##get auxiliary y~unif(0,fx0), that is
###       log.y   <- log.fx0 - rexp(1);
###   
###       ##sampling from slice
###       l1index <- choose(8,1);
###       l2index <- choose(8,2)+l1index;
###       #l3index <- nlengthpar-1;
###   
###       ## identify which lambda to be updated
###       upindex <- 0;
###       flaginx <- FALSE;
###       
###        if(parinx==nlengthpar) { upindex <- 3;
###        }else if(parinx<=l1index) {
###          upindex <- 1;
###         }else {
###          upindex <- 2;
###             ### parinx for lambda2 update
###             indexij <- 1;
###             for(i in 1:(nvcat-1)){
###               for(j in (i+1):nvcat)
###                 { if(indexij==parinx-l1index) {
###                   #cat("equal");
###                   flaginx=TRUE;
###                   break
###                   }else {indexij=indexij+1;}
###                 }
###               if(flaginx==TRUE) break
###             }
###          }
###       #cat("update lambda=", upindex, "\n");
###   
###          n.iter <- 0;
###          flag   <- FALSE;
###          while (!flag & (n.iter<max.iter)) {
###           x1 <- runif(1, bound[1], bound[2]);
###           curpar[parinx] <- x1;
###           if(upindex==1) { lambda1[parinx] <- x1;
###           }else if(upindex==2) { lambda2[i,j] <- x1;
###           }else {eta <- x1; }
###          
###           log.fx1 <- eta.log.fx.slice2(curpar,lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, v.inx.all);
###       
###            #cat("c.log=", log.y, "p.log=", log.fx1, "\n");
###   
###           if (log.y < log.fx1) {
###               flag <- TRUE;
###           } else if (x1 < x0) {
###               bound[1] <- max(init.bound[1],x1);
###           } else if (x1 > x0) {
###               bound[2] <- min(x1,init.bound[2]);
###           }
###           #if(parinx==nlengthpar) { cat("bound for eta=", bound, "\n"); }
###           #cat("bound for par ==", parinx, "=", bound, "\n");
###           n.iter <- n.iter + 1;
###       }
###   
###       if (n.iter >= max.iter) {
###           print("slice failed");
###       }
###       #list(lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, eta=eta, cpar=x1);
###       x1;
###   }
###   
###   




###  ### Updated on 07/02
###  slice.onestep.new0702 <- function(theta, curpar, lambda1, lambda2, eta, inx1, inx2, inx3, nvcat, vcounts, sigma.lambda, v.inx.all, parinx, eta.upbound, max.iter=10000) {
###      nlengthpar <- length(curpar);
###      ### initially only bound for eta is non-negative, lambda's don't have bounds
###      if(parinx<nlengthpar) { bound <- c(-50, 50); 
###      }else { bound <- c(0, eta.upbound);}
###  
###      #get the original value of updating variable
###      if(inx1==1) { # update eta
###          x0=eta;
###       }else if(inx1==2) { ### update lambda1 
###          x0=lambda1[inx2];
###       }else{
###          x0=lambda2[inx2,inx3];
###       }   
###      
###      ##get log f(x0)
###      log.fx0 <- eta.log.fx.slice2(curpar, lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, v.inx.all);
###      ##get auxiliary y~unif(0,fx0), that is
###      log.y   <- log.fx0 - rexp(1);
###  
###         n.iter <- 0;
###         flag   <- FALSE;
###         while (!flag & (n.iter<max.iter)) {
###          x1 <- runif(1, bound[1], bound[2]);
###          if(inx1==1) {
###              eta=x1;
###          }else if(inx1==2) {
###               lambda1[inx2] <- x1;
###          }else{
###               lambda2[inx2, inx3] <- x1;
###          }
###          # the updated parameter list 
###          new.par <- v.par.lst2vec.slice(lambda1, lambda2, eta);
###          log.fx1 <- eta.log.fx.slice2(new.par,lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, v.inx.all);    
###           #cat("c.log=", log.y, "p.log=", log.fx1, "\n");
###  
###          if (log.y < log.fx1) {
###              flag <- TRUE;
###          } else if (x1 < x0) {
###              bound[1] <- max(bound[1],x1);
###          } else if (x1 > x0) {
###              bound[2] <- min(x1,bound[2]);
###          }
###          #if(parinx==nlengthpar) { cat("bound for eta=", bound, "\n"); }
###          #cat("bound for par ==", parinx, "=", bound, "\n");
###          n.iter <- n.iter + 1;
###      }
###  
###      if (n.iter >= max.iter) {
###          print("slice failed");
###      }
###      #list(lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, eta=eta, cpar=x1);
###      x1;
###  }
###  
###  
###      ### Updated on 07/02 -- for case with intercept lambda0 added
###  slice.onestep.new0702.v2 <- function(theta, curpar, lambda0,  lambda1, lambda2, eta, inx1, inx2, inx3, nvcat, vcounts, sigma.lambda, v.inx.all, eta.upbound, max.iter=10000) {
###      nlengthpar <- length(curpar);
###      ### initially only bound for eta is non-negative, lambda's don't have bounds
###      if(inx1>0) { bound <- c(-50, 50); 
###      }else { bound <- c(0, eta.upbound);}
###  
###      #get the original value of updating variable
###      if(inx1==0) { ### update lambda0
###          x0= eta;
###       }else if(inx1==1) { # update eta
###          x0=lambda0;
###       }else if(inx1==2) { ### update lambda1 
###          x0=lambda1[inx2];
###       }else{
###          x0=lambda2[inx2,inx3];
###       }   
###      
###      ##get log f(x0)
###      log.fx0 <- eta.log.fx.slice2.v2(curpar, lambda0, lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, v.inx.all);
###      ##get auxiliary y~unif(0,fx0), that is
###      log.y   <- log.fx0 - rexp(1);
###  
###         n.iter <- 0;
###         flag   <- FALSE;
###         while (!flag & (n.iter<max.iter)) {
###          x1 <- runif(1, bound[1], bound[2]);
###          if(inx1==0){
###              eta <- x1;
###          }else if(inx1==1) {
###              lambda0 <- x1;
###          }else if(inx1==2) {
###               lambda1[inx2] <- x1;
###          }else{
###               lambda2[inx2, inx3] <- x1;
###          }
###          # the updated parameter list 
###          new.par <- v.par.lst2vec.slice.v2(lambda0, lambda1, lambda2, eta);
###          log.fx1 <- eta.log.fx.slice2.v2(new.par,lambda0, lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, v.inx.all);    
###           #cat("c.log=", log.y, "p.log=", log.fx1, "\n");
###  
###          if (log.y < log.fx1) {
###              flag <- TRUE;
###          } else if (x1 < x0) {
###              bound[1] <- max(bound[1],x1);
###          } else if (x1 > x0) {
###              bound[2] <- min(x1,bound[2]);
###          }
###          n.iter <- n.iter + 1;
###      }
###  
###      if (n.iter >= max.iter) {
###          print("slice failed");
###      }
###      #list(lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, eta=eta, cpar=x1);
###      x1;
###  }
###  

###  ### the loglikelihood adjust to new slice
###  eta.log.fx.slice2 <- function(curpar, lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, V.INX.ALL){
###    nttl <- sum(vcounts);
###    npar <- length(curpar);
###    lambdavec <- curpar[1:(npar-1)];
###    
###    ksi <- get.v.prob.slice(lambda1, lambda2, V.INX.ALL);
###    
###    logprioreta <- -2*log(eta*nttl+1);
###    logpriorlambda <- sum(-lambdavec^2/(2*sigma.lambda^2));
###  
###    logprior <- logprioreta+logpriorlambda;
###  
###    ### for theta=0, bound from -Inf, assign a really large number -500
###    logtheta <- log(theta);
###    logtheta[which(theta==0)]= -500;
###  
###    loglike <- lgamma(1/eta)+sum(sapply(1:length(vcounts), function(x) -lgamma(ksi[x]/eta)+logtheta[x]*(ksi[x]/eta+vcounts[x]-1))); 
###  
###    rst <- loglike+logprior;
###    rst
###  }
###  
###  
###   ### updated on 07/02 - for case with intercept lambda0 added
###    ### the loglikelihood adjust to new slice
###  eta.log.fx.slice2.v2 <- function(curpar, lambda0, lambda1, lambda2, eta, theta, nvcat, vcounts, sigma.lambda, V.INX.ALL){
###    nttl <- sum(vcounts);
###    npar <- length(curpar);
###    lambdavec <- curpar[1:(npar-1)];
###    
###    ksi <- get.v.prob.slice.v2(lambda0, lambda1, lambda2, V.INX.ALL);
###    
###    logprioreta <- -2*log(eta*nttl+1);
###    logpriorlambda <- sum(-lambdavec^2/(2*sigma.lambda^2));
###  
###    logprior <- logprioreta+logpriorlambda;
###  
###    ### for theta=0, bound from -Inf, assign a really large number -500
###    logtheta <- log(theta);
###    logtheta[which(theta==0)]= -500;
###  
###    loglike <- lgamma(1/eta)+sum(sapply(1:length(vcounts), function(x) -lgamma(ksi[x]/eta)+logtheta[x]*(ksi[x]/eta+vcounts[x]-1))); 
###  
###    rst <- loglike+logprior;
###    rst
###  }
###  ############################################################################################################
############################################################################################################
############################################################################################################

# Estimation of theta as ksi using log linear model
ksi.loglinear <- function(V.INX.ALL, emp.v.n) {

  v.loglinear.dat <- as.data.frame(cbind(V.INX.ALL, emp.v.n));
  colnames(v.loglinear.dat) <- c("X1","X2","X3","X4","X5","X6","X7","X8","Freq");

  ### To make the sequence of fitted.values of this log linear model is the same with VR.INX.ALL
  ### need to fit the model with factor sequence as: R, X8, X7, X6, ..., X1
  attach(v.loglinear.dat);
  v.loglinear <- glm(Freq ~ factor(X8)+factor(X7)+factor(X6)+factor(X5)+factor(X4)+factor(X3)+factor(X2)+factor(X1), family=poisson, data=v.loglinear.dat);
  ### the fitted values gives the estimated frequence for each category, hence for the probability of each category,
  ### we need to divide by the total number of frequence
  ksi <- (v.loglinear$fitted.values)/sum(emp.v.n);
  
  ksi
}


####### Estimating theta using Dirichlet() method from "Nonparametric Bayes Modeling of Multivariate Categorical Data" by Dunson and Xing


### Update on 06/08
### exactly the same as old.cal.omega when last element of v=1, but much simpler
cal.omega  <- function(n.stick, v) {
 p <- numeric(n.stick);
 p[1] <- v[1];
 if(n.stick>1) {
  p[2:n.stick] <- sapply(2:n.stick, function(x) v[x]*prod(1-v[1:(x-1)]));
 }
 p
} 


#### sampling v from beta(1, alpha) distribution with in the low and up bound 
samp.ind.v <- function(lowbound, upbound, alpha, max.try=80000) {
  try <- 1;
  failed <- 0;
  if(upbound-lowbound<0.01) {
    tmp <- (upbound+lowbound)/2;
  }else{  
         tmp <- rbeta(1,1, alpha);
         while(try<max.try & (tmp<lowbound | tmp>upbound)) {
            tmp <- rbeta(1, 1, alpha);
            try <- try+1;
         }
  }
  if(try==max.try)  { tmp <- (lowbound+upbound)/2; failed=1; }
  list(samp.v=tmp, failed=failed);
}



### The most recent Function of Dirichlet(), based on rewriting the Dir1002()
### For Dir1002(), it takes a hugh memory since it records much more data, especially large matrix
### for estimation of Pr(V) using "Nonparametric Bayes Modeling of Multivariate Categorical Data"
###by Dunson and Xing
xuan.pvdirichlet.fast.1127 <- function(k, data.i, niter, nburn, max.h, simu.v, a.alpha, b.alpha, V.INX.ALL, V.CATS, true.theta){

  n.v <- ncol(simu.v);
  n.obs <- nrow(simu.v);
  pchi.theta <- 0;
  
  ## save for each iteration, 
  rst.thai <- matrix(0, max.h, n.v);
  rst.u <- array(0,n.obs);
  rst.z <- array(0,n.obs);
  rst.v <- array(0,max.h);
  rst.theta <- matrix(0, niter-1, nrow(V.INX.ALL));
  n.hitmax <- 0;
                    
  ###########################  Begin initializing values  ########################################################
  

  ### initialize alpha first based on paper  alpha ~ gamma(a.alpha, b.alpha)
  #rst.alpha <- rgamma(1, a.alpha, b.alpha);
  rst.alpha <- 0.5;
  rst.v <- numeric(max.h);
  rst.v <- rbeta(max.h, 1, rst.alpha);
  rst.v[max.h] <- 1; # to make the sum of sticks equals to 1
  tmp.omega <- cal.omega(max.h, rst.v);

  ### Based on rst.z, rst.v, rst.w (tmp.omega) generate rst.u
  rst.u <- runif(n.obs, 0, 1/(max.h));
  rst.z <- sapply(1:n.obs, function(x) sample(which(rst.u[x]<tmp.omega),1));
  max.z <- max(rst.z);

     z.inx.list <- list();
     for(i in 1:max.z) {
        z.inx.list[[i]] <- which(rst.z==i);
      }

    for(h.iter in 1:max.z) {
      for(j in 1:n.v) {
         # in our case each covariate having 2 level, hence only need to save the 1st probability
         n.vhjeq1 <- sum(simu.v[ z.inx.list[[h.iter]], j ]);
         rst.thai[h.iter,j] <- rdirichlet(1,c(1+length( z.inx.list[[h.iter]] )-n.vhjeq1, 1+n.vhjeq1))[1];
      }
    }
   if(max.h>max.z) {
     s.length <- (max.h-max.z)*n.v;
     tmp.thai <- matrix(rdirichlet(s.length,c(1,1))[,1], (max.h-max.z), n.v);
     rst.thai[(max.z+1):(max.h),1:n.v] <- tmp.thai;
   }
     
       #### get the corresponding v.index number for each subjects
     sub.indx <- apply(simu.v,1, function(x) vtosimvalue(x,V.CATS))+1;

  ############################ Finish Initialize ###################################################################
   
  ##### Begin Iteration
  for(s in 2:niter) {
   max.z <- max(rst.z);
   ### step 1 update u
    tmp.omega <- cal.omega(max.h, rst.v);
    rst.u <- sapply(1:n.obs, function(x) runif(1, 0, tmp.omega[rst.z[x]]));
    #cat("\n finishi (1)");

   ### step 2 update thai
    # make a list such that the h^th element contains the index of subjects with z==h
     z.inx.list <- list();
     for(i in 1:max.z) {
        z.inx.list[[i]] <- which(rst.z==i);
      }
   
   ### step 3 update v  #################################################################  
   if(max.z==1) {
     lowb <- max(rst.u);
     upb <- 1;
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[1] <- rst$samp.v;
     #failed <- rst$failed;
     #cat("for h=1, bound=(", lowb, ",", upb, ")"); 
     rst.v[2:(max.h)] <- rbeta((max.h-1), 1, rst.alpha);
     ############################################################ 
    }else if(max.z==2) {
     # for v1
       # if there is no subject with Z==1
       if(length(z.inx.list[[1]])==0) {
         lowb <- 0;
       }else{
        lowb <- max(rst.u[z.inx.list[[1]]]);
       }
       upb <- 1-max(rst.u[z.inx.list[[2]]])/rst.v[2];
       #cat("for h=1, bound=(", lowb, "," ,upb, ")"); 
       rst <- samp.ind.v(lowb, upb, rst.alpha);
       rst.v[1] <- rst$samp.v;
       #failed <- rst$failed;
       
     # for v2 
     lowb <- min(1,max(rst.u[ z.inx.list[[2]] ])/(1-rst.v[1]));
     upb <- 1;
     #cat("for h=2, bound=(", lowb, "," ,upb, ")");     
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[2] <- rst$samp.v;
     #failed <- failed+rst$failed;
     #if(failed>=1) {   failed <- 1; }
     ### for v3 and beyond 
     rst.v[3:(max.h)] <- rbeta((max.h-2), 1, rst.alpha);
    ############################################################ 
    }else{   ### for (max.z>=3)
    ###(1) for updating rst.v[1]
       ## (1.1) low bound for V_h always come from groups with Z==h
     if(length( z.inx.list[[1]] )==0) {
       lowb <- 0; }else{
       lowb <- max(rst.u[ z.inx.list[[1]] ]);
     }  
       ## (1.2) ## Find upbound for V_(h=1), coming from groups with Z= (h+1), (h+2), ... (max.h)
      # for group of Z equal to 2
     if(length( z.inx.list[[2]] )==0) {
       max1 <- 0;
     }else{
       max1 <- max(rst.u[ z.inx.list[[2]] ])/rst.v[2];
     }
       upb1 <- 1-max1;       
     ## (1.3)
     #  for group of Z great than 3 
     gt.h <- c(3:max.z);
     length.gt.h <- max.z-3+1;
     max.num <- numeric(length.gt.h);
     max.num <- sapply(gt.h, function(x)  max(0, rst.u[ z.inx.list[[x]] ]));  # update on 07/29
     denom.gt.h <- numeric(length.gt.h);
     denom.gt.h[1] <- rst.v[3]*(1-rst.v[2]);  # from group z=3
     # update on 7/29
     if(length.gt.h>=2) { # from group z=4, 5, ...
       denom.gt.h[2:length.gt.h] <- sapply(4:max.z, function(x) rst.v[x]*prod(sapply(2:(x-1), function(y) 1-rst.v[y])));
     } 
     max2 <- min(1,max(0,max.num/denom.gt.h));
     upb2 <- 1-max2;
     #cat("max2==", max2, "\n");
     #upb <- 1-max(max1,max2);
     upb <- min(upb1, upb2);
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[1] <- rst$samp.v;
     failed <- rst$failed;
     #cat("for h=1, bound=(", lowb, ",", upb, ")-max1=", max1, "-ma2=", max2, ", with rst.v[1]==", rst.v[1], "\n");     
   ###  End for updating rst.v[1] with max.z>=3
     
     
   ### (2) for updating rst.v[h], for h=2, ..., (max.z-1) with max.z>=3
     for(h.iter in 2:(max.z-1)) {
       # (2.0) for lower bound
       if( length(z.inx.list[[h.iter]] )==0) {
          lowb <- 0;
        }else{
         lowb <- max(rst.u[ z.inx.list[[h.iter]] ])/prod(sapply(1:(h.iter-1), function(x) 1-rst.v[x]));
        }
       # (2.1) find upper bound max1 for rst.v[h] from group with Z=h+1
       #indxeq.hiterplus1 <- which(rst.z==(h.iter+1));
       if(length( z.inx.list[[h.iter+1]]) ==0) {
         max1 <- 0;
       }else{
          #cat("numer=", max(rst.u[ z.inx.list[[h.iter+1]] ]), "***denom==", (rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x])))), "\n");
          max1 <- max(rst.u[ z.inx.list[[h.iter+1]] ])/(rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x]))));            }
       upb1 <- 1-max1;
       
       #(2.2) find upper bound max2 for rst.v[h] from group with Z=h+2, h+3, ...
       if(h.iter+1>=max.z) { # on the equality with h.iter=max.z-1, the upbound taken care by upb1
          max2 <- 0;
       }else {   
          ### Begin Alternative approach
          whichgt.h <- c((h.iter+2):(max.z));  ## {h*:h*>h.iter}
          lengthgt.h <- length(whichgt.h);
          
          num.gt.h <- numeric(lengthgt.h);
          num.gt.h <- sapply(whichgt.h, function(x) max(0, rst.u[ z.inx.list[[x]]  ]));
          indx.num.gt0 <- which(num.gt.h>0);
          num.gt0 <- num.gt.h[indx.num.gt0];
          denom.gt0 <- numeric(length(indx.num.gt0));
          
          denom.gt0 <- sapply(whichgt.h[indx.num.gt0], function(x) rst.v[x]*prod(sapply(1:(h.iter-1), function(y) 1-rst.v[y]))*prod(sapply((h.iter+1):(x-1), function(z)  1-rst.v[z])));
          
          #at("num.gt.h===", num.gt0, "\n");
          #cat("denom.gt.h===", denom.gt0, "\n");q
          max2 <- min(1,max(0, num.gt0/denom.gt0));
        }
          upb2 <- 1-max2;
          
          upb <- min(upb1, upb2);
          #cat("for h=", h.iter, ", bound=(", lowb, ",", upb, ")--max1=", max1, "--max2=", max2, ",");     
          rst<- samp.ind.v(lowb, upb, rst.alpha);
          rst.v[h.iter] <- rst$samp.v;
          #failed <- failed+rst$failed;
     }
   
     #(3) #### for h.iter==max.z ( hence definitely some subjects fall in this group)
     lowb <- max(rst.u[ z.inx.list[[max.z]] ])/prod(sapply(1:(max.z-1), function(x) 1-rst.v[x]));
     upb <- 1;
     #cat("for h=", max.z, ", bound=(", lowb, ",", upb, ")");    
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[max.z] <- rst$samp.v;
     #failed <- failed+rst$failed;
     #(4) update rst.v[h] for h=max.z+1, ..., max.h
     if((max.z+1)<=max.h){
     rst.v[(max.z+1):(max.h)] <- rbeta((max.h-max.z),1,rst.alpha);
     }
     #if(failed>=1) { failed <- 1; }
   } ## finish for max.z>=3

   ### keept last v=1 for sum length=1
   rst.v[max.h] <- 1;
   tmp.omga <- cal.omega(max.h, rst.v);
   rst.v[which(rst.v[1:(max.h-1)]==1)] <- 1-(10)^(-10); ### to exclude the case that middle part rst.v==1, hence (1-rst.v)=0 in the denominator on Sep17
   #cat(",  (3)");
   # after updating rst.v, need to update tmp.omega   
   ###########  Finish step 3  for updating V   ############### ########### ########### ########### ########### ########### ########### ###########  

  ### Update on 06/08-- to save time, first calculate the p(v|h) and use it in updating Z
         ind.h.prob <- matrix(0, max.h,nrow(V.INX.ALL));
         for(tmp.h in 1:max.h) {
            ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
         }
       
    ### step 4 
     # update omega first, need to calculate for max.h
     min.u <- min(rst.u);
     #cat("min.u==", min.u, ",   ");
     addlength <- sapply(1:(max.h), function(x) sum(tmp.omega[1:x]));
     k.tilt <- min(max.h, min(which(addlength>=(1-min.u))));
     cat(", with k.tilt==", k.tilt, "\n");

       for(i in 1:n.obs) {
         tmp.A <- which(tmp.omega[1:k.tilt]>rst.u[i]);
         #cat("for sub==", i, ", sample from(", tmp.A, ")\n");
         size.A <- length(tmp.A);
         # calculate the multinomial distribution probability
         if(size.A==0) {
          rst.z[i] <- rst.z[i];
         }else if(size.A==1){
          rst.z[i] <- tmp.A[1];
         }else{  
          prob.A <- numeric(size.A);
          prob.A <- ind.h.prob[tmp.A, sub.indx[i]];
          prob.A <- prob.A/sum(prob.A);
          rst.z[i] <- tmp.A[which(rmultinom(1,1,prob.A)==1)];
         }
      }
    #cat(",  (4)");
       #max.z <- max(rst.z);
       #if(max.z==max.h) n.hitmax <- n.hitmax+1;
       #cat("updated with z==", rst.z, " with max.z==", max(rst.z), "\n");
     
   #### step 4.5 -- only update in 0928 function()
   unique.z <- unique(rst.z);
   rank.z <- rank(unique.z);
   update.rstz <- sapply(1:n.obs, function(x) rank.z[which(unique.z==rst.z[x])]);
   update.tmp.omega <- c(tmp.omega[sort(unique.z)], tmp.omega[-unique.z]);
   rst.z <- update.rstz;
   tmp.omega <- update.tmp.omega;
      ### based on each individual length, update rst.v
      update.rstv <- numeric(max.h);
      update.rstv[1] <- tmp.omega[1];
      for(j in 2:max.h) {
         update.rstv[j] <- tmp.omega[j]/prod(sapply(1:(j-1), function(y) (1-update.rstv[y])));
      }   
    rst.v <- update.rstv; # added on 2015/02/12
   
  ### step 2 update thai
     # make a list such that the h^th element contains the index of subjects with z==h
     max.z <- max(rst.z);
     z.inx.list <- list();
     for(i in 1:max.z) {
        z.inx.list[[i]] <- which(rst.z==i);
      }
   
    for(h.iter in 1:max.z) {
      for(j in 1:n.v) {
         # in our case each covariate having 2 level, hence only need to save the 1st probability
         n.vhjeq1 <- sum(simu.v[ z.inx.list[[h.iter]], j ]);
         rst.thai[h.iter,j] <- rdirichlet(1,c(1+length( z.inx.list[[h.iter]] )-n.vhjeq1, 1+n.vhjeq1))[1];
      }
    }
    if(max.h>max.z) {
     s.length <- (max.h-max.z)*n.v;
     tmp.thai <- matrix(rdirichlet(s.length,c(1,1))[,1], (max.h-max.z), n.v);
     rst.thai[(max.z+1):(max.h),1:n.v] <- tmp.thai;
    }
    #cat(",  (2)");

    for(tmp.h in 1:max.h) {
       ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
    } 
    rst.theta[s-1,] <- t(ind.h.prob)%*%tmp.omega;
      
        ### step 5
      # update omega first
    tmp.shape <-  a.alpha+max.z;
     if(max.z==max.h) {
       inx.neq1 <- which(rst.v[1:(max.h)]<1);
       tmp.rate <- Inf;
       rst.alpha <- (10)^(-10);
     }else{
       inx.neq1 <- which(rst.v[1:(max.z)]<1);
       tmp.rate <- b.alpha-log(1-sum(tmp.omega[1:max.z]));
       rst.alpha <- rgamma(1, shape=tmp.shape, rate=tmp.rate);
     }
  
   #cat("alpha (a==", tmp.shape, ", b==", tmp.rate, ")\n");
    #cat(",  (5)");
     ### calculate the pchi after each iteration
    if(s>2) {
      if(s<=(nburn+2)) {
        m <- apply(rst.theta[1:(s-1),],2,mean);
      }else{
        m <- apply(rst.theta[nburn:(s-1),],2,mean);
      }  
    mean.pchi.theta <- probv.pchi(m, true.theta);
    crnt.pchi.theta <- probv.pchi(rst.theta[s-1,], true.theta);
      
    cat(" at End iteration ", s, ", pchi===( ", mean.pchi.theta, ",", crnt.pchi.theta, "), ****max.z===", max(rst.z), ", alpha=", rst.alpha, ", for (K,i)==(", k, ",", data.i, ") with n.hitmax==", n.hitmax, "\n");
   }
   #cat(",  (6)");
   #cat("updated with u==", rst.u, "\n");
   #cat("updated with v==", rst.v[1:max.z], "\n");
   #beyond.length <- c(beyond.length, 1-sum(tmp.omega[1:(max.z)]));
   #if(max.z==max.h) { n.hitmax <- n.hitmax+1; }
   
  } ### end of iteration
   list(theta=rst.theta, mean.theta=m);
  
}



###############################################################################################################


###   ##### It takes hugh memory space since several version of theta matrix (niter*length(theta)), better not use it 
###   xuan.pvdirichlet.fast.1002 <- function(k, data.i, niter, nburn, max.h, simu.v, a.alpha, b.alpha, V.INX.ALL, V.CATS, true.theta){
###   
###     n.v <- ncol(simu.v);
###     n.obs <- nrow(simu.v);
###     pchi.theta <- 0;
###     
###     ## save for each iteration, 
###     rst.thai <- matrix(0, max.h, n.v);
###     rst.u <- array(0,n.obs);
###     rst.z <- array(0,n.obs);
###     rst.v <- array(0,max.h);
###     rst.theta <- matrix(0, niter-1, nrow(V.INX.ALL));
###     rst.theta.rw <- matrix(0, niter-1, nrow(V.INX.ALL));
###       
###     n.hitmax <- 0;
###                       
###     ###########################  Begin initializing values  ########################################################
###     
###   
###     ### initialize alpha first based on paper  alpha ~ gamma(a.alpha, b.alpha)
###     #rst.alpha <- rgamma(1, a.alpha, b.alpha);
###     rst.alpha <- 0.5;
###     rst.v <- numeric(max.h);
###     rst.v <- rbeta(max.h, 1, rst.alpha);
###     rst.v[max.h] <- 1; # to make the sum of sticks equals to 1
###     tmp.omega <- cal.omega(max.h, rst.v);
###   
###     ### Based on rst.z, rst.v, rst.w (tmp.omega) generate rst.u
###     rst.u <- runif(n.obs, 0, 1/(max.h));
###     rst.z <- sapply(1:n.obs, function(x) sample(which(rst.u[x]<tmp.omega),1));
###     max.z <- max(rst.z);
###   
###        z.inx.list <- list();
###        for(i in 1:max.z) {
###           z.inx.list[[i]] <- which(rst.z==i);
###         }
###       for(h.iter in 1:max.z) {
###         for(j in 1:n.v) {
###            # in our case each covariate having 2 level, hence only need to save the 1st probability
###            n.vhjeq1 <- sum(simu.v[ z.inx.list[[h.iter]], j ]);
###            rst.thai[h.iter,j] <- rdirichlet(1,c(1+length( z.inx.list[[h.iter]] )-n.vhjeq1, 1+n.vhjeq1))[1];
###         }
###       }
###      if(max.h>max.z) {
###        s.length <- (max.h-max.z)*n.v;
###        tmp.thai <- matrix(rdirichlet(s.length,c(1,1))[,1], (max.h-max.z), n.v);
###        rst.thai[(max.z+1):(max.h),1:n.v] <- tmp.thai;
###      }
###     
###      all.v <- matrix(0, niter, max.h);
###      all.z <- matrix(0, niter, n.obs);
###      all.u <- matrix(0, niter, n.obs);
###   
###      all.v[1,] <- rst.v;
###      all.z[1,] <- rst.z;
###      all.u[1,] <- rst.u;
###        
###      beyond.length <- 0;
###      all.f <- 0;
###      hit.max <- 0;
###          #### get the corresponding v.index number for each subjects
###        sub.indx <- apply(simu.v,1, function(x) vtosimvalue(x,V.CATS))+1;
###   
###     ############################ Finish Initialize ###################################################################
###      
###     ##### Begin Iteration
###     for(s in 2:niter) {
###      max.z <- max(rst.z);
###      ### step 1 update u
###       tmp.omega <- cal.omega(max.h, rst.v);
###       rst.u <- sapply(1:n.obs, function(x) runif(1, 0, tmp.omega[rst.z[x]]));
###       #cat("\n finishi (1)");
###   
###      ### step 2 update thai
###       # make a list such that the h^th element contains the index of subjects with z==h
###        z.inx.list <- list();
###        for(i in 1:max.z) {
###           z.inx.list[[i]] <- which(rst.z==i);
###         }
###      
###      ### step 3 update v  #################################################################  
###      if(max.z==1) {
###        lowb <- max(rst.u);
###        upb <- 1;
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[1] <- rst$samp.v;
###        failed <- rst$failed;
###        #cat("for h=1, bound=(", lowb, ",", upb, ")"); 
###        rst.v[2:(max.h)] <- rbeta((max.h-1), 1, rst.alpha);
###        ############################################################ 
###       }else if(max.z==2) {
###        # for v1
###          # if there is no subject with Z==1
###          if(length(z.inx.list[[1]])==0) {
###            lowb <- 0;
###          }else{
###           lowb <- max(rst.u[z.inx.list[[1]]]);
###          }
###          upb <- 1-max(rst.u[z.inx.list[[2]]])/rst.v[2];
###          cat("for h=1, bound=(", lowb, "," ,upb, ")"); 
###          rst <- samp.ind.v(lowb, upb, rst.alpha);
###          rst.v[1] <- rst$samp.v;
###          failed <- rst$failed;
###          
###        # for v2 
###        lowb <- min(1,max(rst.u[ z.inx.list[[2]] ])/(1-rst.v[1]));
###        upb <- 1;
###        #cat("for h=2, bound=(", lowb, "," ,upb, ")");     
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[2] <- rst$samp.v;
###        failed <- failed+rst$failed;
###        if(failed>=1) {   failed <- 1; }
###        ### for v3 and beyond 
###        rst.v[3:(max.h)] <- rbeta((max.h-2), 1, rst.alpha);
###       ############################################################ 
###       }else{   ### for (max.z>=3)
###       ###(1) for updating rst.v[1]
###          ## (1.1) low bound for V_h always come from groups with Z==h
###        if(length( z.inx.list[[1]] )==0) {
###          lowb <- 0; }else{
###          lowb <- max(rst.u[ z.inx.list[[1]] ]);
###        }  
###          ## (1.2) ## Find upbound for V_(h=1), coming from groups with Z= (h+1), (h+2), ... (max.h)
###         # for group of Z equal to 2
###        if(length( z.inx.list[[2]] )==0) {
###          max1 <- 0;
###        }else{
###          max1 <- max(rst.u[ z.inx.list[[2]] ])/rst.v[2];
###        }
###          upb1 <- 1-max1;       
###        ## (1.3)
###        #  for group of Z great than 3 
###        gt.h <- c(3:max.z);
###        length.gt.h <- max.z-3+1;
###        max.num <- numeric(length.gt.h);
###        max.num <- sapply(gt.h, function(x)  max(0, rst.u[ z.inx.list[[x]] ]));  # update on 07/29
###        denom.gt.h <- numeric(length.gt.h);
###        denom.gt.h[1] <- rst.v[3]*(1-rst.v[2]);  # from group z=3
###        # update on 7/29
###        if(length.gt.h>=2) { # from group z=4, 5, ...
###          denom.gt.h[2:length.gt.h] <- sapply(4:max.z, function(x) rst.v[x]*prod(sapply(2:(x-1), function(y) 1-rst.v[y])));
###        } 
###        max2 <- min(1,max(0,max.num/denom.gt.h));
###        upb2 <- 1-max2;
###        #cat("max2==", max2, "\n");
###        #upb <- 1-max(max1,max2);
###        upb <- min(upb1, upb2);
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[1] <- rst$samp.v;
###        failed <- rst$failed;
###        #cat("for h=1, bound=(", lowb, ",", upb, ")-max1=", max1, "-ma2=", max2, ", with rst.v[1]==", rst.v[1], "\n");     
###      ###  End for updating rst.v[1] with max.z>=3
###        
###        
###      ### (2) for updating rst.v[h], for h=2, ..., (max.z-1) with max.z>=3
###        for(h.iter in 2:(max.z-1)) {
###          # (2.0) for lower bound
###          if( length(z.inx.list[[h.iter]] )==0) {
###             lowb <- 0;
###           }else{
###            lowb <- max(rst.u[ z.inx.list[[h.iter]] ])/prod(sapply(1:(h.iter-1), function(x) 1-rst.v[x]));
###           }
###          # (2.1) find upper bound max1 for rst.v[h] from group with Z=h+1
###          #indxeq.hiterplus1 <- which(rst.z==(h.iter+1));
###          if(length( z.inx.list[[h.iter+1]]) ==0) {
###            max1 <- 0;
###          }else{
###             #cat("numer=", max(rst.u[ z.inx.list[[h.iter+1]] ]), "***denom==", (rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x])))), "\n");
###             max1 <- max(rst.u[ z.inx.list[[h.iter+1]] ])/(rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x]))));            }
###          upb1 <- 1-max1;
###          
###          #(2.2) find upper bound max2 for rst.v[h] from group with Z=h+2, h+3, ...
###          if(h.iter+1>=max.z) { # on the equality with h.iter=max.z-1, the upbound taken care by upb1
###             max2 <- 0;
###          }else {   
###             ### Begin Alternative approach
###             whichgt.h <- c((h.iter+2):(max.z));  ## {h*:h*>h.iter}
###             lengthgt.h <- length(whichgt.h);
###             
###             num.gt.h <- numeric(lengthgt.h);
###             num.gt.h <- sapply(whichgt.h, function(x) max(0, rst.u[ z.inx.list[[x]]  ]));
###             indx.num.gt0 <- which(num.gt.h>0);
###             num.gt0 <- num.gt.h[indx.num.gt0];
###             denom.gt0 <- numeric(length(indx.num.gt0));
###             
###             denom.gt0 <- sapply(whichgt.h[indx.num.gt0], function(x) rst.v[x]*prod(sapply(1:(h.iter-1), function(y) 1-rst.v[y]))*prod(sapply((h.iter+1):(x-1), function(z)  1-rst.v[z])));
###             
###             #at("num.gt.h===", num.gt0, "\n");
###             #cat("denom.gt.h===", denom.gt0, "\n");q
###             max2 <- min(1,max(0, num.gt0/denom.gt0));
###           }
###             upb2 <- 1-max2;
###             
###             upb <- min(upb1, upb2);
###             #cat("for h=", h.iter, ", bound=(", lowb, ",", upb, ")--max1=", max1, "--max2=", max2, ",");     
###             rst<- samp.ind.v(lowb, upb, rst.alpha);
###             rst.v[h.iter] <- rst$samp.v;
###             failed <- failed+rst$failed;
###        }
###      
###        #(3) #### for h.iter==max.z ( hence definitely some subjects fall in this group)
###        lowb <- max(rst.u[ z.inx.list[[max.z]] ])/prod(sapply(1:(max.z-1), function(x) 1-rst.v[x]));
###        upb <- 1;
###        #cat("for h=", max.z, ", bound=(", lowb, ",", upb, ")");    
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[max.z] <- rst$samp.v;
###        failed <- failed+rst$failed;
###        #(4) update rst.v[h] for h=max.z+1, ..., max.h
###        if((max.z+1)<=max.h){
###        rst.v[(max.z+1):(max.h)] <- rbeta((max.h-max.z),1,rst.alpha);
###        }
###        if(failed>=1) { failed <- 1; }
###      } ## finish for max.z>=3
###   
###      ### keept last v=1 for sum length=1
###      rst.v[max.h] <- 1;
###      tmp.omga <- cal.omega(max.h, rst.v);
###      rst.v[which(rst.v[1:(max.h-1)]==1)] <- 1-(10)^(-10); ### to exclude the case that middle part rst.v==1, hence (1-rst.v)=0 in the denominator on Sep17
###      #cat(",  (3)");
###      # after updating rst.v, need to update tmp.omega   
###      ###########  Finish step 3  for updating V   ############### ########### ########### ########### ########### ########### ########### ###########  
###   
###     ### Update on 06/08-- to save time, first calculate the p(v|h) and use it in updating Z
###            ind.h.prob <- matrix(0, max.h,nrow(V.INX.ALL));
###            for(tmp.h in 1:max.h) {
###               ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
###            }
###          
###       ### step 4 
###        # update omega first, need to calculate for max.h
###        min.u <- min(rst.u);
###        #cat("min.u==", min.u, ",   ");
###        addlength <- sapply(1:(max.h), function(x) sum(tmp.omega[1:x]));
###        k.tilt <- min(max.h, min(which(addlength>=(1-min.u))));
###        cat(", with k.tilt==", k.tilt, "\n");
###   
###          for(i in 1:n.obs) {
###            tmp.A <- which(tmp.omega[1:k.tilt]>rst.u[i]);
###            #cat("for sub==", i, ", sample from(", tmp.A, ")\n");
###            size.A <- length(tmp.A);
###            # calculate the multinomial distribution probability
###            if(size.A==0) {
###             rst.z[i] <- rst.z[i];
###            }else if(size.A==1){
###             rst.z[i] <- tmp.A[1];
###            }else{  
###             prob.A <- numeric(size.A);
###             prob.A <- ind.h.prob[tmp.A, sub.indx[i]];
###             prob.A <- prob.A/sum(prob.A);
###             rst.z[i] <- tmp.A[which(rmultinom(1,1,prob.A)==1)];
###            }
###         }
###       #cat(",  (4)");
###          #max.z <- max(rst.z);
###          #if(max.z==max.h) n.hitmax <- n.hitmax+1;
###          #cat("updated with z==", rst.z, " with max.z==", max(rst.z), "\n");
###        
###      #### step 4.5 -- only update in 0928 function()
###      unique.z <- unique(rst.z);
###      rank.z <- rank(unique.z);
###      update.rstz <- sapply(1:n.obs, function(x) rank.z[which(unique.z==rst.z[x])]);
###      update.tmp.omega <- c(tmp.omega[sort(unique.z)], tmp.omega[-unique.z]);
###      rst.z <- update.rstz;
###      tmp.omega <- update.tmp.omega;
###         ### based on each individual length, update rst.v
###         update.rstv <- numeric(max.h);
###         update.rstv[1] <- tmp.omega[1];
###         for(j in 2:max.h) {
###            update.rstv[j] <- tmp.omega[j]/prod(sapply(1:(j-1), function(y) (1-update.rstv[y])));
###         }   
###   
###        ### step 2 update thai
###        # make a list such that the h^th element contains the index of subjects with z==h
###        z.inx.list <- list();
###        for(i in 1:max.z) {
###           z.inx.list[[i]] <- which(rst.z==i);
###         }
###      
###       for(h.iter in 1:max.z) {
###         for(j in 1:n.v) {
###            # in our case each covariate having 2 level, hence only need to save the 1st probability
###            n.vhjeq1 <- sum(simu.v[ z.inx.list[[h.iter]], j ]);
###            rst.thai[h.iter,j] <- rdirichlet(1,c(1+length( z.inx.list[[h.iter]] )-n.vhjeq1, 1+n.vhjeq1))[1];
###         }
###       }
###       if(max.h>max.z) {
###        s.length <- (max.h-max.z)*n.v;
###        tmp.thai <- matrix(rdirichlet(s.length,c(1,1))[,1], (max.h-max.z), n.v);
###        rst.thai[(max.z+1):(max.h),1:n.v] <- tmp.thai;
###       }
###       #cat(",  (2)");
###   
###       for(tmp.h in 1:max.h) {
###          ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
###       } 
###       rst.theta[s-1,] <- t(ind.h.prob)%*%tmp.omega;
###         
###       if(max.z==1) {
###         rst.theta.rw[s-1,] <- ind.h.prob[1,];
###       }else{
###         rst.theta.rw[s-1,] <- t(ind.h.prob[1:max.z,])%*%(tmp.omega[1:max.z])/(sum(tmp.omega[1:max.z]));
###       }
###   
###           ### step 5
###         # update omega first
###       tmp.shape <-  a.alpha+max.z;
###        if(max.z==max.h) {
###          inx.neq1 <- which(rst.v[1:(max.h)]<1);
###          tmp.rate <- Inf;
###          rst.alpha <- (10)^(-10);
###        }else{
###          inx.neq1 <- which(rst.v[1:(max.z)]<1);
###          tmp.rate <- b.alpha-log(1-sum(tmp.omega[1:max.z]));
###           rst.alpha <- rgamma(1, shape=tmp.shape, rate=tmp.rate);
###        }
###     
###      #cat("alpha (a==", tmp.shape, ", b==", tmp.rate, ")\n");
###       #cat(",  (5)");
###   
###      all.v[s,] <- rst.v;
###      all.z[s,] <- rst.z;
###      all.u[s,] <- rst.u;
###        ### calculate the pchi after each iteration
###       if(s>2) {
###         if(s<=(nburn+2)) {
###           m <- apply(rst.theta[1:(s-1),],2,mean);
###           m.new <- apply(rst.theta.rw[1:(s-1),], 2, mean);
###         }else{
###           m <- apply(rst.theta[nburn:(s-1),],2,mean);
###           m.new <- apply(rst.theta.rw[nburn:(s-1),],2,mean);
###         }  
###       mean.pchi.theta <- probv.pchi(m, true.theta);
###       crnt.pchi.theta <- probv.pchi(rst.theta[s-1,], true.theta);
###         
###       mean.new.pchi.theta <- probv.pchi(m.new, true.theta);
###       crnt.new.pchi.theta <- probv.pchi(rst.theta.rw[s-1,], true.theta);
###       cat(" at End iteration ", s, ", pchi===( ", mean.pchi.theta, ",", crnt.pchi.theta, "),  new.pchi==(", mean.new.pchi.theta, ", ",  crnt.new.pchi.theta, "),  ****max.z===", max(rst.z), ", alpha=", rst.alpha, ", for (K,i)==(", k, ",", data.i, ") with n.hitmax==", n.hitmax, "\n");
###      }
###      #cat(",  (6)");
###      #cat("updated with u==", rst.u, "\n");
###      #cat("updated with v==", rst.v[1:max.z], "\n");
###      beyond.length <- c(beyond.length, 1-sum(tmp.omega[1:(max.z)]));
###      all.f <- c(all.f, failed);
###      hit.max <- c(hit.max, (max.z==max.h));
###      if(max.z==max.h) { n.hitmax <- n.hitmax+1; }
###      
###     } ### end of iteration
###      list(theta=rst.theta, n.hitmax=n.hitmax, beyond.length=beyond.length[-1], all.f=all.f[-1], hit.max=hit.max[-1], theta.rw=rst.theta.rw, all.v=all.v, all.u=all.u, all.z=all.z);
###   }
###   



#posterior sampling of v
posts.samp.v <- function(maxz, maxh, rst.v, rst.u, rst.z, rst.alpha) {

  unique.z <- unique(rst.z);

  cat("unique z==", unique.z, "\n");
  if(maxz==1) {
    lowb <- max(rst.u);
    upb <- 1;
    cat("()=(", lowb,",", upb, ")"); 
    rst.v[1] <- samp.ind.v(lowb,upb,rst.alpha);
    cat(rst.v[1], "\n");
    #sample the rest    
    rst.v[2:maxh] <- rbeta((maxh-maxz),1,rst.alpha);
  }else if(maxz==2) {
    max.u.ind.h <- sapply(1:maxz, function(x) max(0,rst.u[which(rst.z==x)]));
    #for v1
    if(is.element(1,unique.z)){ # if there are subjects with Z=1
      lowb <- max.u.ind.h[1];
      upb <- 1-max.u.ind.h[2]/rst.v[2];
    }else{
      lowb <- 0;
      upb <- 1;
      #rst.v[1] <- rst.v[1]; # to strick constraints (if random in (0,1), the lowbound for v[2] might be over 1, hence larger than upbound
    }
    cat("for 1 ()=(", lowb,",", upb, ")"); 
    rst.v[1] <- samp.ind.v(lowb,upb,rst.alpha);
    cat(rst.v[1], "\n");
    
    #for v2
    if(is.element(2,unique.z)){
      lowb <- min(1,max.u.ind.h[2]/(1-rst.v[1]));
      upb <- 1;
    }else{
      lowb <- 0;
      upb <- 1;
    }
    cat("for 2 ()=(", lowb,",", upb, ")");
    rst.v[2] <- samp.ind.v(lowb,upb,rst.alpha);
    cat(rst.v[2], "\n");
    
    #sample the rest
    rst.v[(maxz+1):maxh] <- rbeta((maxh-maxz),1,rst.alpha);
  }else {
    max.u.ind.h <- sapply(1:maxz, function(x) max(0,rst.u[which(rst.z==x)]));
     #for v1
    lowb <- max.u.ind.h[1];
      exist.h <- unique.z[which(unique.z>1)];
    
      if(is.element(2,exist.h)) {
        # for h=2
        upb1 <- 1-max.u.ind.h[2]/rst.v[2];
        exist.h.no2 <- sort(setdiff(exist.h, 2));
        upb2 <- 1-sapply(exist.h.no2, function(x) max.u.ind.h[x]/(rst.v[x]*prod(1-rst.v[2:(x-1)])));
        upb <- min(upb1, upb2);
        cat("exist.h.no2==", exist.h.no2, "***max.u=", max.u.ind.h[exist.h.no2], ", *** upb1==", upb1, "**** upb2=", upb2, "\n");
       }else{
        upb2 <- 1-sapply(exist.h, function(x) max.u.ind.h[x]/(rst.v[x]*prod(1-rst.v[2:(x-1)])));
        upb <- min(upb2);
        cat("max.u=", max.u.ind.h[exist.h], ", *** upb2===", upb2, "\n");
       }
    cat("for 1 ()=(", lowb,",", upb, ")");  
    rst.v[1] <- samp.ind.v(lowb,upb,rst.alpha);
    cat(rst.v[1], "\n");
   

    #for v2
    lowb <- min(1,max.u.ind.h[2]/(1-rst.v[1]));
      exist.h <- unique.z[which(unique.z>2)];
      if(is.element(3,exist.h)){
       # for h=3
       upb1 <- 1-max.u.ind.h[3]/(rst.v[3]*(1-rst.v[1]));
       exist.h.no3 <- sort(setdiff(exist.h, 3));
       if(length(exist.h.no3)>0){
        # for h=4, ... , max.z
        upb2 <- 1-sapply(exist.h.no3, function(x) max.u.ind.h[x]/(rst.v[x]*(1-rst.v[1])*prod(1-rst.v[3:(x-1)])));
        upb <- min(upb1, upb2);
        cat("max.u==", max.u.ind.h[exist.h.no3], "**** upb1==", upb1, "*** upb2==", upb2, "\n");
       }else{
         upb <- upb1;
       } 
     }else{
       upb2 <- 1-sapply(exist.h, function(x) max.u.ind.h[x]/(rst.v[x]*(1-rst.v[1])*prod(1-rst.v[3:(x-1)])));
       upb <- min(upb2);
     }
     cat("for 2 ()=(", lowb,",", upb, ")\n"); 
     rst.v[2] <- samp.ind.v(lowb,upb,rst.alpha);
     cat(rst.v[2], "\n");

     #sample the rest v3 to v(max.z)
     for(hiter in 3:maxz) {
       if(is.element(hiter, unique.z)){
         lowb <- max.u.ind.h[hiter]/(prod(1-rst.v[1:(hiter-1)]));
         if(hiter==maxz) {
           upb <- 1;
         }else if(hiter==(maxz-1)) {
           upb <- 1-max.u.ind.h[maxz]/(prod(1-rst.v[1:(hiter-1)])*rst.v[maxz]);
         }else{
           exist.h <- unique.z[which(unique.z>hiter)];
           if(length(exist.h)==1) {
             if(is.element(hiter+1,exist.h)) {
               upb <- 1-sapply(exist.h, function(x) max.u.ind.h[x]/(prod(1-rst.v[1:(hiter-1)])*rst.v[x]));
             }else{
               upb <- 1-max.u.ind.h[exist.h]/(prod(1-rst.v[1:(hiter-1)])*prod(1-rst.v[(hiter+1):(exist.h-1)])*rst.v[exist.h]);
             }
           }else{
            upb1 <- 1;
            if(is.element(hiter+1, exist.h)) {
               upb1 <- 1-sapply(exist.h, function(x) max.u.ind.h[x]/(prod(1-rst.v[1:(hiter-1)])*rst.v[x])); }
            exist.h.noplus1 <- setdiff(exist.h, (hiter+1));
            upb2 <- 1-sapply(exist.h.noplus1, function(x) max.u.ind.h[x]/(prod(1-rst.v[1:(hiter-1)])*prod(1-rst.v[(hiter+1):(x-1)])*rst.v[x]));
            upb <- min(upb1,upb2);
           } 
         }
        cat("for ", hiter, " ()=(", lowb, ",", upb, ")");
        rst.v[hiter] <- samp.ind.v(lowb,upb,rst.alpha);
        cat(rst.v[hiter], "\n");
       }else{
         rst.v[hiter] <- 0;
         cat("for ", hiter, "===", rst.v[hiter],"\n");
       }   
        
     }

     rst.v[(maxz+1):maxh] <- rbeta((maxh-maxz),1,rst.alpha);
    }
    rst.v;
  }  


##############################################################################

####     ##multinomial log likelihood for auxillary covariates
####     ##last curpar is eta
####     ##curpar is a vector, need to make it a list first and then
####     ##call get.prob.v
####     ### neta == Sample size
####     v.log.fx0629 <- function(vcount, curtheta, curpar, curomega, nvcats) {
####         ##convert to lst of prob v
####         pv.eta     <- v.par.vec2lst(curpar, nvcats);
####         lst.prob.v <- pv.eta[[1]];
####         eta        <- pv.eta[[2]];
####         prob.v     <- get.prob.v(lst.prob.v, 1);
####         nsize <- sum(vcounts);
####     
####         logomega <- log(curomega);
####         logomega[which(is.infinite(logomega))] <- -1000;
####         logtheta <- log(curtheta);
####         logtheta[which(is.infinite(logtheta))] <- -1000;
####         logvtomega <- log(prob.v*curomega);
####         logvtomega[which(is.infinite(logvtomega))] <- 1000;
####     
####         logOmega <- lgamma(length(prob.v)/eta)-length(prob.v)*lgamma(1/eta)+sum(logomega)*(1/eta-1);
####         logD   <- lgamma(sum(prob.v*curomega)) - sum(logvtomega)+sum(sapply(1:length(prob.v),function(x) logtheta[x]*(vcount[x]+prob.v[x]*curomega[x]-1)));
####         
####         if (!is.null(neta)) {
####             ##uniform shrinkage prior
####             logeta <- -2 * log(nsize * eta + 1);
####         } else {
####             logeta <- 0;
####         }
####     
####         ##bound away from -Inf
####         #logV[which(is.infinite(logV))] <- -1000;
####         rst  <- logD + logOmega + logeta;
####         #cat("logD=", logD, ", logOmega=", logOmega, ", logeta=", logeta, "\n");
####           
####         stopifnot(!is.infinite(rst))
####         rst
####     }
####      
####     slice.onestep0629 <- function(vcount, curtheta, curpar, curomega, parinx, fix.inx, upbound, log.fx, ncats, max.iter=80000) {
####         ##browser();
####         bound <- c(0, upbound);
####         if(fix.inx==999) { ### then update omega
####            x0 <- curomega[parinx];
####         }else{   ### then update pi and eta
####            x0    <- curpar[parinx];
####         }
####     
####         ##get log f(x0)
####         log.fx0 <- log.fx(vcount, curtheta, curpar, curomega, ncats);
####     
####         ##get auxiliary y~unif(0,fx0), that is
####         log.y   <- log.fx0 - rexp(1);
####     
####         ##sampling from slice
####         flag   <- FALSE;
####         n.iter <- 0;
####         while (!flag & (n.iter<max.iter)) {
####             x1             <- runif(1, bound[1], bound[2]);
####             if(fix.inx==999) {
####               curomega[parinx] <- x1;
####             }else{  
####              curpar[parinx] <- x1;
####            }
####             
####             if (!is.null(fix.inx)) {
####               if(fix.inx!=999) {
####                 curpar[fix.inx]<- upbound - x1;
####               } 
####             }
####             log.fx1 <- log.fx(vcount, curtheta, curpar, curomega, ncats);
####     
####             if (log.y < log.fx1) {
####                 flag <- TRUE;
####             } else if (x1 < x0) {
####                 bound[1] <- x1;
####             } else if (x1 > x0) {
####                 bound[2] <- x1;
####             }
####             #cat("bound===", bound, "\n");
####             n.iter <- n.iter + 1;
####         }
####     
####         if (n.iter >= max.iter) {
####             print("slice failed");
####         }
####         x1;
####     }
####     
####     
####     ##vall: the observed/simulated V data
####     ##niter: number of iterations
####     ##log.fx: loglikelihood function
####     ##inits: initial values
####     ## v.inx==V.CATS:categories for all V auxillary covariates
####     ### prior.neta: number of times to do the slice sampling
####     
####     slice.v.par.new0629 <- function(vall, niter, log.fx, v.inx,
####                             prior.neta=NULL, inits=NULL, true.theta) {
####     
####         bound.eta <- 10;
####     
####         ##v categories
####         ncats  <- mapply(length, v.inx);
####         ntheta <- prod(ncats);
####         npv    <- length(v.inx);
####     
####         ##init values
####         if (is.null(inits)) {
####             for (i in 1:npv) {
####                 ##pi_is init at 1/L_i
####                 inits <- c(inits, rep(1/ncats[i], ncats[i]));
####             }
####             ##eta init at 1/N
####             inits <- c(inits, 1/nrow(vall));
####         }
####     
####         ##get count data
####         vcount <- get.v.inx(vall, v.inx, counts=TRUE);
####     
####         ##initialize results
####         rst.theta <- array(NA, dim=c(niter+1, ntheta));
####         rst.omega <- array(NA, dim=c(niter+1, ntheta));
####         rst       <- array(NA, dim=c(niter+1, sum(ncats)+1));
####         rst[1, ]  <- inits;
####         rst.omega[1,] <- rep(1/ntheta, ntheta);
####         rst.theta[1,] <- rep(1/ntheta, ntheta);
####         
####         for (s in 2:(niter+1)) {
####             ##browser();
####             rst[s,] <- rst[s-1,];
####             rst.omega[s,] <- rst.omega[s-1,];
####             ##gibbs theta L
####             #pv.eta      <- v.par.vec2lst(rst[s,], ncats);
####             #lst.prob.v  <- pv.eta[[1]];
####             #eta         <- pv.eta[[2]];
####             #prob.v      <- get.prob.v(lst.prob.v, eta);
####             #cur.dir.par <- vcount + prob.v;
####             #next.theta  <- c(rdirichlet(1, cur.dir.par));
####             #rst.theta[, s-1] <- next.theta;
####     
####             ##Gibbs sampling for pi's
####             cur.start <- 1;
####             for (i in 1:npv) {
####                 cur.end     <- (cur.start-1) + ncats[i];
####                 next.start  <- cur.end + 1;
####     
####                 ##random choose one pi to be 1-sum of the rest pis
####                 cur.fix.inx <- sample(cur.start:cur.end, 1);
####                 for (j in cur.start:cur.end) {
####                     ##fixed pi will not be sampled
####                     if (cur.fix.inx == j) next;
####     
####                     ##find bound
####                     up.bound <- rst[s,j] + rst[s,cur.fix.inx];
####     
####                     ##slice sampling
####                     slc.pi <- slice.onestep0629(vcount, next.theta, rst[s,], rst.omega[s,], j, cur.fix.inx, up.bound, v.log.fx0629, ncats, prior.neta);
####     
####                     ##update
####                     rst[s, j] <- slc.pi;
####                     rst[s, cur.fix.inx] <- up.bound - slc.pi;
####                 }
####     
####                 cur.start <- next.start;
####             }
####     
####             ###Gibbs sampling for omega
####             for(j in 1:(ntheta-1)) {
####                if(j==1) {
####                  up.bound <- 1;
####                }else{
####                  up.bound <- 1-sum(rst.omega[s,1:(j-1)]);
####                }
####                if(up.bound==0) {
####                  rst.omega[s,j:(ntheta-1)] <- 0;
####                  break;
####                }  
####                slc.omega <- slice.onestep0629(vcount, next.theta, rst[s,], rst.omega[s,], j, 999, up.bound, v.log.fx0629, ncats, prior.neta);
####                # update
####                rst.omega[s,j] <- slc.omega;
####                #cat("for j==", j, ", update with ", slc.omega, "\n");
####             }
####             rst.omega[s, ntheta] <- 1-sum(rst.omega[s,1:(ntheta-1)]);
####     
####             ##Gibbs/Slice sampling for eta
####             slc.eta <- slice.onestep0629(vcount, next.theta, rst[s,], rst.omega[s,], ncol(rst), 1, bound.eta, v.log.fx0629, ncats, prior.neta);
####             rst[s, ncol(rst)] <- slc.eta;
####     
####             ### update theta
####             pv.eta <- v.par.vec2lst(rst[s,], ncats);
####             lst.prob.pi <- pv.eta[[1]];
####             probpi <- get.prob.v(lst.prob.pi,1);
####             omegapi <- rst.omega[s,]*probpi;
####             rst.theta[s,] <- rdirichlet(1, vcounts+omegapi);
####             cat("at iteration S=", s, ", the c.pchi==", probv.pchi(rst.theta[s,], true.theta), ", the a.pchi==", probv.pchi(apply(rst.theta[1:s,],2,mean), true.theta), "\n");
####         }
####     
####         list(post.pieta=rst, post.omega=rst.omega, post.theta=rst.theta);
####     }
####     
####     
####     prv.xuan0629 <- function(niter, nburn, v, v.log.fx, V.CATS, ttln, true.theta)
####     {
####     
####       post.v.par <- slice.v.par.new0629(v, niter, v.log.fx, V.CATS, ttln);
####       prob.v.dist <- post.v.par$post.prob.v;
####     
####       mean.theta.dan <- apply(prob.v.dist[,(nburn+1):niter], 1, mean);
####       mse.dan <- t(true.theta-mean.theta.dan)%*%(true.theta-mean.theta.dan);
####       pchi.probv.rst <- probv.pchi(mean.theta.dan, true.theta);
####     
####       list(all.theta=prob.v.dist, eta=post.v.par$post.pieta[,17], mean.theta=mean.theta.dan, mse.theta=mse.dan, pchi.theta=pchi.probv.rst);
####     }
####     



##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################


####### Standard Uncongenial Analysis

### obtain the estimates of Delta, Alpha, Gamma in the uncongenial imputation model
uncong.impt <- function(SIMU.V.RST, SIMU.R.RST.SUM, MISS.Y.RST) {

  YV.CATS <- list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1));
  ### Y1VR.INX.ALL gives all the possible combination of Y1VR
  Y1V.INX.ALL <- get.v.all(YV.CATS);

  ### prepare for calculating frequency for each combination of (Y1, V, R)
  SIMU.Y1V <- cbind(MISS.Y.RST[,1], SIMU.V.RST);
  ### calculate in SIMU.Y1VR, the frequency for each combination of Y1VR
  emp.Y1V.n <- get.v.inx(SIMU.Y1V, YV.CATS, counts=TRUE);

  ### save the data preparing for logistic regression of P(Y1|V,R),
  ### need for each category the number of success and number of total
  cmpt1 <- as.data.frame(cbind(Y1V.INX.ALL, emp.Y1V.n));
  colnames(cmpt1) <- c("Y1", "V1", "V2","V3","V4","V5","V6","V7","V8", "Freq");

  ### calculate for each category combination, the total number of frequency
  ### since put the outcome Y1 in the first, finding the corresponding sucuess and failure is very easy
  totalforV <- rep(NA,length(emp.Y1V.n)/length(YV.CATS[[1]]));
  for (i in 1:length(totalforV))
  {  totalforV[i] <- cmpt1[i,10]+cmpt1[i+256,10]; }

  ### the cmpt1.ana is the data used for logistic regression
  cmpt1.ana <-as.data.frame(cbind(1,0, cmpt1[1:256,-1], totalforV));
  colnames(cmpt1.ana) <- c("T", "prevy", "V1","V2","V3","V4","V5","V6","V7","V8","FreqY=0","Total");

  ###(1.2) based on the logistic regression for Y2
  ### Find the subjects with R greater or equal than 2, hence have Y2 observed
  indx.rgeq2 <- which(SIMU.R.RST.SUM>=2);

  SIMU.Y2Y1V <- cbind(MISS.Y.RST[indx.rgeq2,2], MISS.Y.RST[indx.rgeq2,1], SIMU.V.RST[indx.rgeq2,]);

  ### for Y2 imputation model is free of missing pattern
  Y2Y1V.CATS <-list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1));
  Y2Y1V.INX.ALL <- get.v.all(Y2Y1V.CATS);
  ### calculate the number of subjects in every one of all possible combination category of (Y2,Y1,V)
  emp.Y2Y1V.n <- get.v.inx(SIMU.Y2Y1V, Y2Y1V.CATS, counts=TRUE);
  cmpt2 <- as.data.frame(cbind(Y2Y1V.INX.ALL, emp.Y2Y1V.n));
  colnames(cmpt2) <- c("Y2", "Y1", "V1","V2","V3","V4","V5","V6","V7","V8","Freq");

  ### get the independent category combination based on (Y1, V)
  totalforV <- rep(NA, nrow(Y2Y1V.INX.ALL)/length(Y2Y1V.CATS[[1]]));
  templength <- length(totalforV);
  for (i in 1:templength)
   {  totalforV[i] <- cmpt2[i,11]+cmpt2[i+templength,11]; }

  cmpt2.ana <-as.data.frame(cbind(2, cmpt2[1:templength,-1], totalforV));
  colnames(cmpt2.ana) <- c("T", "prevy", "V1","V2","V3","V4","V5","V6","V7","V8","FreqY=0","Total");

  ### find the subjects with R=3
  indx.r3 <- which(SIMU.R.RST.SUM==3);
  SIMU.Y3Y2V <- cbind(MISS.Y.RST[indx.r3,3], MISS.Y.RST[indx.r3,2], SIMU.V.RST[indx.r3,]);
  
  Y3Y2V.CATS <-list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1));
  emp.Y3Y2V.n <- get.v.inx(SIMU.Y3Y2V, Y3Y2V.CATS, counts=TRUE);    
  Y3Y2V.INX.ALL <- get.v.all(Y3Y2V.CATS);
  cmpt3 <- as.data.frame(cbind(Y3Y2V.INX.ALL, emp.Y3Y2V.n));
  colnames(cmpt3) <- c("Y3","Y2","V1","V2","V3","V4","V5","V6","V7","V8","Freq");

  totalforVR <- rep(NA,nrow(Y3Y2V.INX.ALL)/length(Y3Y2V.CATS[[1]]));
  templength <- length(totalforVR);
  for (i in 1:templength)
   {  totalforVR[i] <- cmpt3[i,11]+cmpt3[i+templength,11]; }
  sum(totalforVR);

  cmpt3.ana <-as.data.frame(cbind(3, cmpt3[1:templength,-1], totalforVR));
  colnames(cmpt3.ana) <- c("T", "prevy", "V1","V2","V3","V4","V5","V6","V7","V8","FreqY=0","Total");

  ### based on the available observed data, estimate parameters in imputation model
  un.ana <- rbind(cmpt1.ana, cmpt2.ana, cmpt3.ana);
  colnames(un.ana) <- c("T", "prevy", "V1","V2","V3","V4","V5","V6","V7","V8","FreqY=0","Total");

  fit1.uncong <- glm(cbind(un.ana[,12]-un.ana[,11],un.ana[,11]) ~ as.factor(T)+prevy+as.factor(V1)+as.factor(V2)+as.factor(V3)+as.factor(V4)+as.factor(V5)+as.factor(V6)+as.factor(V7)+as.factor(V8), family=binomial(logit), data=un.ana);

  ##### Estimation of parameters in imputation model based on Uncongenial analysis##############################################
  est.Delta <- fit1.uncong$coefficient[1:3];
  est.Gamma <- fit1.uncong$coefficient[4];
  est.Alpha <- fit1.uncong$coefficient[5:12];

  list(est.Delta=est.Delta, est.Gamma=est.Gamma, est.Alpha=est.Alpha);
}




###### For standard Uncongenial Analysis, USE this Function 
uncong.trans.0630<- function(nimpt, theta, alpha, gamma, delta, y, v, r, BETA) {
     
   est.beta.Uncong <- matrix(0,nimpt,3);
   est.beta.Uncong.sd <- matrix(0, nimpt, 3);
   est.p.uncong.optim <- matrix(0, nimpt, 5);
   est.beta.uncong.optim <- matrix(0, nimpt, 3);
   est.beta.uncong.optim.sd <- matrix(0, nimpt, 3);
   prob.var.Uncong.var <- matrix(0, nimpt, 5);

   origin.y <- y;
   ttime <- c(1,1,2,2,3,3);
   ttime <- ttime-1;
   
   for (i in 1:nimpt){
     #impute missing Y in Uncongenial Analysis based on imputation model, Alpha only vector as same throught all T
     y <- SIMU.Y.AMAR.xuan(delta, alpha, gamma, v, r, origin.y);  ### model is problematic...in impt with intercept, in imputation, the intercept is not included...
     cat("+++++++++++++++++++++++++++++++++ dim(y)=", dim(y),"\n");
     
     t10 <- c(0, 0, length(which(y[,1]==1)), length(which(y[,1]==0)));
     t20 <- c(1, 0, length(intersect(which(y[,1]==0), which(y[,2]==1))), length(intersect(which(y[,1]==0), which(y[,2]==0))));
     t21 <- c(1, 1, length(intersect(which(y[,1]==1), which(y[,2]==1))), length(intersect(which(y[,1]==1), which(y[,2]==0))));
     t30 <- c(2, 0, length(intersect(which(y[,2]==0), which(y[,3]==1))), length(intersect(which(y[,2]==0), which(y[,3]==0))));
     t31 <- c(2, 1, length(intersect(which(y[,2]==1), which(y[,3]==1))), length(intersect(which(y[,2]==1), which(y[,3]==0))));
     
     data.inf <- as.data.frame(rbind(t10,t20,t21,t30,t31));
     colnames(data.inf) <- c("T", "prevy", "NY1", "NY0");
     data.inf.rst <- glm(cbind(NY1, NY0) ~ T + prevy, family=binomial(logit), data=data.inf);
     est.beta.Uncong[i,] <- data.inf.rst$coefficient;
     # the standard error for each parameter sqrt(U_j)
     est.beta.Uncong.sd[i,] <- summary(data.inf.rst)$coefficients[,2];
     
     freqy <- freq.y(y);
     rst <- optim(data.inf.rst$coefficient, impt.fullloglikelihood, freq=freqy, hessian=TRUE, control=list(maxit=10000, fnscale=-1));
     est.beta.uncong.optim[i,] <- rst$par;
     varcov.beta <- solve(-1*(rst$hessian));
     est.beta.uncong.optim.sd[i,] <- sqrt(diag(varcov.beta));
     
     est.p.uncong.optim[i,] <- get.prob.5(rst$par);
     grad.p <- jacobian(get.prob.5, rst$par);
     prob.var.Uncong.var[i,] <-  diag(grad.p%*%varcov.beta%*%t(grad.p));
     cat("finish Imputated data============ ", i, "\n");
   }
   #write.table(est.beta.Uncong, file="uncongN100_AMAR.txt");
   ### get the average mean based on all imputated data sets
   est.Beta <- apply(est.beta.Uncong,2,mean);
   ### obtain the between imputation variance
   est.Beta.var.between <- apply(est.beta.Uncong,2,var); ### B
   est.beta.Uncong.var <- est.beta.Uncong.sd*est.beta.Uncong.sd; ### U_j
   est.Beta.var.in <- apply(est.beta.Uncong.var, 2, mean); ### mean(U_j)
   est.Beta.var <- sapply(1:3, function(x) (est.Beta.var.in[x]+(1+1/nimpt)*est.Beta.var.between[x]));

   est.Beta.optim <- apply(est.beta.uncong.optim,2,mean);
   ### obtain the between imputation variance
   est.Beta.var.between <- apply(est.beta.uncong.optim,2,var); ### B
   est.beta.Uncong.var <- est.beta.uncong.optim.sd*est.beta.uncong.optim.sd; ### U_j
   est.Beta.var.in <- apply(est.beta.Uncong.var, 2, mean); ### mean(U_j)
   est.Beta.optim.var <- sapply(1:3, function(x) (est.Beta.var.in[x]+(1+1/nimpt)*est.Beta.var.between[x]));
 
   true.prob <- get.prob.5(BETA);
   est.prob <- apply(est.p.uncong.optim, 2, mean);

   prob.var.between <- apply(est.p.uncong.optim, 2, var);
   prob.var.in <- apply(prob.var.Uncong.var, 2, mean);
   prob.var <- sapply(1:5, function(x) (prob.var.in[x]+(1+1/nimpt)*prob.var.between[x]));
   
   ## obtain the mean of within imputation variance  
   rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   rst;
}   




###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

######## Standard congenial analysis (two method since two approaches for gradient calculation)


### posterior log likelihood as function of (beta, alpha, gamma), return the (-1)*posterior loglikelihood and Delta
x.U <- function(u.q, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
  u.Beta <- u.q[1:3];
  u.Alpha <- u.q[4:11];
  u.Gamma <- u.q[12];

  n.T <- ncol(u.y);
  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
  logPrior <- sum(-u.q^2/(2*u.sigma^2));
  Eta.matrix <- new.obs.likelihood(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
  logLike = sum(Eta.matrix);

  U= -(logLike + logPrior);
  #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
  list(U=U, Delta=u.Delta);
}


##### update on 04/20/2015 by using .C 
newx.U <- function(u.q, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
  u.Beta <- u.q[1:3];
  u.Alpha <- u.q[4:11];
  u.Gamma <- u.q[12];
  rst <- 0;
  n.T <- ncol(u.y);
  nsub <- nrow(u.y);
  ncolv <- ncol(u.v);

  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
  #logPrior <- sum(-u.q^2/(2*u.sigma^2));
   
  vvec <- as.vector(t(u.v));
  vec.alpha <- as.vector(u.fullAlpha);
  yvec <- as.vector(t(u.y));
  yvec[which(is.na(yvec))] <- 0;
  delta <- as.vector(t(u.Delta));

  ### log.Like gives the log(likelihood) of the observed outcome 
  logLike <- .C("ttl_obsloglikelihood", as.integer(n.T), as.integer(nsub), as.double(vvec), as.integer(u.r), as.integer(ncolv), as.integer(yvec), as.double(delta), as.double(vec.alpha), as.double(u.Gamma), as.double(rst))[[10]];
  
  #Eta.matrix <- new.obs.likelihood(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
  #logLike = sum(Eta.matrix);

  #U= -(logLike + logPrior);
  U= -(logLike);
                                        #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
  U;
}



#### update on 04/20/2015
#### realize that the gradient of Gamma is almost 0 for matching imputation model and transition model
#### get rid of Gamma gradient in HMC to check whether better sampling obtained  
newx.U2 <- function(u.q, u.Gamma, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
  u.Beta <- u.q[1:3];
  u.Alpha <- u.q[4:11];
  #u.Gamma <- u.q[12];
  rst <- 0;
  n.T <- ncol(u.y);
  nsub <- nrow(u.y);
  ncolv <- ncol(u.v);

  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
  logPrior <- sum(-u.q^2/(2*u.sigma^2));
   
  vvec <- as.vector(t(u.v));
  vec.alpha <- as.vector(u.fullAlpha);
  yvec <- as.vector(t(u.y));
  yvec[which(is.na(yvec))] <- 0;
  delta <- as.vector(t(u.Delta));

  ### log.Like gives the log(likelihood) of the observed outcome 
  logLike <- .C("ttl_obsloglikelihood", as.integer(n.T), as.integer(nsub), as.double(vvec), as.integer(u.r), as.integer(ncolv), as.integer(yvec), as.double(delta), as.double(vec.alpha), as.double(u.Gamma), as.double(rst))[[10]];
  
  #Eta.matrix <- new.obs.likelihood(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
  #logLike = sum(Eta.matrix);
  
  U= -(logLike + logPrior);
  #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
  list(U=U, Delta=u.Delta);
  #U;
}




###   x.U.parallel <- function(u.q, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
###     u.Beta <- u.q[1:3];
###     u.Alpha <- u.q[4:11];
###     u.Gamma <- u.q[12];
###   
###     n.T <- ncol(u.y);
###     u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
###     u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
###     logPrior <- sum(-u.q^2/(2*u.sigma^2));
###     Eta.matrix <- new.obs.likelihood.parallel(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
###     logLike = sum(Eta.matrix);
###   
###     U= -(logLike + logPrior);
###     #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
###     list(U=U, Delta=u.Delta);
###   }
###   

###   ### Sam as x.U, but it returns only the (-1)*posterior loglikelihood
###   x.U.Uonly.parallel <- function(u.q, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
###     u.Beta <- u.q[1:3];
###     u.Alpha <- u.q[4:11];
###     u.Gamma <- u.q[12];
###   
###     n.T <- ncol(u.y);
###     u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
###     u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
###     logPrior <- sum(-u.q^2/(2*u.sigma^2));
###     Eta.matrix <- new.obs.likelihood.parallel(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
###     logLike = sum(Eta.matrix);
###   
###     U= -(logLike + logPrior);
###     #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
###     U;
###   }
###   


### Sam as x.U, but it returns only the (-1)*posterior loglikelihood
x.U.Uonly <- function(u.q, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
  u.Beta <- u.q[1:3];
  u.Alpha <- u.q[4:11];
  u.Gamma <- u.q[12];

  n.T <- ncol(u.y);
  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
  logPrior <- sum(-u.q^2/(2*u.sigma^2));
  Eta.matrix <- new.obs.likelihood(u.v, u.r, u.y, u.Delta, u.fullAlpha, u.Gamma);
  logLike = sum(Eta.matrix);

  U= -(logLike + logPrior);
  #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
  U;
}






#### (1) The standard congenial analysis using grad() for gradient calculation


### hmc function for posterior sampling of parameter in imputation and transition model
new.hmcStandard <- function(headername, true.u, v,r, y, q, sigma, n.iter, epsilon, L, V.ALL, theta, Marginal.method, diag.sigma.p){

#x is the matrix of the covariates with first column=1 as intercept
n.obs = nrow(v); #get number of observations	
n.param = length(q);
n.T = ncol(y);

#post.samp used to store the parameters value at every iteration
post.samp <- matrix(NA, n.iter, n.param);
post.samp.U <- rep(0, n.iter);

### save the current value of position variable
current.q <- q;
true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
### n.accept used to record how many times the proposed is accepted
n.accept <- 0;
acceptrate <- 0;
### calculate for each iteration and save the value of position variable for each iteration
for (i in 1:n.iter){
#cat("\n Begin Iter=", i, "#################################################################\n");
  if(i==1) {
    #current.beta <- current.q[1:3];    
    #current.alpha <- current.q[4:11];
    #current.gamma <- current.q[12];    
    #c.fullalpha <- cbind(current.alpha, current.alpha, current.alpha);
    #once theta is updated, Delta needs to be recalculated, and U, grad.U
    #current.Delta <- get.Deltat(n.T, c.fullalpha, current.beta, current.gamma, V.ALL, theta, Marginal.method);  
   
    rst.xU <- x.U(current.q, v, r, y, sigma, theta, V.ALL, Marginal.method);
    current.U <- rst.xU$U;
    current.Delta <- rst.xU$Delta;
    current.g <- grad(x.U.Uonly, x=current.q, method="simple",u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
   }

    samp <- new.getSamp(v,r, y, current.q, current.U, current.g, sigma, n.param, epsilon, L, current.Delta, theta, diag.sigma.p, V.ALL, acceptrate, Marginal.method);

  ### update the position value, U, and gradient of U
     current.q <- samp$q;
     current.Delta <- samp$Delta;
     current.U <- samp$U;
     current.g <- samp$g;
     current.update <- samp$Update; ### 1 for accept the prosed; 0 for reject the propose
         
  ### save the result of iteration i
     post.samp[i,] <- current.q;
     post.samp.U[i] <- current.U;

     n.accept <- n.accept+1*current.update;
     acceptrate <- (n.accept)/(i);
   
        ### in case the acceptance rate is small, quit to save time 
        if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
        {
          cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
          list(samp=post.samp, acceptrate=n.accept/n.iter)
          break
        }
         
         if (current.update==0)
         {
          cat("###Reject at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         } else {
          cat("###Accept at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         }
       }
  rst <- cbind(post.samp, post.samp.U);
  rst;
}


### getSamp is the function to complete one update of Hamiltonian Monte Carlo
new.getSamp <- function(v,r,y, gs.q, gs.U, gs.g, sigma, n.param, gs.et, L, gs.Delta, gs.theta, diag.sigma.p, V.ALL, acptrate, Marginal.method){
### first save the current value of position variable q
### generate random momentum variable p with the same length with q with normal distribution (0,1)
### save the current value of gradient of du/dq

current.q = gs.q;
current.g = gs.g;
current.Delta= gs.Delta;
current.U = gs.U;

n.T = ncol(y);

var.p <- matrix(0, n.param, n.param);
diag(var.p) <- diag.sigma.p^2;
mean.p <- rep(0,n.param);
p <- rmvnorm(1, mean.p, var.p);
q = current.q;
### save the newly generated momentum variable as current.p
### then the current value of K(p)=(p^T%*%p)/2 as current.K based on current.p
### then the current value of H=K(p)+U(q), save the value as current.H
current.p = p;
current.K = sum(current.p^2/(2*diag.sigma.p^2));
current.H = current.U + current.K;

gs.et <- runif(1, 0.8*gs.et, 1.2*gs.et);
gamma.lowbound <- (-3)*diag.sigma.p[12];
gamma.upbound <- (-1)*gamma.lowbound;

for(leap in 1:L){
    p = p - gs.et*gs.g/2;
    # q = q + gs.et*sapply(1:n.param, function(x) p[x]/(diag.sigma.p[x])^2);
    q = q + gs.et*p/(diag.sigma.p^2);
    while(q[12]< gamma.lowbound|| q[12]>gamma.upbound){
          if(q[12]< gamma.lowbound) {
              q[12]=2*gamma.lowbound-q[12];
              p[12]=-p[12];
          }
          if(q[12]> gamma.upbound) {
              q[12]=2*gamma.upbound-q[12];
              p[12]=-p[12];
          }
      }
    #### since q is new, need update gs.Delta based on this q first
        #tmp.beta <- q[1:3];
        #tmp.alpha <- q[4:11];
        #tmp.gamma <- q[12];
        #tmp.full.alpha <- cbind(tmp.alpha, tmp.alpha, tmp.alpha);
        #tmp.Delta <- get.Deltat(n.T, tmp.full.alpha, tmp.beta, tmp.gamma, V.ALL, gs.theta, Marginal.method);    
    ### true gradient based on all subjects (1st choice)
    gs.g = grad(x.U.Uonly, x=q, method="simple", u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=gs.theta, vall=V.ALL, marginal.method=Marginal.method);
    p = p - gs.et*gs.g/2;   
}

 
### calculate the U, K, and H based on proposed q after L steps of HMC
proposed.U = x.U.Uonly(q,v,r,y,sigma,gs.theta,V.ALL,Marginal.method);
proposed.K = sum(p^2/(2*diag.sigma.p^2));
proposed.H = proposed.U+proposed.K;

### calculate the acceptance probability
acceptProb = min(1, exp(current.H - proposed.H));

sqr.diff.q <- sum((current.q-q)^2);

#cat("square diff q is ", sqr.diff.q, ". c.U-p.U=", current.U-proposed.U, ", c.K-p.K=", current.K-proposed.K, "the accept Prob=", acceptProb ,"\n, c.H=", current.H, ", p.H=", proposed.H);
#cat("\n et=", gs.et, ", diff(q)^2=", sqr.diff.q, "\n");

   ### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
   if(runif(1,0,1) < acceptProb){
       tmp.beta <- q[1:3];
       tmp.alpha <- q[4:11];
       tmp.gamma <- q[12];
       tmp.full.alpha <- cbind(tmp.alpha, tmp.alpha, tmp.alpha);
       tmp.Delta <- get.Deltat(n.T, tmp.full.alpha, tmp.beta, tmp.gamma, V.ALL, gs.theta, Marginal.method);
       list(q = q, Update=1, delta=tmp.Delta, g=gs.g, U=proposed.U);  # accept
   }else{
       list(q = current.q, Update=0, delta=current.Delta, g=current.g, U=current.U);  # reject
   }
   
}


########################################################################################################################
#### (2) The standard congenial analysis using new.gradient.1127.post.loglike() for gradient calculation

###### udpate on 11/27/2014
### gradient of observed loglikelihood only (the log prior is not included) w.r.t (beta, alpha, gamma)
new.gradient.1127.loglike <- function(v, r, y.obs, delta, q, theta, NTIME, SIMU.MARGIN){
    beta <- q[1:3];
    alpha <- q[4:11];
    full.alpha <- cbind(alpha, alpha, alpha);
    gamma <- q[12];
    prob.impt <- new.all.obs.imputatemdl(v, r, y.obs, delta, full.alpha, gamma);
    
    inxreq2 <- which(r==2);
    inxreq3 <- which(r==3);
    inxrgeq2 <- union(inxreq2, inxreq3);
    inxrgeq3 <- which(r>=3);
    
    inxr2y0 <- intersect(inxrgeq2, which(y.obs[,1]==0))
    inxr2y1 <- intersect(inxrgeq2, which(y.obs[,1]==1))
    inxr3y0 <- intersect(inxrgeq3, which(y.obs[,2]==0))
    inxr3y1 <- intersect(inxrgeq3, which(y.obs[,2]==1))

    jac <- jacobian.full(y=q, func=get.Deltat.grad, parms=c(NTIME, theta));
    jac <- jac[1:5,1:12]; ### obtain d(delta)/d(beta, alpha, gamma)
    
    dftoddelta <- numeric(5);
    dftoddelta[1] <- sum(y.obs[,1]-prob.impt[,1]);
    dftoddelta[2] <- sum(y.obs[inxr2y0,2]-prob.impt[inxr2y0,2]);
    dftoddelta[3] <- sum(y.obs[inxr2y1,2]-prob.impt[inxr2y1,2]);
    dftoddelta[4] <- sum(y.obs[inxr3y0,3]-prob.impt[inxr3y0,3]);
    dftoddelta[5] <- sum(y.obs[inxr3y1,3]-prob.impt[inxr3y1,3]);
    dftodalpha <- sapply(1:8, function(x) sum((y.obs[,1]-prob.impt[,1])*v[,x])+sum((y.obs[inxrgeq2,2]-prob.impt[inxrgeq2,2])*v[inxrgeq2,x])+sum((y.obs[inxrgeq3,3]-prob.impt[inxrgeq3,3])*v[inxrgeq3,x]));
    dftogamma <- sum((y.obs[inxrgeq2,2]-prob.impt[inxrgeq2,2])*y.obs[inxrgeq2,1])+sum((y.obs[inxrgeq3,3]-prob.impt[inxrgeq3,3])*y.obs[inxrgeq3,2]);
        
    f.dftobeta <- sapply(1:3, function(x) dftoddelta%*%jac[1:5,x]);
    f.dftodalpha <- sapply(1:8, function(x) dftodalpha[x]+dftoddelta%*%jac[1:5,3+x]);
    f.dftogamma <- dftogamma+dftoddelta%*%jac[1:5,12];
    
    f.grad <- c(f.dftobeta, f.dftodalpha, f.dftogamma);
    f.grad;
}    


###### udpate on 11/27/2014
### gradient of the posterior loglikelihood (= observed loglikelihood+ prior loglikelihood) w.r.t. (Beta, Alpha, Gamma)
new.gradient.1127.post.loglike <- function(v, r, y.obs, delta, q, sigma, theta, NTIME, SIMU.MARGIN){
    beta <- q[1:3];
    alpha <- q[4:11];
    full.alpha <- cbind(alpha,alpha,alpha);
    gamma <- q[12];
    
    f.grad.loglike <- new.gradient.1127.loglike(v,r,y.obs,delta, q, theta, NTIME, SIMU.MARGIN);
    f.grad.prior <- -c(beta,alpha,gamma)/sigma^2;
    f.grad <- -(f.grad.loglike+f.grad.prior);
    f.grad;
}    


######### Used in new approach of gradient calculation. get ddelta/(dbeta, dalpha, dgamma)
get.Deltat.grad <- function(t=0, q, pars) {
   with (as.list(c(q, pars)), {
    ##all delta t
    rst.delta <- matrix(0, NTIME, 2);
    nvall<- nrow(V.INX.ALL);
 
    ##Pr(Y0=1|V)=0 as Y0=1 with prob 1
    prev.pyt.gvn.v <- rep(0, nvall);
    beta <- q[1:3];
    alpha <- q[4:11];
    alpha <- cbind(alpha, alpha, alpha);
    gamma <- q[12];

    ##Pr(V|Y0=y) for (all v & y=0) as col 1, and (all v & y=1) as col 2
    ##hence the 1st column=probv, and the 2nd column=0 as Y0=1 with prob 0
    v.gvn.prevy0 <- est.theta;
    v.gvn.prevy1 <- rep(0, nvall);
    
    for (tt in 1:NTIME) {
     alpha.t <- alpha[,tt];

     if(tt==1) 
     { 
       rst1 <- nleqslv.fn.t1(tt, alpha.t, beta, gamma, V.INX.ALL, v.gvn.prevy0, v.gvn.prevy1, SIMU.MARGIN);
       rst.delta[1,1] <- rst1; }
     else
     {  rst2 <- nleqslv.fn.t2(tt, alpha.t, beta, gamma, V.INX.ALL, v.gvn.prevy0, v.gvn.prevy1, SIMU.MARGIN);
        rst.delta[tt,] <- rst2;
     }
          
      ### Pr(Yt=1|V) for all v
        prev.pyt.gvn.v <- AMAR.eygvnvonimputate.fast(rst.delta[tt, ], tt, V.INX.ALL, prev.pyt.gvn.v, alpha.t, gamma);
  
        ### Pr(V|Yt) for all v and Yt in (0,1)
        v.gvn.prevy <- AMAR.vgvnprevy(est.theta, prev.pyt.gvn.v);
        v.gvn.prevy0 <- v.gvn.prevy[,1];
        v.gvn.prevy1 <- v.gvn.prevy[,2];
        #cat("For T=",tt, ", Delta=", rst.delta[tt,], ".\n");
     } 
    return(as.list(c(rst.delta[1,1], rst.delta[2,1], rst.delta[2,2], rst.delta[3,1], rst.delta[3,2])));
  })
}


#####################################################################################################################################
####### Compare two different approach for calculating the gradient
cmpr.gradient.time <- function(test){

      pmpt <- proc.time();
      true.g <-grad(x.U.Uonly, x=q, method="simple",u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=true.theta, vall=V.INX.ALL, marginal.method=3);
      apptime <- proc.time()-pmpt;
      apptime;
      cat(apptime, "\n");
      
      #############################################
      
      pmpt <- proc.time();
      dftoddelta <- numeric(5)
      prob.impt <- new.all.obs.imputatemdl(v, r, y, DELTA.ALL, ALPHA, GAMMA)
      inxrgeq2 <- which(r>=2);
      inxrgeq3 <- which(r==3);
      inxr2y0 <- intersect(inxrgeq2, which(y[,1]==0))
      inxr2y1 <- intersect(inxrgeq2, which(y[,1]==1))
      inxr3y0 <- intersect(inxrgeq3, which(y[,2]==0))
      inxr3y1 <- intersect(inxrgeq3, which(y[,2]==1))
      dftoddelta[1] <- sum(y[,1]-prob.impt[,1])
      dftoddelta[2] <- sum(y[inxr2y0,2]-prob.impt[inxr2y0,2])
      dftoddelta[3] <- sum(y[inxr2y1,2]-prob.impt[inxr2y1,2])
      dftoddelta[4] <- sum(y[inxr3y0,3]-prob.impt[inxr3y0,3])
      dftoddelta[5] <- sum(y[inxr3y1,3]-prob.impt[inxr3y1,3])
      
      jac <- jacobian.full(y=c(q), func=get.Deltat.grad, parms=c(NTIME, true.theta, SIMU.MARGIN));
      jac <- jac[1:5,1:12];
      
      dftodalpha <- numeric(8);
      dftodalpha <- sapply(1:8, function(x) sum((y[,1]-prob.impt[,1])*v[,x])+sum((y[inxrgeq2,2]-prob.impt[inxrgeq2,2])*v[inxrgeq2,x])+sum((y[inxrgeq3,3]-prob.impt[inxrgeq3,3])*v[inxrgeq3,x]));
      
      inxreq2 <- which(r==2);
      inxreq3 <- inxrgeq3;
      dftogamma <- sum((y[inxreq2,2]-prob.impt[inxreq2,2])*y[inxreq2,1])+sum((y[inxreq3,2]-prob.impt[inxreq3,2])*y[inxreq3,1]+(y[inxreq3,3]-prob.impt[inxreq3,3])*y[inxreq3,2]);
      
      f.dftobeta <- sapply(1:3, function(x) dftoddelta%*%jac[1:5,x]);
      f.dftodalpha <- sapply(1:8, function(x) dftodalpha[x]+dftoddelta%*%jac[1:5,3+x]);
      f.dftogamma <- dftogamma+dftoddelta%*%jac[1:5,12];
      
      f.grad <- c(f.dftobeta, f.dftodalpha, f.dftogamma);
      apptime <- proc.time()-pmpt;
      apptime;
      cat(apptime, "\n");
}



#####################################################################################################################################

### update from new.hmcstandard with the grad() calculation replaced grad(x.U.Uonly) by new.gradient.1127.post.loglike() 
### hmc function for posterior sampling of parameter in imputation and transition model
new.hmcStandard.f <- function(true.u, v, r, y.obs, q, sigma, n.iter, epsilon, L, V.ALL, theta, Marginal.method, diag.sigma.p){

#x is the matrix of the covariates with first column=1 as intercept
n.obs = nrow(v); #get number of observations	
n.param = length(q);
n.T = ncol(y.obs);

#post.samp used to store the parameters value at every iteration
post.samp <- matrix(NA, n.iter, n.param+1);

### save the current value of position variable
current.q <- q;
true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
### n.accept used to record how many times the proposed is accepted
n.accept <- 0;
acceptrate <- 0;
### calculate for each iteration and save the value of position variable for each iteration
for (i in 1:n.iter){
#cat("\n Begin Iter=", i, "#################################################################\n");
  if(i==1) {
    rst.xU <- x.U(current.q, v, r, y.obs, sigma, theta, V.ALL, Marginal.method);
    current.U <- rst.xU$U;
    current.Delta <- rst.xU$Delta;
    #current.g2 <- grad(x.U.Uonly, x=current.q, method="simple",u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
    current.g <- new.gradient.1127.post.loglike(v, r, y.obs, current.Delta, current.q, sigma, theta, n.T, Marginal.method);
    cat("c.g==", round(current.g, digit=3), "\n");
    }

    samp <- new.getSamp.f(v,r, y.obs, current.q, current.U, current.g, sigma, n.param, epsilon, L, current.Delta, theta, diag.sigma.p, V.ALL, acceptrate, Marginal.method);

  ### update the position value, U, and gradient of U
     current.q <- samp$q;
     current.Delta <- samp$Delta;
     current.U <- samp$U;
     current.g <- samp$g;
     current.update <- samp$Update; ### 1 for accept the prosed; 0 for reject the propose
          
  ### save the result of iteration i
     post.samp[i,] <- c(current.q, current.U);

     n.accept <- n.accept+1*current.update;
     acceptrate <- (n.accept)/(i);
   

        ### in case the acceptance rate is small, quit to save time 
        if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
        {
          cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
          list(samp=post.samp, acceptrate=n.accept/n.iter)
          break
        }
         
         if (current.update==0)
         {
          cat("###Reject at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         } else {
          cat("###Accept at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         }
       }

    post.samp;
}

################################################################################################################
### update from new.getSamp.f with the grad() calculation replaced grad(x.U.Uonly) by new.gradient.1127.post.loglike() 
### getSamp is the function to complete one update of Hamiltonian Monte Carlo
new.getSamp.f <- function(v,r,y.obs, gs.q, gs.U, gs.g, sigma, n.param, gs.et, L, gs.Delta, gs.theta, diag.sigma.p, V.ALL, acptrate, Marginal.method){
### first save the current value of position variable q
### generate random momentum variable p with the same length with q with normal distribution (0,1)
### save the current value of gradient of du/dq

current.q = gs.q;
current.g = gs.g;
current.Delta= gs.Delta;
current.U = gs.U;

n.T = ncol(y.obs);

var.p <- matrix(0, n.param, n.param);
diag(var.p) <- diag.sigma.p^2;
mean.p <- rep(0,n.param);
p <- rmvnorm(1, mean.p, var.p);
#cat("current p is", p, ".\n");

q = current.q;
### save the newly generated momentum variable as current.p
### then the current value of K(p)=(p^T%*%p)/2 as current.K based on current.p
### then the current value of H=K(p)+U(q), save the value as current.H
current.p = p;
current.K = sum(current.p^2/(2*diag.sigma.p^2));
current.H = current.U + current.K;

gs.et <- runif(1, 0.8*gs.et, 1.2*gs.et);

gamma.lowbound <- -3*diag.sigma.p[12];
gamma.upbound <- (-1)*gamma.lowbound;

for(leap in 1:L){
    p = p - gs.et*gs.g/2;
    # q = q + gs.et*sapply(1:n.param, function(x) p[x]/(diag.sigma.p[x])^2);
    q = q + gs.et*p/(diag.sigma.p^2);
    while(q[12]< gamma.lowbound|| q[12]>gamma.upbound){
          if(q[12]< gamma.lowbound) {
              q[12]=2*gamma.lowbound-q[12];
              p[12]=-p[12];
          }
          if(q[12]> gamma.upbound) {
              q[12]=2*gamma.upbound-q[12];
              p[12]=-p[12];
          }
      }
    #### since q is new, need update gs.Delta based on this q first
        tmp.beta <- q[1:3];
        tmp.alpha <- q[4:11];
        tmp.gamma <- q[12];
        tmp.full.alpha <- cbind(tmp.alpha, tmp.alpha, tmp.alpha);
        tmp.Delta <- get.Deltat(n.T, tmp.full.alpha, tmp.beta, tmp.gamma, V.ALL, gs.theta, Marginal.method);
    
    #gs.g = grad(x.U.Uonly, x=q, method="simple", u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=gs.theta, vall=V.ALL, marginal.method=Marginal.method);
    gs.g <- new.gradient.1127.post.loglike(v, r, y.obs, tmp.Delta, q, sigma, gs.theta, n.T, Marginal.method); ### new approach for gradient calculation, save time for large sample size
    p = p - gs.et*gs.g/2;   
}

### calculate the U, K, and H based on proposed q after L steps of HMC
proposed.U = x.U.Uonly(q,v,r,y.obs,sigma,gs.theta,V.ALL,Marginal.method);
proposed.K = sum(p^2/(2*diag.sigma.p^2));
proposed.H = proposed.U+proposed.K;

### calculate the acceptance probability
acceptProb = min(1, exp(current.H - proposed.H));

sqr.diff.q <- sum((current.q-q)^2);
cat("\n et=", gs.et, ", diff(q)^2=", sqr.diff.q, "\n");

### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
if(runif(1,0,1) < acceptProb){
    list(q = q, Update=1,  delta=tmp.Delta, g=gs.g, U=proposed.U);  # accept
}else{
    list(q = current.q, Update=0, delta=current.Delta, g=current.g, U=current.U);  # reject
}

}





##########################################################################################
##########################################################################################
##########################################################################################

###### For using Dirichelt method for p(v) estimation, use this method, check on 07/07/2015
xuan.pvdirichlet.multilevel <- function(k, data.i, niter, nburn, max.h, simu.v, a.alpha, b.alpha, V.INX.ALL, V.CATS, true.theta){

  n.v <- ncol(simu.v);
  n.obs <- nrow(simu.v);
  pchi.theta <- 0;
  
  ## save for each iteration,
  ind.v.length <- sapply(1:length(V.CATS), function(x) length(V.CATS[[x]]));
  ttl.level <- sum(ind.v.length);
  prev.length.ind.v <- c(0, sapply(1:length(V.CATS), function(x) sum(ind.v.length[1:x])));
  rst.thai <- matrix(0, max.h, ttl.level);
  
  rst.u <- array(0,n.obs);
  rst.z <- array(0,n.obs);
  rst.v <- array(0,max.h);
  rst.theta <- matrix(0, niter-1, nrow(V.INX.ALL));
  n.hitmax <- 0;
                    
  ###########################  Begin initializing values  ########################################################
  

  ### initialize alpha first based on paper  alpha ~ gamma(a.alpha, b.alpha)
  #rst.alpha <- rgamma(1, a.alpha, b.alpha);
  rst.alpha <- 0.5;
  rst.v <- numeric(max.h);
  rst.v <- rbeta(max.h, 1, rst.alpha);
  rst.v[max.h] <- 1; # to make the sum of sticks equals to 1
  tmp.omega <- cal.omega(max.h, rst.v);

  ### Based on rst.z, rst.v, rst.w (tmp.omega) generate rst.u
  rst.u <- runif(n.obs, 0, 1/(max.h));
  rst.z <- sapply(1:n.obs, function(x) sample(which(rst.u[x]<tmp.omega),1));
  unique.z <- unique(rst.z);
  rank.z <- rank(unique.z);
  update.rstz <- sapply(1:n.obs, function(x) rank.z[which(unique.z==rst.z[x])]);
  update.tmp.omega <- c(tmp.omega[sort(unique.z)], tmp.omega[-unique.z]);
  rst.z <- update.rstz;
  tmp.omega <- update.tmp.omega;
  ### based on each individual length, update rst.v
     update.rstv <- numeric(max.h);
     update.rstv[1] <- tmp.omega[1];
     for(j in 2:max.h) {
        update.rstv[j] <- tmp.omega[j]/prod(sapply(1:(j-1), function(y) (1-update.rstv[y])));
     }   
  rst.v <- update.rstv; # added on 2015/02/12
  max.z <- max(rst.z);
  
   z.inx.list <- list();
     for(i in 1:max.z) {
        z.inx.list[[i]] <- which(rst.z==i);
      }

    for(h.iter in 1:max.z) {
      tmp.thai <- NULL;
      for(j in 1:n.v) { # for each categorical variable
         # in our case each covariate having at least 2 levels
         freq.v <- sapply(1:ind.v.length[j], function(x) length(which(simu.v[ z.inx.list[[h.iter]], j ]==V.CATS[[j]][x])));
         prob.vj <- freq.v/sum(freq.v);
         tmp.thai <- c(tmp.thai, rdirichlet(1, prob.vj));
      }
      #cat("tmp.thai==", tmp.thai, "\n");
      rst.thai[h.iter,] <- tmp.thai;
    }
  
   if(max.h>max.z) {
     for(h.iter in (max.z+1):max.h){
        tmp.thai <- NULL;
        for(j in 1:n.v) {
            tmp.thai <- c(tmp.thai, rdirichlet(1, rep(1, ind.v.length[j])));
        }
     rst.thai[h.iter, ] <- tmp.thai;
     }
   }  
     
       #### get the corresponding v.index number for each subjects
     sub.indx <- apply(simu.v,1, function(x) vtosimvalue(x,V.CATS))+1;

  ############################ Finish Initialize ###################################################################
   
  ##### Begin Iteration
  for(s in 2:niter) {
   max.z <- max(rst.z);
   ### step 1 update u
    tmp.omega <- cal.omega(max.h, rst.v);
    rst.u <- sapply(1:n.obs, function(x) runif(1, 0, tmp.omega[rst.z[x]]));
    #cat("\n finishi (1)");

   ### step 2 update thai
    # make a list such that the h^th element contains the index of subjects with z==h
     z.inx.list <- list();
     for(i in 1:max.z) {
        z.inx.list[[i]] <- which(rst.z==i);
      }
   
   ### step 3 update v  #################################################################  
   if(max.z==1) {
     lowb <- max(rst.u);
     upb <- 1;
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[1] <- rst$samp.v;
     failed <- rst$failed;
     #cat("for h=1, bound=(", lowb, ",", upb, ")"); 
     rst.v[2:(max.h)] <- rbeta((max.h-1), 1, rst.alpha);
     ############################################################ 
    }else if(max.z==2) {
     # for v1
       # if there is no subject with Z==1
       if(length(z.inx.list[[1]])==0) {
         lowb <- 0;
       }else{
        lowb <- max(rst.u[z.inx.list[[1]]]);
       }
       upb <- 1-max(rst.u[z.inx.list[[2]]])/rst.v[2];
       #cat("for h=1, bound=(", lowb, "," ,upb, ")"); 
       rst <- samp.ind.v(lowb, upb, rst.alpha);
       rst.v[1] <- rst$samp.v;
       failed <- rst$failed;
       
     # for v2 
     lowb <- min(1,max(rst.u[ z.inx.list[[2]] ])/(1-rst.v[1]));
     upb <- 1;
     #cat("for h=2, bound=(", lowb, "," ,upb, ")");     
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[2] <- rst$samp.v;
     failed <- failed+rst$failed;
     if(failed>=1) {   failed <- 1; }
     ### for v3 and beyond 
     rst.v[3:(max.h)] <- rbeta((max.h-2), 1, rst.alpha);
    ############################################################ 
    }else{   ### for (max.z>=3)
    ###(1) for updating rst.v[1]
       ## (1.1) low bound for V_h always come from groups with Z==h
     if(length( z.inx.list[[1]] )==0) {
       lowb <- 0; }else{
       lowb <- max(rst.u[ z.inx.list[[1]] ]);
     }  
       ## (1.2) ## Find upbound for V_(h=1), coming from groups with Z= (h+1), (h+2), ... (max.h)
      # for group of Z equal to 2
     if(length( z.inx.list[[2]] )==0) {
       max1 <- 0;
     }else{
       max1 <- max(rst.u[ z.inx.list[[2]] ])/rst.v[2];
     }
       upb1 <- 1-max1;       
     ## (1.3)
     #  for group of Z great than 3 
     gt.h <- c(3:max.z);
     length.gt.h <- max.z-3+1;
     max.num <- numeric(length.gt.h);
     max.num <- sapply(gt.h, function(x)  max(0, rst.u[ z.inx.list[[x]] ]));  # update on 07/29
     denom.gt.h <- numeric(length.gt.h);
     denom.gt.h[1] <- rst.v[3]*(1-rst.v[2]);  # from group z=3
     # update on 7/29
     if(length.gt.h>=2) { # from group z=4, 5, ...
       denom.gt.h[2:length.gt.h] <- sapply(4:max.z, function(x) rst.v[x]*prod(sapply(2:(x-1), function(y) 1-rst.v[y])));
     } 
     max2 <- min(1,max(0,max.num/denom.gt.h));
     upb2 <- 1-max2;
     #cat("max2==", max2, "\n");
     #upb <- 1-max(max1,max2);
     upb <- min(upb1, upb2);
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[1] <- rst$samp.v;
     failed <- rst$failed;
     #cat("for h=1, bound=(", lowb, ",", upb, ")-max1=", max1, "-ma2=", max2, ", with rst.v[1]==", rst.v[1], "\n");     
   ###  End for updating rst.v[1] with max.z>=3
     
     
   ### (2) for updating rst.v[h], for h=2, ..., (max.z-1) with max.z>=3
     for(h.iter in 2:(max.z-1)) {
       # (2.0) for lower bound
       if( length(z.inx.list[[h.iter]] )==0) {
          lowb <- 0;
        }else{
         lowb <- max(rst.u[ z.inx.list[[h.iter]] ])/prod(sapply(1:(h.iter-1), function(x) 1-rst.v[x]));
        }
       # (2.1) find upper bound max1 for rst.v[h] from group with Z=h+1
       #indxeq.hiterplus1 <- which(rst.z==(h.iter+1));
       if(length( z.inx.list[[h.iter+1]]) ==0) {
         max1 <- 0;
       }else{
          #cat("numer=", max(rst.u[ z.inx.list[[h.iter+1]] ]), "***denom==", (rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x])))), "\n");
          max1 <- max(rst.u[ z.inx.list[[h.iter+1]] ])/(rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x]))));            }
       upb1 <- 1-max1;
       
       #(2.2) find upper bound max2 for rst.v[h] from group with Z=h+2, h+3, ...
       if(h.iter+1>=max.z) { # on the equality with h.iter=max.z-1, the upbound taken care by upb1
          max2 <- 0;
       }else {   
          ### Begin Alternative approach
          whichgt.h <- c((h.iter+2):(max.z));  ## {h*:h*>h.iter}
          lengthgt.h <- length(whichgt.h);
          
          num.gt.h <- numeric(lengthgt.h);
          num.gt.h <- sapply(whichgt.h, function(x) max(0, rst.u[ z.inx.list[[x]]  ]));
          indx.num.gt0 <- which(num.gt.h>0);
          num.gt0 <- num.gt.h[indx.num.gt0];
          denom.gt0 <- numeric(length(indx.num.gt0));
          
          denom.gt0 <- sapply(whichgt.h[indx.num.gt0], function(x) rst.v[x]*prod(sapply(1:(h.iter-1), function(y) 1-rst.v[y]))*prod(sapply((h.iter+1):(x-1), function(z)  1-rst.v[z])));
          
          #at("num.gt.h===", num.gt0, "\n");
          #cat("denom.gt.h===", denom.gt0, "\n");q
          max2 <- min(1,max(0, num.gt0/denom.gt0));
        }
          upb2 <- 1-max2;
          
          upb <- min(upb1, upb2);
          #cat("for h=", h.iter, ", bound=(", lowb, ",", upb, ")--max1=", max1, "--max2=", max2, ",");     
          rst<- samp.ind.v(lowb, upb, rst.alpha);
          rst.v[h.iter] <- rst$samp.v;
          failed <- failed+rst$failed;
     }
   
     #(3) #### for h.iter==max.z ( hence definitely some subjects fall in this group)
     lowb <- max(rst.u[ z.inx.list[[max.z]] ])/prod(sapply(1:(max.z-1), function(x) 1-rst.v[x]));
     upb <- 1;
     #cat("for h=", max.z, ", bound=(", lowb, ",", upb, ")");    
     rst <- samp.ind.v(lowb, upb, rst.alpha);
     rst.v[max.z] <- rst$samp.v;
     failed <- failed+rst$failed;
     #(4) update rst.v[h] for h=max.z+1, ..., max.h
     if((max.z+1)<=max.h){
     rst.v[(max.z+1):(max.h)] <- rbeta((max.h-max.z),1,rst.alpha);
     }
     if(failed>=1) { failed <- 1; }
   } ## finish for max.z>=3

   ### keept last v=1 for sum length=1
   rst.v[max.h] <- 1;
   tmp.omga <- cal.omega(max.h, rst.v);
   rst.v[which(rst.v[1:(max.h-1)]==1)] <- 1-(10)^(-10); ### to exclude the case that middle part rst.v==1, hence (1-rst.v)=0 in the denominator on Sep17
   #cat(",  (3)");
   # after updating rst.v, need to update tmp.omega   
   ###########  Finish step 3  for updating V   ############### ########### ########### ########### ########### ########### ########### ###########  

  ### Update on 06/08-- to save time, first calculate the p(v|h) and use it in updating Z
         ind.h.prob <- matrix(0, max.h,nrow(V.INX.ALL));
         for(tmp.h in 1:max.h) {
            tmp.thai <- rst.thai[tmp.h,];
            ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(v) prod(tmp.thai[sapply(1:n.v, function(x) 1+length(which(V.CATS[[x]] < v[x]))+prev.length.ind.v[x])]));
            #ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
         }
         #cat("summary of ind.h.prob==", summary(apply(ind.h.prob,1,sum)), "\n");
       
    ### step 4 
     # update omega first, need to calculate for max.h
     min.u <- min(rst.u);
     #cat("min.u==", min.u, ",   ");
     addlength <- sapply(1:(max.h), function(x) sum(tmp.omega[1:x]));
     k.tilt <- min(max.h, min(which(addlength>=(1-min.u))));
     #cat(", with k.tilt==", k.tilt, "\n");

       for(i in 1:n.obs) {
         tmp.A <- which(tmp.omega[1:k.tilt]>rst.u[i]);
         #cat("for sub==", i, ", sample from(", tmp.A, ")\n");
         size.A <- length(tmp.A);
         # calculate the multinomial distribution probability
         if(size.A==0) {
          rst.z[i] <- rst.z[i];
         }else if(size.A==1){
          rst.z[i] <- tmp.A[1];
         }else{  
          prob.A <- numeric(size.A);
          prob.A <- ind.h.prob[tmp.A, sub.indx[i]];
          prob.A <- prob.A/sum(prob.A);
          rst.z[i] <- tmp.A[which(rmultinom(1,1,prob.A)==1)];
         }
      }
    #cat(",  (4)");
       #max.z <- max(rst.z);
       #if(max.z==max.h) n.hitmax <- n.hitmax+1;
       #cat("updated with z==", rst.z, " with max.z==", max(rst.z), "\n");
     
   #### step 4.5 -- only update in 0928 function()
   unique.z <- unique(rst.z);
   rank.z <- rank(unique.z);
   update.rstz <- sapply(1:n.obs, function(x) rank.z[which(unique.z==rst.z[x])]);
   update.tmp.omega <- c(tmp.omega[sort(unique.z)], tmp.omega[-unique.z]);
   rst.z <- update.rstz;
   tmp.omega <- update.tmp.omega;
      ### based on each individual length, update rst.v
      update.rstv <- numeric(max.h);
      update.rstv[1] <- tmp.omega[1];
      for(j in 2:max.h) {
         update.rstv[j] <- tmp.omega[j]/prod(sapply(1:(j-1), function(y) (1-update.rstv[y])));
      }
    rst.v <- update.rstv;  # added on 2015/02/12
    max.z <- max(rst.z);

  ### step 2 update thai
     # make a list such that the h^th element contains the index of subjects with z==h
     z.inx.list <- list();
     for(i in 1:max.z) {
        z.inx.list[[i]] <- which(rst.z==i);
      }

        for(h.iter in 1:max.z) {
      tmp.thai <- NULL;
      for(j in 1:n.v) { # for each categorical variable
         # in our case each covariate having 2 level, hence only need to save the 1st probability
         freq.v <- sapply(1:ind.v.length[j], function(x) length(which(simu.v[ z.inx.list[[h.iter]], j ]==V.CATS[[j]][x])));
         prob.vj <- freq.v/sum(freq.v);
         tmp.thai <- c(tmp.thai, rdirichlet(1, prob.vj));
      }
      #cat("tmp.thai==", tmp.thai, "\n");
      rst.thai[h.iter,] <- tmp.thai;
    }
  
   if(max.h>max.z) {
     for(h.iter in (max.z+1):max.h){
        tmp.thai <- NULL;
        for(j in 1:n.v) {
            tmp.thai <- c(tmp.thai, rdirichlet(1, rep(1, ind.v.length[j])));
        }
     rst.thai[h.iter, ] <- tmp.thai;
     }
   }  

      ind.h.prob <- matrix(0, max.h, nrow(V.INX.ALL));
      for(tmp.h in 1:max.h) {
            tmp.thai <- rst.thai[tmp.h,];
            ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(v) prod(tmp.thai[sapply(1:n.v, function(x) 1+length(which(V.CATS[[x]] < v[x]))+prev.length.ind.v[x])]));
            #ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
         }
       rst.theta[s-1,] <- t(ind.h.prob)%*%tmp.omega;
       #cat("sum of theta==", sum(rst.theta[s-1,]), "\n");
      
        ### step 5
      # update omega first
    tmp.shape <-  a.alpha+max.z;
     if(max.z==max.h) {
       inx.neq1 <- which(rst.v[1:(max.h)]<1);
       tmp.rate <- Inf;
       rst.alpha <- (10)^(-10);
     }else{
       inx.neq1 <- which(rst.v[1:(max.z)]<1);
       tmp.rate <- b.alpha-log(1-sum(tmp.omega[1:max.z]));
       rst.alpha <- rgamma(1, shape=tmp.shape, rate=tmp.rate);
     }
  
   #cat("alpha (a==", tmp.shape, ", b==", tmp.rate, ")\n");
    #cat(",  (5)");
     ### calculate the pchi after each iteration
    if(s>2) {
      if(s<=(nburn+2)) {
        m <- apply(rst.theta[1:(s-1),],2,mean);
      }else{
        m <- apply(rst.theta[nburn:(s-1),],2,mean);
      }  
    mean.pchi.theta <- probv.pchi(m, true.theta);
    crnt.pchi.theta <- probv.pchi(rst.theta[s-1,], true.theta);
      
    cat(" at End iteration ", s, ", pchi===( ", mean.pchi.theta, ",", crnt.pchi.theta, "), ****max.z===", max(rst.z), ", alpha=", rst.alpha, ", for (K,i)==(", k, ",", data.i, ") with n.hitmax==", n.hitmax, "\n");
   }
   #cat(",  (6)");
   #cat("updated with u==", rst.u, "\n");
   #cat("updated with v==", rst.v[1:max.z], "\n");
   #beyond.length <- c(beyond.length, 1-sum(tmp.omega[1:(max.z)]));
   if(max.z==max.h) { n.hitmax <- n.hitmax+1; }
   
  } ### end of iteration
   list(theta=rst.theta, n.hitmax=n.hitmax);
}




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################









###   new.hmcStandard.parallel <- function(true.u, v,r, y, q, sigma, n.iter, epsilon, L, V.ALL, theta, Marginal.method, diag.sigma.p){
###   
###   #x is the matrix of the covariates with first column=1 as intercept
###   n.obs = nrow(v); #get number of observations	
###   n.param = length(q);
###   n.T = ncol(y);
###   
###   #post.samp used to store the parameters value at every iteration
###   post.samp <- matrix(NA, n.iter, n.param);
###   post.samp.U <- rep(0, n.iter);
###   
###   ### save the current value of position variable
###   current.q <- q;
###   true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
###   ### n.accept used to record how many times the proposed is accepted
###   n.accept <- 0;
###   acceptrate <- 0;
###   ### calculate for each iteration and save the value of position variable for each iteration
###   for (i in 1:n.iter){
###   #cat("\n Begin Iter=", i, "#################################################################\n");
###     if(i==1) {
###       #current.beta <- current.q[1:3];    
###       #current.alpha <- current.q[4:11];
###       #current.gamma <- current.q[12];    
###       #c.fullalpha <- cbind(current.alpha, current.alpha, current.alpha);
###       #once theta is updated, Delta needs to be recalculated, and U, grad.U
###       #current.Delta <- get.Deltat(n.T, c.fullalpha, current.beta, current.gamma, V.ALL, theta, Marginal.method);  
###      
###       rst.xU <- x.U.parallel(current.q, v, r, y, sigma, theta, V.ALL, Marginal.method);
###       current.U <- rst.xU$U;
###       current.Delta <- rst.xU$Delta;
###       current.g <- grad(x.U.Uonly.parallel, x=current.q, method="simple",u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
###      }
###   
###       samp <- new.getSamp.parallel(v,r, y, current.q, current.U, current.g, sigma, n.param, epsilon, L, current.Delta, theta, diag.sigma.p, V.ALL, acceptrate, Marginal.method);
###   
###     ### update the position value, U, and gradient of U
###        current.q <- samp$q;
###        current.Delta <- samp$Delta;
###        current.U <- samp$U;
###        current.g <- samp$g;
###        current.update <- samp$Update; ### 1 for accept the prosed; 0 for reject the propose
###            
###     ### save the result of iteration i
###        post.samp[i,] <- current.q;
###        post.samp.U[i] <- current.U;
###   
###        n.accept <- n.accept+1*current.update;
###        acceptrate <- (n.accept)/(i);
###      
###           ### in case the acceptance rate is small, quit to save time 
###           if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
###           {
###             cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
###             list(samp=post.samp, acceptrate=n.accept/n.iter)
###             break
###           }
###            
###            if (current.update==0)
###            {
###             cat("###Reject at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
###            } else {
###             cat("###Accept at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
###            }
###          }
###     rst <- cbind(post.samp, post.samp.U);
###     rst;
###   }
###   
###   
###   ### getSamp is the function to complete one update of Hamiltonian Monte Carlo
###   new.getSamp.parallel <- function(v,r,y, gs.q, gs.U, gs.g, sigma, n.param, gs.et, L, gs.Delta, gs.theta, diag.sigma.p, V.ALL, acptrate, Marginal.method){
###   ### first save the current value of position variable q
###   ### generate random momentum variable p with the same length with q with normal distribution (0,1)
###   ### save the current value of gradient of du/dq
###   
###   current.q = gs.q;
###   current.g = gs.g;
###   current.Delta= gs.Delta;
###   current.U = gs.U;
###   
###   n.T = ncol(y);
###   
###   var.p <- matrix(0, n.param, n.param);
###   diag(var.p) <- diag.sigma.p^2;
###   mean.p <- rep(0,n.param);
###   p <- rmvnorm(1, mean.p, var.p);
###   q = current.q;
###   ### save the newly generated momentum variable as current.p
###   ### then the current value of K(p)=(p^T%*%p)/2 as current.K based on current.p
###   ### then the current value of H=K(p)+U(q), save the value as current.H
###   current.p = p;
###   current.K = sum(current.p^2/(2*diag.sigma.p^2));
###   current.H = current.U + current.K;
###   
###   gs.et <- runif(1, 0.8*gs.et, 1.2*gs.et);
###   gamma.lowbound <- (-3)*diag.sigma.p[12];
###   gamma.upbound <- (-1)*gamma.lowbound;
###   
###   for(leap in 1:L){
###       p = p - gs.et*gs.g/2;
###       # q = q + gs.et*sapply(1:n.param, function(x) p[x]/(diag.sigma.p[x])^2);
###       q = q + gs.et*p/(diag.sigma.p^2);
###       while(q[12]< gamma.lowbound|| q[12]>gamma.upbound){
###             if(q[12]< gamma.lowbound) {
###                 q[12]=2*gamma.lowbound-q[12];
###                 p[12]=-p[12];
###             }
###             if(q[12]> gamma.upbound) {
###                 q[12]=2*gamma.upbound-q[12];
###                 p[12]=-p[12];
###             }
###         }
###       #### since q is new, need update gs.Delta based on this q first
###           #tmp.beta <- q[1:3];
###           #tmp.alpha <- q[4:11];
###           #tmp.gamma <- q[12];
###           #tmp.full.alpha <- cbind(tmp.alpha, tmp.alpha, tmp.alpha);
###           #tmp.Delta <- get.Deltat(n.T, tmp.full.alpha, tmp.beta, tmp.gamma, V.ALL, gs.theta, Marginal.method);    
###       ### true gradient based on all subjects (1st choice)
###       gs.g = grad(x.U.Uonly.parallel, x=q, method="simple", u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=gs.theta, vall=V.ALL, marginal.method=Marginal.method);
###       p = p - gs.et*gs.g/2;   
###   }
###   
###    
###   ### calculate the U, K, and H based on proposed q after L steps of HMC
###   proposed.U = x.U.Uonly.parallel(q,v,r,y,sigma,gs.theta,V.ALL,Marginal.method);
###   proposed.K = sum(p^2/(2*diag.sigma.p^2));
###   proposed.H = proposed.U+proposed.K;
###   
###   ### calculate the acceptance probability
###   acceptProb = min(1, exp(current.H - proposed.H));
###   
###   sqr.diff.q <- sum((current.q-q)^2);
###   
###   #cat("square diff q is ", sqr.diff.q, ". c.U-p.U=", current.U-proposed.U, ", c.K-p.K=", current.K-proposed.K, "the accept Prob=", acceptProb ,"\n, c.H=", current.H, ", p.H=", proposed.H);
###   #cat("\n et=", gs.et, ", diff(q)^2=", sqr.diff.q, "\n");
###   
###      ### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
###      if(runif(1,0,1) < acceptProb){
###          tmp.beta <- q[1:3];
###          tmp.alpha <- q[4:11];
###          tmp.gamma <- q[12];
###          tmp.full.alpha <- cbind(tmp.alpha, tmp.alpha, tmp.alpha);
###          tmp.Delta <- get.Deltat(n.T, tmp.full.alpha, tmp.beta, tmp.gamma, V.ALL, gs.theta, Marginal.method);
###          list(q = q, Update=1, delta=tmp.Delta, g=gs.g, U=proposed.U);  # accept
###      }else{
###          list(q = current.q, Update=0, delta=current.Delta, g=current.g, U=current.U);  # reject
###      }
###      
###   }
###   





###   new.hmcStandard.f.parallel <- function(true.u, v, r, y.obs, q, sigma, n.iter, epsilon, L, V.ALL, theta, Marginal.method, diag.sigma.p){
###   
###   #x is the matrix of the covariates with first column=1 as intercept
###   n.obs = nrow(v); #get number of observations	
###   n.param = length(q);
###   n.T = ncol(y.obs);
###   
###   #post.samp used to store the parameters value at every iteration
###   post.samp <- matrix(NA, n.iter, n.param+1);
###   
###   ### save the current value of position variable
###   current.q <- q;
###   true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
###   ### n.accept used to record how many times the proposed is accepted
###   n.accept <- 0;
###   acceptrate <- 0;
###   ### calculate for each iteration and save the value of position variable for each iteration
###   for (i in 1:n.iter){
###   #cat("\n Begin Iter=", i, "#################################################################\n");
###     if(i==1) {
###       rst.xU <- x.U.parallel(current.q, v, r, y.obs, sigma, theta, V.ALL, Marginal.method);
###       current.U <- rst.xU$U;
###       current.Delta <- rst.xU$Delta;
###       #current.g2 <- grad(x.U.Uonly, x=current.q, method="simple",u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
###       current.g <- new.gradient.1127.post.loglike(v, r, y.obs, current.Delta, current.q, sigma, theta, n.T, Marginal.method);
###       cat("c.g==", round(current.g, digit=3), "\n");
###       }
###   
###       samp <- new.getSamp.f.parallel(v,r, y.obs, current.q, current.U, current.g, sigma, n.param, epsilon, L, current.Delta, theta, diag.sigma.p, V.ALL, acceptrate, Marginal.method);
###   
###     ### update the position value, U, and gradient of U
###        current.q <- samp$q;
###        current.Delta <- samp$Delta;
###        current.U <- samp$U;
###        current.g <- samp$g;
###        current.update <- samp$Update; ### 1 for accept the prosed; 0 for reject the propose
###             
###     ### save the result of iteration i
###        post.samp[i,] <- c(current.q, current.U);
###   
###        n.accept <- n.accept+1*current.update;
###        acceptrate <- (n.accept)/(i);
###      
###   
###           ### in case the acceptance rate is small, quit to save time 
###           if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
###           {
###             cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
###             list(samp=post.samp, acceptrate=n.accept/n.iter)
###             break
###           }
###            
###            if (current.update==0)
###            {
###             cat("###Reject at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
###            } else {
###             cat("###Accept at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
###            }
###          }
###   
###       post.samp;
###   }
###   
################################################################################################################
###   ### update from new.getSamp.f with the grad() calculation replaced grad(x.U.Uonly) by new.gradient.1127.post.loglike() 
###   ### getSamp is the function to complete one update of Hamiltonian Monte Carlo
###   new.getSamp.f.parallel<- function(v,r,y.obs, gs.q, gs.U, gs.g, sigma, n.param, gs.et, L, gs.Delta, gs.theta, diag.sigma.p, V.ALL, acptrate, Marginal.method){
###   ### first save the current value of position variable q
###   ### generate random momentum variable p with the same length with q with normal distribution (0,1)
###   ### save the current value of gradient of du/dq
###   
###   current.q = gs.q;
###   current.g = gs.g;
###   current.Delta= gs.Delta;
###   current.U = gs.U;
###   
###   n.T = ncol(y.obs);
###   
###   var.p <- matrix(0, n.param, n.param);
###   diag(var.p) <- diag.sigma.p^2;
###   mean.p <- rep(0,n.param);
###   p <- rmvnorm(1, mean.p, var.p);
###   #cat("current p is", p, ".\n");
###   
###   q = current.q;
###   ### save the newly generated momentum variable as current.p
###   ### then the current value of K(p)=(p^T%*%p)/2 as current.K based on current.p
###   ### then the current value of H=K(p)+U(q), save the value as current.H
###   current.p = p;
###   current.K = sum(current.p^2/(2*diag.sigma.p^2));
###   current.H = current.U + current.K;
###   
###   gs.et <- runif(1, 0.8*gs.et, 1.2*gs.et);
###   
###   gamma.lowbound <- -3*diag.sigma.p[12];
###   gamma.upbound <- (-1)*gamma.lowbound;
###   
###   for(leap in 1:L){
###       p = p - gs.et*gs.g/2;
###       # q = q + gs.et*sapply(1:n.param, function(x) p[x]/(diag.sigma.p[x])^2);
###       q = q + gs.et*p/(diag.sigma.p^2);
###       while(q[12]< gamma.lowbound|| q[12]>gamma.upbound){
###             if(q[12]< gamma.lowbound) {
###                 q[12]=2*gamma.lowbound-q[12];
###                 p[12]=-p[12];
###             }
###             if(q[12]> gamma.upbound) {
###                 q[12]=2*gamma.upbound-q[12];
###                 p[12]=-p[12];
###             }
###         }
###       #### since q is new, need update gs.Delta based on this q first
###           tmp.beta <- q[1:3];
###           tmp.alpha <- q[4:11];
###           tmp.gamma <- q[12];
###           tmp.full.alpha <- cbind(tmp.alpha, tmp.alpha, tmp.alpha);
###           tmp.Delta <- get.Deltat(n.T, tmp.full.alpha, tmp.beta, tmp.gamma, V.ALL, gs.theta, Marginal.method);
###       
###       #gs.g = grad(x.U.Uonly, x=q, method="simple", u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=gs.theta, vall=V.ALL, marginal.method=Marginal.method);
###       gs.g <- new.gradient.1127.post.loglike(v, r, y.obs, tmp.Delta, q, sigma, gs.theta, n.T, Marginal.method); ### new approach for gradient calculation, save time for large sample size
###       p = p - gs.et*gs.g/2;   
###   }
###   
###   ### calculate the U, K, and H based on proposed q after L steps of HMC
###   proposed.U = x.U.Uonly.parallel(q,v,r,y.obs,sigma,gs.theta,V.ALL,Marginal.method);
###   proposed.K = sum(p^2/(2*diag.sigma.p^2));
###   proposed.H = proposed.U+proposed.K;
###   
###   ### calculate the acceptance probability
###   acceptProb = min(1, exp(current.H - proposed.H));
###   
###   sqr.diff.q <- sum((current.q-q)^2);
###   cat("\n et=", gs.et, ", diff(q)^2=", sqr.diff.q, "\n");
###   
###   ### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
###   if(runif(1,0,1) < acceptProb){
###       list(q = q, Update=1,  delta=tmp.Delta, g=gs.g, U=proposed.U);  # accept
###   }else{
###       list(q = current.q, Update=0, delta=current.Delta, g=current.g, U=current.U);  # reject
###   }
###   
###   }
###   











### hmc function for posterior sampling of parameter in imputation and transition model
new.hmcStandard.noGamma <- function(true.u, v,r, y, q, u.gamma, sigma, n.iter, epsilon, L, V.ALL, theta, Marginal.method, diag.sigma.p){

#x is the matrix of the covariates with first column=1 as intercept
n.obs = nrow(v); #get number of observations	
n.param = length(q);
n.T = ncol(y);

#post.samp used to store the parameters value at every iteration
post.samp <- matrix(NA, n.iter, n.param);
post.samp.U <- rep(0, n.iter);

### save the current value of position variable
current.q <- q;
#true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7);
### n.accept used to record how many times the proposed is accepted
n.accept <- 0;
acceptrate <- 0;
### calculate for each iteration and save the value of position variable for each iteration
for (i in 1:n.iter){
#cat("\n Begin Iter=", i, "#################################################################\n");
  if(i==1) {   
    current.U <- newx.U2(current.q, u.gamma, v, r, y, sigma, theta, V.ALL, Marginal.method);
    current.g <- grad(newx.U2, x=current.q, method="simple", u.Gamma=u.gamma, u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
    cat("current.U==", current.U, "**** current.g===", current.g);
   }

    samp <- new.getSamp.noGamma(v,r, y, current.q, current.U, current.g, u.gamma, sigma, n.param, epsilon, L, current.Delta, theta, diag.sigma.p, V.ALL, acceptrate, Marginal.method);

  ### update the position value, U, and gradient of U
     current.q <- samp$q;
     current.U <- samp$U;
     current.g <- samp$g;
     current.update <- samp$Update; ### 1 for accept the prosed; 0 for reject the propose
         
  ### save the result of iteration i
     post.samp[i,] <- current.q;
     post.samp.U[i] <- current.U;

     n.accept <- n.accept+1*current.update;
     acceptrate <- (n.accept)/(i);
   
        ### in case the acceptance rate is small, quit to save time 
        if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
        {
          cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
          list(samp=post.samp, acceptrate=n.accept/n.iter)
          break
        }
         
         if (current.update==0)
         {
          cat("###Reject at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         } else {
          cat("###Accept at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         }
       }
  rst <- cbind(post.samp, post.samp.U);
  rst;
}


   

### getSamp is the function to complete one update of Hamiltonian Monte Carlo
new.getSamp.noGamma <- function(v,r,y, gs.q, gs.U, gs.g, gs.gamma, sigma, n.param, gs.et, L, gs.Delta, gs.theta, diag.sigma.p, V.ALL, acptrate, Marginal.method){
### first save the current value of position variable q
### generate random momentum variable p with the same length with q with normal distribution (0,1)
### save the current value of gradient of du/dq

current.q = gs.q;
current.g = gs.g;
current.U = gs.U;

n.T = ncol(y);

var.p <- matrix(0, n.param, n.param);
diag(var.p) <- diag.sigma.p^2;
mean.p <- rep(0,n.param);
p <- rmvnorm(1, mean.p, var.p);
q = current.q;
### save the newly generated momentum variable as current.p
### then the current value of K(p)=(p^T%*%p)/2 as current.K based on current.p
### then the current value of H=K(p)+U(q), save the value as current.H
current.p = p;
current.K = sum(current.p^2/(2*diag.sigma.p^2));
current.H = current.U + current.K;

gs.et <- runif(1, 0.8*gs.et, 1.2*gs.et);
for(leap in 1:L){
    p = p - gs.et*gs.g/2;
    # q = q + gs.et*sapply(1:n.param, function(x) p[x]/(diag.sigma.p[x])^2);
    q = q + gs.et*p/(diag.sigma.p^2);
    gs.g = grad(newx.U2, x=q, method="simple", u.Gamma=gs.gamma, u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=gs.theta, vall=V.ALL, marginal.method=Marginal.method);
    p = p - gs.et*gs.g/2;   
}

 
### calculate the U, K, and H based on proposed q after L steps of HMC
proposed.U = newx.U2(q, gs.gamma, v,r,y,sigma,gs.theta,V.ALL,Marginal.method);
proposed.K = sum(p^2/(2*diag.sigma.p^2));
proposed.H = proposed.U+proposed.K;

### calculate the acceptance probability
acceptProb = min(1, exp(current.H - proposed.H));

sqr.diff.q <- sum((current.q-q)^2);

#cat("square diff q is ", sqr.diff.q, ". c.U-p.U=", current.U-proposed.U, ", c.K-p.K=", current.K-proposed.K, "the accept Prob=", acceptProb ,"\n, c.H=", current.H, ", p.H=", proposed.H);
#cat("\n et=", gs.et, ", diff(q)^2=", sqr.diff.q, "\n");

   ### generate a uniform(0,1) variable to determin whether accept or reject the proposed state
   if(runif(1,0,1) < acceptProb){
       list(q = q, Update=1, g=gs.g, U=proposed.U);  # accept
   }else{
       list(q = current.q, Update=0, g=current.g, U=current.U);  # reject
   }
   
}


#### update on 04/20/
#### since Gamma gradient is almost 0
newx.U2 <- function(u.q, u.Gamma, u.v, u.r, u.y, u.sigma, theta, vall, marginal.method) {
  u.Beta <- u.q[1:3];
  u.Alpha <- u.q[4:11];
  #u.Gamma <- u.q[12];
  rst <- 0;
  n.T <- ncol(u.y);
  nsub <- nrow(u.y);
  ncolv <- ncol(u.v);

  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
  
  logPrior <- sum(-u.q^2/(2*u.sigma^2));
   
  vvec <- as.vector(t(u.v));
  vec.alpha <- as.vector(u.fullAlpha);
  yvec <- as.vector(t(u.y));
  yvec[which(is.na(yvec))] <- 0;
  delta <- as.vector(t(u.Delta));

  ### log.Like gives the log(likelihood) of the observed outcome 
  logLike <- .C("ttl_obsloglikelihood", as.integer(n.T), as.integer(nsub), as.double(vvec), as.integer(u.r), as.integer(ncolv), as.integer(yvec), as.double(delta), as.double(vec.alpha), as.double(u.Gamma), as.double(rst))[[10]];
  
  U= -(logLike + logPrior);
  #cat("logPrior=", logPrior, ", logLike=", logLike, "U=", U, ".\n");
  U;
}



### obs.loglikelihood value is invariant of Gamma value
obs.loglikelihood <- function(x, u.Gamma, u.v, u.r, u.y, theta, vall, marginal.method) {
  
  u.Beta <- x[1:3];
  u.Alpha <- x[4:11];
  #u.Gamma <- x[12];
  rst <- 0;
  n.T <- ncol(u.y);
  nsub <- nrow(u.y);
  ncolv <- ncol(u.v);

  u.fullAlpha <- cbind(u.Alpha, u.Alpha, u.Alpha);
  u.Delta <- get.Deltat(n.T, u.fullAlpha, u.Beta, u.Gamma, vall, theta, marginal.method);
  vvec <- as.vector(t(u.v));
  vec.alpha <- as.vector(u.fullAlpha);
  yvec <- as.vector(t(u.y));
  yvec[which(is.na(yvec))] <- 0;
  delta <- as.vector(t(u.Delta));

  ### log.Like gives the log(likelihood) of the observed outcome 
  logLike <- .C("ttl_obsloglikelihood", as.integer(n.T), as.integer(nsub), as.double(vvec), as.integer(u.r), as.integer(ncolv), as.integer(yvec), as.double(delta), as.double(vec.alpha), as.double(u.Gamma), as.double(rst))[[10]];
  logLike;
}


### update 06/20/2015
runMLEcong <- function(inx) {
    SIMU.MARGIN  <- 3;
    V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
    V.INX.ALL     <- get.v.all(V.CATS);
    BETA <- c(0.5, 0.25, 0.3);  
    headername <- "~/Documents/Aux/hmc/Simu8Vfor20140213/size200/";

K.inx <- c(rep(1,30), rep(2,30), rep(3,30), rep(4,10), rep(20,30), rep(21,30), rep(22,30), rep(23,30), rep(24,30), rep(25,30), rep(26,30), rep(27,30), rep(28,30), rep(29,30));
I.inx <- c(rep(c(1:30),3), c(1:10), rep(c(1:30),10) );
KI.inx <- cbind(K.inx, I.inx);
    
   #qrst <- read.table(file=paste(headername, "Size200_ttl400_mlecong_parcl_rst.txt", sep=""), row.names=NULL);
   #qrst <- qrst[,-1];

   alltheta <- read.table(file=paste(headername, "Size200_ttl400_all_mean_dir_theta.txt",sep=""));
   alltheta <- as.matrix(alltheta);

   k <- KI.inx[inx,1];
   i <- KI.inx[inx,2];
   if(k<20) {
     vry <- read.table(file=paste(headername, "/data/Size200_K",k, "_data_",i, "_vry.txt",sep="")); 
   }else{
     vry <- read.table(file=paste(headername, "/add300/data/Size200_K",k, "_data",i,"_vry.txt",sep=""));
   }
   vry <- as.matrix(vry);
   v.obs <- vry[,1:8];
   r.obs <- vry[,9];
   y.obs <- vry[,10:12];
   theta <- alltheta[inx,];

   q.init <- rep(1,11);
   #rst <- optim(qrst[inx,1:11], obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=alltheta[inx,], vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
    rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=alltheta[inx,], vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1)); 
    c.value <- rst$value;
   diff <- 1;
   while(diff>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=alltheta[inx,], vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       diff <- rst$value-c.value;
       c.value <- rst$value;
   } 
   rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
   var <- (solve(-1*(rst$hessian)));
   se.beta <- sqrt(diag(var[1:3, 1:3]));   
   sig1 <- 0.90;
   low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
   up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
   inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
   sig2 <- 0.95;
   low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
   up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
   inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

   ### the standard error for p
   true.p <- get.prob.5(BETA);
   p <- get.prob.5(rst$par[1:3]);
   se.p <- get.se.p(rst$par[1:3], rst$hessian);
     #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));


   low.p <- p-qnorm((1+sig1)/2)*se.p;
   up.p <- p+qnorm((1+sig1)/2)*se.p;
   inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));

   low.p <- p-qnorm((1+sig2)/2)*se.p;
   up.p <- p+qnorm((1+sig2)/2)*se.p;
   inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));

   #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value), file=paste(headername, "0619/Size200_data", inx, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "0619/Size200_data", inx, "_mlecong_CI_rst.txt", sep=""));
}




runMLEcong20000 <- function(inx) {
    SIMU.MARGIN  <- 3;
    V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
    V.INX.ALL     <- get.v.all(V.CATS);
    BETA <- c(0.5, 0.25, 0.3);  
    headername <- "~/Documents/Aux/hmc/Simu8Vfor20140213/size20000/";
    K.inx <- c(rep(1,7), rep(2,7), rep(3,7), rep(4,5), rep(5,4), sapply(6:12, function(x) rep(x,10)), sapply(20:29, function(x) rep(x, 30)));
    I.inx <- c(rep(c(1:7),3), c(1:5), c(1:4), rep(c(1:10), 7), rep(c(1:30),10) );
    KI.inx <- cbind(K.inx, I.inx);

   #qrst <- read.table(file=paste(headername, "Size200_ttl400_mlecong_CG_rst.txt", sep=""));
   #qrst <- as.matrix(qrst);

   alltheta <- read.table(file=paste(headername, "Size20000_ttl400_all_mean_slice_theta.txt",sep=""));
   alltheta <- as.matrix(alltheta);
  
   k <- KI.inx[inx,1];
   i <- KI.inx[inx,2];
    if(k<8) {
     vry <- read.table(file=paste(headername, "data/K",k, "_data_",i, "_vry.txt",sep="")); 
   }else if(8<=k && k<20){
     vry <- read.table(file=paste(headername, "data/Size20000_K",k, "_data_",i,"_vry.txt",sep=""));
   }else{
     vry <- read.table(file=paste(headername, "add300/data/Size20000_K", k, "_data", i, "_vry.txt", sep=""));
   } 
   vry <- as.matrix(vry);
   v.obs <- vry[,1:8];
   r.obs <- vry[,9];
   y.obs <- vry[,10:12];
   theta <- alltheta[inx,];

   q.init <- rep(1,11);
   rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
   dif <- 1;
   c.value <- rst$value;
    
   while(dif>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- rst$value-c.value;
       c.value <- rst$value;
   }
   rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v.obs, u.r=r.obs, u.y=y.obs, theta=theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
   var <- (solve(-1*(rst$hessian)));
   se.beta <- sqrt(diag(var[1:3, 1:3]));   
   sig1 <- 0.90;
   low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
   up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
   inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
   sig2 <- 0.95;
   low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
   up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
   inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

   ### the standard error for p
   true.p <- get.prob.5(BETA);
   p <- get.prob.5(rst$par[1:3]);
   se.p <- get.se.p(rst$par[1:3], rst$hessian);
      #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));

   low.p <- p-qnorm((1+sig1)/2)*se.p;
   up.p <- p+qnorm((1+sig1)/2)*se.p;
   inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));

   low.p <- p-qnorm((1+sig2)/2)*se.p;
   up.p <- p+qnorm((1+sig2)/2)*se.p;
   inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));

   #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value), file=paste(headername, "0619/Size20000_data", inx, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "0619/Size20000_data", inx, "_mlecong_CI_rst.txt", sep=""));
}



#### summarize the bias and MSE for Beta and conditional probability
### for input of rst, the first 3 columns should be estimates for BETA on 03/12/2015
obtain.rst.summary <- function(rst, BETA) {
    bias.beta <- apply(rst[,1:3],2,mean)-BETA;
    sd.bias.beta <- apply(rst[,1:3], 2, sd)/sqrt(nrow(rst))
    a <- cbind(bias.beta, sd.bias.beta);
        
    true.prob <- get.prob.5(BETA);
    prob <- t(apply(rst[,1:3], 1, function(x) get.prob.5(x)));
    bias.prob <- apply(prob,2,mean)-true.prob;
    sd.bias.prob <- apply(prob, 2, sd)/sqrt(nrow(rst));
    b <- cbind(bias.prob, sd.bias.prob);
        
    mse.beta <- sapply(1:3, function(x) 10^3*mean((rst[,x]-BETA[x])^2));
    sd.mse.beta <- sapply(1:3, function(x) 10^3*sd((rst[,x]-BETA[x])^2)/sqrt(nrow(rst)));
    c <- cbind(mse.beta, sd.mse.beta)
        
    mse.p <- sapply(1:5, function(x) 10^3*mean((prob[,x]-true.prob[x])^2));
    sd.mse.p <- sapply(1:5, function(x) 10^3*sd((prob[,x]-true.prob[x])^2/sqrt(nrow(prob))));
    d <- cbind(mse.p, sd.mse.p);
    
    smry <- rbind(a, c, b, d);
    smry;
}





### update 06/22/2015
std.uncong <- function(inx, nImptData) {
    SIMU.MARGIN  <- 3;
    V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
    V.INX.ALL     <- get.v.all(V.CATS);
    BETA <- c(0.5, 0.25, 0.3);  
    headername <- "~/Documents/Aux/hmc/Simu8Vfor20140213/size200/";

K.inx <- c(rep(1,30), rep(2,30), rep(3,30), rep(4,10), rep(20,30), rep(21,30), rep(22,30), rep(23,30), rep(24,30), rep(25,30), rep(26,30), rep(27,30), rep(28,30), rep(29,30));
I.inx <- c(rep(c(1:30),3), c(1:10), rep(c(1:30),10) );
KI.inx <- cbind(K.inx, I.inx);
   
   alltheta <- read.table(file=paste(headername, "Size200_ttl400_all_mean_dir_theta.txt",sep=""));
   alltheta <- as.matrix(alltheta);

   k <- KI.inx[inx,1];
   i <- KI.inx[inx,2];
   if(k<20) {
     vry <- read.table(file=paste(headername, "/data/Size200_K",k, "_data_",i, "_vry.txt",sep="")); 
   }else{
     vry <- read.table(file=paste(headername, "/add300/data/Size200_K",k, "_data",i,"_vry.txt",sep=""));
   }
   vry <- as.matrix(vry);
   v.obs <- vry[,1:8];
   r.obs <- vry[,9];
   y.obs <- vry[,10:12];
   theta <- alltheta[inx,];

     uncong.impt.rst <- uncong.impt(v.obs, r.obs, y.obs);
      est.Alpha <- uncong.impt.rst$est.Alpha;
      est.Delta <- uncong.impt.rst$est.Delta;
      est.Gamma <- uncong.impt.rst$est.Gamma;
      rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v.obs, r.obs, BETA);
 
      betap <- rst[c(7:9, 13:17)];
      sd.betap <- sqrt(rst[c(10:12, 18:22)]);
      
      #write.table(c(rst, est.Alpha, est.Gamma), file=paste(headername, "0622/Size200_K", k, "_data", i, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));
      write.table(rst, file=paste(headername, "0630/Size200_K", k, "_data", i, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));

      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "0630/Size200_K", k, "_data", i, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));
    #write.table(c(inx90, inx95), file=paste(headername, "0622/Size200_data", inx, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""))
 }



### update 06/22/2015
std.uncong.20000 <- function(inx, nImptData) {

    SIMU.MARGIN  <- 3;
    V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
    V.INX.ALL     <- get.v.all(V.CATS);
    BETA <- c(0.5, 0.25, 0.3);  
    headername <- "~/Documents/Aux/hmc/Simu8Vfor20140213/size20000/";
    K.inx <- c(rep(1,7), rep(2,7), rep(3,7), rep(4,5), rep(5,4), sapply(6:12, function(x) rep(x,10)), sapply(20:29, function(x) rep(x, 30)));
    I.inx <- c(rep(c(1:7),3), c(1:5), c(1:4), rep(c(1:10), 7), rep(c(1:30),10) );
    KI.inx <- cbind(K.inx, I.inx);

   alltheta <- read.table(file=paste(headername, "Size20000_ttl400_all_mean_slice_theta.txt",sep=""));
   alltheta <- as.matrix(alltheta);
  
   k <- KI.inx[inx,1];
   i <- KI.inx[inx,2];
    if(k<8) {
     vry <- read.table(file=paste(headername, "data/K",k, "_data_",i, "_vry.txt",sep="")); 
   }else if(8<=k && k<20){
     vry <- read.table(file=paste(headername, "data/Size20000_K",k, "_data_",i,"_vry.txt",sep=""));
   }else{
     vry <- read.table(file=paste(headername, "add300/data/Size20000_K", k, "_data", i, "_vry.txt", sep=""));
   } 
   vry <- as.matrix(vry);
   v.obs <- vry[,1:8];
   r.obs <- vry[,9];
   y.obs <- vry[,10:12];
   theta <- alltheta[inx,];

     uncong.impt.rst <- uncong.impt(v.obs, r.obs, y.obs);
      est.Alpha <- uncong.impt.rst$est.Alpha;
      est.Delta <- uncong.impt.rst$est.Delta;
      est.Gamma <- uncong.impt.rst$est.Gamma;
      rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v.obs, r.obs, BETA);

      betap <- rst[c(7:9, 13:17)];
      sd.betap <- sqrt(rst[c(10:12, 18:22)]);

       write.table(rst, file=paste(headername, "0630/Size20000_K", k, "_data", i, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));
    #write.table(c(rst, est.Alpha, est.Gamma), file=paste(headername, "0622/Size20000_data", inx, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));

      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "0630/Size20000_K", k, "_data", i, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));
     #write.table(c(inx90, inx95), file=paste(headername, "0622/Size20000_data", inx, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""))
 }

  

### on 06/30
### obs.loglikelihood value is invariant of Gamma value
freq.y <- function(y) {
    nttl <- nrow(y);
    freqn <- numeric(10);

    freqn[1] <- length(which(y[,1]==0));
    freqn[2] <- nttl-freqn[1];

    freqn[3] <- length(intersect(which(y[,1]==0), which(y[,2]==0))); ### (y2=0, y1=0)
    freqn[4] <- length(intersect(which(y[,1]==1), which(y[,2]==0))); ### (y2=0, y1=1)
    freqn[5] <- length(intersect(which(y[,1]==0), which(y[,2]==1))); ### (y2=1, y1=0)
    freqn[6] <- length(intersect(which(y[,1]==1), which(y[,2]==1))); ### (y2=1, y1=1)

    freqn[7] <- length(intersect(which(y[,3]==0), which(y[,2]==0))); ### (y3=0, y2=0)
    freqn[8] <- length(intersect(which(y[,3]==0), which(y[,2]==1))); ### (y3=0, y2=1)
    freqn[9] <- length(intersect(which(y[,3]==1), which(y[,2]==0))); ### (y3=1, y2=0)
    freqn[10] <- length(intersect(which(y[,3]==1), which(y[,2]==1))); ### (y3=1, y2=1)
  
    freqn;
}

### on 06/30
### the loglikelihood for inference of beta in uncongenial analysis for 'complete' imputed dataset
impt.fullloglikelihood<- function(x, freq) {
  tmp <- numeric(10);

  tmp[1] <- 1-expit(x[1]);
  tmp[2] <- expit(x[1]);
  tmp[3] <- 1-expit(x[1]+x[2]);
  tmp[4] <- 1-expit(x[1]+x[2]+x[3]);
  tmp[5] <- 1-tmp[3];
  tmp[6] <- 1-tmp[4];
  tmp[7] <- 1-expit(x[1]+2*x[2]);
  tmp[8] <- 1-expit(x[1]+2*x[2]+x[3]);
  tmp[9] <- 1-tmp[7];
  tmp[10] <- 1-tmp[8];

  rst <- sum(sapply(1:10, function(m) log(tmp[m])*freq[m]));
  rst;    
}


### the mle of p based on mle of beta for congenial analysis
get.se.p <- function(beta, hes) {
 
     var.beta <- (-1)*solve(hes);
     
     eta <- numeric(5);
     eta[1] <- beta[1];
     eta[2] <- sum(beta[1:2]);
     eta[3] <- sum(beta[1:3]);
     eta[4] <- beta[1]+2*beta[2];
     eta[5] <- beta[1]+2*beta[2]+beta[3];    

     var.eta <- numeric(5);
     var.eta[1] <- var.beta[1,1];
     var.eta[2] <- var.beta[1,1]+var.beta[2,2]+2*var.beta[1,2];
     var.eta[3] <- var.beta[1,1]+var.beta[2,2]+var.beta[3,3]+2*var.beta[1,2]+2*var.beta[1,3]+2*var.beta[2,3];
     var.eta[4] <- var.beta[1,1]+4*var.beta[2,2]+4*var.beta[1,2];
     var.eta[5] <- var.beta[1,1]+4*var.beta[2,2]+var.beta[3,3]+2*var.beta[1,3]+4*var.beta[1,2]+4*var.beta[2,3];

     var.p <- numeric(5);
     var.p <- sapply(1:5, function(x) var.eta[x]*(exp(eta[x]/2)+exp(-1*eta[x]/2))^(-4));
     se.p <- sqrt(var.p);

     se.p
}





###   #### Prv.dirichlet for multilevel
###   xuan.pvdirichlet.multilevel.simple <- function(k, data.i, niter, nburn, max.h, simu.v, a.alpha, b.alpha, V.INX.ALL, V.CATS, true.theta){
###   
###     n.v <- ncol(simu.v);
###     n.obs <- nrow(simu.v);
###     pchi.theta <- 0;
###     
###     ## save for each iteration,
###    
###     ###
###     ind.v.length <- sapply(1:length(V.CATS), function(x) length(V.CATS[[x]]));
###     ttl.level <- sum(ind.v.length);
###     prev.length.ind.v <- c(0, sapply(1:length(V.CATS), function(x) sum(ind.v.length[1:x])));
###     rst.thai <- matrix(0, max.h, ttl.level);
###     
###     rst.u <- array(0,n.obs);
###     rst.z <- array(0,n.obs);
###     rst.v <- array(0,max.h);
###     rst.theta <- matrix(0, niter-1, nrow(V.INX.ALL));
###     n.hitmax <- 0;
###                       
###     ###########################  Begin initializing values  ######################################################
###     
###   
###     ### initialize alpha first based on paper  alpha ~ gamma(a.alpha, b.alpha)
###     #rst.alpha <- rgamma(1, a.alpha, b.alpha);
###     rst.alpha <- 0.5;
###     rst.v <- numeric(max.h);
###     rst.v <- rbeta(max.h, 1, rst.alpha);
###     rst.v[max.h] <- 1; # to make the sum of sticks equals to 1
###     tmp.omega <- cal.omega(max.h, rst.v);
###   
###     ### Based on rst.z, rst.v, rst.w (tmp.omega) generate rst.u
###     rst.u <- runif(n.obs, 0, 1/(max.h));
###     rst.z <- sapply(1:n.obs, function(x) sample(which(rst.u[x]<tmp.omega),1));
###     max.z <- max(rst.z);
###   
###        z.inx.list <- list();
###        for(i in 1:max.z) {
###           z.inx.list[[i]] <- which(rst.z==i);
###         }
###        cat("max.z=", max.z, "\n");
###     
###       for(h.iter in 1:max.z) {
###         tmp.thai <- NULL;
###         #cat("h=", h.iter, ", size=", length(z.inx.list[[h.iter]]), "\n");
###         #for(j in 1:n.v) { # for each categorical variable
###         #   prob.vj <- sapply(1:ind.v.length[j], function(x) length(which(simu.v[ z.inx.list[[h.iter]],j]==V.CATS[[j]][x])))/length(z.inx.list[[h.iter]]);
###         #   tmp.thai <- c(tmp.thai, rdirichlet(1, prob.vj));
###         #}
###         tmp.thai.lst <- sapply(1:n.v, function(j) sapply(1:ind.v.length[j], function(x) length(which(simu.v[ z.inx.list[[h.iter]],j]==V.CATS[[j]][x])))/length(z.inx.list[[h.iter]]));
###         for(i in 1:n.v) { tmp.thai <- c(tmp.thai, tmp.thai.lst[[i]]); }
###         rst.thai[h.iter,] <- tmp.thai;
###       }
###     
###       if(max.h>max.z){
###          for(h.iter in (max.z+1):max.h) {
###            tmp.thai.lst <- sapply(1:n.v, function(x) rdirichlet(1, c(rep(1,length(V.CATS[[x]])))));
###            tmp.thai <- NULL;
###            for(i in 1:n.v) { tmp.thai <- c(tmp.thai, tmp.thai.lst[[i]]); }
###            rst.thai[h.iter,] <- tmp.thai;
###            #s.length <- (max.h-max.z)*ttl.level;
###            #tmp.thai <- matrix(rdirichlet(s.length,c(1,1))[,1], (max.h-max.z), ttl.level);
###            #rst.thai[(max.z+1):(max.h),1:ttl.level] <- tmp.thai;
###           }
###        }  
###        
###          #### get the corresponding v.index number for each subjects
###        sub.indx <- apply(simu.v,1, function(x) vtosimvalue(x,V.CATS))+1;
###   
###     ############################ Finish Initialize ###################################################################
###      
###     ##### Begin Iteration
###     for(s in 2:niter) {
###      max.z <- max(rst.z);
###      ### step 1 update u
###       tmp.omega <- cal.omega(max.h, rst.v);
###       rst.u <- sapply(1:n.obs, function(x) runif(1, 0, tmp.omega[rst.z[x]]));
###       #cat("\n finishi (1)");
###   
###      ### step 2 update thai
###       # make a list such that the h^th element contains the index of subjects with z==h
###        z.inx.list <- list();
###        for(i in 1:max.z) {
###           z.inx.list[[i]] <- which(rst.z==i);
###         }
###      
###      ### step 3 update v  #################################################################  
###      if(max.z==1) {
###        lowb <- max(rst.u);
###        upb <- 1;
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[1] <- rst$samp.v;
###        failed <- rst$failed;
###        #cat("for h=1, bound=(", lowb, ",", upb, ")"); 
###        rst.v[2:(max.h)] <- rbeta((max.h-1), 1, rst.alpha);
###        ############################################################ 
###       }else if(max.z==2) {
###        # for v1
###          # if there is no subject with Z==1
###          if(length(z.inx.list[[1]])==0) {
###            lowb <- 0;
###          }else{
###           lowb <- max(rst.u[z.inx.list[[1]]]);
###          }
###          upb <- 1-max(rst.u[z.inx.list[[2]]])/rst.v[2];
###          cat("for h=1, bound=(", lowb, "," ,upb, ")"); 
###          rst <- samp.ind.v(lowb, upb, rst.alpha);
###          rst.v[1] <- rst$samp.v;
###          failed <- rst$failed;
###          
###        # for v2 
###        lowb <- min(1,max(rst.u[ z.inx.list[[2]] ])/(1-rst.v[1]));
###        upb <- 1;
###        #cat("for h=2, bound=(", lowb, "," ,upb, ")");     
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[2] <- rst$samp.v;
###        failed <- failed+rst$failed;
###        if(failed>=1) {   failed <- 1; }
###        ### for v3 and beyond 
###        rst.v[3:(max.h)] <- rbeta((max.h-2), 1, rst.alpha);
###       ############################################################ 
###       }else{   ### for (max.z>=3)
###       ###(1) for updating rst.v[1]
###          ## (1.1) low bound for V_h always come from groups with Z==h
###        if(length( z.inx.list[[1]] )==0) {
###          lowb <- 0; }else{
###          lowb <- max(rst.u[ z.inx.list[[1]] ]);
###        }  
###          ## (1.2) ## Find upbound for V_(h=1), coming from groups with Z= (h+1), (h+2), ... (max.h)
###         # for group of Z equal to 2
###        if(length( z.inx.list[[2]] )==0) {
###          max1 <- 0;
###        }else{
###          max1 <- max(rst.u[ z.inx.list[[2]] ])/rst.v[2];
###        }
###          upb1 <- 1-max1;       
###        ## (1.3)
###        #  for group of Z great than 3 
###        gt.h <- c(3:max.z);
###        length.gt.h <- max.z-3+1;
###        max.num <- numeric(length.gt.h);
###        max.num <- sapply(gt.h, function(x)  max(0, rst.u[ z.inx.list[[x]] ]));  # update on 07/29
###        denom.gt.h <- numeric(length.gt.h);
###        denom.gt.h[1] <- rst.v[3]*(1-rst.v[2]);  # from group z=3
###        # update on 7/29
###        if(length.gt.h>=2) { # from group z=4, 5, ...
###          denom.gt.h[2:length.gt.h] <- sapply(4:max.z, function(x) rst.v[x]*prod(sapply(2:(x-1), function(y) 1-rst.v[y])));
###        } 
###        max2 <- min(1,max(0,max.num/denom.gt.h));
###        upb2 <- 1-max2;
###        #cat("max2==", max2, "\n");
###        #upb <- 1-max(max1,max2);
###        upb <- min(upb1, upb2);
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[1] <- rst$samp.v;
###        failed <- rst$failed;
###        #cat("for h=1, bound=(", lowb, ",", upb, ")-max1=", max1, "-ma2=", max2, ", with rst.v[1]==", rst.v[1], "\n");     
###      ###  End for updating rst.v[1] with max.z>=3
###        
###        
###      ### (2) for updating rst.v[h], for h=2, ..., (max.z-1) with max.z>=3
###        for(h.iter in 2:(max.z-1)) {
###          # (2.0) for lower bound
###          if( length(z.inx.list[[h.iter]] )==0) {
###             lowb <- 0;
###           }else{
###            lowb <- max(rst.u[ z.inx.list[[h.iter]] ])/prod(sapply(1:(h.iter-1), function(x) 1-rst.v[x]));
###           }
###          # (2.1) find upper bound max1 for rst.v[h] from group with Z=h+1
###          #indxeq.hiterplus1 <- which(rst.z==(h.iter+1));
###          if(length( z.inx.list[[h.iter+1]]) ==0) {
###            max1 <- 0;
###          }else{
###             #cat("numer=", max(rst.u[ z.inx.list[[h.iter+1]] ]), "***denom==", (rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x])))), "\n");
###             max1 <- max(rst.u[ z.inx.list[[h.iter+1]] ])/(rst.v[(h.iter+1)]*prod(sapply(1:(h.iter-1), function(x) (1-rst.v[x]))));            }
###          upb1 <- 1-max1;
###          
###          #(2.2) find upper bound max2 for rst.v[h] from group with Z=h+2, h+3, ...
###          if(h.iter+1>=max.z) { # on the equality with h.iter=max.z-1, the upbound taken care by upb1
###             max2 <- 0;
###          }else {   
###             ### Begin Alternative approach
###             whichgt.h <- c((h.iter+2):(max.z));  ## {h*:h*>h.iter}
###             lengthgt.h <- length(whichgt.h);
###             
###             num.gt.h <- numeric(lengthgt.h);
###             num.gt.h <- sapply(whichgt.h, function(x) max(0, rst.u[ z.inx.list[[x]]  ]));
###             indx.num.gt0 <- which(num.gt.h>0);
###             num.gt0 <- num.gt.h[indx.num.gt0];
###             denom.gt0 <- numeric(length(indx.num.gt0));
###             
###             denom.gt0 <- sapply(whichgt.h[indx.num.gt0], function(x) rst.v[x]*prod(sapply(1:(h.iter-1), function(y) 1-rst.v[y]))*prod(sapply((h.iter+1):(x-1), function(z)  1-rst.v[z])));
###             
###             #at("num.gt.h===", num.gt0, "\n");
###             #cat("denom.gt.h===", denom.gt0, "\n");q
###             max2 <- min(1,max(0, num.gt0/denom.gt0));
###           }
###             upb2 <- 1-max2;
###             
###             upb <- min(upb1, upb2);
###             #cat("for h=", h.iter, ", bound=(", lowb, ",", upb, ")--max1=", max1, "--max2=", max2, ",");     
###             rst<- samp.ind.v(lowb, upb, rst.alpha);
###             rst.v[h.iter] <- rst$samp.v;
###             failed <- failed+rst$failed;
###        }
###      
###        #(3) #### for h.iter==max.z ( hence definitely some subjects fall in this group)
###        lowb <- max(rst.u[ z.inx.list[[max.z]] ])/prod(sapply(1:(max.z-1), function(x) 1-rst.v[x]));
###        upb <- 1;
###        #cat("for h=", max.z, ", bound=(", lowb, ",", upb, ")");    
###        rst <- samp.ind.v(lowb, upb, rst.alpha);
###        rst.v[max.z] <- rst$samp.v;
###        failed <- failed+rst$failed;
###        #(4) update rst.v[h] for h=max.z+1, ..., max.h
###        if((max.z+1)<=max.h){
###        rst.v[(max.z+1):(max.h)] <- rbeta((max.h-max.z),1,rst.alpha);
###        }
###        if(failed>=1) { failed <- 1; }
###      } ## finish for max.z>=3
###   
###      ### keept last v=1 for sum length=1
###      rst.v[max.h] <- 1;
###      tmp.omga <- cal.omega(max.h, rst.v);
###      rst.v[which(rst.v[1:(max.h-1)]==1)] <- 1-(10)^(-10); ### to exclude the case that middle part rst.v==1, hence (1-rst.v)=0 in the denominator on Sep17
###      #cat(",  (3)");
###      # after updating rst.v, need to update tmp.omega   
###      ###########  Finish step 3  for updating V   ############### ########### ########### ########### ########### ########### ########### ###########  
###   
###     ### Update on 06/08-- to save time, first calculate the p(v|h) and use it in updating Z
###            ind.h.prob <- matrix(0, max.h, nrow(V.INX.ALL));
###            for(tmp.h in 1:max.h) {
###               tmp.thai <- rst.thai[tmp.h,];
###                ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1, function(v) prod(tmp.thai[sapply(1:n.v, function(x) 1+length(which(V.CATS[[x]] < v[x]))+prev.length.ind.v[x])]));
###               #ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
###            }
###          #cat("dim of ind.h.prob==", dim(ind.h.prob), ", s of ind.h.prob==", apply(ind.h.prob, 1, sum), "\n");
###       ### step 4 
###        # update omega first, need to calculate for max.h
###        min.u <- min(rst.u);
###        #cat("min.u==", min.u, ",   ");
###        addlength <- sapply(1:(max.h), function(x) sum(tmp.omega[1:x]));
###        k.tilt <- min(max.h, min(which(addlength>=(1-min.u))));
###        cat(", with k.tilt==", k.tilt, "\n");
###   
###          for(i in 1:n.obs) {
###            tmp.A <- which(tmp.omega[1:k.tilt]>rst.u[i]);
###            #cat("for sub==", i, ", sample from(", tmp.A, ")\n");
###            size.A <- length(tmp.A);
###            # calculate the multinomial distribution probability
###            if(size.A==0) {
###             rst.z[i] <- rst.z[i];
###            }else if(size.A==1){
###             rst.z[i] <- tmp.A[1];
###            }else{
###             # cat("i==", i, ", tmp.A==", tmp.A, ", sub.indx[i]=", sub.indx[i], "***");   
###             prob.A <- numeric(size.A);
###             prob.A <- ind.h.prob[tmp.A, sub.indx[i]];
###               # cat("prob.A==", prob.A, "***");
###             prob.A <- prob.A/sum(prob.A);
###             rst.z[i] <- tmp.A[which(rmultinom(1,1,prob.A)==1)];
###            }
###         }
###       #cat(",  (4)");
###          #max.z <- max(rst.z);
###          #if(max.z==max.h) n.hitmax <- n.hitmax+1;
###          #cat("updated with z==", rst.z, " with max.z==", max(rst.z), "\n");
###        
###      #### step 4.5 -- only update in 0928 function()
###      unique.z <- unique(rst.z);
###      rank.z <- rank(unique.z);
###      update.rstz <- sapply(1:n.obs, function(x) rank.z[which(unique.z==rst.z[x])]);
###      update.tmp.omega <- c(tmp.omega[sort(unique.z)], tmp.omega[-unique.z]);
###      rst.z <- update.rstz;
###      tmp.omega <- update.tmp.omega;
###         ### based on each individual length, update rst.v
###         update.rstv <- numeric(max.h);
###         update.rstv[1] <- tmp.omega[1];
###         for(j in 2:max.h) {
###            update.rstv[j] <- tmp.omega[j]/prod(sapply(1:(j-1), function(y) (1-update.rstv[y])));
###         }   
###      rst.v <- update.rstv;
###      max.z <- max(rst.z);
###      
###     ### step 2 update thai
###        # make a list such that the h^th element contains the index of subjects with z==h
###        
###        z.inx.list <- list();
###        for(i in 1:max.z) {
###           z.inx.list[[i]] <- which(rst.z==i);
###         }
###   
###        for(h.iter in 1:max.z) {
###         tmp.thai <- NULL;
###           #cat("h=", h.iter, ", size=", length(z.inx.list[[h.iter]]), "\n");
###       
###         for(j in 1:n.v) { # for each categorical variable
###            # in our case each covariate having 2 level, hence only need to save the 1st probability
###            #n.vhjeq1 <- sum(simu.v[ z.inx.list[[h.iter]], j ]);
###            #prob.vj <- table(simu.v[ z.inx.list[[h.iter]],j])/length(z.inx.list[[h.iter]]);
###             prob.vj <- sapply(1:ind.v.length[j], function(x) length(which(simu.v[ z.inx.list[[h.iter]],j]==V.CATS[[j]][x])))/length(z.inx.list[[h.iter]]);
###            #prob.vj[which(prob.vj==0)] <- (10)^(-3);
###            #cat("sum p==", sum(prob.vj), ", dir h==", h.iter, ", j=", j, "***");
###            tmp.thai <- c(tmp.thai, rdirichlet(1, prob.vj));
###         }
###           #cat("tmp.thai==", tmp.thai, "\n");
###         rst.thai[h.iter,] <- tmp.thai;
###       }
###      
###      if(max.h>max.z){
###       for(h.iter in (max.z+1):max.h) {
###        #cat("max.dir h=", h.iter, "***");
###        tmp.thai.lst <- sapply(1:n.v, function(x) rdirichlet(1, c(rep(1,length(V.CATS[[x]])))));
###        tmp.thai <- NULL;
###        for(i in 1:n.v) { tmp.thai <- c(tmp.thai, tmp.thai.lst[[i]]); }
###        rst.thai[h.iter,] <- tmp.thai;
###        #s.length <- (max.h-max.z)*ttl.level;
###        #tmp.thai <- matrix(rdirichlet(s.length,c(1,1))[,1], (max.h-max.z), ttl.level);
###        #rst.thai[(max.z+1):(max.h),1:ttl.level] <- tmp.thai;
###        }
###      }
###      
###            for(tmp.h in 1:max.h) {
###               tmp.thai <- rst.thai[tmp.h,];
###               ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1, function(v) prod(tmp.thai[sapply(1:n.v, function(x) 1+length(which(V.CATS[[x]] < v[x]))+prev.length.ind.v[x])]));
###               #ind.h.prob[tmp.h,] <- apply(V.INX.ALL,1,function(x) prod(rst.thai[tmp.h,]*(1-x)+(1-rst.thai[tmp.h,])*x));
###            }
###             rst.theta[s-1,] <- t(ind.h.prob)%*%tmp.omega;
###         
###           ### step 5
###         # update omega first
###       tmp.shape <-  a.alpha+max.z;
###        if(max.z==max.h) {
###          inx.neq1 <- which(rst.v[1:(max.h)]<1);
###          tmp.rate <- Inf;
###          rst.alpha <- (10)^(-10);
###        }else{
###          inx.neq1 <- which(rst.v[1:(max.z)]<1);
###          tmp.rate <- b.alpha-log(1-sum(tmp.omega[1:max.z]));
###          rst.alpha <- rgamma(1, shape=tmp.shape, rate=tmp.rate);
###        }
###     
###      #cat("alpha (a==", tmp.shape, ", b==", tmp.rate, ")\n");
###       #cat(",  (5)");
###        ### calculate the pchi after each iteration
###       if(s>2) {
###         if(s<=(nburn+2)) {
###           m <- apply(rst.theta[1:(s-1),],2,mean);
###         }else{
###           m <- apply(rst.theta[nburn:(s-1),],2,mean);
###         }  
###       mean.pchi.theta <- probv.pchi(m, true.theta);
###       crnt.pchi.theta <- probv.pchi(rst.theta[s-1,], true.theta);
###         
###       cat(" at End iteration ", s, ", pchi===( ", mean.pchi.theta, ",", crnt.pchi.theta, "), ****max.z===", max(rst.z), ", alpha=", rst.alpha, ", for (K,i)==(", k, ",", data.i, ") with n.hitmax==", n.hitmax, "\n");
###      }
###      #cat(",  (6)");
###      #cat("updated with u==", rst.u, "\n");
###      #cat("updated with v==", rst.v[1:max.z], "\n");
###      #beyond.length <- c(beyond.length, 1-sum(tmp.omega[1:(max.z)]));
###      if(max.z==max.h) { n.hitmax <- n.hitmax+1; }
###      
###     } ### end of iteration
###      list(theta=rst.theta, n.hitmax=n.hitmax);
###   }
###   






##### Use NUTS for HMC of Congenial analysis on 07/07/2015
nuts <- function(headername, true.u, v,r, y, q, sigma, n.iter, epsilon, V.ALL, theta, Marginal.method, diag.sigma.p){

#x is the matrix of the covariates with first column=1 as intercept
n.obs = nrow(v); #get number of observations	
n.param = length(q);
n.T = ncol(y);

#post.samp used to store the parameters value at every iteration
post.samp <- matrix(NA, n.iter, n.param);
post.samp.U <- rep(0, n.iter);

### save the current value of position variable
current.q <- q;
true.q <- c(0.5, 0.25, 0.3, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.3);
### n.accept used to record how many times the proposed is accepted
n.accept <- 0;
acceptrate <- 0;
### calculate for each iteration and save the value of position variable for each iteration
for (i in 1:n.iter){
#cat("\n Begin Iter=", i, "#################################################################\n");
  if(i==1) {
    current.U <- newx.U2(current.q, 0, v, r, y, sigma, theta, V.ALL, Marginal.method);
    current.g <- grad(newx.U2, x=current.q, method="simple", u.Gamma=0, u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
   }

    samp <- nuts.onestep(v,r, y, current.q, current.U, current.g, sigma, n.param, epsilon, theta, diag.sigma.p, V.ALL, Marginal.method);

  ### update the position value, U, and gradient of U
     if(t(samp-current.q)%*%(samp-current.q) ==0) {
         current.update <- 0;
     }else{
         current.update <- 1;
     }
     ### 1 for accept the prosed; 0 for reject the propose
     current.q <- samp;
     current.U <- newx.U2(current.q, 0, v, r, y, sigma, theta, V.ALL, Marginal.method);
     current.g <- grad(newx.U2, x=current.q, method="simple", u.Gamma=0, u.v=v, u.r=r, u.y=y, u.sigma=sigma, theta=theta, vall=V.ALL, marginal.method=Marginal.method);
         
  ### save the result of iteration i
     post.samp[i,] <- current.q;
     post.samp.U[i] <- current.U;

     n.accept <- n.accept+1*current.update;
     acceptrate <- (n.accept)/(i);
   
        ### in case the acceptance rate is small, quit to save time 
        if((acceptrate < 0.1 && i>6)||(acceptrate < 0.1 && i>n.iter/2))
        {
          cat("Low acceptance rate=", acceptrate, "at iteration", i, ", quit\n");
          list(samp=post.samp, acceptrate=n.accept/n.iter)
          break
        }
         
         if (current.update==0)
         {
          cat("###Reject at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         } else {
          cat("###Accept at Iter", i, ", c.acceptrate= ", acceptrate, ",\n with c.q= ", round(current.q,digit=3), ".\n c.U=", current.U, "true.U=", true.u);
         }
       }
  rst <- cbind(post.samp, post.samp.U);
  rst;
}



 

### one step of NUTS for HMC in AMAR congenial analysis
nuts.onestep <- function(v.obs, r.obs, y.obs, gs.q, gs.U, gs.g, sigma, n.param, gs.et, gs.theta, diag.sigma.p, V.INX.ALL, SIMU.MARGIN){

   #### resample p    
   var.p <- matrix(0, n.param, n.param);
   diag(var.p) <- diag.sigma.p^2;
   mean.p <- rep(0,n.param);
   p <- as.numeric(rmvnorm(1, mean.p, var.p));

   ### for AMAR case the Delta value changes with gamma, making obs likelihood invariant of gamma, hence make it 0
   gs.gamma <- 0;
   
   ### save the newly generated momentum variable as current.p
   ### then the current value of K(p)=(p^T%*%p)/2 as current.K based on current.p
   ### then the current value of H=K(p)+U(q), save the value as current.H
   current.q = gs.q;
   current.U = gs.U;
   current.p = p;
   current.K = sum(current.p^2/(2*diag.sigma.p^2));
   current.H = current.U + current.K;
   current.exp = exp(-current.H);
   mu <- runif(1, 0, current.exp);
   
   q.neg <- gs.q;
   q.pos <- gs.q;
   p.neg <- p;
   p.pos <- p;
   j <- 0;
   proposed.q <- current.q;
   n <- 1;
   s <- 1;
   
   while(s==1) {
     u.v <- sample(c(-1,1), size=1, prob=c(0.5, 0.5)); ### choose a direction vj
     if(u.v== -1){
        rst <- recurse(q.neg, p.neg, mu, u.v, j, gs.et, v.obs, r.obs, y.obs, gs.gamma, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN);
        q.neg <- rst$q.neg;
        p.neg <- rst$p.neg;
        u.q <- rst$u.q;
        u.n <- rst$u.n;
        u.s <- rst$u.s;
     }else{
        rst <- recurse(q.pos, p.pos, mu, u.v, j, gs.et, v.obs, r.obs, y.obs, gs.gamma, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN);
        q.pos <- rst$q.pos;
        p.pos <- rst$p.pos;
        u.q <- rst$u.q;
        u.n <- rst$u.n;
        u.s <- rst$u.s;
     }
   
     if(u.s==1) {
        prob <- min(1, u.n/n);
        cat("prob==", prob, "****");
        if(runif(1,0,1) <= prob) {  proposed.q <- u.q; }
     }###end if(u.s==1)
   
     n <- n+u.n;
     s <- u.s*(I(t(q.pos-q.neg)%*%p.neg >= 0)*I(t(q.pos-q.neg)%*%p.pos >= 0));
     j <- j+1;
   }### end while    
   cat("n==", n, ", j==", j, "\n");
   proposed.q;

}



recurse <- function(q, p, mu, u.v, j, gs.et, v.obs, r.obs, y.obs, gs.gamma, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN){
   delta.max <- 1000;
   if(j == 0){ ### base case take one leapfrog step in the direction u.v
       gs.g = grad(newx.U2, x=q, method="simple", u.Gamma=gs.gamma, u.v=v.obs, u.r=r.obs, u.y=y.obs, u.sigma=sigma, theta=gs.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN);
       u.p <- p+u.v*gs.et*gs.g/2;
       u.q <- q+u.v*gs.et*u.p;
       gs.g = grad(newx.U2, x=u.q, method="simple", u.Gamma=gs.gamma, u.v=v.obs, u.r=r.obs, u.y=y.obs, u.sigma=sigma, theta=gs.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN);
       u.p <- u.p+u.v*gs.et*gs.g/2;
       proposed.U = newx.U2(u.q, gs.gamma, v.obs, r.obs, y.obs, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN);
       proposed.K = sum(u.p^2/(2*diag.sigma.p^2));
       proposed.H = proposed.U+proposed.K;
       proposed.exp <- exp(-proposed.H);
       u.n <- I(mu <= proposed.exp);
       u.s <- I(log(proposed.exp) > log(mu)-delta.max);
       cat("j==", j, "**** u.n==", u.n, "\n");
       q.neg <- u.q;
       p.neg <- u.p;
       q.pos <- u.q;
       p.pos <- u.p;
       #return(q.neg=u.q, p.neg=u.p, q.pos=u.q, p.pos=u.p, u.q=u.q, u.n=u.n, u.s=u.s); 
    }else{  ### Recursion, implicitly build the left and right subtrees
       rst <- recurse(q, p, mu, u.v, j-1, gs.et, v.obs, r.obs, y.obs, gs.gamma, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN);
       u.s <- rst$u.s;
       q.neg <- rst$q.neg;
       p.neg <- rst$p.neg;
       q.pos <- rst$q.pos;
       p.pos <- rst$p.pos;
       u.q <- rst$u.q;
       u.n <- rst$u.n;
       cat("j==", j, "**** u.n==", u.n, "\n");
       if(u.s== 1) {
           if(u.v== -1) {
                rst <-  recurse(q.neg, p.neg, mu, u.v, j-1, gs.et, v.obs, r.obs, y.obs, gs.gamma, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN);
                q.neg <- rst$q.neg;
                p.neg <- rst$p.neg;
                q.prime2 <- rst$u.q;
                n.prime2 <- rst$u.n;
                s.prime2 <- rst$u.s;
           }else{
                rst <-  recurse(q.pos, p.pos, mu, u.v, j-1, gs.et, v.obs, r.obs, y.obs, gs.gamma, sigma, gs.theta, V.INX.ALL, SIMU.MARGIN);
                q.pos <- rst$q.pos;
                p.pos <- rst$p.pos;
                q.prime2 <- rst$u.q;
                n.prime2 <- rst$u.n;
                s.prime2 <- rst$u.s;
           }
           prob <- n.prime2/(u.n+n.prime2);
           cat("u.n==", u.n, "*** n2==", n.prime2, "***prob==", prob, "****");
           if(runif(1,0,1) <= prob) {  u.q <- q.prime2; }
           u.s <- s.prime2*I(t(q.pos-q.neg)%*%p.neg >=0)*I(t(q.pos-q.neg)%*%p.pos >=0);
           u.n <- u.n+n.prime2;
       }### end if
    }
   list(q.neg=q.neg, p.neg=p.neg, q.pos=q.pos, p.pos=p.pos, u.q=u.q, u.n=u.n, u.s=u.s);  
}









#### on 07/14/2015
#### AMAR case from (simulate data, estimate theta, run uncong, and run mle cong)
runAMARfull <- function(headername, size, data.i) {
   source("~/Documents/Aux/hmc/test_Xuan_AMAR_noGamma.R");
   dyn.load("~/Documents/Aux/hmc/Xuan_AMAR_bictool_backup_v2.so");

   SIMU.MARGIN  <- 3;
   V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
   V.INX.ALL     <- get.v.all(V.CATS);
   true.theta <- read.table(file="~/Documents/Aux/hmc/K1_truetheta.txt");
   true.theta <- true.theta[[1]];
   NTIME <- 3; #not including t=0 baseline    
   BETA  <- c(0.5, 0.25, 0.3);
   ALPHA <- matrix(c(0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7,0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7),8, NTIME);
   GAMMA <- 0.3;
   phi.v <- c(0.6, 0.7, 0.5, 0.4, 0, 0, 0, 0);
   PHI <- c(-3.5, phi.v, 0.5, 0);
   DELTA.ALL<- get.Deltat(NTIME, ALPHA, BETA, GAMMA, V.INX.ALL, true.theta, SIMU.MARGIN);

   
    v     <- simu.V(size, V.INX.ALL, true.theta);
    y.full     <- simu.Y(DELTA.ALL, ALPHA, GAMMA, v);
    r.ind <- simu.R(y.full, v, PHI);
    r <- NTIME-sapply(1:size, function(x) sum(r.ind[x,]));
    y.obs  <- set.missing.y(y.full, r.ind);
    emp.v.n <- get.v.inx(v, V.CATS, counts=TRUE);
    write.table(cbind(v, r, y.obs), file=paste(headername, "alldata/Size", size, "_data", data.i, "_vry.txt", sep=""));

    probv.NITER <- 5000;
    probv.NBURN <- 1000;
    eta.init <- 1;
    eta.lowbound <- 0;
    eta.upbound <- 10;
    prvslice.rst <- prv.slice(V.INX.ALL, emp.v.n, probv.NITER, probv.NBURN, eta.init, eta.lowbound, eta.upbound, true.theta);
    write.table(prvslice.rst$mean.theta, file=paste(headername, "theta/Size", size, "_data", data.i, "_mean_slice_theta.txt",sep=""));         
    est.theta <- prvslice.rst$mean.theta;

    uncong.impt.rst <- uncong.impt(v, r, y.obs);
    est.Alpha <- uncong.impt.rst$est.Alpha;
    est.Delta <- uncong.impt.rst$est.Delta;
    est.Gamma <- uncong.impt.rst$est.Gamma;
    nImptData <- 1000;
    rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v, r, BETA);
    #rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   
    betap <- rst[c(7:9, 13:17)];
    sd.betap <- sqrt(rst[c(10:12, 18:22)]);
    write.table(rst, file=paste(headername, "alluncong/Size", size, "_data", data.i, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));

      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "alluncong/Size", size, "_data", data.i, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));

     
       q.init <- rep(1,11);
       rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- 1;
       c.value <- rst$value;
    
       while(dif>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- rst$value-c.value;
       c.value <- rst$value;
   }
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
    var <- (solve(-1*(rst$hessian)));
    se.beta <- sqrt(diag(var[1:3, 1:3]));   
    sig1 <- 0.90;
    low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
    inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
    sig2 <- 0.95;
    low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
    inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

    ### the standard error for p
    true.p <- get.prob.5(BETA);
    p <- get.prob.5(rst$par[1:3]);
    se.p <- get.se.p(rst$par[1:3], rst$hessian);
      #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));
    low.p <- p-qnorm((1+sig1)/2)*se.p;
    up.p <- p+qnorm((1+sig1)/2)*se.p;
    inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    low.p <- p-qnorm((1+sig2)/2)*se.p;
    up.p <- p+qnorm((1+sig2)/2)*se.p;
    inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value, se.beta, p, se.p), file=paste(headername, "allcong/Size", size, "_data", data.i, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "allcong/Size", size, "_data", data.i, "_mlecong_CI_rst.txt", sep=""));
}  
   







###### update on 07/15
###### originally the missing y setting is wrong, resimulate data only keep original v to save computation time of p(v) estimation
revise.AMAR <- function(headername, size, data.i) {

   source("~/Documents/Aux/hmc/test_Xuan_AMAR_noGamma.R");
   dyn.load("~/Documents/Aux/hmc/Xuan_AMAR_bictool_backup_v2.so");

   SIMU.MARGIN  <- 3;
   V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
   V.INX.ALL     <- get.v.all(V.CATS);
   true.theta <- read.table(file="~/Documents/Aux/hmc/K1_truetheta.txt");
   true.theta <- true.theta[[1]];
   NTIME <- 3; #not including t=0 baseline    
   BETA  <- c(0.5, 0.25, 0.3);
   ALPHA <- matrix(c(0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7,0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7),8, NTIME);
   GAMMA <- 0.3;
   phi.v <- c(0.6, 0.7, 0.5, 0.4, 0, 0, 0, 0);
   PHI <- c(-3.5, phi.v, 0.5, 0);
   DELTA.ALL<- get.Deltat(NTIME, ALPHA, BETA, GAMMA, V.INX.ALL, true.theta, SIMU.MARGIN);

   vry <- read.table(file=paste(headername, "data/Size", size, "_data", data.i, "_vry.txt", sep=""));
   vry <- as.matrix(vry);   
   v     <- vry[,1:8];
   y.full     <- simu.Y(DELTA.ALL, ALPHA, GAMMA, v);
   r.ind <- simu.R(y.full, v, PHI);
   r <- NTIME-sapply(1:size, function(x) sum(r.ind[x,]));
   y.obs  <- set.missing.y(y.full, r.ind);
   emp.v.n <- get.v.inx(v, V.CATS, counts=TRUE);
   write.table(cbind(v, r, y.obs), file=paste(headername, "data/Size", size, "_data", data.i, "_vry.txt", sep=""));

    est.theta <- read.table(file=paste(headername, "theta/Size", size, "_data", data.i, "_mean_slice_theta.txt",sep="")); 
    est.theta <- est.theta[[1]];        

    uncong.impt.rst <- uncong.impt(v, r, y.obs);
    est.Alpha <- uncong.impt.rst$est.Alpha;
    est.Delta <- uncong.impt.rst$est.Delta;
    est.Gamma <- uncong.impt.rst$est.Gamma;
    nImptData <- 1000;
    rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v, r, BETA);
    #rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   
    betap <- rst[c(7:9, 13:17)];
    sd.betap <- sqrt(rst[c(10:12, 18:22)]);
    write.table(rst, file=paste(headername, "uncong/Size", size, "_data", data.i, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));
    
      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "uncong/Size", size, "_data", data.i, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));

     
       q.init <- rep(1,11);
       rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- 1;
       c.value <- rst$value;
    
       while(dif>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- rst$value-c.value;
       c.value <- rst$value;
   }
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
    var <- (solve(-1*(rst$hessian)));
    se.beta <- sqrt(diag(var[1:3, 1:3]));   
    sig1 <- 0.90;
    low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
    inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
    sig2 <- 0.95;
    low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
    inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

    ### the standard error for p
    true.p <- get.prob.5(BETA);
    p <- get.prob.5(rst$par[1:3]);
    se.p <- get.se.p(rst$par[1:3], rst$hessian);
      #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));
    low.p <- p-qnorm((1+sig1)/2)*se.p;
    up.p <- p+qnorm((1+sig1)/2)*se.p;
    inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    low.p <- p-qnorm((1+sig2)/2)*se.p;
    up.p <- p+qnorm((1+sig2)/2)*se.p;
    inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value, se.beta, p, se.p), file=paste(headername, "cong/Size", size, "_data", data.i, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "cong/Size", size, "_data", data.i, "_mlecong_CI_rst.txt", sep=""));


}






#### write on 07/16/2015 
#### For AMAR size200 with 8V --- reuse the original v (after finding bug in simulating y.obs) 
revise.AMAR.size200<- function(headername, size, inx) {

   source("~/Documents/Aux/hmc/test_Xuan_AMAR_noGamma.R");
   dyn.load("~/Documents/Aux/hmc/Xuan_AMAR_bictool_backup_v2.so");

   SIMU.MARGIN  <- 3;
   V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
   V.INX.ALL     <- get.v.all(V.CATS);
   true.theta <- read.table(file="~/Documents/Aux/hmc/K1_truetheta.txt");
   true.theta <- true.theta[[1]];
   NTIME <- 3; #not including t=0 baseline    
   BETA  <- c(0.5, 0.25, 0.3);
   ALPHA <- matrix(c(0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7,0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7),8, NTIME);
   GAMMA <- 0.3;
   phi.v <- c(0.6, 0.7, 0.5, 0.4, 0, 0, 0, 0);
   PHI <- c(-3.5, phi.v, 0.5, 0);
   DELTA.ALL<- get.Deltat(NTIME, ALPHA, BETA, GAMMA, V.INX.ALL, true.theta, SIMU.MARGIN);

  KIindex <- read.table(file="~/Documents/Aux/hmc/Simu8Vfor20140213/size200/KIindex_size200.txt");
  KIindex <- as.matrix(KIindex);
  k <- KIindex[inx,1];
  i <- KIindex[inx,2];
  
   if(k < 10) {
      vry <- read.table(file=paste(headername, "data/Size200_K",k,"_data_",i, "_vry.txt",sep=""));
   }else{
      vry <- read.table(file=paste(headername, "add300/data/Size200_K", k, "_data", i, "_vry.txt", sep=""));
   }
   vry <- as.matrix(vry);   
   v     <- vry[,1:8];
   y.full     <- simu.Y(DELTA.ALL, ALPHA, GAMMA, v);
   r.ind <- simu.R(y.full, v, PHI);
   r <- NTIME-sapply(1:size, function(x) sum(r.ind[x,]));
   y.obs  <- set.missing.y(y.full, r.ind);
   #emp.v.n <- get.v.inx(v, V.CATS, counts=TRUE);
   write.table(cbind(v, r, y.obs), file=paste(headername, "alldata/Size", size, "_data", inx, "_vry.txt", sep=""));
###   
   ###  vry <- read.table(file=paste(headername, "alldata/Size200_data", inx, "_vry.txt", sep=""));
   ###  vry <- as.matrix(vry);   
   ###  v     <- vry[,1:8];
   ###  r <- vry[,9];
   ###  y.obs <- vry[,10:12];
   ###
   
   all.theta <- read.table(file=paste(headername, "Size200_ttl400_all_mean_dan_theta.txt",sep=""));
   all.theta <- as.matrix(all.theta);
   est.theta <- all.theta[inx,];        

    uncong.impt.rst <- uncong.impt(v, r, y.obs);
    est.Alpha <- uncong.impt.rst$est.Alpha;
    est.Delta <- uncong.impt.rst$est.Delta;
    est.Gamma <- uncong.impt.rst$est.Gamma;
    nImptData <- 1000;
    rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v, r, BETA);
    #rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   
    betap <- rst[c(7:9, 13:17)];
    sd.betap <- sqrt(rst[c(10:12, 18:22)]);
    write.table(rst, file=paste(headername, "alluncong/Size", size, "_data", inx, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));
    
      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "alluncong/Size", size, "_data", inx, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));

     
       q.init <- rep(1,11);
       rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- 1;
       c.value <- rst$value;
    
       while(dif>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- rst$value-c.value;
       c.value <- rst$value;
   }
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
    var <- (solve(-1*(rst$hessian)));
    se.beta <- sqrt(diag(var[1:3, 1:3]));   
    sig1 <- 0.90;
    low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
    inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
    sig2 <- 0.95;
    low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
    inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

    ### the standard error for p
    true.p <- get.prob.5(BETA);
    p <- get.prob.5(rst$par[1:3]);
    se.p <- get.se.p(rst$par[1:3], rst$hessian);
      #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));
    low.p <- p-qnorm((1+sig1)/2)*se.p;
    up.p <- p+qnorm((1+sig1)/2)*se.p;
    inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    low.p <- p-qnorm((1+sig2)/2)*se.p;
    up.p <- p+qnorm((1+sig2)/2)*se.p;
    inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value, se.beta, p, se.p), file=paste(headername, "allcong/Size", size, "_data", inx, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "allcong/Size", size, "_data", inx, "_mlecong_CI_rst.txt", sep=""));


}







#### write on 07/16/2015 
#### For AMAR size200 with 8V --- from beginning to end
all.AMAR.size200<- function(headername, size, inx) {

   source("~/Documents/Aux/hmc/test_Xuan_AMAR_noGamma.R");
   dyn.load("~/Documents/Aux/hmc/Xuan_AMAR_bictool_backup_v2.so");

   SIMU.MARGIN  <- 3;
   V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
   V.INX.ALL     <- get.v.all(V.CATS);
   true.theta <- read.table(file="~/Documents/Aux/hmc/K1_truetheta.txt");
   true.theta <- true.theta[[1]];
   NTIME <- 3; #not including t=0 baseline    
   BETA  <- c(0.5, 0.25, 0.3);
   ALPHA <- matrix(c(0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7,0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7),8, NTIME);
   GAMMA <- 0.3;
   phi.v <- c(0.6, 0.7, 0.5, 0.4, 0, 0, 0, 0);
   PHI <- c(-3.5, phi.v, 0.5, 0);
   DELTA.ALL<- get.Deltat(NTIME, ALPHA, BETA, GAMMA, V.INX.ALL, true.theta, SIMU.MARGIN);


   v  <- simu.V(size, V.INX.ALL, true.theta);
   y.full     <- simu.Y(DELTA.ALL, ALPHA, GAMMA, v);
   r.ind <- simu.R(y.full, v, PHI);
   r <- NTIME-sapply(1:size, function(x) sum(r.ind[x,]));
   y.obs  <- set.missing.y(y.full, r.ind);
   #emp.v.n <- get.v.inx(v, V.CATS, counts=TRUE);
   write.table(cbind(v, r, y.obs), file=paste(headername, "alldata/Size", size, "_data", inx, "_vry.txt", sep=""));

   probv.NITER <- 5000;
   probv.NBURN <- 1000;
   max.h <- 80;
   pv.dir <- xuan.pvdirichlet.fast.1127(1,1, probv.NITER, probv.NBURN, max.h, v, a.alpha, b.alpha, V.INX.ALL, V.CATS, true.theta);
   write.table(pv.dir$mean.theta, file=paste(headername, "alltheta/Size", size, "_data", inx, "_mean_dir_theta.txt",sep=""));
   est.theta <- pv.dir$mean.theta;
   
   #all.theta <- read.table(file=paste(headername, "Size200_ttl400_all_mean_dan_theta.txt",sep=""));
   #all.theta <- as.matrix(all.theta);
   #est.theta <- all.theta[inx,];        

    uncong.impt.rst <- uncong.impt(v, r, y.obs);
    est.Alpha <- uncong.impt.rst$est.Alpha;
    est.Delta <- uncong.impt.rst$est.Delta;
    est.Gamma <- uncong.impt.rst$est.Gamma;
    nImptData <- 1000;
    rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v, r, BETA);
    #rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   
    betap <- rst[c(7:9, 13:17)];
    sd.betap <- sqrt(rst[c(10:12, 18:22)]);
    write.table(c(rst, betap, sd.betap), file=paste(headername, "alluncong/Size", size, "_data", inx, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));
    
      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "alluncong/Size", size, "_data", inx, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));

     
       q.init <- rep(1,11);
       rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- 1;
       c.value <- rst$value;
    
       while(dif>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- rst$value-c.value;
       c.value <- rst$value;
   }
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
    var <- (solve(-1*(rst$hessian)));
    se.beta <- sqrt(diag(var[1:3, 1:3]));   
    sig1 <- 0.90;
    low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
    inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
    sig2 <- 0.95;
    low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
    inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

    ### the standard error for p
    true.p <- get.prob.5(BETA);
    p <- get.prob.5(rst$par[1:3]);
    se.p <- get.se.p(rst$par[1:3], rst$hessian);
      #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));
    low.p <- p-qnorm((1+sig1)/2)*se.p;
    up.p <- p+qnorm((1+sig1)/2)*se.p;
    inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    low.p <- p-qnorm((1+sig2)/2)*se.p;
    up.p <- p+qnorm((1+sig2)/2)*se.p;
    inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value, se.beta, p, se.p), file=paste(headername, "allcong/Size", size, "_data", inx, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "allcong/Size", size, "_data", inx, "_mlecong_CI_rst.txt", sep=""));


}






revise.AMAR.size20000<- function(headername, size, inx) {

   source("~/Documents/Aux/hmc/test_Xuan_AMAR_noGamma.R");
   dyn.load("~/Documents/Aux/hmc/Xuan_AMAR_bictool_backup_v2.so");

   SIMU.MARGIN  <- 3;
   V.CATS      <- list(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1), c(0,1));
   V.INX.ALL     <- get.v.all(V.CATS);
   true.theta <- read.table(file="~/Documents/Aux/hmc/K1_truetheta.txt");
   true.theta <- true.theta[[1]];
   NTIME <- 3; #not including t=0 baseline    
   BETA  <- c(0.5, 0.25, 0.3);
   ALPHA <- matrix(c(0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7, 0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7,0.4, 0.3, 0.5, 0.9,0.8, 0.6, 0.3, 0.7),8, NTIME);
   GAMMA <- 0.3;
   phi.v <- c(0.6, 0.7, 0.5, 0.4, 0, 0, 0, 0);
   PHI <- c(-3.5, phi.v, 0.5, 0);
   DELTA.ALL<- get.Deltat(NTIME, ALPHA, BETA, GAMMA, V.INX.ALL, true.theta, SIMU.MARGIN);

   KIindex <- read.table(file="~/Documents/Aux/hmc/Simu8Vfor20140213/size20000/KIindex_size20000.txt");
   KIindex <- as.matrix(KIindex);
   k <- KIindex[inx,1];
   i <- KIindex[inx,2];

   if(k <= 7) {
      vry <- read.table(file=paste(headername, "data/K",k,"_data_",i, "_vry.txt",sep=""));
   }else if(k >7 && k <=12) {
      vry <- read.table(file=paste(headername, "data/Size20000_K", k, "_data_", i, "_vry.txt", sep=""));
   }else{
      vry <- read.table(file=paste(headername, "add300/data/Size20000_K", k, "_data", i, "_vry.txt", sep="")); 
   }
   vry <- as.matrix(vry);   
   v     <- vry[,1:8];
   y.full     <- simu.Y(DELTA.ALL, ALPHA, GAMMA, v);
   r.ind <- simu.R(y.full, v, PHI);
   r <- NTIME-sapply(1:size, function(x) sum(r.ind[x,]));
   y.obs  <- set.missing.y(y.full, r.ind);
   #emp.v.n <- get.v.inx(v, V.CATS, counts=TRUE);
   write.table(cbind(v, r, y.obs), file=paste(headername, "alldata/Size", size, "_data", inx, "_vry.txt", sep=""));

   all.theta <- read.table(file=paste(headername, "Size20000_ttl400_all_mean_slice_theta.txt",sep=""));
   all.theta <- as.matrix(all.theta);
   est.theta <- all.theta[inx,];        

    uncong.impt.rst <- uncong.impt(v, r, y.obs);
    est.Alpha <- uncong.impt.rst$est.Alpha;
    est.Delta <- uncong.impt.rst$est.Delta;
    est.Gamma <- uncong.impt.rst$est.Gamma;
    nImptData <- 1000;
    rst <- uncong.trans.0630(nImptData, est.theta, est.Alpha, est.Gamma, est.Delta, y.obs, v, r, BETA);
    #rst <- c(est.Beta, est.Beta.var, est.Beta.optim, est.Beta.optim.var, est.prob, prob.var);
   
    betap <- rst[c(7:9, 13:17)];
    sd.betap <- sqrt(rst[c(10:12, 18:22)]);
    write.table(rst, file=paste(headername, "alluncong/Size", size, "_data", inx, "_stduncong_par_nimpt", nImptData, ".txt", sep=""));
    
      sig1 <- 0.90;
      sig2 <- 0.95;
      low.b <- betap-qnorm((1+sig1)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig1)/2)*sd.betap;
      par <- c(BETA, get.prob.5(BETA));
      inx90 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
     
      low.b <- betap-qnorm((1+sig2)/2)*sd.betap;
      up.b <- betap+qnorm((1+sig2)/2)*sd.betap;
      inx95 <- sapply(1:length(par), function(x) (low.b[x]<=par[x])&&(up.b[x]>=par[x]));   
      write.table(c(inx90, inx95), file=paste(headername, "alluncong/Size", size, "_data", inx, "_stduncong_CI_nimpt", nImptData, ".txt", sep=""));

     
       q.init <- rep(1,11);
       rst <- optim(q.init, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- 1;
       c.value <- rst$value;
    
       while(dif>(10)^(-2)) {
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1));
       dif <- rst$value-c.value;
       c.value <- rst$value;
   }
       rst <- optim(rst$par, obs.loglikelihood, hessian=TRUE, u.Gamma=0, u.v=v, u.r=r, u.y=y.obs, theta=est.theta, vall=V.INX.ALL, marginal.method=SIMU.MARGIN, control=list(maxit=10000, fnscale= -1), method="CG");

    ### the standard error for beta
    var <- (solve(-1*(rst$hessian)));
    se.beta <- sqrt(diag(var[1:3, 1:3]));   
    sig1 <- 0.90;
    low.b <- rst$par[1:3]- qnorm((1+sig1)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig1)/2)*se.beta;
    inx90.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));   
    sig2 <- 0.95;
    low.b <- rst$par[1:3]- qnorm((1+sig2)/2)*se.beta;
    up.b <- rst$par[1:3]+ qnorm((1+sig2)/2)*se.beta;
    inx95.beta <- sapply(1:3, function(x) (low.b[x]<=BETA[x])&&(up.b[x]>=BETA[x]));

    ### the standard error for p
    true.p <- get.prob.5(BETA);
    p <- get.prob.5(rst$par[1:3]);
    se.p <- get.se.p(rst$par[1:3], rst$hessian);
      #### same se.p with the following method as used in uncongenial by Delta method
      ### grad.p <- jacobian(get.prob.5, rst$par);
      ### se.p <- sqrt(diag(grad.p%*%var%*%t(grad.p)));
    low.p <- p-qnorm((1+sig1)/2)*se.p;
    up.p <- p+qnorm((1+sig1)/2)*se.p;
    inx90.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    low.p <- p-qnorm((1+sig2)/2)*se.p;
    up.p <- p+qnorm((1+sig2)/2)*se.p;
    inx95.p <- sapply(1:5, function(x) (low.p[x]<=true.p[x])&&(up.p[x]>=true.p[x]));
    #CI.ind[inx,] <- c(inx90.beta, inx95.beta, inx90.p, inx95.p);

   write.table(c(rst$par, rst$convergence, rst$value, se.beta, p, se.p), file=paste(headername, "allcong/Size", size, "_data", inx, "_mlecong_CG_parcl_rst.txt", sep=""));
   write.table(c(inx90.beta, inx95.beta, inx90.p, inx95.p), file=paste(headername, "allcong/Size", size, "_data", inx, "_mlecong_CI_rst.txt", sep=""));


}



##### update on 07/16/2015 for summarizing the result of AMAR
AMAR.summary <- function(headername, size, ttlndata, BETA) {
    rst <- NULL;
    for(data.i in 1:ttlndata) {
    tmp <-read.table(file=paste(headername, "alluncong/Size", size, "_data", data.i, "_stduncong_CI_nimpt1000.txt",sep=""));
    tmp <- tmp[[1]];
    rst <- rbind(rst, tmp);
    }
    uncong.ci <- rst;
    uncong.mean.ci <- apply(uncong.ci,2,mean)
    uncong.ci <- cbind(uncong.mean.ci[1:8], uncong.mean.ci[9:16])  ### 90%, and 95% rate for uncong

    rst <- NULL;
    for(data.i in 1:ttlndata) {
        tmp <-read.table(file=paste(headername, "allcong/Size", size, "_data", data.i, "_mlecong_CI_rst.txt",sep=""));
        tmp <- tmp[[1]];di
        rst <- rbind(rst, tmp);
    }
    cong.ci <- rst;
    cong.mean.ci <- apply(cong.ci,2,mean);
    cong.ci <- cbind(cong.mean.ci[c(1:3, 7:11)], cong.mean.ci[c(4:6, 12:16)]) ### 90%, and 95% rate for cong
    
    ci.summary <- cbind(uncong.ci, cong.ci);
    ci.summary <- as.data.frame(ci.summary);
    colnames(ci.summary) <- c("90% uncong", "95% uncong", "90% cong", "95% cong");
    write.table(ci.summary, file=paste(headername, "Size", size, "_ttl", ttlndata, "_ci_summary.txt",sep=""));
    
    #### obtain the bias and MSE results for uncong and cong in AMAR case
    rst <- NULL;
    for(data.i in 1:ttlndata) {
        tmp <-read.table(file=paste(headername, "allcong/Size", size, "_data", data.i, "_mlecong_CG_parcl_rst.txt",sep=""));
        tmp <- tmp[[1]];
        rst <- rbind(rst, tmp);
    }
    cong.beta <- rst[,1:3];
    cong.summary <- obtain.rst.summary(cong.beta, BETA);
    
    rst <- NULL;
    for(data.i in 1:ttlndata) {
        tmp <-read.table(file=paste(headername, "alluncong/Size", size, "_data", data.i, "_stduncong_par_nimpt1000.txt",sep=""));
        tmp <- tmp[[1]];
        rst <- rbind(rst, tmp);
    }
    uncong.beta <- rst[,1:3];
    uncong.summary <- obtain.rst.summary(uncong.beta, BETA);
    
    par.summary <- cbind(uncong.summary, cong.summary);
    par.summary <- as.data.frame(par.summary);
    colnames(par.summary) <- c("uncong.bias", "uncong.se", "cong.bias", "cong.se");
    write.table(par.summary, file=paste(headername, "Size", size, "_ttl", ttlndata, "_par_compare_rst.txt",sep=""));

    list(ci.summary=ci.summary, par.summary=par.summary)
}



