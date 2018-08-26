
#include <R.h>
#include <Rmath.h>


void XuanAMAR_eytgvnprevyv1(double *delta, double *alpha, double *gamma,
	       double *v, int *nalpha, int *t, double *rst) {
  int    i;
  double e0, e1;
  
  //printf("delta= %f \n", delta); 
  e0 = *delta;

  for (i=0; i < *nalpha; i++) {
    e0 = e0 + alpha[i] * v[i];
  }

  if (1 == *t) {
    //for prevy=0, that is t=1 the baseline
    e1 = e0;
    rst[0] = exp(e0)/(1+exp(e0));
    rst[1] = exp(e1)/(1+exp(e1));
  } else {
    e1 = e0 + (*gamma);
    //printf("e0= %f, e1= %f", e0, e1);
    rst[0] = exp(e0)/(1+exp(e0)); //if prevy=0
    rst[1] = exp(e1)/(1+exp(e1)); //prevy=1
  }
    //Xuan: for the base outcome, get rst[0] and rst[1]
    //Xuan: it gives P(Y_it=1|Y_i,(t-1)=prevy, V_i)
  //printf("Pr(Yt=1|Y_(t-1)=y,v)=%f, %f \n", rst[0],rst[1]);
}


//Working for AMAR case, based on imputation model get E(Y_t|V, Y_(t-1)=y) for y=prevy
void XuanAMAR_eytgvnprevyv(double *delta, double *alpha, double *gamma,
	       double *v, int *nalpha, int *t, double *rst) {
  int    i;
  double e0, e1;
  
  // if(1 < *t)
  //  { printf("inside XuanAMAR_eytgvnprevyv/n"); }
  e0 = *delta;

  for (i=0; i < *nalpha; i++) {
    e0 = e0 + alpha[i] * v[i];
  }

  if (1 == *t) {
    //for prevy=0, that is t=1 the baseline
    e1 = e0;
    rst[0] = exp(e0)/(1+exp(e0));
    rst[1] = exp(e1)/(1+exp(e1));
  } else {
    e1 = e0 + (*gamma);
    //printf("e0= %f, e1= %f", e0, e1);
    rst[0] = exp(e0)/(1+exp(e0)); //if prevy=0
    rst[1] = exp(e1)/(1+exp(e1)); //prevy=1
  }
    //Xuan: for the base outcome, get rst[0] and rst[1]
    //Xuan: it gives P(Y_it=1|Y_i,(t-1)=prevy, V_i)
  // printf("Pr(Yt=1|Y_(t-1)=y,v)=%f, %f \n", rst[0],rst[1]);
}




//Xuan: some modification -- working for AMAR case -- transitonal model
//For T=1, as only 1 delta corresponding to Y0=0 (since Y0=0 with prob=1)
void XuanAMAR_transitional_getdiff_t1(double *delta, double *prevy, double *alphat, double *beta, double *gamma, int *t, double *vecvall, double *vgvnprevy0, double *vgvnprevy1, int *ncolv, int *nrowv, int *margin, double *rst) {

    int i,j;
    double lhs=0, rhs=0;  //left hand side probability
    double curv[*ncolv];
    double curpyt[2];
 
    //get lhs. i.e. transition 
    if (1 == *margin) {
	lhs = beta[0];
    } else if (2 == *margin) {
	lhs = beta[0] + beta[*t+1];
    } else if (3 == *margin) {
      lhs = beta[0] + beta[1]*((*t)-1)+beta[2]*(*prevy);
    } else if (4 == *margin) {
        if (0 == *t) {
	    lhs = beta[0];
	} else {
	    lhs = beta[1];
	}
    }
   // lhs 1 for prevy=0
    lhs = exp(lhs)/(1+exp(lhs));
 
    //get rhs
    //Xuan: nrowv= total possible realization of v
    for (i=0; i < *nrowv; i++) {
      for (j=0; j<(*ncolv); j++)
      	{ curv[j] = vecvall[(i)*(*ncolv)+j]; }
      //curpyt gives Pr(Yt=1|Y_(t-1)=y, v) for y=0, and 1
      XuanAMAR_eytgvnprevyv(delta, alphat, gamma, curv, ncolv, t, curpyt);
	//  printf("%f, %f /n", curpyt[0], curpyt[1]);
        // sum_y Pr(Yt=1|Y_(t-1)=y, v)*Pr(v|Y_(t-1)=y) 
        if(*prevy==0) { rhs = rhs + curpyt[0]*vgvnprevy0[i];}
        else { rhs = rhs + curpyt[1]*vgvnprevy1[i]; }
    }

    // printf("prob are %f, %f /n", lhs, rhs);
    
   *rst = log(rhs)-log(lhs);   
}



//Xuan: some modification -- working for AMAR case -- transitonal model
// For T>=2, two deltas corresponds to previous Y=0 and 1
void XuanAMAR_transitional_getdiff_t2(double *deltas, double *alphat, double *beta, double *gamma, int *t, double *vecvall, double *vgvnprevy0, double *vgvnprevy1, int *ncolv, int *nrowv, int *margin, double *rst) {

    int i,j;
    double lhs1=0, rhs1=0, lhs2=0, rhs2=0, delta=0;  //left hand side probability
    double curv[*ncolv];
    double curpyt[2];
   
    // printf("deltas are %f, %f  \n", (deltas[0]), (deltas[1]));
    //get lhs. i.e. transition 
    if (3 == *margin) {
	lhs1 = beta[0] + beta[1]*((*t)-1);
        lhs2 = lhs1 + beta[2];
    } else {
	    lhs1 = beta[1];
            lhs2 = beta[1];
    }
    
   // lhs 1 for prevy=0
    lhs1 = exp(lhs1)/(1+exp(lhs1));
    lhs2 = exp(lhs2)/(1+exp(lhs2));
    //  printf("lhs are %f, %f \n", lhs1, lhs2);
    //get rhs
    //Xuan: nrowv= total possible realization of v
    for (i=0; i < *nrowv; i++) {
      for (j=0; j<(*ncolv); j++)
      	{ curv[j] = vecvall[(i)*(*ncolv)+j]; }
      //curpyt gives Pr(Yt=1|Y_(t-1)=y, v) for y=0, and 1
       
      delta=deltas[0];	
       XuanAMAR_eytgvnprevyv1(&delta, alphat, gamma, curv, ncolv, t, curpyt);
       // sum_y Pr(Yt=1|Y_(t-1)=y, v)*Pr(v|Y_(t-1)=y) 
       // printf("Pr(Yt=1|Y_(t-1)=y,v)=%f, %f/n", curpyt[0], curpyt[1]);
         rhs1 += curpyt[0]*vgvnprevy0[i];
	 delta=deltas[1];
	XuanAMAR_eytgvnprevyv1(&delta, alphat, gamma, curv, ncolv, t, curpyt);
         rhs2 += curpyt[1]*vgvnprevy1[i]; 
       
    }

    //  printf("rhs2 are %f, %f \n", rhs1, rhs2);
    
   rst[0] = log(rhs1)-log(lhs1);   
   rst[1] = log(rhs2)-log(lhs2);
}








/////////////////////////


void eytgvnprevyvr(int *t, double *v, int *r, int *prevy, double *delta, double *alpha, int *ncolv, double *gamma, double *rst) {
   
  int i;
  double e0=0;
  e0=*delta + (*gamma) * (*prevy);
  //printf("d= %f *** y-1== %d  gamma==%f   e0==%f", *delta, *prevy, *gamma, e0);
  for (i=0; i< *ncolv; i++) {
    e0=e0+ alpha[i] * v[i];  // e0=delta + V*alpha1
   }
  //printf("*** e0== %f", e0);
  *rst = 1-1/(1+exp(e0));
}



void ttl_obsloglikelihood(int *nT, int *nsub, double *v, int *r, int *ncolv, int *yobs, double *delta, double *alpha, double *gamma, double *rst){

  int i, j, t, curr, prevy=0;
    double curpyt=0, indobs=0,  curdelta=0; 
    double curv[*ncolv];
    int cury[*nT];
    double alphat[*ncolv];
   *rst = 0;
 
    for(i=1; i<=(*nsub); i++) {
      //printf("\n i= %d -----", i);
      for(j=0; j<(*ncolv); j++) {
         curv[j] = v[(i-1) * (*ncolv)+j];   
       }
      for(j=0; j<(*nT); j++){
         cury[j] = yobs[(i-1) * (*nT)+j];
       }  
      curr = r[i-1];
      prevy = 0;
      indobs = 0;
    
      for(t=1; t<= curr; t++) { // for Y1, Y2, ... Y(r) 
      //calculate the E(Yt|Y(t-1),V,R)
      for(j=0; j<(*ncolv); j++) {
    	alphat[j] = alpha[(t-1)* (*ncolv)+j];
      }
      curdelta = delta[(t-1) * 2+prevy];
      eytgvnprevyvr(&t, curv, &curr, &prevy, &curdelta, alphat, ncolv, gamma, &curpyt); /// for a specific (V,R)
      indobs += log(curpyt * (cury[t-1]) + (1-curpyt) * (1-cury[t-1]));
      prevy = cury[t-1]; 
      }
      //printf("i th ll== %f", indobs);
    *rst += indobs;
    //printf("\n");
  }
}


