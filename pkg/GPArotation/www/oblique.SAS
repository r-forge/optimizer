/*
This file contains the SAS PROC IML code that appears in the paper:

Gradient Projection Algorithms and Software for Arbitrary
Rotation Criteria in Factor Analysis.

by

Coen A. Bernaards and Robert I. Jennrich

Website: http://www.stat.ucla.edu/research
*/

proc iml;
 
start GPFoblq;
 al=1;
 RESET NONAME;
 print ,"ITER         Q    LOG10      ALPHA";
 Ti = inv(T);
 L=A * Ti`;
 run vgQ;
 ft=Q;
 G=-(L` * Gq * Ti)`;
 do iter = 0 to 500;
  TG=T # G;
  Gp=G - T * diag(TG[+,]);
  s=sqrt(trace(Gp` * Gp));
  sl = log(s)/log(10);
  print ,iter Q sl al;
  if (s < 0.00001) then goto skip;
  al=2*al;
  do i = 0 to 10;
    X=T-al * Gp;
	X2=X#X;
	v=1/sqrt(X2[+,]);
	Tt=X * diag(v);
	Ti=inv(Tt);
	L=A * Ti`;
	run vgQ;
	if (Q < ft-.5*s*s*al) then
     i=10;
    else
     al=al/2;
  end;
  T=Tt;
  ft=Q;
  G=-(L` * Gq * Ti)`;
  skip: if (s < 0.00001) then iter=500;
  end;
 RESET NAME;
 Th=t;
 Lh=L;
 phi=T` * T;
finish GPFoblq;

/* Quartimin rotation */
start vgQ;
  L2 = L#L;
  k = ncol(L);
  M = J(k)-I(k);
  Q = sum(L2 # (L2 * M))/4;
  Gq = L # (L2 * M);
finish vgQ; 

/* Bentler's criterion */
/*start vgQ;
  L2 = L#L;
  M = L2` * L2;
  D = diag(M);
  Q = -(log(det(M))-log(det(D)))/4;
  Gq = -L # (L2 * (inv(M)-inv(D)));
finish vgQ; */

/* Crawford-Ferguson.    
  kappa = 0 Quartimax
  kappa = 1/p  Varimax   
  kappa = k/(2*p) Equamax
  kappa = (k-1)/(p+k-2) Parsimax
  kappa = 1  Factor parsimony */  
/*start vgQ;
  kappa = 1;
  k = ncol(L);
  p = nrow(L);
  N = J(k)-I(k);
  M = J(p)-I(p);
  L2 = L#L; 
  f1 = (1-kappa)* trace(L2` * (L2 * N))/4;
  f2 = kappa * trace(L2` * (M * L2))/4;
  Q = f1 + f2;
  Gq = (1-kappa) * (L # (L2 * N)) + kappa * (L # (M * L2));
finish vgQ; */

/* Oblimin family.
  gamma = 0 Quartimin
  gamma =.5 Bi-quartimin
  gamma = 1 Covarimin */
/*start vgQ;
  gamma = 0;
  k = ncol(L);
  p = nrow(L);
  N = J(k)-I(k); 
  L2 = L#L;
  Q = sum(L2 # ((I(p)-gamma # J(p,p,1/p)) * L2 * N))/4;
  Gq = L # ((I(p)-gamma # J(p,p,1/p)) * L2 * N);
finish vgQ; */

/* Partially specified target rotation.
  Needs weight matrix W with 1's at specified values, 0 otherwise
  e.g. W = {1 0,1 0,1 0,1 0, 0 1,0 1,0 1,0 1};
  When W has only 1's this is procrustes rotation
  Needs a Target matrix Target with hypothesized factor loadings.
  e.g. Target = J(8,2,0); (i.e. only zeroes) */
/*start vgQ;
  W = {1 0,
       1 0,
	   1 0,
	   1 0,
	   0 1,
	   0 1,
	   0 1,
	   0 1};
  Target = J(8,2,0);
  Btilde = W # Target;
  Q = sum((W#L-Btilde)##2);
  Gq = 2*(W#L-Btilde);
finish vgQ; */

/* Simplimax.
  k: Number of close to zero loadings */ 
/*start vgQ; 
  k = 8;
  L2=L#L;
  B =L2;
  B[rank(B)]=L2;
  tr=B[k];
  Imat = L2 <= tr;
  Q = sum(Imat # L2);
  Gq = 2*Imat # L;
finish vgQ; */

/* Oblimax */
/*start vgQ;
  Q = -(log(sum(L#L#L#L))-2*log(sum(L#L)));
  Gq = -(4*(L#L#L)/sum(L#L#L#L)-4*L/sum(L#L));
finish vgQ; */

/* Geomin 
  epsilon: needed to avoid numerical problems.*/
/*start vgQ;
  k = ncol(L);
  p = nrow(L);
  eps = 0.01;
  L2 = L#L+eps;
  LL2= log(L2);
  SL = LL2[,+];
  pro = exp(SL/k);
  Q = sum(pro);
  RM = repeat(pro,1,k);
  Gq = (2/k)*(L/L2)#RM;
finish vgQ;*/

/* Target rotation.
  Needs Target with hypothesized factor loadings.
  e.g. Target = J(8,2,0); (i.e. only zeroes) */
/*start vgQ;
  Target = J(8,2,1);
  Btilde = Target;
  Q = sum((L-Btilde)##2);
  Gq = 2*(L-Btilde);
finish vgQ; */

/* Infomax. McKeon (1968) */
/*start vgQ;
  L2 = L#L;
  k = ncol(L);
  p = nrow(L);
  SS = trace(L` * L);
  S1 = L2[,+];
  S2 = L2[+,];
  E0 = L2/SS;
  E1 = S1/SS;
  E2 = S2/SS;
  Q0 = sum(-E0 # log(E0));
  Q1 = sum(-E1 # log(E1));
  Q2 = sum(-E2 # log(E2));
  Q = log(k) + Q0 - Q1 - Q2;
  H = -log(E0) + 1;
  alpha0 = sum(L2 # H)/(SS * SS);
  G0 = H/SS - alpha0 * J(p, k, 1);
  H1 = -log(E1) + 1;
  alpha1 = S1` * H1/(SS * SS);
  G1 = repeat(H1, 1, k)/SS - alpha1 * J(p, k, 1);
  H2 = -log(E2) + 1;
  alpha2 = S2 * H2`/(SS * SS);
  G2 = repeat(H2, p, 1)/SS - alpha2 * J(p, k, 1);
  Gq = 2 # L # (G0 - G1 - G2);
finish vgQ;*/
    
do;
  /* 8 physical variables from Harman */
  A ={.830 -.396,
      .818 -.469,
      .777 -.470,
      .798 -.401,
      .786  .500,
      .672  .458,
      .594  .444,
      .647  .333};
  T = I(2);
  run GPFoblq;
  print ,T Lh phi;
end;
quit;
