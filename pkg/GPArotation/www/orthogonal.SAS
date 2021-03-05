/*
This file contains the SAS PROC IML code that appears in the paper:

Gradient Projection Algorithms and Software for Arbitrary
Rotation Criteria in Factor Analysis.

by

Coen A. Bernaards and Robert I. Jennrich

Website: http://www.stat.ucla.edu/research
*/

proc iml;
 
start GPForth;
 al=1;
 RESET NONAME;
 print ,"ITER         Q    LOG10      ALPHA";
 L=A * T;
 run vgQ;
 ft=Q;
 G=A`*Gq;
 do iter = 0 to 500;
  M=T` * G;
  S=(M + M`)/2;
  Gp=G - T*S;
  s=sqrt(trace(Gp` * Gp));
  sl = log(s)/log(10);
  print ,iter ft sl al;
  if (s < 0.00001) then goto skip;
  al=2*al;
  do i = 0 to 10;
    X=T-al*Gp;
    call SVD (u,d,v,X);
    Tt=U*V`;
    L=A*Tt;
    run vgQ;
	if (Q < ft-.5*s*s*al) then
     i=10;
    else
     al=al/2;
  end;
  T=Tt;
  ft=Q;
  G=A`*Gq;
  skip: if (s < 0.00001) then iter=500;
  end;
 RESET NAME;
 Lh=A*T;
finish GPForth;

/* Quartimax subroutine */
start vgQ;
  L2 = L#L;
  Q = -(trace(L2`*L2)/4);
  Gq = -L#L#L;
finish vgQ;

/* Varimax subroutine */
/*start vgQ;
  L2 = L#L;
  n  = nrow(L2);
  cm = L2[+,]/n;
  QL = L2 - repeat(cm,n,1);
  Q = -(trace(QL`*QL)/4);
  Gq = -L#QL;
finish vgQ; */

/* Bentler's criterion */
/*start vgQ;
  L2 = L#L;
  M = L2` * L2;
  D = diag(M);
  Q = -(log(det(M))-log(det(D)))/4;
  Gq = -L # (L2 * (inv(M)-inv(D)));
finish vgQ; */

/* Minimum entropy */
/*start vgQ;
  L2 = L#L;
  Q = -sum(L2 # log(L2))/2;
  Gq = -(L # log(L2) + L);
finish vgQ; */

/* Crawford-Ferguson  
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
  gamma =0 Quartimin
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

/* Target rotation.
  Needs Target with hypothesized factor loadings.
  e.g. Target = J(8,2,0); (i.e. only zeroes) */
/*start vgQ;
  Target = J(8,2,1);
  Btilde = Target;
  Q = sum((L-Btilde)##2);
  Gq = 2*(L-Btilde);
finish vgQ; */

/* Partially specified target rotation.
  Needs weight matrix W with 1's at specified values, 0 otherwise
  e.g. W = {1 0,1 0,1 0,1 0, 0 1,0 1,0 1,0 1};
  When W has only 1's this is procrustes rotation
  Needs a Target matrix Target with hypothesized factor loadings.
  e.g. Target = J(8,2,0); (i.e. only zeroes)*/
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
 Needs k: Number of close to zero loadings */ 
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

/* Tandem Principle I. Comrey (1967) */
/*start vgQ;
  LL = L * L`;
  LL2 = LL # LL;
  L2 = L#L;
  Q = -trace(L2` * LL2 * L2);
  Gq1 = 4 * L # (LL2 * L2);
  Gq2 = 4 * (LL # (L2 * L2`)) * L;
  Gq  = -Gq1 - Gq2;
finish vgQ; */

/* Tandem Principle II. Comrey (1967) */
/*start vgQ;
  LL = L * L`;
  LL2 = LL # LL;
  L2 = L#L;
  Q = trace(L2` * (1-LL2) * L2);
  Gq1 = 4 * L # ((1-LL2) * L2);
  Gq2 = 4 * (LL # (L2 * L2`)) * L;
  Gq  = Gq1 - Gq2;
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

/* McCammon minimum entropy ratio. */
/*start vgQ;
  L2 = L#L;
  k = ncol(L);
  pn = nrow(L);
  SS = trace(L` * L);
  M1 = J(pn, pn, 1);
  S2 = L2[+,];
  den = repeat(s2, pn, 1);
  P = L2/den;
  Q1 = -sum(P # log(P));
  H = -(log(P) + 1);
  R = M1 * L2;
  G1 = H/R - M1 * (L2 # H/(R#R));
  p2 = S2/SS;
  Q2 = -sum(p2 # log(p2));
  h = -(log(p2) + 1);
  alpha = h * p2`;
  G2 = repeat(h,pn,1)/SS - alpha * J(pn, k, 1);
  Gq = 2 * L # (G1/Q1 - G2/Q2);
  Q = log(Q1) - log(Q2);
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
  run GPForth;
  print ,T Lh; 
end;
quit;
