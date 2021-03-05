/*
This file contains the SAS PROC IML code to obtain a random matrix
that can be used for a starting value. Deatils are in the paper:

Gradient Projection Algorithms and Software for Arbitrary.
Rotation Criteria in Factor Analysis.

by

Coen A. Bernaards and Robert I. Jennrich

input required:
k = Number of dimensions requested. 
output:
rm = orthogonal random matrix.
*/

proc iml;

k = 3;

rn = I(k);
do i = 1 to k;
 do j = 1 to k;
  rn[i,j] = normal(0);
 end;
end;

call qr(rm,r,piv,lindep,rn);
print rm;

quit;

