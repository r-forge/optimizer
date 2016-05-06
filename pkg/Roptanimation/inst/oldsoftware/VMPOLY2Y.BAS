5 PRINT "DRIVER -- GENERAL PARAMETER ESTIMATION DRIVER - 851017,851128"
10 REM INPUT "FILE FOR CONSOLE IMAGE = ";G$
15 REM IF G$="" THEN LET G$="NUL"
20 REM OPEN G$ FOR APPEND AS 3: REM !! may be non-standard
25 REM PRINT #3,"DRIVER -- GENERAL PARAMETER ESTIMATION DRIVER - 851017,851128"
30 DIM B(25),X(25),O(25,3),Y(91,4),PP(25,2),SS(25,2)
35 REM in case of extra line
40 DIM G(25),T(25),C(25),H(25,25): REM VM DIMs
44 DIM A(25,25): REM POSTVM DIM to save Hessian
45 REM in case of extra line
50 REM PRINT #3,TIME$,DATE$: REM !!
55 PRINT TIME$,DATE$: REM !!
60 REM CALLS:
65 REM    FUNCTION MINIMIZER  -- line 1000
70 REM    SETUP               -- line 3000
75 REM    POST-SOL'N ANALYSIS -- line 6000
80 REM    RESULT ANALYSIS     -- line 4500
85 REM
90 REM X$, T$ and U$ are used for responses and time
95 GOSUB 3000: REM setup of user-specified problem
100 LET I8=0: REM count of iterations or gradient evaluations
105 LET I9=0: REM count of function or sum of squares evaluations
110 PRINT P$: REM display problem title (user supplied)
115 REM PRINT #3,P$
120 GOSUB 500: REM display bounds
125 PRINT "ENTER INITIAL VALUES FOR PARAMETERS ( [cr] = Y ) ";
130 INPUT X$
135 REM PRINT #3,"ENTER INITIAL VALUES FOR PARAMETERS ( [cr] = Y ) ";X$
140 IF X$="N" THEN 175
145 IF X$="n" THEN 175
150 FOR J=1 TO N
155 PRINT "B(";J;")=";
160 INPUT B(J)
165 REM PRINT #3,"B(";J;")=";B(J)
170 NEXT J
175 INPUT " are masks or bounds to be set or altered ([cr] = no) ";X$
180 REM PRINT #3," are masks or bounds to be set or altered ([cr] = no) ";X$
185 IF X$="y" OR X$="Y" THEN 195
190 GOTO 285: REM note that already masked parameters stay masked
195 FOR J=1 TO N
200 PRINT "b(";J;")=";B(J);"  mask ([cr] = no) ";
205 INPUT X$
210 REM PRINT #3,"b(";J;")=";B(J);"  mask ([cr] = no) ";X$
215 LET O(J,3)=0: REM masked parameter
220 IF X$="y" OR X$="Y" THEN 280: REM if masked, no need for bounds
225 LET O(J,3)=1: REM free parameter
230 PRINT "lower bound =";O(J,1);"  ([cr] = no change) ";
235 INPUT X$
240 REM PRINT #3,"lower bound =";O(J,1);"  ([cr] = no change) ";X$
245 IF X$="" THEN GOTO 255: REM no change
250 LET O(J,1)=VAL(X$): REM new lower bound
255 PRINT "upper bound =";O(J,2);"  ([cr] = no change) ";
260 INPUT X$
265 REM PRINT #3,"upper bound =";O(J,2);"  ([cr] = no change) ";X$
270 IF X$="" THEN GOTO 280: REM no change
275 LET O(J,2)=VAL(X$): REM new upper bound
280 NEXT J
285 LET S1=0: REM set steplength parameters within minimizer
290 LET S2=0: REM or change these lines of code
295 GOSUB 560: REM mask and bound check
300 LET T$=TIME$: REM !! SYSTEM DEPENDENT
305 GOSUB 1000: REM function minimization routine
310 LET U$=TIME$: REM !! system dependent
315 GOSUB 435: REM !! system dependent time display routine
320 REM PRINT "CALCULATED FUNCTION MINIMUM = ";F0
325 REM PRINT #3,"CALCULATED FUNCTION MINIMUM = ";F0
330 REM PRINT "PARAMETER ESTIMATES"
335 REM PRINT #3,"PARAMETER ESTIMATES"
340 FOR J=1 TO N
345 REM PRINT "B(";J;") = ";X(J)
350 REM PRINT #3,"B(";J;") = ";X(J)
355 NEXT J
360 GOSUB 6000: REM post-solution analysis
365 IF F0<=F1 THEN 375 ELSE 305: REM no lower point found
367 REM INPUT "continue from lower point found ([cr] = Y) ";X$
369 REM PRINT #3,"continue from lower point found ([cr] = Y) ";X$
371 REM IF X$<>"N" AND X$<>"n" THEN 110
373 REM finished
375 INPUT " call problem-dependent results-analysis ( [cr] = N ) ";X$
385 REM PRINT #3," call problem-dependent results-analysis ( [cr] = N ) ";X$
390 IF X$="Y" THEN GOSUB 4500
395 IF X$="y" THEN GOSUB 4500
400 PRINT " another set of starting parameters ( [cr] = N ) ";
405 INPUT X$
410 REM PRINT #3," another set of starting parameters ( [cr] = N ) ";X$
415 IF X$="Y" THEN 100
420 IF X$="y" THEN 100
425 REM CLOSE #3
430 STOP
435 LET T1=VAL(LEFT$(T$,2)): REM !!
440 LET T2=VAL(MID$(T$,4,2)): REM !!
445 LET T3=VAL(RIGHT$(T$,2)): REM !!
450 LET U1=VAL(LEFT$(U$,2)): REM !!
455 LET U2=VAL(MID$(U$,4,2)): REM !!
460 LET U3=VAL(RIGHT$(U$,2)): REM !!
465 LET T4=T3+60*T2+3600*T1: REM !!
470 LET U4=U3+60*U2+3600*U1: REM !!
475 LET T5=U4-T4: REM !!
480 IF T5<0 THEN LET T5=T5+24*3600!: REM !! ASSUME NEVER OVER 24 HRS
485 REM PRINT "  ELAPSED SECS=";T5;" AFTER";I8;" GRAD &";I9;" FN EVAL": REM !!
490 REM PRINT #3,"  ELAPSED SECS=";T5;" AFTER";I8;" GRAD &";I9;" FN EVAL": REM !!
495 RETURN: REM !!
500 PRINT "bounds on parameters for problem"
505 REM PRINT #3,"bounds on parameters for problem"
510 FOR J=1 TO N
515 IF O(J,3)=0 THEN 535
520 PRINT O(J,1);" <=   b(";J;")  <=  ";O(J,2)
525 REM PRINT #3,O(J,1);" <=   b(";J;")  <=  ";O(J,2)
530 GOTO 545
535 PRINT " ****  b(";J;") WILL BE MASKED (FIXED) ****"
540 REM PRINT #3," ****  b(";J;") WILL BE MASKED (FIXED) ****"
545 NEXT J
550 RETURN
555 REM check bounds and set up indicator
560 LET I5=N: REM count masks set
565 FOR J=1 TO N
570 IF O(J,3)=0 THEN  665
575 LET O(J,3)=1: REM reset free or active bound (below)
580 LET I5=I5-1: REM not a masked parameter
585 IF B(J)>O(J,1) THEN  610
590 PRINT "parameter B(";J;") reset from ";B(J);"  to  ";O(J,1)
595 REM PRINT #3,"parameter B(";J;") reset from ";B(J);"  to  ";O(J,1)
600 LET B(J)=O(J,1)
605 LET O(J,3)=-2 : REM lower bound active
610 IF B(J)<O(J,2) THEN  665
615 PRINT "parameter B(";J;") reset from ";B(J);"  to  ";O(J,2)
620 REM PRINT #3,"parameter B(";J;") reset from ";B(J);"  to  ";O(J,2)
625 IF O(J,3)>-2 THEN 655:REM check if lower, upper bounds equal (masked)
630 LET O(J,3)=0: REM set mask on
635 LET I5=I5+1: REM adjust mask count
640 PRINT " mask imposed on parameter B(";J;") at upper bound ";O(J,2)
645 REM PRINT #3," mask imposed on parameter B(";J;") at upper bound ";O(J,2)
650 GOTO 665
655 LET B(J)=O(J,2)
660 LET O(J,3)=-1: REM upper bound active
665 LET X(J)=B(J): REM to ensure "best" parameters are saved
670 NEXT J
675 RETURN
1000 REM PRINT "VM -- VARIABLE METRIC WITH BOUNDS CONSTRAINTS -- 851118"
1004 REM PRINT #3,"VM -- VARIABLE METRIC WITH BOUNDS CONSTRAINTS -- 851118"
1008 REM CALLS:  
1012 REM     FUNCTION F(B)  -- line 2000
1016 REM     GRADIENT G(B)  -- line 2500
1020 REM     ENVRON (computing environment) -- line 7120
1024 REM
1028 REM INPUTS TO THE ROUTINE:
1032 REM    B() -- a vector of initial parameter estimates
1036 REM    N   -- the number of parameters in the function F(B)
1040 REM    O( , ) -- masks and bounds
1044 REM    I5  -- the number of masked parameters. Must be zero if
1048 REM           masks are not used.
1052 REM OUTPUT FROM THE ROUTINE:
1056 REM    X() -- a vector of final parameter estimates for the
1060 REM           values of the parameters which minimize the function.
1064 REM    F0  -- value of the function at the minimum
1068 REM    G() -- value of the gradient at the minimum
1072 REM    A() -- value of the Hessian at the minimum
1076 REM    I8  -- the number of gradient evaluations
1080 REM    I9  -- the number of function evaluations
1084 REM
1088 GOSUB 7120: REM computing environment
1092 REM DIM B(N),X(N),G(N),T(N),H(N,N),C(N) needed, A(N,N) in POSTVM
1096 LET J8=1: REM no of fn eval before parameter display
1100 LET J7=1: REM counter for parameter display
1104 GOSUB 1648: REM function in F
1108 IF I3=0 THEN 1124
1112 PRINT "FUNCTION NOT COMPUTABLE AT INITIAL POINT"
1116 REM PRINT #3,"FUNCTION NOT COMPUTABLE AT INITIAL POINT"
1120 STOP
1124 LET I8=I8+1: REM count gradient evaluation
1128 LET I7=I8: REM records last gradient=steepest descent count
1132 GOSUB 2500: REM gradient in G( )
1136 LET F0=F: REM save lowest value so far
1140 FOR I=1 TO N
1144 FOR J=1 TO N
1148 IF I8>I7 THEN LET A(I,J)=H(I,J): REM save dispersion information
1152 LET H(I,J)=0
1156 NEXT J
1160 LET H(I,I)=1: REM initialize Hessian inverse to unit matrix
1164 NEXT I
1168 LET I7=I8: REM save itn count at last steepest desc. dirn.
1172 LET I3=0: REM ensure "non-computable function" flag not set
1176 REM PRINT I8;" GRADIENTS, ";I9;" FNS, LOSS FN=";F0
1180 REM PRINT #3,I8;" GRADIENTS, ";I9;" FNS, LOSS FN=";F0
1184 GOSUB 1664: REM display parameters
1188 REM top of iteration
1192 FOR I=1 TO N
1196 LET X(I)=B(I): REM save "best" parameters
1200 LET C(I)=G(I): REM save last gradient
1204 NEXT I
1208 FOR J=1 TO N: REM adjust gradient for constraints
1212 IF O(J,3)=0 THEN 1240: REM masked
1216 IF O(J,3)=1 THEN 1244: REM free
1220 IF (O(J,3)+1.5)*G(J)<0 THEN 1240: REM L.M. non-negative
1224 LET O(J,3)=1: REM free parameter
1228 REM PRINT "freeing parameter ";J
1232 REM PRINT #3,"freeing parameter ";J
1236 GOTO 1244
1240 LET G(J)=0: REM active constraint
1244 NEXT J
1248 REM form  H * G, and GT * H * G
1252 LET G8=0: REM gradient times step
1256 FOR I=1 TO N
1260 LET S9=0: REM accumulator for H * G row I
1264 FOR J=1 TO N
1268 LET S9=S9-H(I,J)*G(J): REM note negative sign
1272 NEXT J
1276 LET T(I)=S9: REM so that T() solves Newton eqns.
1280 IF O(I,3)=1 THEN 1288
1284 LET T(I)=0: REM constraint forces search component to zero
1288 LET G8=G8+T(I)*G(I): REM gradient projection on search direction
1292 NEXT I
1293 REM IF I8=I7 THEN PRINT " GRADIENT PROJECTION NORM=";-G8
1294 REM IF I8=I7 THEN PRINT #3," GRADIENT PROJECTION NORM=";-G8
1296 IF G8<0 THEN 1316
1300 REM PRINT "UPHILL SEARCH DIRECTION"
1304 REM PRINT #3,"UPHILL SEARCH DIRECTION"
1308 IF I7=I8 THEN 1644: REM if dirn. is sd, can't continue
1312 GOTO 1140: REM otherwise reset to steepest descent (SD)
1316 LET D1=-G8
1320 LET S1=1: REM step reduction line search
1324 LET S8=S1: REM find max-step and check new constraint
1328 FOR J=1 TO N
1332 IF O(J,3)<=0 THEN 1368: REM masked or active constraint
1336 IF T(J)=0 THEN 1368: REM no sense in checking
1340 IF T(J)>0 THEN 1356: REM check upper bound
1344 LET S7=(O(J,1)-X(J))/T(J): REM distance to lower bd.
1348 IF S7>S8 THEN 1368: REM smaller step in effect
1352 GOTO 1364
1356 LET S7=(O(J,2)-X(J))/T(J): REM distance to upper bd.
1360 IF S8<=S7 THEN 1368: REM already smaller step in effect (<= not <)
1364 LET S8=S7
1368 NEXT J
1372 LET I6=0: REM counter for parameters unchanged in step
1376 FOR J=1 TO N
1380 IF O(J,3)=0 THEN 1400: REM mask check (count as unchanged)
1384 LET B(J)=X(J)+S8*T(J)
1388 IF E5+X(J)<>E5+B(J) THEN 1400
1392 LET I6=I6+1
1396 LET B(J)=X(J): REM to avoid "drift"
1400 NEXT J: REM end of step
1404 IF I6<N-I5 THEN 1432: REM check if no unchanged parameters
1408 IF I8=I7 THEN 1644: REM done if already steepest descent
1412 REM PRINT "Steepest Descent ";
1416 REM PRINT #3,"Steepest Descent ";
1420 GOTO 1140
1424 REM note that this convergence test is sufficient if bounded 
1428 REM constraints are released at start of each iteration
1432 GOSUB 1648: REM evaluate function
1436 IF I3=1 THEN 1444: REM check if function computable
1440 IF F<F0-D1*S8*.0001 THEN 1460: REM acceptable point
1444 LET S1=.2*S8: REM reduce stepsize
1448 REM PRINT "*";: REM symbol * means stepsize is being reduced
1452 REM PRINT #3,"*";
1456 GOTO 1324
1460 LET F0=F: REM save new lowest value
1464 FOR J=1 TO N: REM check bounds
1468 IF O(J,3)<=0 THEN 1520: REM masked or constrained parameter
1472 IF (B(J)-O(J,1))>E9*(ABS(O(J,1))+1) THEN 1492: REM check lower bound
1476 REM PRINT "lower ";
1480 REM PRINT #3,"lower ";
1484 LET O(J,3)=-2: REM lower bound active
1488 GOTO 1512
1492 IF (O(J,2)-B(J))<=E9*(ABS(O(J,2))+1) THEN 1500: REM check upper bound
1496 GOTO 1520: REM no constraint
1500 LET O(J,3)=-1: REM upper bound active
1504 REM PRINT "upper ";
1508 REM PRINT #3,"upper ";
1512 REM PRINT "bound activated on parameter ";J
1516 REM PRINT #3,"bound activated on parameter ";J
1520 NEXT J
1524 LET I8=I8+1: REM prepare update -- get new gradient
1528 GOSUB 2500
1532 LET I3=0: REM & reset flag to indicate update failure
1536 LET D1=0: REM BFGS update of inverse Hessian approximation
1540 LET D2=0
1544 FOR I=1 TO N
1548 LET T(I)=T(I)*S8
1552 LET C(I)=G(I)-C(I): REM change in gradient
1553 IF O(I,3)=0 THEN LET C(I)=0: REM FOR MASK
1556 LET D1=D1+T(I)*C(I)
1560 NEXT I
1564 FOR I=1 TO N
1568 LET S9=0
1572 FOR J=1 TO N
1576 LET S9=S9+H(I,J)*C(J)
1580 NEXT J
1584 LET X(I)=S9
1588 LET D2=D2+S9*C(I)
1592 NEXT I
1596 IF D1<=0 THEN 1628: REM test of positive definite-ness of H
1600 LET D2=1+D2/D1
1604 FOR I=1 TO N
1608 FOR J=1 TO N
1612 LET H(I,J)=H(I,J)-(T(I)*X(J)+X(I)*T(J)-D2*T(I)*T(J))/D1
1616 NEXT J: REM above is main update equation
1620 NEXT I
1624 GOTO 1176
1628 REM PRINT "UPDATE NOT POSSIBLE"
1632 REM PRINT #3,"UPDATE NOT POSSIBLE"
1636 LET I3=1: REM &
1640 GOTO 1140: REM return to try steepest descent
1644 RETURN: REM end of routine
1648 LET I3=0: REM loss function evaluation -- reset flag
1652 GOSUB 2000: REM user function subroutine
1656 LET I9=I9+1
1660 RETURN
1664 IF J8=0 THEN RETURN: REM no parameter display
1668 IF I9<J7*J8 THEN RETURN: REM check if to be printed
1672 LET J7=INT(I9/J8)+1: REM increment parameter display control
1676 REM PRINT "parameters";
1680 REM PRINT #3,"parameters";
1684 FOR J=1 TO N
1688 LET Q$=""
1692 IF B(J)-O(J,1)<=E9*(ABS(O(J,1))+E9) THEN LET Q$="L": REM lower bd
1696 IF O(J,2)-B(J)<=E9*(ABS(O(J,2))+E9) THEN LET Q$="U": REM upper bd
1700 IF O(J,3)=0 THEN LET Q$="M": REM masked 
1704 REM PRINT "  ";B(J);Q$;
1708 REM PRINT #3,"  ";B(J);Q$;
1712 IF 5*INT(J/5)<>J THEN 1732
1716 REM PRINT
1720 REM PRINT "          ";
1724 REM PRINT #3,
1728 REM PRINT #3,"          ";
1732 NEXT J
1736 REM PRINT
1740 REM PRINT #3,
1744 RETURN
2000 LET I3=0: REM polygon.fn
2010 LET F=0: REM will contain area and the constraints via penalty
2020 FOR L=3 TO L9: REM compute areas one by one
2030 LET Q1=B(L-2): LET Q2=B(L-1): REM two sides
2040 LET Q3=B(L+L8): REM angle since l8=l9-3
2050 LET Q4=Q1*Q2*SIN(Q3): REM area component
2060 LET F=F-Q4
2070 NEXT L
2080 LET Q8=-.5*F: REM records area
2090 REM now form the constraints into Y() before adding in penalty
2100 REM lengths from origin don't need constraint, so ignore these
2110 REM also use square of distances to avoid sqrt calculation
2115 LET L=0: REM constraint counter
2116 LET Q5=0: REM penalty term
2117 LET Q6=0: REM set 1 for index of highest violated constraint
2120 FOR L1=2 TO L9-1: REM first counter
2130 FOR L2=L1+1 TO L9: REM second counter
2140 LET L=L+1
2150 LET Q1=B(L1-1): REM ra
2160 LET Q2=B(L2-1): REM rb
2170 LET Q3=0: REM angle
2180 FOR L3=L1+1 TO L2
2190 LET Q3=Q3+B(L3+L8): REM add up angles
2200 NEXT L3
2210 LET Q4=Q1*Q1+Q2*Q2-2*Q1*Q2*COS(Q3): REM square of distance
2220 LET Y(L,1)=1#-Q4: REM slack --> 0 as constraint approached
2230 IF Y(L,1)>0 THEN 2240
2232 LET Q6=L: REM save index if violated
2234 LET F=B9: REM BIG
2236 GOTO 2280
2240 LET Q5=Q5-Q9*LOG(Y(L,1))
2250 NEXT L2
2260 NEXT L1
2270 LET F=F+Q5
2280 IF Q6>0 OR Q8<=QQ0 THEN 2310: REM don't replot ?? May want to anyway
2285 LET QQ0=Q8: REM reset area
2290 GOTO 9000
2300 REM ELSE WANT TO JUST DISPLAY NUMBERS
2310 GOTO 9144
2500 FOR L=1 TO N: REM rem zero gradient components
2510 LET G(L)=0
2530 NEXT L
2532 REM ??? how to ensure Y exists before we use it here
2540 IF Q6>0 THEN STOP: REM don't evaluate gradient if infeasible
2541 REM at least for now
2550 FOR L=3 TO L9: REM compute areas one by one
2560 LET Q1=B(L-2): LET Q2=B(L-1): REM two sides
2565 LET Q3=B(L+L8): REM angle since l8=l9-3
2570 REM LET Q4=Q1*Q2*SIN(Q3): REM area component
2580 LET G(L-2)=G(L-2)-Q2*SIN(Q3)
2584 LET G(L-1)=G(L-1)-Q1*SIN(Q3)
2586 LET G(L+L8)=G(L+L8)-Q1*Q2*COS(Q3)
2588 REM could probably make this more efficient??
2590 NEXT L
2615 LET L=0: REM constraint counter
2616 LET Q5=0: REM penalty term
2617 LET Q6=0: REM set 1 for index of highest violated constraint
2620 FOR L1=2 TO L9-1: REM first counter
2630 FOR L2=L1+1 TO L9: REM second counter
2640 LET L=L+1
2650 LET Q1=B(L1-1): REM ra
2660 LET Q2=B(L2-1): REM rb
2670 LET Q3=0: REM angle
2680 FOR L3=L1+1 TO L2
2690 LET Q3=Q3+B(L3+L8): REM add up angles
2700 NEXT L3
2710 REM LET Q4=Q1*Q1+Q2*Q2-2*Q1*Q2*COS(Q3): REM square of distance
2711 REM LET Y(L,1)=1#-Q4: REM slack --> 0 as constraint approached
2720 IF Y(L,1)<=0 THEN STOP: REM for safety???
2722 LET G(L1-1)=G(L1-1)+2#*Q9*(Q1-Q2*COS(Q3))/Y(L,1)
2723 REM NOTE TERM IS -Q9*LOG(Y) AND Y=1-DIST SO +VE CONTRIBUTION
2724 LET G(L2-1)=G(L2-1)+2#*Q9*(Q2-Q1*COS(Q3))/Y(L,1)
2725 FOR L3=L1+1 TO L2
2726 LET G(L3+L8)=G(L3+L8)+2#*Q9*Q1*Q2*SIN(Q3)/Y(L,1)
2727 NEXT L3: REM count 1 for each included angle
2730 NEXT L2
2740 NEXT L1
2750 RETURN
3000 PRINT "polygon.fn -- area inside plane polygon -- 890413"
3010 REM PRINT #3,"polygon.fn -- area inside plane polygon -- 890413"
3020 LET P$="polygon.fn -- area inside plane polygon"
3030 INPUT "number of vertices = ";L9
3040 IF L9>14 THEN 3030: REM limit of 14 points --> 28-3=25 parameters
3050 LET M=(L9-2)*(L9-1)/2: REM number of distance constraints
3060 PRINT "there are ";M;" distance constraints"
3070 REM PRINT #3,"there are ";M;" distance constraints"
3080 LET N=2*L9-3
3100 PRINT "First point is always at the origin (0,0) cartesian"
3110 PRINT "Second point is at (B(1),0) cartesian or polar"
3120 PRINT "B(1)..B(L9-1) give radial coordinates of points 2..L9"
3130 PRINT "B(L9)..B(2*L9-3) give angle coordinates of points 3..L9"
3135 LET L8=L9-3: REM so l+l8 indexes angles as l=3..l9
3140 PRINT "  since point 2 fixed onto x axis (angular coordinate 0)"
3150 PRINT "Distances between points can be worked out by cosine rule for"
3160 PRINT "triangles i.e. D = sqrt(ra^2 + rb^2 - 2 ra rb cos(angle)"
3170 PRINT
3200 REM PRINT #3, "First point is always at the origin (0,0) cartesian"
3210 REM PRINT #3, "Second point is at (B(1),0) cartesian or polar"
3220 REM PRINT #3, "B(1)..B(L9-1) give radial coordinates of points 2..L9"
3230 REM PRINT #3, "B(L9)..B(2*L9-3) give angle coordinates of points 3..L9"
3240 REM PRINT #3, "  since point 2 fixed onto x axis (angular coordinate 0)"
3250 REM PRINT #3, "Distances between points can be worked out by cosine rule for"
3260 REM PRINT #3, "triangles i.e. D = sqrt(ra^2 + rb^2 - 2 ra rb cos(angle)"
3270 REM PRINT #3,
3300 INPUT "penalty function parameter (zero for no penalty):";Q9
3310 REM PRINT #3,"penalty function parameter (zero for no penalty):";Q9
3320 FOR L=1 TO L9-1
3330 LET O(L,1)=0: LET O(L,2)=1: LET O(L,3)=1: REM radii
3340 NEXT L
3350 FOR L=L9 TO N: REM loop over angles for points 2..l9
3360 LET O(L,1)=0: REM positive
3365 REM??  impose monotonicity via I3 flag -- ?? not necessary -- angles >0
3370 LET O(L,2)=4*ATN(1): REM pi max or it is a reflection
3380 LET O(L,3)=1
3390 NEXT L
3400 REM set parameters to a regular polygon of size .98
3402 LET Q7=.999999
3404 LET Q6=4*ATN(1)/L9: REM pi/no. of vertices
3406 LET Q5=Q7*SIN(Q6): REM q6/l9 = alpha
3408 LET B(1)=Q5: REM radial posn of point 2
3410 LET Q1=Q5: REM x2
3412 LET Q2=0:  REM y2
3414 FOR L=3 TO L9: REM loop over remaining points
3416 LET B(L+L8)=Q6: REM alpha -- polar angle is alpha * (point index - 2)
3418 LET Q1=Q1+Q5*COS(2*(L-2)*Q6)
3420 LET Q2=Q2+Q5*SIN(2*(L-2)*Q6): REM could simplify, but not repeated
3422 LET B(L-1)=SQR(Q1*Q1+Q2*Q2):  REM compute distance from origin
3430 NEXT L
3436 LET QQ0=0: REM Max area found yet
3438 SCREEN 2: REM for graphics
3450 RETURN
3455 REM ??? needs
3456 REM 1) gradient*****
3457 REM 2) more and different penalty fns
3458 REM 3) possible Lagrangian formulation
3459 REM 4) graphics
3460 REM 5) title and blurb on problem
3461 REM 6) NIS promos
3462 REM 7) sets of initial conditions e.g. regular polygons, points on arc
3499 REM         **** done
4500 CLS: REM to reset problem
4505 PRINT "reset penalty parameter from :";Q9;"  to  ";
4507 INPUT Q9
4510 REM PRINT #3,"reset penalty parameter to :";Q9
4515 LET QQ0=0: REM reset area to force a replot
4520 RETURN
6000 REM PRINT "POST-SOLUTION ANALYSIS FOR VARIABLE METRIC"
6005 REM PRINT #3,"POST-SOLUTION ANALYSIS FOR VARIABLE METRIC"
6010 REM POSTVM.BAS 860405
6015 REM CALLS
6020 REM     FUNCTION F  -- line 2000
6025 REM     RESIDUAL R1 -- line 3500 (only if J6=1)
6030 REM
6035 REM INPUTS TO THE ROUTINE
6040 REM    X() -- the vector of final parameter estimates
6045 REM    N   -- the number of parameters in the function F(B)
6050 REM    M   -- the number of data points (0 for a .FN)
6055 REM    O( , ) -- masks and bounds
6060 REM    Y( , ) -- the data array for the problem
6065 REM    A( , ) -- the Hessian inverse approximation (from VM)
6070 REM    G() -- the gradient at X() (from VM)
6075 REM    I5  -- the number of masked (fixed) parameters
6080 REM    J6  -- flag which is 1 if residuals can be computed
6085 REM    F0  -- value of the function at the supposed minimum
6090 REM
6095 REM    B9, E5, E9 -- constants provided by ENVRON
6100 REM
6105 REM OUTPUT FROM THE ROUTINE
6110 REM    F1 -- lower function value than F0 (if found, else a
6115 REM          higher value to prevent a restart in the DRIVER)
6120 REM    X() -- new parameter values if a lower point is found
6125 REM
6130 REM need array T(N) to store tilt angles
6135 IF M>N-I5 THEN 6150
6140 LET F9=-1
6145 GOTO 6155
6150 LET F9=F0/(M-N+I5)
6155 REM PRINT "LOSS FUNCTION PER DATA POINT (SIGMA^2) =";F9
6160 REM PRINT #3,"LOSS FUNCTION PER DATA POINT (SIGMA^2) =";F9
6165 REM PRINT
6170 REM PRINT #3,
6175 PRINT "SOLUTION WITH DISPERSION MEASURES AND GRADIENT OF LOSS FUNCTION"
6180 REM PRINT #3,"SOLUTION WITH DISPERSION MEASURES AND GRADIENT OF LOSS FUNCTION"
6185 REM PRINT
6190 REM PRINT #3,
6195 FOR J=1 TO N: REM includes masked parameters
6200 LET T(J)=A(J,J)*F9
6205 REM PRINT "B(";J;") = ";B(J);" DISPERSION = ";T(J);"  GRAD(";J;") = ";G(J)
6210 REM PRINT #3,"B(";J;") = ";B(J);" DISPERSION = ";T(J);"  GRAD(";J;") = ";G(J)
6215 NEXT J
6220 REM PRINT
6225 REM PRINT #3,
6230 LET T5=1/E9: REM used for active bounds ('large' F(B))
6235 LET I6=0: REM flag for active bounds (-2=l.b.,-1=u.b.,0=inactive)
6240 FOR J=1 TO N
6245 LET E8=E5*SQR(E9): REM now compute points left and right
6250 LET E8=E8*(ABS(X(J))+E8): REM stepsize
6255 REM PRINT "B(";J;")=";
6260 REM PRINT #3,"B(";J;")=";
6265 IF O(J,3)=0 THEN 6480: REM masked
6270 LET B(J)=X(J)-E8
6275 IF B(J)<O(J,1) THEN 6305
6280 LET I3=0
6285 GOSUB 2000: REM function value
6290 IF I3=1 THEN 6310
6295 IF F<F0 THEN 6650: REM lower function value
6300 GOTO 6315
6305 LET I6=-2: REM active lower bound
6310 LET F=T5: REM failure
6315 LET F3=F
6320 LET B(J)=X(J)+E8
6325 IF B(J)>O(J,2) THEN 6355
6330 LET I3=0
6335 GOSUB 2000: REM function value
6340 IF I3=1 THEN 6360
6345 IF F<F0 THEN 6650: REM lower function value
6350 GOTO 6365
6355 LET I6=-1: REM active upper bound
6360 LET F=T5: REM failure
6365 LET F4=F
6370 LET B(J)=X(J): REM reset parameter value
6375 REM LET C1=.5*(F4-F3)/E8: REM linear term -- should be zero
6380 REM LET C2=(F4+F3-2*F0)/(2*E8*E8): REM quadratic term
6385 REM IF I6=-2 THEN LET C2=(F4-F0)/(E8*E8): REM symmetrized curvature
6390 REM IF I6=-1 THEN LET C2=(F3-F0)/(E8*E8): REM when bounds active
6395 REM LET C0=1+C1*C1: REM denominator for curvature
6400 REM LET C2=C2/(C0*SQR(C0))
6405 REM LET C0=B9: REM set large in case of singularity
6410 REM IF C2<>0 THEN LET C0=1/C2: REM radius of curvature for problem
6415 REM LET T(J)=-45*ATN(C1)/ATN(1): REM tilt -- DIMension T in HJ, NM
6420 REM PRINT USING "##.######^^^^";B(J);: REM !!
6425 REM PRINT #3,USING "##.######^^^^";B(J);: REM !!
6430 REM PRINT "  step=";
6435 REM PRINT #3,"  step=";
6440 REM PRINT USING "##.##^^^^";E8;: REM !!
6445 REM PRINT #3,USING "##.##^^^^";E8;: REM !!
6450 REM PRINT "  f-, f+";
6455 REM PRINT #3,"  f-, f+";
6460 REM PRINT USING " ##.####^^^^";F3;F4: REM !!
6465 REM PRINT #3,USING " ##.####^^^^";F3;F4: REM !!
6470 REM LET X(J)=C0: REM save radius of curvature in x()
6475 GOTO 6500
6480 REM PRINT USING "##.######^^^^";B(J);: REM !!
6485 REM PRINT " MASKED"
6490 REM PRINT #3,USING "##.######^^^^";B(J);: REM !!
6495 REM PRINT #3," MASKED"
6500 NEXT J
6505 REM PRINT "        best function value found is            ";
6510 REM PRINT #3,"        best function value found is            ";
6515 REM PRINT USING " ##.####^^^^";F0: REM !!
6520 REM PRINT #3,USING " ##.####^^^^";F0: REM !!
6525 REM PRINT
6530 REM PRINT #3,
6535 REM PRINT "radii of curvature for surface along axial directions"
6540 REM PRINT #3,"radii of curvature for surface along axial directions"
6545 REM PRINT "  & tilt angle in degrees"
6550 REM PRINT #3,"  & tilt angle in degrees"
6555 REM PRINT
6560 REM PRINT #3,
6565 FOR J=1 TO N
6570 IF O(J,3)=0 THEN 6615
6575 REM PRINT "for B(";J;")  R. OF CURV. = ";
6580 REM PRINT #3,"for B(";J;")  R. OF CURV. = ";
6585 REM PRINT USING "##.###^^^^";X(J);: REM !!
6590 REM PRINT #3,USING "##.###^^^^";X(J);: REM !!
6595 REM PRINT "   tilt = ";
6600 REM PRINT #3,"   tilt = ";
6605 REM PRINT USING "####.#####";T(J): REM !!
6610 REM PRINT #3,USING "####.#####";T(J): REM !!
6615 NEXT J
6625 FOR J=1 TO N
6630 LET X(J)=B(J): REM reset parameter for restart
6635 NEXT J
6640 LET F1=F0+ABS(F0)+1: REM to ensure no restart
6645 RETURN: REM to DRIVER
6650 PRINT
6655 REM PRINT #3 ,
6660 PRINT "lower function value found"
6665 REM PRINT #3,"lower function value found"
6675 FOR K=1 TO J
6680 LET X(K)=B(K): REM reset parameter values -- note, j included
6685 NEXT K
6690 LET F1=F: REM to indicate lower function value
6695 RETURN
7000 REM ENVRON -- SUBROUTINE TO DETERMINE MACHINE PRECISION
7010 REM                AND SET SOME STANDARD VARIABLES
7020 REM INPUTS:  NONE
7030 REM
7040 REM OUTPUTS:
7050 REM    B9 -- a large number (currently 1E+35)
7060 REM    E5 -- a number for relative equality tests (currently 10)
7070 REM    E9 -- the machine precision,  E9=MIN (X : 1+X>1)
7080 REM    E6 -- the radix of arithmetic
7090 REM    J1 -- the number of radis digits in mantissa of
7100 REM          floating-point numbers
7110 REM
7120 LET D1=1: REM use a variable to avoid constants when double
7130 REM       precision is invoked
7140 LET E5=10: REM arbitrary scaling for additive equality tests
7150 LET B9=1E+35: REM big number, not necessarily biggest possible
7160 LET E6=1: REM initial value for radix 
7170 LET E9=1: REM initial value for machine precision
7180 LET E9=E9/2: REM start of loop to decrease estimated machine precision
7190 LET D0=E6+E9: REM force storage of sum into a floating-point scalar
7200 IF D0>E6 THEN 7180: REM repeat reduction while (1+E9) > 1
7210 LET E9=E9*2: REM restore smallest E9 which gives (1+E9) > 1
7220 LET E6=E6+1: REM try different radix values
7230 LET D0=E6+E9
7240 IF D0>E6 THEN 7220: REM until a shift is observed
7250 LET J1=1: REM initial count of radix digits in mantissa
7260 LET E9=1: REM use radix for exact machine precision
7270 LET E9=E9/E6: REM loop while dividing by radix
7280 LET J1=J1+1: REM increment counter
7290 LET D0=D1+E9: REM add tp 1
7300 IF D0>D1 THEN 7270: REM test and repeat until equality
7310 LET E9=E9*E6: REM recover last value of machine precision
7320 LET J1=J1-1: REM and adjust the number of digits
7330 RETURN
9000 CLS: REM clear screen before plotting polygon
9004 REM do we need to set the mode (or can we??)
9008 REM set up the actual points in space
9012 LET PP(1,1)=0: LET PP(1,2)=0: REM point 1
9016 LET PP(2,1)=B(1): LET PP(2,2)=0: REM point 2
9020 LET AA=0
9024 FOR LL1=3 TO L9
9028 LET RR=B(LL1-1): REM radius
9032 LET AA=AA+B(L8+LL1): REM angle
9036 REM PRINT "point ";LL1;" radius & angle =";RR;",";45*AA/ATN(1)
9040 LET PP(LL1,1)=RR*COS(AA): LET PP(LL1,2)=RR*SIN(AA)
9044 NEXT LL1: REM points now drawn
9048 REM scaling factors for plot
9049 REM VV, WW are vertical and horizontal nos. of pixels for full screen
9050 REM VS, HS are final scaling factors for horizontal and vertical scaling
9051 REM VW is the aspect ratio (varies with the screen)
9052 LET VW=.69 : REM width/height ratio for particular screen
9053 REM VW is .46 for L Miotti AT, .69 for Zenith
9054 REM Increase VW to make polygon "taller", decrease to make "fatter"
9057 LET VV=200: LET WW=640: REM cga
9058 REM LET VV=348: LET WW=720: REM herc
9059 LET VS=VV*140/200: LET HS=WW*140!/(VW*640!): REM adjustable to screen
9060 LET VO=INT(.7*VV): REM from top
9064 LET WO=INT(.2*WW): REM from left
9072 LET CC=7: REM line colour (default is green for OK)
9076 FOR LL1=1 TO L9
9080 LET SS(LL1,1)=INT(WO+HS*PP(LL1,1))
9084 LET SS(LL1,2)=INT(VO-VS*PP(LL1,2))
9088 REM PRINT P(LL1,1);",";PP(LL1,2);TAB(40);SS(LL1,1);",";SS(LL1,2)
9092 NEXT LL1: REM end scaling
9096 REM LET Z$=INKEY$
9100 REM IF Z$<>"P" AND Z$<>"p" THEN 9092
9104 REM CLS: SCREEN 2
9108 LINE (SS(1,1),SS(1,2))-(SS(2,1),SS(2,2)),CC
9112 FOR LL1=3 TO L9
9116 LINE -(SS(LL1,1),SS(LL1,2)),CC
9120 NEXT LL1
9124 LINE -(SS(1,1),SS(1,2)),CC
9128 LET CC=CC+1: IF CC<16 THEN 9104
9132 REM LET Z$=INKEY$
9136 REM IF Z$<>"P" AND Z$<>"p" THEN 9124
9140 REM SCREEN 2: SCREEN 0: CLS
9144 LOCATE 19,1: PRINT "AREA=";QQ0
9148 LOCATE 20,1: PRINT "evaluations: fn = ";I9;"  grad =";I8;"  fn=";F0;
9149 PRINT "       "
9152 RETURN
