#Example from: Gradient Projection Algorithms and Software for 
#  Arbitrary Rotation Criteria in Factor Analysis.
#  by Coen A. Bernaards and Robert I. Jennrich
#  Website: http://www.stat.ucla.edu/research

 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()


fuzz <- 1e-6 
all.ok <- TRUE  

#  quartimax (orthogonal) rotation of Harman's 8 physical variables.

data("Harman", package="GPArotation")

qHarman  <- GPForth(Harman8, Tmat=diag(2), method="quartimax")
qHarman2 <- quartimax(Harman8) 

 if( fuzz < max(abs(qHarman$Lh - qHarman2$loadings))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qHarman$Lh, digits=18)
    cat("difference:\n")
    print(qHarman$Lh - qHarman2$loadings, digits=18)
    all.ok <- FALSE  
    } 

#qHarman$Th - qHarman2$Th

 tst <- t(matrix(c(
  0.898754567491920398, 0.194823580226859222,
  0.933943406208487592, 0.129748657024604030,
  0.902131483644799892, 0.103864268239045668,
  0.876508251941102934, 0.171284220753554678,
  0.315572019798302239, 0.876476069451083251,
  0.251123191235179066, 0.773488941629975613,
  0.198007116064591759, 0.714678376605717203,
  0.307857241091366252, 0.659334451631046314
  ), 2, 8))
 
 if( fuzz < max(abs(qHarman$Lh - tst ))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qHarman$Lh, digits=18)
    cat("difference:\n")
    print(qHarman$Lh - tst, digits=18)
    all.ok <- FALSE  
    } 

 tst <- t(matrix(c(
   0.790828307905322436, 0.612038060430562525,
  -0.612038060430562525, 0.790828307905322214
  ), 2, 2))
 
 if( fuzz < max(abs(qHarman$Th - tst ))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qHarman$Th, digits=18)
    cat("difference:\n")
    print(qHarman$Th - tst, digits=18)
    all.ok <- FALSE  
    } 





#  quartimin (oblique) rotation of Harman's 8 physical variables.

qminHarman  <- GPFoblq(Harman8, Tmat=diag(2), method="quartimin")
qminHarman2 <- quartimin(Harman8) 

 if( fuzz < max(abs(qminHarman$Lh - qminHarman2$loadings))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qminHarman$Lh, digits=18)
    cat("difference:\n")
    print(qminHarman$Lh - qminHarman2$loadings, digits=18)
    all.ok <- FALSE  
    } 


 tst <- t(matrix(c(
   0.8918217697289939627,  0.0560146456758183961,
   0.9536799985772628219, -0.0232460005406671701,
   0.9291498623396581280, -0.0465027396531852502,
   0.8766828510822184395,  0.0336582451338717017,
   0.0136988312985193428,  0.9250013826349388069,
  -0.0172668087945964319,  0.8212535444941218010,
  -0.0524468998178311899,  0.7649536381341245361,
   0.0858880630098148856,  0.6831160953442911854
   ),2, 8))				       
					       
 if( fuzz < max(abs(qminHarman$Lh - tst ))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qminHarman$Lh, digits=18)
    cat("difference:\n")
    print(qminHarman$Lh - tst, digits=18)
    all.ok <- FALSE  
    } 


 tst <- t(matrix(c(
  1.000000000000000000, 0.472747617396915065,
  0.472747617396915065, 1.000000000000000000
   ),2, 2))				       
					       
 if( fuzz < max(abs(qminHarman$Phi - tst ))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qminHarman$Phi, digits=18)
    cat("difference:\n")
    print(qminHarman$Phi - tst, digits=18)
    all.ok <- FALSE  
    } 


 tst <- t(matrix(c(
   0.878125245495924522, 0.836723841642554422,
  -0.478430823863515542, 0.547625065922776710
   ),2, 2))				       
					       
 if( fuzz < max(abs(qminHarman$Th - tst ))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(qminHarman$Th, digits=18)
    cat("difference:\n")
    print(qminHarman$Th - tst, digits=18)
    all.ok <- FALSE  
    } 


cat("tests completed.\n")

if (! all.ok) stop("some tests FAILED")

