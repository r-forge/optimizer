## 0README.txt

This file simply describes a number of examples and test scripts
related to the optimx package. 

Note that even methods that are generally quite good still show
failures in some of these examples.

John C. Nash 2023-06-13

The files are:

3Rosen.R -- this illustrates three different extensions of the famous
            Rosenbrock banana-shaped valley problem which has just 2
            parameters. The three extensions to more parameters give
            different function, gradient and hessian values.

argclash.R -- Illustration of some issues when the objective function
            uses information from dot arguments and the name is the
            same as one that is an argument of a called function, 
            e.g., gradient approximation

axsearch-Ex.R -- a simple illustration of the use and output of the
            axsearch() function to get "tilt" and "radius of curvature"
            (roc) in each parameter direction at a point on a function 
            surface

dropnmbk.R -- A little script to show that method nmbk fails when 
            started with a parameter on a bound.

fhop.R -- To run many methods on some or all of the problems in the
            funconstrain package with upper, lower, both, or no 
            bounds.

hessapptest.R -- to show how to change to origin of the Hessian and
           Jacobian approximations used.

hessian-used.R -- Illustration of use of hessian in optimization with
            four methods from optimx.

HessQuality.R -- tests on funconstrain set to check symmetry of 
           returned Hessian

hobbs.R -- Various applications of methods to the Hobbs Weed problem
           from Nash, Compact numerical methods for computers, 
           Adam Hilger:Bristol, 1979 and Second Edition, IOP: Bristol,
           1990. This is a problem that has some nasty features yet
           is very simple to state.

jonesrun.R -- tests using a function from Owen Jones, Robert Maillardet, 
           and Andrew Robinson, 
           Scientific Programming and Simulation Using R, Second Edition, 
           Chapman & Hall/CRC, Boca Raton Page 230

ncgtests.R -- various tests of the revised R conjugate gradients minimizer
	   ncg by J C Nash

optimrgrapprox.R -- Examples of calls to different gradient approximations

rosenbrock.R -- diverse examples using the genrose extended Rosenbrock
           test function

scalechk-Ex.R -- examples of calls to scalechk() function for checking the
           parameter scaling in optimization problems

simplefuntst.R -- Extremely simple example used to show syntax of calls.

snewtonbtest.R -- bounded tests using snewton methods

specctrlhobbs.R -- tests passing of special controls to bobyqa using Hobbs.

ssqbtest.R -- a simple sum of squares test showing differences between
              opm() and optimx()

testbmstep.R -- Examination of bounded step within optimizers.

trig1507.R -- More' Garbow and Hillstrom, 1981, problem 26 tests.
           Note that this function is also part of the funconstrain
           set and fhop.R will use it also.

trystarttests.R -- testing the starttests option in older optimx()
	   function (deprecated).

woodfn.R -- well-known 4 parameter optimization test problem

woodtest.R -- various tests on the wood function.


