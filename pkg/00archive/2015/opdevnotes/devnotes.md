      o lbfgsb3 explanation poor. No convcode at the moment. 140901

      o Rvmmin gives warning when called with no gradient. 
        Use control "usenumDeriv=TRUE" to overcome this. 130411

      o tn and tnbc do not have proper controls yet for iterations and
	trace

      o The bdmsk vector is set clumsily for Rvmmin.

      o Where should input checks be done? We should be able to leave them
	out more easily to make code more efficient.

      o Get some routines out of the "all" set. Create a new "ALL+" category
	that is really all the routines, and leave "all" as all reasonable ones.
        See 2013-08-02 note about UOBYQA.

      o May want to tweak bobyqa, newuoa default controls.

      o Sets of control defaults so we can get different examples to show up
	properly, and also to allow "soft" vs. "hard" tries of problems e.g.,
	with tighter or looser iteration and tolerance settings.
        Note changes 130806 to 140504 in maxfeval lines (more in 13 than 14??)
	Need to standardize, then test. May not be able to fully reproduce 
	results.

      o Simplify, simplify, simplify

VERSION 2014-??-??

      o  Rvmmin failure in some cases traced to mcontrol having extra elements.
	This applies to other methods too, and will need constant attention.

      o   trace(optimx, exit = quote(print(output)))
          now works for optimx. Thanks to Gabor Grothendieck.

      o Need to add tn and tnbc (Rtnmin choice)

      o L-BFGS-B version 3 (Fortran) implemented and available as choice
	"lbfgsb3". Note that there is an apparent glitch in the previous
	versions, according to J.L. Morales and J. Nocedal.
        L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN routines 
        for large scale bound constrained optimization (2011, November), 
        ACM Transactions on Mathematical Software, volume 38, 
        pages 7:1 to 7:4.
