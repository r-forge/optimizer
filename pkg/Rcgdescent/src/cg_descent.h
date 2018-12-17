#ifndef CG_DESCENT_H_
#define CG_DESCENT_H_

#include <cmath>
#include <climits>
#include <cfloat>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cstdio>

const double ZERO = 0.0;
const double ONE = 1.0;

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

typedef struct cg_com_struct /* common variables */
{
    /* parameters computed by the code */
    long              n ; /* problem dimension, saved for reference */
    long             nf ; /* number of function evaluations */
    long             ng ; /* number of gradient evaluations */
    int         QuadOK ; /* T (quadratic step successful) */
    int       UseCubic ; /* T (use cubic step) F (use secant step) */
    int           neps ; /* number of time eps updated */
    int       PertRule ; /* T => estimated error in function value is eps*Ck,
                            F => estimated error in function value is eps */
    int          QuadF ; /* T => function appears to be quadratic */
    double   SmallCost ; /* |f| <= SmallCost => set PertRule = F */
    double       alpha ; /* stepsize along search direction */
    double           f ; /* function value for step alpha */
    double          df ; /* function derivative for step alpha */
    double       fpert ; /* perturbation is eps*|f| if PertRule is T */
    double         eps ; /* current value of eps */
    double         tol ; /* computing tolerance */
    double          f0 ; /* old function value */
    double         df0 ; /* old derivative */
    double          Ck ; /* average cost as given by the rule:
                            Qk = Qdecay*Qk + 1, Ck += (fabs (f) - Ck)/Qk */
    double    wolfe_hi ; /* upper bound for slope in Wolfe test */
    double    wolfe_lo ; /* lower bound for slope in Wolfe test */
    double   awolfe_hi ; /* upper bound for slope, approximate Wolfe test */
    int         AWolfe ; /* F (use Wolfe line search)
                                T (use approximate Wolfe line search)
                                do not change user's AWolfe, this value can be
                                changed based on AWolfeFac */
    int          Wolfe ; /* T (means code reached the Wolfe part of cg_line */
    double         rho ; /* either Parm->rho or Parm->nan_rho */
    double    alphaold ; /* previous value for stepsize alpha */
    double          *x ; /* current iterate */
    double      *xtemp ; /* x + alpha*d */
    double          *d ; /* current search direction */
    double          *g ; /* gradient at x */
    double      *gtemp ; /* gradient at x + alpha*d */
  std::function<double(double *, long)> cg_value;
  //    double   (*cg_value) (double *, long) ; /* f = cg_value (x, n) */
  std::function<void(double *, double *, long)> cg_grad;
  //    void      (*cg_grad) (double *, double *, long) ; /* cg_grad (g, x, n) */
  std::function<double(double *, double *, long)> cg_valgrad;
  //    double (*cg_valgrad) (double *, double *, long) ; /* f = cg_valgrad (g,x,n)*/
    cg_parameter *Parm ; /* user parameters */
} cg_com ;

/* prototypes */

static int cg_Wolfe
(
    double   alpha, /* stepsize */
    double       f, /* function value associated with stepsize alpha */
    double    dphi, /* derivative value associated with stepsize alpha */
    cg_com    *Com  /* cg com */
) ;

static int cg_tol
(
    double     gnorm, /* gradient sup-norm */
    cg_com    *Com    /* cg com */
) ;

static int cg_line
(
    cg_com   *Com  /* cg com structure */
) ;

static int cg_contract
(
    double    *A, /* left side of bracketing interval */
    double   *fA, /* function value at a */
    double   *dA, /* derivative at a */
    double    *B, /* right side of bracketing interval */
    double   *fB, /* function value at b */
    double   *dB, /* derivative at b */
    cg_com  *Com  /* cg com structure */
) ;

static int cg_evaluate
(
    const char    *what, /* fg = evaluate func and grad, g = grad only,f = func only*/
    const char     *nan, /* y means check function/derivative values for nan */
    cg_com   *Com
) ;

static double cg_cubic
(
    double  a,
    double fa, /* function value at a */
    double da, /* derivative at a */
    double  b,
    double fb, /* function value at b */
    double db  /* derivative at b */
) ;

static void cg_matvec
(
    double *y, /* product vector */
    double *A, /* dense matrix */
    double *x, /* input vector */
    int     n, /* number of columns of A */
    long     m, /* number of rows of A */
    int     w  /* T => y = A*x, F => y = A'*x */
) ;

static void cg_trisolve
(
    double *x, /* right side on input, solution on output */
    double *R, /* dense matrix */
    int     m, /* leading dimension of R */
    int     n, /* dimension of triangular system */
    int     w  /* T => Rx = y, F => R'x = y */
) ;

static double cg_inf
(
    double *x, /* vector */
    long     n /* length of vector */
) ;

static void cg_scale0
(
    double *y, /* output vector */
    double *x, /* input vector */
    double  s, /* scalar */
    int     n /* length of vector */
) ;

static void cg_scale
(
    double *y, /* output vector */
    double *x, /* input vector */
    double  s, /* scalar */
    long     n /* length of vector */
) ;

static void cg_daxpy0
(
    double     *x, /* input and output vector */
    double     *d, /* direction */
    double  alpha, /* stepsize */
    int         n  /* length of the vectors */
) ;

static void cg_daxpy
(
    double     *x, /* input and output vector */
    double     *d, /* direction */
    double  alpha, /* stepsize */
    long         n  /* length of the vectors */
) ;

static double cg_dot0
(
    double *x, /* first vector */
    double *y, /* second vector */
    int     n /* length of vectors */
) ;

static double cg_dot
(
    double *x, /* first vector */
    double *y, /* second vector */
    long     n /* length of vectors */
) ;

static void cg_copy0
(
    double *y, /* output of copy */
    double *x, /* input of copy */
    int     n  /* length of vectors */
) ;

static void cg_copy
(
    double *y, /* output of copy */
    double *x, /* input of copy */
    long     n  /* length of vectors */
) ;

static void cg_step
(
    double *xtemp, /*output vector */
    double     *x, /* initial vector */
    double     *d, /* search direction */
    double  alpha, /* stepsize */
    long         n  /* length of the vectors */
) ;

static void cg_init
(
    double *x, /* input and output vector */
    double  s, /* scalar */
    long     n /* length of vector */
) ;

static double cg_update_2
(
    double *gold, /* old g */
    double *gnew, /* new g */
    double    *d, /* d */
    long        n /* length of vectors */
) ;

static double cg_update_inf
(
    double *gold, /* old g */
    double *gnew, /* new g */
    double    *d, /* d */
    long        n /* length of vectors */
) ;

static double cg_update_ykyk
(
    double *gold, /* old g */
    double *gnew, /* new g */
    double *Ykyk,
    double *Ykgk,
    long        n /* length of vectors */
) ;

static double cg_update_inf2
(
    double   *gold, /* old g */
    double   *gnew, /* new g */
    double      *d, /* d */
    double *gnorm2, /* 2-norm of g */
    long          n /* length of vectors */
) ;

static double cg_update_d
(
    double      *d,
    double      *g,
    double    beta,
    double *gnorm2, /* 2-norm of g */
    long          n /* length of vectors */
) ;

static void cg_Yk
(
    double    *y, /*output vector */
    double *gold, /* initial vector */
    double *gnew, /* search direction */
    double  *yty, /* y'y */
    long        n  /* length of the vectors */
) ;

/* static void cg_printParms */
/* ( */
/*     cg_parameter  *Parm */
/* ) ; */

#endif // CG_DESCENT_H_
