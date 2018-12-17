#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif  // R_NO_REMAP

#include <Rinternals.h>
#include <R_ext/Parse.h>

#include "cg_user.h"
#include <ctime>
#include <string>
#include <sstream>

using std::endl;

namespace CGMIN {

  // This class was taken from the BOOM library header boom_r_tools.hpp.
  // Adapted under the LGPL license.
  // 
  // A class to handle protection of R objects in an exception safe
  // way.  Define a single RMemoryProtector at the start of any
  // function where R memory needs to be allocated.  Then instead of
  // writing
  // PROTECT(my_r_object);
  // /* do stuff */
  // UNPROTECT(1);
  //
  // you write
  // RMemoryProtector protector;
  // protector.protect(my_r_object);
  //  /* do stuff */
  //
  // There is no need to call UNPROTECT, which is handled by the class
  // destructor.
  class RMemoryProtector {
   public:
    RMemoryProtector() : protection_count_(0) {}

    ~RMemoryProtector() {
      UNPROTECT(protection_count_);
    }

    // Args:
    //   r_object: An R object that needs protecting for the life of
    //     this RMemoryProtector object.
    // Returns:
    //   The protected r_object.
    SEXP protect(SEXP r_object) {
      PROTECT(r_object);
      ++protection_count_;
      return r_object;
    }

   private:
    int protection_count_;
  };

  //===========================================================================
  // Borrowed from the Boom library under the LGPL license.
  // Get an item from a list.
  //
  // Args:
  //   list:  An R list containing named elements.
  //   name:  The name of the element in 'list' to be returned.
  //   expect_answer: If true, then an exception will be thrown if list[[name]]
  //     is not present.
  //
  // Returns:
  //   If the requested element is present, it is returned.  Otherwise
  //   R_NilValue is returned if expect_answer is false.  An exception is raised
  //   if the requested element is not present and expect_answer is true.
  SEXP getListElement(SEXP list, const std::string &name,
                      bool expect_answer = false);
  
  SEXP getListElement(SEXP list, const std::string &name, bool expect_answer) {
    SEXP elmt = R_NilValue;
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if(Rf_isNull(names)){
      std::ostringstream err;
      err << "Attempt to use getListElement in a list with"
          << " no 'names' attribute." << endl
          << "You were searching for the name: " << name << endl;
      throw std::runtime_error(err.str());
    }
    for(int i = 0; i < Rf_length(list); i++)
      if(name == CHAR(STRING_ELT(names, i))){
        elmt = VECTOR_ELT(list, i);
        break;
      }
    if (expect_answer && elmt == R_NilValue) {
      std::ostringstream warning;
      warning << "Could not find list element named: " << name << endl;
      Rf_warning(warning.str().c_str());
    }
    return elmt;
  }

  //===========================================================================
  // Borrowed from the Boom library under the LGPL license.
  //
  // Convert a C++ string into an R string (a character vector of length 1).
  SEXP ToRString(const std::string &s) {
    SEXP ans;
    RMemoryProtector protector;
    protector.protect(ans = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(ans, 0, Rf_mkChar(s.c_str()));
    return ans;
  }

  // Borrowed from the Boom library under the LGPL license.
  std::string ToString(SEXP r_string) {
    if (TYPEOF(r_string) == CHARSXP) {
      return CHAR(r_string);
    } else if(Rf_isString(r_string)){
      return CHAR(STRING_ELT(r_string, 0));
    } else {
      throw std::runtime_error("ToString could not convert its "
                               "argument to a string");
    }
    return "";
  }

  // Adapted from the Boom library under the LGPL license.
  SEXP ToRVector(const double *x, long n) {
    RMemoryProtector protector;
    SEXP ans = protector.protect(Rf_allocVector(REALSXP, n));
    double *data = REAL(ans);
    for(int i = 0; i < n; ++i) data[i] = x[i];
    return ans;
  }
  
  //===========================================================================
  // A representation of an R vector created from C++.  This is not simply the
  // vector elements, but a named object created in a specific environment.
  class RVector {
   public:
    // Args:
    //   r_environment:  The environment in which to create the vector.
    RVector(SEXP r_environment)
        : r_env_(r_environment),
          name_(""),
          r_symbol_(R_NilValue),
          r_vector_(R_NilValue),
          vector_data_(nullptr),
          size_(0)
    {}

    // Set the R vector to the given values.  If the vector has not been
    // previously set, then calling this method creates a named vector object in
    // the R environment held by this object.  The object will be removed by R's
    // garbage collector when this RVector goes out of scope.
    //
    // If the vector has already been set, then an exception will be raised if
    // the supplied size differs from the size of the object being stored.
    void set(double *x, long n) {
      if (size_ == 0) {
        std::ostringstream argname;
        argname << "x_" << clock();
        name_ = argname.str();
        size_ = n;
        r_vector_ = protector_.protect(Rf_allocVector(REALSXP, size_));
        vector_data_ = REAL(r_vector_);

        protector_.protect(r_symbol_ = Rf_install(name_.c_str()));
        Rf_defineVar(r_symbol_, r_vector_, r_env_);
      }
      if (size_ != n) {
        std::ostringstream err;
        err << "Wrong size assignment to non-empty vector adapter.";
        throw std::runtime_error(err.str());
      }
      for (int i = 0; i < size_; ++i) {
        vector_data_[i] = x[i];
      }
    }

    // A handle to the actual R object.
    SEXP value() {return r_vector_;}

    // The name of the vector as seen by R.
    const std::string &name() const {return name_;}
    
   private:
    // Protects the SEXP's allocated by this object.
    RMemoryProtector protector_;

    // The environment in which the vector is to be created.
    SEXP r_env_;

    // The name of the vector in R.
    std::string name_;

    // The symbol for the vector used by R's internal symbol table.
    SEXP r_symbol_;

    // The vector.
    SEXP r_vector_;

    // A pointer to the entries in r_vector_.
    double *vector_data_;

    // The size of r_vector_.  Equal to Rf_length(r_vector_).
    long size_;
  };

  //===========================================================================
  // Wrap an R function, with the signature 'scalar = f(vector_x)'.  The
  // function is represented by its name and its environment.  The actual
  // function pointer only manipulated through R, so it is not explicitly
  // represented here.
  class RFunctionAdapter {
   public:
    // Args:
    //   r_function_info: An R list containing
    //   - The name of a function.
    //   - The environment in which the function is to be evaluated.
    RFunctionAdapter(SEXP r_function_info)
        : function_name_(ToString(getListElement(
              r_function_info, "function.name", true))),
          r_env_(getListElement(r_function_info, "env")),
          function_argument_(r_env_)
    {}

    // Args:
    //   x, dim: The C representation of the input vector to the function.
    //
    // Returns:
    //   An R function call for evaluating the function at the given values.
    //   Ready to be passed to Rf_eval.
    SEXP function_call(double *x, long dim) {
      RMemoryProtector function_call_protector;
      function_argument_.set(x, dim);
      ParseStatus parse_status = PARSE_NULL;
      SEXP r_call = function_call_protector.protect(R_ParseVector(
          ToRString(call_string()),
          1,
          &parse_status,
          R_NilValue));
      if (parse_status != PARSE_OK) {
        std::ostringstream err;
        err << "Could not parse expression: " << call_string_;
        throw std::runtime_error(err.str());
      }
      return r_call;
    }

    // The function call as seen by R, e.g. 'f(x)'.
    const std::string &call_string() {
      if (call_string_.empty()) {
        std::ostringstream callstring;
        callstring << function_name_ << "("
                   << function_argument_.name()
                   << ")";
        call_string_ = callstring.str();
      }
      return call_string_;
    }

    // The R environment in which the function call takes place.
    SEXP r_env() {return r_env_;}
    
   private:
    // The protector handles the PROTECT / UNPROTECT stuff for the SEXP's
    // defined by this object.
    RMemoryProtector protector;

    // The name of the R function, from R.
    std::string function_name_;

    // The environment in which the function is to be evaluated.
    SEXP r_env_;

    // The function call to be evaluated, as a string.  E.g. 'f(x)'
    std::string call_string_;

    // Provides a handle to the numerical value of the function argument, and
    // its name.
    RVector function_argument_;
  };

  //---------------------------------------------------------------------------
  // A function adapter for the case of scalar.y <- f(vector.x)
  // Suitable for the target function in the min problem.
  class RScalarFunctionAdapter : public RFunctionAdapter {
   public:
    RScalarFunctionAdapter(SEXP r_scalar_valued_function_info)
        : RFunctionAdapter(r_scalar_valued_function_info) {}
    
    double operator()(double *x, long dim) {
      return Rf_asReal(Rf_eval(VECTOR_ELT(function_call(x, dim), 0),
                               r_env()));
    }
  };

  //---------------------------------------------------------------------------
  // A function adapter for the case of vector.y <- f(vector.x).
  // Suitable for the gradient in the min problem.
  class RVectorFunctionAdapter : public RFunctionAdapter {
   public:
    RVectorFunctionAdapter(SEXP r_vector_valued_function_info)
        : RFunctionAdapter(r_vector_valued_function_info) {}
    void operator()(double *gradient, double *x, long dim) {
      RMemoryProtector protector;
      SEXP r_gradient = protector.protect(Rf_eval(
          VECTOR_ELT(function_call(x, dim), 0),
          r_env()));
      if (gradient) {
        double *rgrad_data = REAL(r_gradient);
        for (long i = 0; i < dim; ++i) {
          gradient[i] = rgrad_data[i];
        }
      }
    }
  };

  //===========================================================================
  // Translate the list of control parameters supplied by R to the cg_parameter
  // struct.
  void FillControlParams(SEXP r_cg_param_list, cg_parameter *param) {
    param->PrintFinal = false;
    param->PrintLevel = 0;
    param->PrintParms = false;
    param->LBFGS = Rf_asInteger(getListElement(r_cg_param_list, "LBFGS", true));
    param->memory = Rf_asInteger(getListElement(r_cg_param_list, "memory", true));
    param->SubCheck = Rf_asInteger(getListElement(r_cg_param_list, "SubCheck", true));
    param->SubSkip = Rf_asInteger(getListElement(r_cg_param_list, "SubSkip", true));
    param->eta0 = Rf_asReal(getListElement(r_cg_param_list, "eta0", true));
    param->eta1 = Rf_asReal(getListElement(r_cg_param_list, "eta1", true));
    param->eta2 = Rf_asReal(getListElement(r_cg_param_list, "eta2", true));
    param->AWolfe = Rf_asInteger(getListElement(r_cg_param_list, "AWolfe", true));
    param->AWolfeFac = Rf_asReal(getListElement(r_cg_param_list, "AWolfeFac", true));
    param->Qdecay = Rf_asReal(getListElement(r_cg_param_list, "Qdecay", true));
    param->nslow = Rf_asInteger(getListElement(r_cg_param_list, "nslow", true));
    param->StopRule = Rf_asLogical(getListElement(r_cg_param_list, "StopRule", true));
    param->StopFac = Rf_asReal(getListElement(r_cg_param_list, "StopFac", true));
    param->PertRule = Rf_asLogical(getListElement(r_cg_param_list, "PertRule", true));
    param->eps = Rf_asReal(getListElement(r_cg_param_list, "eps", true));
    param->egrow = Rf_asReal(getListElement(r_cg_param_list, "egrow", true));
    param->QuadStep = Rf_asLogical(getListElement(r_cg_param_list, "QuadStep", true));
    param->QuadCutOff = Rf_asReal(getListElement(r_cg_param_list, "QuadCutOff", true));
    param->QuadSafe = Rf_asReal(getListElement(r_cg_param_list, "QuadSafe", true));
    param->UseCubic = Rf_asLogical(getListElement(r_cg_param_list, "UseCubic", true));
    param->CubicCutOff = Rf_asReal(getListElement(r_cg_param_list, "CubicCutOff", true));
    param->SmallCost = Rf_asReal(getListElement(r_cg_param_list, "SmallCost", true));
    param->debug = Rf_asLogical(getListElement(r_cg_param_list, "debug", true));
    param->debugtol = Rf_asReal(getListElement(r_cg_param_list, "debugtol", true));
    param->step = Rf_asReal(getListElement(r_cg_param_list, "step", true));
    param->maxit = Rf_asInteger(getListElement(r_cg_param_list, "maxit", true));
    param->ntries = Rf_asInteger(getListElement(r_cg_param_list, "ntries", true));
    param->ExpandSafe = Rf_asReal(getListElement(r_cg_param_list, "ExpandSafe", true));
    param->SecantAmp = Rf_asReal(getListElement(r_cg_param_list, "SecantAmp", true));
    param->RhoGrow = Rf_asReal(getListElement(r_cg_param_list, "RhoGrow", true));
    param->neps = Rf_asInteger(getListElement(r_cg_param_list, "neps", true));
    param->nshrink = Rf_asInteger(getListElement(r_cg_param_list, "nshrink", true));
    param->nline = Rf_asInteger(getListElement(r_cg_param_list, "nline", true));
    param->restart_fac = Rf_asReal(getListElement(r_cg_param_list, "restart_fac", true));
    param->nan_rho = Rf_asReal(getListElement(r_cg_param_list, "nan_rho", true));
    param->nan_decay = Rf_asReal(getListElement(r_cg_param_list, "nan_decay", true));
  }

  // The return value of cg_descent is an integer code.  Translate it to a human
  // readable string.  
  std::string TranslateReturnStatus(int code) {
    std::string return_status;
    switch (code) {
      case -2:  
        return_status = " (function value became nan)";
        break;
      case -1:
        return_status = "(starting function value is nan)";
        break;
      case 0:
        return_status = "(convergence tolerance satisfied)";
        break;
      case 1:
        return_status = "(change in func <= feps*|f|)";
        break;
      case 2:
        return_status = "(total iterations exceeded maxit)";
        break;
      case 3:
        return_status = "(slope always negative in line search)";
        break;
      case 4:
        return_status = "(number secant iterations exceed nsecant)";
        break;
      case 5:
        return_status = "(search direction not a descent direction)";
        break;
      case 6:
        return_status = "(line search fails in initial interval)";
        break;
      case 7:
        return_status = "(line search fails during bisection)";
        break;
      case 8:
        return_status = "(line search fails during interval update)";
        break;
      case 9:
        return_status = "(debugger is on and the function value increases)";
        break;
      case 10:
        return_status = "(out of memory)";
        break;
      default:
        return_status = "unrecognized return status";
    }
    return return_status;
  }
}  // namespace CGMIN


extern "C" {

  // Args:
  //   Rfun: The function to be optimized.  The function's signature should be
  //     Rfun(x, ...), where x is a numeric vector, and it should return a
  //     single numeric value.
  //   Rgradient: An R function providing the gradient of Rfun.  The function's
  //     signature should be grad <- Rgradient(x, ...).  The return value is the
  //     gradient of Rfun at (x, ...).
  //   Rinitial_values: A vector of initial values where the optimization should
  //     begin.
  //
  // Returns:
  //   A list containing the optimal vector, the function value that it
  //   produces, and a collection of measurements describing the optimization
  //   process.
  SEXP cgminu_wrapper(SEXP Rfun,
                      SEXP Rgradient,
                      SEXP Rinitial_values,
                      SEXP Rcontrol) {
    CGMIN::RMemoryProtector protector;
    try {
      CGMIN::RScalarFunctionAdapter target_fun(Rfun);
      CGMIN::RVectorFunctionAdapter gradient_fun(Rgradient);
      // Passing an empty valgrad to cgmin tells the code that we don't have a
      // separate valgrad function.  target_fun and gradient_fun are required, but
      // valgrad is optional.
      std::function<double(double *, double*, long)> valgrad;
    
      int dimension = Rf_length(Rinitial_values);
      double *data = REAL(Rinitial_values);
      cg_stats stats;
      cg_parameter params;

      // It might make sense to put work and grad_tol in the control params, but
      // they are hard coded for now.
      double *work = NULL;
      double grad_tol = 1e-8;

      // Some control params are not to be touched even by experts.  Calling
      // cg_default sets all the parameters to their default values, even
      // (especially) the expert ones.
      cg_default(&params);

      // Then the parameters suitable for user fiddling can be set from the
      // RControl list.
      CGMIN::FillControlParams(Rcontrol, &params);

      // This is where the work gets done.
      int return_code = cg_descent(
          data,
          dimension,
          &stats,
          &params,
          grad_tol,
          target_fun,
          gradient_fun,
          valgrad,
          work);
      // Translate the return status code to a human readable string.
      std::string return_status = CGMIN::TranslateReturnStatus(return_code);

      // Package everything up and return to the user.
      //
      // IMPORTANT: The names for this list are set in the R function that calls
      // this code.  If you add return values here, be sure to adjust the names
      // in that R function accordingly.
      SEXP r_optimal_vector = protector.protect(CGMIN::ToRVector(data, dimension));
      SEXP r_function_optimum = protector.protect(Rf_ScalarReal(stats.f));
      SEXP r_function_evals = protector.protect(Rf_ScalarInteger(stats.nfunc));
      SEXP r_gradient_evals = protector.protect(Rf_ScalarInteger(stats.ngrad));
      SEXP r_return_status = protector.protect(CGMIN::ToRString(return_status.c_str()));
      SEXP r_converged = protector.protect(Rf_ScalarLogical(return_code == 0));
      SEXP ans = protector.protect(Rf_allocVector(VECSXP, 6));
      SET_VECTOR_ELT(ans, 0, r_optimal_vector);
      SET_VECTOR_ELT(ans, 1, r_function_optimum);
      SET_VECTOR_ELT(ans, 2, r_function_evals);
      SET_VECTOR_ELT(ans, 3, r_gradient_evals);
      SET_VECTOR_ELT(ans, 4, r_return_status);
      SET_VECTOR_ELT(ans, 5, r_converged);
      return ans;
    } catch (std::exception &e) {
      Rf_error("Caught an exception with the following error message: \n%s",
               e.what());
    } catch (...) {
      Rf_error("Caught an unknown exception");
    }
    return R_NilValue;
  }
  
}  // extern "C"
