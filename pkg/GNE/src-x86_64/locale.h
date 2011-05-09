/** 
 * @file  locale.h
 * @brief header file for error messages
 *
 * @author Christophe Dutang
 *
 *
 * Copyright (C) 2009, Christophe Dutang. 
 * All rights reserved.
 *
 */



/* Localization */
#include <R.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("GNE", String)
#else
#define _(String) (String)
#endif

