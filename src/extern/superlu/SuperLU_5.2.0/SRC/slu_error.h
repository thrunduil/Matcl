/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/** @file slu_util.h
 * \brief Utility header file 
 *
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 *
 */

#ifndef __SUPERLU_ERROR /* allow multiple inclusions */
#define __SUPERLU_ERROR

#include "superlu_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
#ifndef __STDC__
#include <malloc.h>
#endif
*/
#include <assert.h>
#include "superlu_enum_consts.h"

/***********************************************************************
 * Macros
 ***********************************************************************/

#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#endif /* __SUPERLU_UTIL */
