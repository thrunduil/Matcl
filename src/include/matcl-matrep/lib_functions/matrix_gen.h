/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#pragma once

#include "matcl-scalar/lib_functions/scalar_gen.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/matrix/permvec.h"

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 
#pragma warning(disable:4275)	//non dll-interface class used as base for dll-interface class

namespace matcl
{

//-----------------------------------------------------------------
//                      sequences
//-----------------------------------------------------------------
// create sequence of Real numbers from s to e with step 1.0
MATCL_MATREP_EXPORT Matrix  range(Real s, Real e);

// create sequence of Real numbers from s to e with step i
MATCL_MATREP_EXPORT Matrix  range(Real s, Real i, Real e);

// create sequence of Float numbers from s to e with step 1.0
MATCL_MATREP_EXPORT Matrix  frange(Float s, Float e);

// create sequence of Float numbers from s to e with step i
MATCL_MATREP_EXPORT Matrix  frange(Float s, Float i, Float e);

// create sequence of Integer numbers from s to e with step 1
MATCL_MATREP_EXPORT Matrix  irange(Integer s, Integer e);

// create sequence of Integer numbers from s to e with step i
MATCL_MATREP_EXPORT Matrix  irange(Integer s, Integer i, Integer e);

// create sequence of length n of equally spaced Real numbers from s to e
MATCL_MATREP_EXPORT Matrix  linspace(Real s, Real e, Integer n);

// create sequence of length n of equally spaced Float numbers from s to e
MATCL_MATREP_EXPORT Matrix  flinspace(Float s, Float e, Integer n);

// create sequence of length n of logarithmically equally spaced Real numbers 
// from 10^s to 10^e
MATCL_MATREP_EXPORT Matrix  logspace(Real s, Real e, Integer n);

// create sequence of length n of logarithmically equally spaced Float numbers 
// from 10^s to 10^e
MATCL_MATREP_EXPORT Matrix  flogspace(Float s, Float e, Integer n);

//-----------------------------------------------------------------
//                      zeroes
//-----------------------------------------------------------------
// create dense matrix with Real zeroes of size rxc
MATCL_MATREP_EXPORT Matrix  zeros(Integer r, Integer c);

// create dense matrix with Integer zeroes of size rxc
MATCL_MATREP_EXPORT Matrix  izeros(Integer r, Integer c);

// create dense matrix with Float zeroes of size rxc
MATCL_MATREP_EXPORT Matrix  fzeros(Integer r, Integer c);

// create dense matrix with Complex zeroes of size rxc
MATCL_MATREP_EXPORT Matrix  czeros(Integer r, Integer c);

// create dense matrix with Float_complex zeroes of size rxc
MATCL_MATREP_EXPORT Matrix  fczeros(Integer r, Integer c);

// create dense matrix with zeroes of size rxc of type given by value_code vt
MATCL_MATREP_EXPORT Matrix  zeros(Integer r, Integer c, matcl::value_code vt);

// create dense matrix with zeroes of size rxc of type Object with type info ti
MATCL_MATREP_EXPORT Matrix  zeros(ti::ti_object ti,Integer r, Integer c);

// create sparse matrix with Real zeroes of size rxc, with allocated space for nnz 
// elements
MATCL_MATREP_EXPORT Matrix  spzeros(Integer r, Integer c,Integer nnz = 0);

// create sparse matrix with Integer zeroes of size rxc, with allocated space for nnz 
// elements
MATCL_MATREP_EXPORT Matrix  ispzeros(Integer r, Integer c,Integer nnz = 0);

// create sparse matrix with Float zeroes of size rxc, with allocated space for nnz 
// elements
MATCL_MATREP_EXPORT Matrix  fspzeros(Integer r, Integer c,Integer nnz = 0);

// create sparse matrix with Complex zeroes of size rxc, with allocated space for nnz 
// elements
MATCL_MATREP_EXPORT Matrix  cspzeros(Integer r, Integer c,Integer nnz = 0);

// create sparse matrix with Float_complex zeroes of size rxc, with allocated space 
// for nnz elements
MATCL_MATREP_EXPORT Matrix  fcspzeros(Integer r, Integer c,Integer nnz = 0);

// create sparse matrix with zeroes of size rxc of type given by value_code vt, with
// allocated space for nnz elements
MATCL_MATREP_EXPORT Matrix  spzeros(Integer r, Integer c,Integer nnz, matcl::value_code);

// create object matrix with zeroes of size rxc of type Object with type info ti, with
// allocated space for nnz elements
MATCL_MATREP_EXPORT Matrix  spzeros(ti::ti_object ti,Integer r, Integer c,Integer nnz = 0);

// create band matrix with Real zeroes of size rxc, with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  bzeros(Integer r, Integer c,Integer fd, Integer ld);

// create band matrix with Integer zeroes of size rxc, with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  ibzeros(Integer r, Integer c,Integer fd, Integer ld);

// create band matrix with Float zeroes of size rxc, with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  fbzeros(Integer r, Integer c,Integer fd, Integer ld);

// create band matrix with Complex zeroes of size rxc, with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  cbzeros(Integer r, Integer c,Integer fd, Integer ld);

// create band matrix with Float_complex zeroes of size rxc, with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  fcbzeros(Integer r, Integer c,Integer fd, Integer ld);

// create band matrix with zeroes of size rxc of type given by value code, with first 
// diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  bzeros(Integer r, Integer c,Integer fd, Integer ld, matcl::value_code);

// create band matrix with zeroes of size rxc of type Object with type info ti, with 
// first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  bzeros(ti::ti_object ti,Integer r, Integer c,Integer fd, Integer ld);

//-----------------------------------------------------------------
//                      ones
//-----------------------------------------------------------------
// create dense matrix with Real ones of size rxc
MATCL_MATREP_EXPORT Matrix  ones(Integer r, Integer c);

// create dense matrix with Integer ones of size rxc
MATCL_MATREP_EXPORT Matrix  iones(Integer r, Integer c);

// create dense matrix with Float ones of size rxc
MATCL_MATREP_EXPORT Matrix  fones(Integer r, Integer c);

// create dense matrix with Float_complex ones of size rxc
MATCL_MATREP_EXPORT Matrix  fcones(Integer r, Integer c);

// create dense matrix with Complex ones of size rxc
MATCL_MATREP_EXPORT Matrix  cones(Integer r, Integer c);

// create dense matrix with ones of size rxc of type given by value_code vt
MATCL_MATREP_EXPORT Matrix  ones(Integer r, Integer c, matcl::value_code vt);

// create dense matrix with ones of size rxc of type Object with type info ti
MATCL_MATREP_EXPORT Matrix  ones(ti::ti_object ti,Integer r, Integer c);

// create sparse matrix with Real ones of size rxc
MATCL_MATREP_EXPORT Matrix  spones(Integer r, Integer c);

// create sparse matrix with Integer ones of size rxc
MATCL_MATREP_EXPORT Matrix  ispones(Integer r, Integer c);

// create sparse matrix with Float ones of size rxc
MATCL_MATREP_EXPORT Matrix  fspones(Integer r, Integer c);

// create sparse matrix with Complex ones of size rxc
MATCL_MATREP_EXPORT Matrix  cspones(Integer r, Integer c);

// create sparse matrix with Float_complex ones of size rxc
MATCL_MATREP_EXPORT Matrix  fcspones(Integer r, Integer c);

// create sparse matrix with ones of size rxc of type given by value_code vt
MATCL_MATREP_EXPORT Matrix  spones(Integer r, Integer c, matcl::value_code vt);

// create sparse matrix with ones of size rxc of type Object with type info ti
MATCL_MATREP_EXPORT Matrix  spones(ti::ti_object ti,Integer r, Integer c);

// create band matrix with Real ones of size rxc
MATCL_MATREP_EXPORT Matrix  bones(Integer r, Integer c);

// create band matrix with Integer ones of size rxc
MATCL_MATREP_EXPORT Matrix  ibones(Integer r, Integer c);

// create band matrix with Float ones of size rxc
MATCL_MATREP_EXPORT Matrix  fbones(Integer r, Integer c);

// create band matrix with Complex ones of size rxc
MATCL_MATREP_EXPORT Matrix  cbones(Integer r, Integer c);

// create band matrix with Float_complex ones of size rxc
MATCL_MATREP_EXPORT Matrix  fcbones(Integer r, Integer c);

// create band matrix with ones of size rxc of type given by value_code vt
MATCL_MATREP_EXPORT Matrix  bones(Integer r, Integer c, matcl::value_code vt);

// create band matrix with ones of size rxc of type Object with type info ti
MATCL_MATREP_EXPORT Matrix  bones(ti::ti_object ti,Integer r, Integer c);

//-----------------------------------------------------------------
//                      eye
//-----------------------------------------------------------------
// create dense matrix with Real ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  eye(Integer r, Integer c);

// create dense matrix with Real ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  eye(Integer r);

// create dense matrix with Integer ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  ieye(Integer r, Integer c);

// create dense matrix with Integer ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  ieye(Integer r);

// create dense matrix with Float ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  feye(Integer r, Integer c);

// create dense matrix with Float ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  feye(Integer r);

// create dense matrix with Complex ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  ceye(Integer r, Integer c);

// create dense matrix with Complex ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  ceye(Integer r);

// create dense matrix with Float_complex ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  fceye(Integer r, Integer c);

// create dense matrix with Float_complex ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  fceye(Integer r);

// create dense matrix with ones on diagonal of size rxc of type given by 
// value_code vt
MATCL_MATREP_EXPORT Matrix  eye(Integer r, Integer c, matcl::value_code vt);

// create dense matrix with ones on diagonal of size rxr of type given by 
// value_code vt
MATCL_MATREP_EXPORT Matrix  eye(Integer r, matcl::value_code vt);

// create dense matrix with ones on diagonal of size rxc of type Object with 
// type info ti
MATCL_MATREP_EXPORT Matrix  eye(ti::ti_object ti,Integer r, Integer c);

// create dense matrix with ones on diagonal of size rxr of type Object with 
// type info ti
MATCL_MATREP_EXPORT Matrix  eye(ti::ti_object ti,Integer r);

// create sparse matrix with Real ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  speye(Integer r, Integer c);

// create sparse matrix with Real ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  speye(Integer r);

// create sparse matrix with Integer ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  ispeye(Integer r, Integer c);

// create sparse matrix with Integer ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  ispeye(Integer r);

// create sparse matrix with Float ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  fspeye(Integer r, Integer c);

// create sparse matrix with Float ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  fspeye(Integer r);

// create sparse matrix with Complex ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  cspeye(Integer r, Integer c);

// create sparse matrix with Complex ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  cspeye(Integer r);

// create sparse matrix with Float_complex ones on diagonal of size rxc
MATCL_MATREP_EXPORT Matrix  fcspeye(Integer r, Integer c);

// create sparse matrix with Float_complex ones on diagonal of size rxr
MATCL_MATREP_EXPORT Matrix  fcspeye(Integer r);

// create sparse matrix with ones on diagonal of size rxc of type given by 
// value_code vt
MATCL_MATREP_EXPORT Matrix  speye(Integer r, Integer c, matcl::value_code vt);

// create sparse matrix with ones on diagonal of size rxr of type given by 
// value_code vt
MATCL_MATREP_EXPORT Matrix  speye(Integer r, matcl::value_code vt);

// create sparse matrix with ones on diagonal of size rxc of type Object with 
// type info ti
MATCL_MATREP_EXPORT Matrix  speye(ti::ti_object ti,Integer r, Integer c);

// create sparse matrix with ones on diagonal of size rxr of type Object with 
// type info ti
MATCL_MATREP_EXPORT Matrix  speye(ti::ti_object ti,Integer r);

// create band matrix with Real ones on diagonal of size rxc with first diag fd
// and last diag ld
MATCL_MATREP_EXPORT Matrix  beye(Integer r, Integer c, Integer fd, Integer ld);

// create band matrix with Real ones on diagonal of size rxr with first diag fd
// and last diag ld
MATCL_MATREP_EXPORT Matrix  beye(Integer r, Integer fd, Integer ld);

// create band matrix with Integer ones on diagonal of size rxc with first diag fd
// and last diag ld
MATCL_MATREP_EXPORT Matrix  ibeye(Integer r, Integer c, Integer fd, Integer ld);

// create band matrix with Integer ones on diagonal of size rxr with first diag fd
// and last diag ld
MATCL_MATREP_EXPORT Matrix  ibeye(Integer r, Integer fd, Integer ld);

// create band matrix with Float ones on diagonal of size rxc with first diag fd 
// and last diag ld
MATCL_MATREP_EXPORT Matrix  fbeye(Integer r, Integer c, Integer fd, Integer ld);

// create band matrix with Float ones on diagonal of size rxr with first diag fd 
// and last diag ld
MATCL_MATREP_EXPORT Matrix  fbeye(Integer r, Integer fd, Integer ld);

// create band matrix with Complex ones on diagonal of size rxc with first diag fd
// and last diag ld
MATCL_MATREP_EXPORT Matrix  cbeye(Integer r, Integer c, Integer fd, Integer ld);

// create band matrix with Complex ones on diagonal of size rxr with first diag fd 
// and last diag ld
MATCL_MATREP_EXPORT Matrix  cbeye(Integer r, Integer fd, Integer ld);

// create band matrix with Float_complex ones on diagonal of size rxc with first diag fd 
// and last diag ld
MATCL_MATREP_EXPORT Matrix  fcbeye(Integer r, Integer c, Integer fd, Integer ld);

// create band matrix with Float_complex ones on diagonal of size rxr with first diag fd
// and last diag ld
MATCL_MATREP_EXPORT Matrix  fcbeye(Integer r, Integer fd, Integer ld);

// create band matrix with ones on diagonal of size rxc of type given by 
// value_code vt with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  beye(Integer r, Integer c, Integer fd, Integer ld, matcl::value_code vt);

// create band matrix with ones on diagonal of size rxr of type given by 
// value_code vt with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  beye(Integer r, Integer fd, Integer ld, matcl::value_code vt);

// create sparse matrix with ones on diagonal of size rxc of type Object with 
// type info ti with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  beye(ti::ti_object ti,Integer r, Integer c, Integer fd, Integer ld);

// create sparse matrix with ones on diagonal of size rxr of type Object with 
// type info ti with first diag fd and last diag ld
MATCL_MATREP_EXPORT Matrix  beye(ti::ti_object ti,Integer r, Integer fd, Integer ld);

//-----------------------------------------------------------------
//                      diag
//-----------------------------------------------------------------
// create dense matrix with elements on diagonal d given by v of size kxk, where 
// k = length(v) + abs(d)
MATCL_MATREP_EXPORT Matrix  diag(const Matrix& v, Integer d= 0);

// create band matrix with elements on diagonal d given by v of size kxk, where 
// k = length(v) + abs(d)
MATCL_MATREP_EXPORT Matrix  bdiag(const Matrix &v, Integer d = 0);

// create sparse matrix with elements on diagonal d given by v of size kxk, where 
// k = length(v) + abs(d)
MATCL_MATREP_EXPORT Matrix  spdiag(const Matrix &v, Integer d = 0);

// create rxc dense matrix from the columns of A and places them along the diagonals 
// specified by d 
MATCL_MATREP_EXPORT Matrix  diags(const Matrix &A, const Matrix &d, Integer r, Integer c);

// create rxc band matrix from the columns of A and places them along the diagonals 
// specified by d 
MATCL_MATREP_EXPORT Matrix  bdiags(const Matrix &A, const Matrix &d, Integer r, Integer c);

// create rxc sparse matrix from the columns of A and places them along the diagonals 
// specified by d 
MATCL_MATREP_EXPORT Matrix  spdiags(const Matrix &A, const Matrix &d, Integer r, Integer c);

//-----------------------------------------------------------------
//                  random numbers generator state
//-----------------------------------------------------------------
// type of random number generator state
// random number generator is thread local; class rand_state is defined
// in matcl-scalar
class rand_state;

// return independent copy of current state of random number generator
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT 
rand_state                  get_rand_state();

// set current state of random number generator, random number generator is
// thread local, therefore this function affects only current thread
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT 
void                        set_rand_state(rand_state st);

// get state of random number generator in current thread
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT 
const rand_state&           global_rand_state();

// initialize state of random number generator pointed by rand_ptr with seed s;
// this affects only current thread;  
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
void                        init_genrand(unsigned long s, const rand_state& rand_ptr);

// create local state of random number generator, restore previous state
// after exiting from current scope
// redeclaration of function defined in matcl-scalar
MATCL_SCALAR_EXPORT 
details::local_rand_state_raii
                            local_rand_state(unsigned long s);  

//-----------------------------------------------------------------
//                      random numbers
//-----------------------------------------------------------------
// generate uniformly distributed Real random number on (0, 1)
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Real                        rand(const rand_state& rand_ptr);

// generate uniformly distributed Float random number on (0, 1)
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Float                       frand(const rand_state& rand_ptr);

// generate uniformly distributed Integer random number on [min,max], 
// where min, max is the range of all possible values of Integer scalars
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Integer                     irand(const rand_state& rand_ptr);

// generate uniformly distributed Complex random number as Complex(rand(),rand())
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Complex                     crand(const rand_state& rand_ptr);

// generate uniformly distributed Float_complex random number as 
// Float_complex(frand(),frand())
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Float_complex               fcrand(const rand_state& rand_ptr);

// generate normally distributed Real random number
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Real                        randn(const rand_state& rand_ptr);

// generate normally distributed Float random number
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Float                       frandn(const rand_state& rand_ptr);

// generate normally distributed Complex random number as Complex(randn(),randn())
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Complex                     crandn(const rand_state& rand_ptr);

// generate normally distributed Float_complex random number as 
// Float_complex(frandn(),frandn())
// redeclaration of function defined in matcl-scalar; default value of rand_ptr
// is global_rand_state();
MATCL_SCALAR_EXPORT 
Float_complex               fcrandn(const rand_state& rand_ptr);

// create dense matrix of size rxc of uniformly distributed Real random numbers on (0, 1)
MATCL_MATREP_EXPORT Matrix  rand(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of uniformly distributed Integer random numbers on [min,max], 
// where min, max is the range of all possible values of Integer scalars
MATCL_MATREP_EXPORT Matrix  irand(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of uniformly distributed Float random numbers on (0, 1)
MATCL_MATREP_EXPORT Matrix  frand(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of uniformly distributed Complex random numbers, i.e. 
// Complex(rand(),rand())
MATCL_MATREP_EXPORT Matrix  crand(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of uniformly distributed Float_complex random numbers, i.e. 
// Float_complex(frand(),frand())
MATCL_MATREP_EXPORT Matrix  fcrand(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of uniformly distributed random numbers of type given by
// value_code
MATCL_MATREP_EXPORT Matrix  rand(Integer r, Integer c, matcl::value_code vt, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of normally distributed Real random numbers
MATCL_MATREP_EXPORT Matrix  randn(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of normally distributed Float random numbers
MATCL_MATREP_EXPORT Matrix  frandn(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of normally distributed Complex random numbers, i.e. 
// Complex(randn(),randn())
MATCL_MATREP_EXPORT Matrix  crandn(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of normally distributed Float_complex random numbers, i.e. 
// Float_complex(frandn(),frandn())
MATCL_MATREP_EXPORT Matrix  fcrandn(Integer r, Integer c, const rand_state& rand_ptr 
                                = global_rand_state());

// create dense matrix of size rxc of normally distributed random numbers of type given by
// value_code
MATCL_MATREP_EXPORT Matrix  randn(Integer r, Integer c, matcl::value_code vt, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
// distributed Real random numbers on (0, 1)
MATCL_MATREP_EXPORT Matrix  sprand(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
// distributed Integer random numbers on [min,max],  where min, max is the range of all 
// possible values of Integer scalars
MATCL_MATREP_EXPORT Matrix  isprand(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
// distributed Float random numbers on (0, 1)
MATCL_MATREP_EXPORT Matrix  fsprand(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
// distributed Complex random numbers, i.e. Complex(rand(),rand())
MATCL_MATREP_EXPORT Matrix  csprand(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
// distributed Float_complex random numbers, i.e. Float_complex(frand(),frand())
MATCL_MATREP_EXPORT Matrix  fcsprand(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of uniformly
// distributed numbers of type given by value_code
MATCL_MATREP_EXPORT Matrix  sprand(Integer r, Integer c, Real d,matcl::value_code vt, 
                                const rand_state& rand_ptr = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of normally
// distributed Real random numbers
MATCL_MATREP_EXPORT Matrix  sprandn(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of normally
// distributed Float random numbers
MATCL_MATREP_EXPORT Matrix  fsprandn(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of normally
// distributed Complex random numbers, i.e. Complex(randn(),randn())
MATCL_MATREP_EXPORT Matrix  csprandn(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of normally
// distributed Float_complex random numbers, i.e. Float_complex(frandn(),frandn())
MATCL_MATREP_EXPORT Matrix  fcsprandn(Integer r, Integer c, Real d, const rand_state& rand_ptr 
                                = global_rand_state());

// create sparse matrix of size rxc with approximatelly d*r*c nonzero elements of normally
// distributed numbers of type given by value_code
MATCL_MATREP_EXPORT Matrix  sprandn(Integer r, Integer c, Real d,matcl::value_code vt,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
// distributed Real random numbers on (0, 1)
MATCL_MATREP_EXPORT Matrix  rand_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
// distributed Integer random numbers on [min,max],  where min, max is the range of all 
// possible values of Integer scalars
MATCL_MATREP_EXPORT Matrix  irand_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
// distributed Float random numbers on (0, 1)
MATCL_MATREP_EXPORT Matrix  frand_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
// distributed Complex random numbers, i.e. Complex(rand(),rand())
MATCL_MATREP_EXPORT Matrix  crand_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
// distributed Float_complex random numbers, i.e. Float_complex(frand(),frand())
MATCL_MATREP_EXPORT Matrix  fcrand_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of uniformly
// distributed random numbers of type given by value_code
MATCL_MATREP_EXPORT Matrix  rand_band(Integer r, Integer c, Integer fd, Integer ld,matcl::value_code vt,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of normally
// distributed Real random numbers
MATCL_MATREP_EXPORT Matrix  randn_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of normally
// distributed Float random numbers
MATCL_MATREP_EXPORT Matrix  frandn_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of normally
// distributed Complex random numbers, i.e. Complex(randn(),randn())
MATCL_MATREP_EXPORT Matrix  crandn_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of normally
// distributed Float_complex random numbers, i.e. Float_complex(frandn(),frandn())
MATCL_MATREP_EXPORT Matrix  fcrandn_band(Integer r, Integer c, Integer fd, Integer ld,
                                const rand_state& rand_ptr = global_rand_state());

// create band matrix of size rxc with first diagonal fd and last diagonal ld of normally
// distributed random numbers of type given by value_code
MATCL_MATREP_EXPORT Matrix  randn_band(Integer r, Integer c, Integer fd, Integer ld, matcl::value_code vt,
                                const rand_state& rand_ptr = global_rand_state());

// create random permutation of integers 1:n
MATCL_MATREP_EXPORT permvec    randperm(Integer n, const rand_state& rand_ptr = global_rand_state());

//-----------------------------------------------------------------
//              fine grained constructors
//-----------------------------------------------------------------
// create dense matrix of size rows x cols with Integer zeroes
MATCL_MATREP_EXPORT Matrix  make_integer_dense(Integer rows,Integer cols);

// create dense matrix of size rows x cols with Real zeroes
MATCL_MATREP_EXPORT Matrix  make_real_dense(Integer rows,Integer cols);

// create dense matrix of size rows x cols with Float zeroes
MATCL_MATREP_EXPORT Matrix  make_float_dense(Integer rows,Integer cols);

// create dense matrix of size rows x cols with Complex zeroes
MATCL_MATREP_EXPORT Matrix  make_complex_dense(Integer rows,Integer cols);

// create dense matrix of size rows x cols with Float_complex zeroes
MATCL_MATREP_EXPORT Matrix  make_float_complex_dense(Integer rows,Integer cols);

// create dense matrix of size rows x cols with Object zeroes with type_info ti
MATCL_MATREP_EXPORT Matrix  make_object_dense(ti::ti_object ti,Integer rows, Integer cols);

// create dense matrix of size rows x cols with elements copied from Integer array
// arr (column oriented storage assumed); optional argument set leading dimension 
// (distance between consecutive columns) of the array arr
MATCL_MATREP_EXPORT Matrix  make_integer_dense(Integer rows,Integer cols, const Integer *arr);
MATCL_MATREP_EXPORT Matrix  make_integer_dense(Integer rows,Integer cols, const Integer *arr, Integer ld);

// create dense matrix of size rows x cols with elements copied from Real array
// arr (column oriented storage assumed); optional argument set leading dimension 
// (distance between consecutive columns) of the array arr
MATCL_MATREP_EXPORT Matrix  make_real_dense(Integer rows,Integer cols,const Real *arr);
MATCL_MATREP_EXPORT Matrix  make_real_dense(Integer rows,Integer cols,const Real *arr, Integer ld);

// create dense matrix of size rows x cols with elements copied from Float array
// arr (column oriented storage assumed); optional argument set leading dimension 
// (distance between consecutive columns) of the array arr
MATCL_MATREP_EXPORT Matrix  make_float_dense(Integer rows,Integer cols,const Float *arr);
MATCL_MATREP_EXPORT Matrix  make_float_dense(Integer rows,Integer cols,const Float *arr, Integer ld);

// create dense matrix of size rows x cols with elements copied from Complex array
// arr (column oriented storage assumed); optional argument set leading dimension 
// (distance between consecutive columns) of the array arr
MATCL_MATREP_EXPORT Matrix  make_complex_dense(Integer rows,Integer cols,const Complex *arr);
MATCL_MATREP_EXPORT Matrix  make_complex_dense(Integer rows,Integer cols,const Complex *arr, Integer ld);

// create dense matrix of size rows x cols with elements copied from Float_complex array
// arr (column oriented storage assumed); optional argument set leading dimension 
// (distance between consecutive columns) of the array arr
MATCL_MATREP_EXPORT Matrix  make_float_complex_dense(Integer rows,Integer cols,
                                const Float_complex *arr);
MATCL_MATREP_EXPORT Matrix  make_float_complex_dense(Integer rows,Integer cols,
                                const Float_complex *arr, Integer ld);

// create dense matrix of size rows x cols with elements copied from Real arrays
// ar_r, ar_i containing real and complex part (column oriented storage assumed); 
MATCL_MATREP_EXPORT Matrix  make_complex_dense(Integer rows,Integer cols,const Real *ar_r,
                                const Real *ar_i);

// create dense matrix of size rows x cols with elements copied from Float arrays
// ar_r, ar_i containing real and complex part (column oriented storage assumed)
MATCL_MATREP_EXPORT Matrix  make_float_complex_dense(Integer rows,Integer cols,const Float *ar_r,
                                const Float *ar_i);

// create dense matrix of size rows x cols with elements copied from Object array
// arr (column oriented storage assumed) and convert elements to Object with type_info ti
// optional argument set leading dimension (distance between consecutive columns)
// of the array arr
MATCL_MATREP_EXPORT Matrix  make_object_dense(ti::ti_object ti,Integer rows,Integer cols,
                                const Object *arr);
MATCL_MATREP_EXPORT Matrix  make_object_dense(ti::ti_object ti,Integer rows,Integer cols,
                                const Object *arr, Integer ld);

// create dense matrix of size rows x cols with Integer elements val
MATCL_MATREP_EXPORT Matrix  make_integer_dense(Integer val,Integer rows,Integer cols);

// create dense matrix of size rows x cols with Real elements val
MATCL_MATREP_EXPORT Matrix  make_real_dense(Real val,Integer rows,Integer cols);

// create dense matrix of size rows x cols with Float elements val
MATCL_MATREP_EXPORT Matrix  make_float_dense(Float val,Integer rows,Integer cols);

// create dense matrix of size rows x cols with Complex elements val
MATCL_MATREP_EXPORT Matrix  make_complex_dense(Complex val,Integer rows,Integer cols);

// create dense matrix of size rows x cols with Float_complex elements val
MATCL_MATREP_EXPORT Matrix  make_float_complex_dense(Float_complex val,Integer rows,Integer cols);

// create dense matrix of size rows x cols with Object elements val
MATCL_MATREP_EXPORT Matrix  make_object_dense(const Object& val,Integer rows,Integer cols);

// create dense matrix of size rows x cols with elements stored in the array arr 
// (column oriented storage assumed) with leading dimension ld (distance between
// consecutive columns) array arr is not internally copied; the array arr cannot
// be constant, condition std::decay<T> == T must be satisfied
//
// this is unsafe function; one should ensure that this matrix and all matrices 
// created from this matrix directly or indirectly are destroyed before destroying
// the array arr
template<class T, class Enable = typename details::enable_if_matcl_scalar_not_object<T, void>::type>
MATCL_MATREP_EXPORT Matrix  make_dense_foreign(Integer rows,Integer cols, T* arr, Integer ld);
MATCL_MATREP_EXPORT Matrix  make_dense_foreign(ti::ti_object ti, Integer rows,Integer cols, Object* arr, 
                                Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Integer zeroes
MATCL_MATREP_EXPORT Matrix  make_integer_band(Integer rows,Integer cols, Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Real zeroes
MATCL_MATREP_EXPORT Matrix  make_real_band(Integer rows,Integer cols, Integer fd, Integer ld);

// create band matrix of size rows x cols with with first diagonal fd and last diagonal ud,
// with Float zeroes
MATCL_MATREP_EXPORT Matrix  make_float_band(Integer rows,Integer cols, Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Complex zeroes
MATCL_MATREP_EXPORT Matrix  make_complex_band(Integer rows,Integer cols, Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Float_complex zeroes
MATCL_MATREP_EXPORT Matrix  make_float_complex_band(Integer rows,Integer cols, Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Object zeroes with type_info ti
MATCL_MATREP_EXPORT Matrix  make_object_band(ti::ti_object ti,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Integer elements val
MATCL_MATREP_EXPORT Matrix  make_integer_band(Integer val,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Real elements val
MATCL_MATREP_EXPORT Matrix  make_real_band(Real val,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Float elements val
MATCL_MATREP_EXPORT Matrix  make_float_band(Float val,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Complex elements val
MATCL_MATREP_EXPORT Matrix  make_complex_band(Complex val,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Float_complex elements val
MATCL_MATREP_EXPORT Matrix  make_float_complex_band(Float_complex val,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ud,
// with Object elements val
MATCL_MATREP_EXPORT Matrix  make_object_band(const Object& val,Integer rows,Integer cols, 
                                Integer fd, Integer ld);

// create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
// with Integer elements
MATCL_MATREP_EXPORT Matrix  make_integer_sparse(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
// with Real elements
MATCL_MATREP_EXPORT Matrix  make_real_sparse(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
// with Float elements
MATCL_MATREP_EXPORT Matrix  make_float_sparse(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
// with Complex elements
MATCL_MATREP_EXPORT Matrix  make_complex_sparse(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
// with Float_complex elements
MATCL_MATREP_EXPORT Matrix  make_float_complex_sparse(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with nnz = 0, and capacity nzmax, 
// with Object elements of type given by ti
MATCL_MATREP_EXPORT Matrix  make_object_sparse(ti::ti_object ti,Integer rows,Integer cols, 
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are arrays of length at least nnz; allocate space for at least nzmax Integer elements
MATCL_MATREP_EXPORT Matrix  make_integer_sparse(const Integer *trip_r, const Integer *trip_c, 
                                const Integer *trip_x, Integer r, Integer c, Integer nnz, 
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are arrays of length at least nnz; allocate space for at least  nzmax Real elements
MATCL_MATREP_EXPORT Matrix  make_real_sparse(const Integer *trip_r, const Integer *trip_c, 
                                const Real *trip_x, Integer r, Integer c, Integer nnz,
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are arrays of length at least nnz; allocate space for at least  nzmax Float elements
MATCL_MATREP_EXPORT Matrix  make_float_sparse(const Integer *trip_r, const Integer *trip_c, 
                                const Float *trip_x, Integer r, Integer c, Integer nnz,
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are arrays of length at least nnz; allocate space for at least  nzmax Complex elements
MATCL_MATREP_EXPORT Matrix  make_complex_sparse(const Integer *trip_r, const Integer *trip_c,
                                const Complex *trip_x, Integer r, Integer c, Integer nnz,
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are arrays of length at least nnz; allocate space for at least nzmax Float_complex elements
MATCL_MATREP_EXPORT Matrix  make_float_complex_sparse(const Integer *trip_r, const Integer *trip_c,
                                const Float_complex *trip_x, Integer r, Integer c, Integer nnz,
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_re
// trip_im are arrays of length at least nnz; real and imaginary elements are stored in separate 
// arrays trip_re and trip_im; allocate space for at least nzmax Complex elements
MATCL_MATREP_EXPORT Matrix  make_complex_sparse(const Integer *trip_r, const Integer *trip_c, 
                                const Real *trip_re, const Real *trip_im, Integer r, Integer c, 
                                Integer nnz, Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_re
// trip_im are arrays of length at least nnz; real and imaginary elements are stored in separate 
// arrays trip_re and trip_im; allocate space for at least nzmax Float_complex elements
MATCL_MATREP_EXPORT Matrix  make_float_complex_sparse(const Integer *trip_r, const Integer *trip_c, 
                                const Float *trip_re, const Float *trip_im, Integer r, Integer c, 
                                Integer nnz, Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are arrays of length at least nnz; allocate space for at least nzmax Object elements of type
// given by type_info ti
MATCL_MATREP_EXPORT Matrix  make_object_sparse(ti::ti_object ti, const Integer *trip_r, const Integer 
                                *trip_c, const Object *trip_x, Integer r, Integer c, Integer nnz,
                                Integer nzmax = 0);

// create sparse matrix of size r x c from triplet representation; trip_c, trip_r, trip_x
// are vectors of equal length nnz; allocate space for at least nzmax elements
MATCL_MATREP_EXPORT Matrix  make_sparse_matrix(const Matrix& trip_r, const Matrix& trip_c, 
                                const Matrix& trip_x, Integer r, Integer c, Integer nzmax = 0);

// create dense matrix of size rows x cols with zeroes of type given by value_code
MATCL_MATREP_EXPORT Matrix  make_dense_matrix(Integer rows,Integer cols, matcl::value_code);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld 
// with zeroes of type given by value_code
MATCL_MATREP_EXPORT Matrix  make_band_matrix(Integer rows,Integer cols, Integer fd, Integer ld,
                                matcl::value_code);

// create sparse matrix of size rows x cols with capacity for nzmax nonzeroest elements
// of type given by value_code
MATCL_MATREP_EXPORT Matrix  make_sparse_matrix(Integer rows,Integer cols, Integer nzmax, 
                                matcl::value_code);

//-----------------------------------------------------------------
//              uninitialized matrices
//-----------------------------------------------------------------
// create dense matrix of size rows x cols; elements are not initialized; 
// on return ptr_data is the pointer to stored Integer data
MATCL_MATREP_EXPORT Matrix  make_integer_dense_noinit(Integer rows,Integer cols, Integer*& ptr_data);

// create dense matrix of size rows x cols; elements are not initialized; 
// on return ptr_data is the pointer to stored Real data
MATCL_MATREP_EXPORT Matrix  make_real_dense_noinit(Integer rows,Integer cols, Real*& ptr_data);

// create dense matrix of size rows x cols; elements are not initialized; 
// on return ptr_data is the pointer to stored Float data
MATCL_MATREP_EXPORT Matrix  make_float_dense_noinit(Integer rows,Integer cols, Float*& ptr_data);

// create dense matrix of size rows x cols; elements are not initialized; 
// on return ptr_data is the pointer to stored Complex data
MATCL_MATREP_EXPORT Matrix  make_complex_dense_noinit(Integer rows,Integer cols, Complex*& ptr_data);

// create dense matrix of size rows x cols; elements are not initialized; 
// on return ptr_data is the pointer to stored Float_complex data
MATCL_MATREP_EXPORT Matrix  make_float_complex_dense_noinit(Integer rows,Integer cols, 
                                Float_complex*& ptr_data);

// create dense matrix of size rows x cols; elements are not initialized; 
// on return ptr_data is the pointer to stored Object data with type_info ti
MATCL_MATREP_EXPORT Matrix  make_object_dense_noinit(ti::ti_object,Integer rows,Integer cols,
                                Object*& ptr_data);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
// elements are not initialized; matrix stores Integer elements
MATCL_MATREP_EXPORT Matrix  make_integer_band_noinit(Integer rows,Integer cols, Integer fd, 
                                Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
// elements are not initialized; matrix stores Real elements
MATCL_MATREP_EXPORT Matrix  make_real_band_noinit(Integer rows,Integer cols, Integer fd, 
                                Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
// elements are not initialized; matrix stores Float elements
MATCL_MATREP_EXPORT Matrix  make_float_band_noinit(Integer rows,Integer cols, Integer fd, 
                                Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
// elements are not initialized; matrix stores Complex elements
MATCL_MATREP_EXPORT Matrix  make_complex_band_noinit(Integer rows,Integer cols, Integer fd, 
                                Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
// elements are not initialized; matrix stores Float_complex elements
MATCL_MATREP_EXPORT Matrix  make_float_complex_band_noinit(Integer rows,Integer cols, Integer fd, 
                                Integer ld);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld; 
// elements are not initialized; matrix stores Object elements with type_info ti
MATCL_MATREP_EXPORT Matrix  make_object_band_noinit(ti::ti_object ti,Integer rows, Integer cols, 
                                Integer fd, Integer ld);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
// elements are not initialized; matrix stores Integer elements
MATCL_MATREP_EXPORT Matrix  make_integer_sparse_noinit(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
// elements are not initialized; matrix stores Real elements
MATCL_MATREP_EXPORT Matrix  make_real_sparse_noinit(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
// elements are not initialized; matrix stores Float elements
MATCL_MATREP_EXPORT Matrix  make_float_sparse_noinit(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
// elements are not initialized; matrix stores Complex elements
MATCL_MATREP_EXPORT Matrix  make_complex_sparse_noinit(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
// elements are not initialized; matrix stores Float_complex elements
MATCL_MATREP_EXPORT Matrix  make_float_complex_sparse_noinit(Integer rows,Integer cols, Integer nzmax = 0);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements; 
// elements are not initialized; matrix stores Object elements with type_info ti
MATCL_MATREP_EXPORT Matrix  make_object_sparse_noinit(ti::ti_object ti, Integer rows, Integer cols, 
                                Integer nzmax = 0);

// create dense matrix of size rows x cols with uninitialized elements of type given by value_code
MATCL_MATREP_EXPORT Matrix  make_dense_noinit(Integer rows,Integer cols, matcl::value_code);

// create band matrix of size rows x cols with first diagonal fd and last diagonal ld
// uninitialized elements of type given by value_code
MATCL_MATREP_EXPORT Matrix  make_band_noinit(Integer rows,Integer cols, Integer fd,
                                Integer ld, matcl::value_code);

// create sparse matrix of size rows x cols with capacity for nzmax nonzero elements with
// uninitialized elements of type given by value_code
MATCL_MATREP_EXPORT Matrix  make_sparse_noinit(Integer rows,Integer cols, Integer nzmax, matcl::value_code);

};

#pragma warning(pop)