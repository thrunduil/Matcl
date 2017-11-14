/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-simd/simd.h"
#include "matcl-core/float/twofold.h"

namespace matcl { namespace simd
{

template<class Float_type, class ... Args>
Float_type  short_horner(Float_type x, Args ... coef);

template<int Poly_size, class Arg_type, class Coef_type>
Arg_type    horner(Arg_type x, const Coef_type* poly);

template<class Arg_type, class Coef_type>
Arg_type    horner(Arg_type x, int poly_size, const Coef_type* poly);

template<class Arg_type, class Coef_type>
Arg_type    horner_abs(Arg_type x, int poly_size, const Coef_type* poly);

template<class Arg_type, class Coef_type>
Arg_type    twofold_horner(Arg_type x, int poly_size, const Coef_type* poly);

template<class Arg_type, class Coef_type>
Arg_type    horner_cond(Arg_type x, int poly_size, const Coef_type* poly);

template<int Poly_size, class Arg_type, class Coef_type>
Arg_type    estrin(Arg_type x, const Coef_type* poly);

template<class Arg_type, class Coef_type>
Arg_type    estrin(Arg_type x, int poly_size, const Coef_type* poly);

#if 0
// evaluate a polynomial at point x using the Horner's scheme
//
//      res = sum_{i = 0}^{N - 1} a_i * x^i
//
// polynomial is represented as an array
//      poly = {a_0, a_1, ..., a_{N-1}}
//
// Template arguments:
//  Poly_size   - number of coefficients (N) of the polynomial
//  Float_type  - floating point type (float, double, or a simd type)
template<int Poly_size, class Float_type>
Float_type  horner(const Float_type* poly, Float_type x);

// evaluate a polynomial at point x using the Horner's scheme
//
//      res = sum_{i = 0}^{N - 1} |a_i| * x^i
//
// coefficients of the polynomial are constructed from absolute values
// of given polynomial, represented as an array:
//      poly = {a_0, a_1, ..., a_{N-1}}
//
// Template arguments:
//  Poly_size   - number of coefficients (N) of the polynomial
//  Float_type  - floating point type (float, double, or a simd type)
template<int Poly_size, class Float_type>
Float_type  horner_abs(const Float_type* poly, Float_type x);

template<int Poly_size, class Float_type>
Float_type  twofold_horner(const Float_type* poly, Float_type x,
                twofold<Float_type>* err);

// return the condition number of a polynomial at point x; 
// The relative accuracy of computed polynomial p_app(x) is given by
//      |p(x) - p_app(x)| / |p(x)| <= cond(p, x) + O(u^2)
// where p(x) is the true value and cond(p, x) is returned condition number
template<int Poly_size, class Float_type>
Float_type  poly_cond(const Float_type* poly, Float_type x);

#endif
}};

#include "matcl-simd/details/poly/poly_eval.inl"
