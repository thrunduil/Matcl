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

#include "matcl-core/lib_functions/constants.h"
#include <limits>

namespace matcl { namespace details
{

//internal use only
#define matcl_CONST_LOG2E       1.442695040888963407359925
#define matcl_CONST_LOG10E      0.4342944819032518276511289
#define matcl_CONST_E           2.718281828459045235360287471352662497e+00
#define matcl_CONST_PI          3.141592653589793238462643383279502884e+00
#define matcl_CONST_HALF_PI     1.570796326794896619231321691639751442e+00
#define matcl_CONST_LN_2        6.931471805599453094172321214581765680e-01
#define matcl_CONST_LN_10       2.302585092994045684017991454684364207e+00
#define matcl_CONST_SQRT_2      1.414213562373095048801688724209698078e+00
#define matcl_CONST_SQRT_1_2    7.071067811865475244008443621048490392e-01

}};

namespace matcl 
{

inline Integer constants::max_int()   { return std::numeric_limits<Integer>::max();};
inline Integer constants::min_int()   { return std::numeric_limits<Integer>::min();};

inline Complex constants::i()         { return Complex(0, 1);};
inline Float_complex constants::f_i() { return Float_complex(0, 1);};

inline Real constants::eps()          { return std::numeric_limits< double >::epsilon();};
inline Real constants::min_real()     { return std::numeric_limits< double >::min();};
inline Real constants::max_real()     { return std::numeric_limits< double >::max();};
inline Real constants::inf()          { return std::numeric_limits< double >::infinity();};
inline Real constants::nan()          { return std::numeric_limits< double >::quiet_NaN();};
inline Real constants::e()            { return matcl_CONST_E;};
inline Real constants::pi()           { return matcl_CONST_PI;};
inline Real constants::pi_2()         { return matcl_CONST_HALF_PI;};
inline Real constants::log2e()        { return matcl_CONST_LOG2E;};
inline Real constants::log10e()       { return matcl_CONST_LOG10E;};
inline Real constants::ln2()          { return matcl_CONST_LN_2;};
inline Real constants::ln10()         { return matcl_CONST_LN_10;};
inline Real constants::sqrt2()        { return matcl_CONST_SQRT_2;};
inline Real constants::sqrt1_2()      { return matcl_CONST_SQRT_1_2;};

inline Float constants::f_eps()       { return std::numeric_limits< Float >::epsilon();};
inline Float constants::f_min_real()  { return std::numeric_limits< Float >::min();};
inline Float constants::f_max_real()  { return std::numeric_limits< Float >::max();};
inline Float constants::f_inf()       { return std::numeric_limits< Float >::infinity();};
inline Float constants::f_nan()       { return std::numeric_limits< Float >::quiet_NaN();};
inline Float constants::f_e()         { return static_cast<Float>(matcl_CONST_E);};
inline Float constants::f_pi()        { return static_cast<Float>(matcl_CONST_PI);};
inline Float constants::f_pi_2()      { return static_cast<Float>(matcl_CONST_HALF_PI);};
inline Float constants::f_log2e()     { return static_cast<Float>(matcl_CONST_LOG2E);};
inline Float constants::f_log10e()    { return static_cast<Float>(matcl_CONST_LOG10E);};
inline Float constants::f_ln2()       { return static_cast<Float>(matcl_CONST_LN_2);};
inline Float constants::f_ln10()      { return static_cast<Float>(matcl_CONST_LN_10);};
inline Float constants::f_sqrt2()     { return static_cast<Float>(matcl_CONST_SQRT_2);};
inline Float constants::f_sqrt1_2()   { return static_cast<Float>(matcl_CONST_SQRT_1_2);};

template<> inline 
Real constants::eps<>()               { return constants::eps(); };

template<> inline 
Real constants::min_real<>()          { return constants::min_real(); };

template<> inline 
Real constants::max_real<>()          { return constants::max_real(); };

template<> inline 
Real constants::inf<>()               { return constants::inf(); };

template<> inline 
Real constants::nan<>()               { return constants::nan(); };

template<> inline 
Real constants::e<>()                 { return constants::e(); };

template<> inline 
Real constants::pi<>()                { return constants::pi(); };

template<> inline 
Real constants::pi_2<>()              { return constants::pi_2(); };

template<> inline 
Real constants::log2e<>()             { return constants::log2e(); };

template<> inline 
Real constants::log10e<>()            { return constants::log10e(); };

template<> inline 
Real constants::ln2<>()               { return constants::ln2(); };

template<> inline 
Real constants::ln10<>()              { return constants::ln10(); };

template<> inline 
Real constants::sqrt2<>()             { return constants::sqrt2(); };

template<> inline 
Real constants::sqrt1_2<>()           { return constants::sqrt1_2(); };

template<> inline 
complex<Real> constants::i<>()        { return constants::i(); };

template<> inline 
Float constants::eps<>()              { return constants::f_eps(); };

template<> inline 
Float constants::min_real<>()         { return constants::f_min_real(); };

template<> inline 
Float constants::max_real<>()         { return constants::f_max_real(); };

template<> inline
Float constants::inf<>()              { return constants::f_inf(); };

template<> inline 
Float constants::nan<>()              { return constants::f_nan(); };

template<> inline 
Float constants::e<>()                { return constants::f_e(); };

template<> inline 
Float constants::pi<>()               { return constants::f_pi(); };

template<> inline 
Float constants::pi_2<>()             { return constants::f_pi_2(); };

template<> inline 
Float constants::log2e<>()            { return constants::f_log2e(); };

template<> inline 
Float constants::log10e<>()           { return constants::f_log10e(); };

template<> inline 
Float constants::ln2<>()              { return constants::f_ln2(); };

template<> inline 
Float constants::ln10<>()             { return constants::f_ln10(); };

template<> inline 
Float constants::sqrt2<>()            { return constants::f_sqrt2(); };

template<> inline 
Float constants::sqrt1_2<>()          { return constants::f_sqrt1_2(); };

template<> inline 
complex<Float> constants::i<>()       { return constants::f_i(); };

};
