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

#include "matcl-scalar/lib_functions/manip.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-scalar/lib_functions/func_unary.h"

#pragma warning(push)
#pragma warning(disable: 4244) //conversion from 'matcl::Real' to 'matcl::Float', possible loss of data

namespace matcl { namespace details
{

template<class To, class From> struct convert_scalar_helper;

//Integer
template<> struct convert_scalar_helper<Integer,Integer>    
{ 
    static Integer eval(Integer v) { return v; }
};
template<> struct convert_scalar_helper<Integer,Float>    
{ 
    static Integer eval(Float v) { return (Integer)v; }
};
template<> struct convert_scalar_helper<Integer,Real>    
{ 
    static Integer eval(Real v) { return (Integer)v; }
};
template<> struct convert_scalar_helper<Integer,Float_complex>    
{ 
    static Integer eval(const Float_complex& v) { return (Integer)real(v); }
};
template<> struct convert_scalar_helper<Integer,Complex>    
{ 
    static Integer eval(const Complex& v) { return (Integer)real(v); }
};
template<> struct convert_scalar_helper<Integer,Object>    
{ 
    static Integer eval(const Object& v) { return cast_integer(v); }
};

//Float
template<> struct convert_scalar_helper<Float,Integer>    
{ 
    static Float eval(Integer v) { return (Float)v; }
};
template<> struct convert_scalar_helper<Float,Float>    
{ 
    static Float eval(Float v) { return v; }
};
template<> struct convert_scalar_helper<Float,Real>    
{ 
    static Float eval(Real v) { return (Float)v; }
};
template<> struct convert_scalar_helper<Float,Float_complex>    
{ 
    static Float eval(const Float_complex& v) { return real(v); }
};
template<> struct convert_scalar_helper<Float,Complex>    
{ 
    static Float eval(const Complex& v) { return (Float)real(v); }
};
template<> struct convert_scalar_helper<Float,Object>    
{ 
    static Float eval(const Object& v) { return cast_float(v); }
};

//Real
template<> struct convert_scalar_helper<Real,Integer>    
{ 
    static Real eval(Integer v) { return (Real)v; }
};
template<> struct convert_scalar_helper<Real,Float>    
{ 
    static Real eval(Float v) { return (Real)v; }
};
template<> struct convert_scalar_helper<Real,Real>    
{ 
    static Real eval(Real v) { return v; }
};
template<> struct convert_scalar_helper<Real,Float_complex>    
{ 
    static Real eval(const Float_complex& v) { return (Real)real(v); }
};
template<> struct convert_scalar_helper<Real,Complex>    
{ 
    static Real eval(const Complex& v) { return real(v); }
};
template<> struct convert_scalar_helper<Real,Object>    
{ 
    static Real eval(const Object& v) { return cast_real(v); }
};

//Float_complex
template<> struct convert_scalar_helper<Float_complex,Integer>    
{ 
    static Float_complex eval(Integer v) { return (Float_complex)(Real)v; }
};
template<> struct convert_scalar_helper<Float_complex,Float>    
{ 
    static Float_complex eval(Float v) { return Float_complex(v); }
};
template<> struct convert_scalar_helper<Float_complex,Real>    
{ 
    static Float_complex eval(Real v) { return (Float_complex)v; }
};
template<> struct convert_scalar_helper<Float_complex,Float_complex>    
{ 
    static Float_complex eval(const Float_complex& v) { return v; }
};
template<> struct convert_scalar_helper<Float_complex,Complex>    
{ 
    static Float_complex eval(const Complex& v) { return Float_complex(v); }
};
template<> struct convert_scalar_helper<Float_complex,Object>    
{ 
    static Float_complex eval(const Object& v) { return cast_float_complex(v); }
};

//Complex
template<> struct convert_scalar_helper<Complex,Integer>    
{ 
    static Complex eval(Integer v) { return (Complex)(Real)v; }
};
template<> struct convert_scalar_helper<Complex,Float>    
{ 
    static Complex eval(Float v) { return (Complex)(Real)v; }
};
template<> struct convert_scalar_helper<Complex,Real>    
{ 
    static Complex eval(Real v) { return (Complex)v; }
};
template<> struct convert_scalar_helper<Complex,Float_complex>    
{ 
    static Complex eval(const Float_complex& v) { return Complex(v); }
};
template<> struct convert_scalar_helper<Complex,Complex>    
{ 
    static Complex eval(const Complex& v) { return v; }
};
template<> struct convert_scalar_helper<Complex,Object>    
{ 
    static Complex eval(const Object& v) { return cast_complex(v); }
};

//Object
template<> struct convert_scalar_helper<Object,Integer>    
{ 
    static Object eval(Integer v) { return Object(v); }
};
template<> struct convert_scalar_helper<Object,Float>    
{ 
    static Object eval(Float v) { return Object(v); }
};
template<> struct convert_scalar_helper<Object,Real>    
{ 
    static Object eval(Real v) { return Object(v); }
};
template<> struct convert_scalar_helper<Object,Float_complex>    
{ 
    static Object eval(const Float_complex& v) { return Object(v); }
};
template<> struct convert_scalar_helper<Object,Complex>    
{ 
    static Object eval(const Complex& v) { return Object(v); }
};
template<> struct convert_scalar_helper<Object,Object>    
{ 
    static Object eval(const Object& v) { return v; }
};

}}

#pragma warning(pop)

namespace matcl
{

template<class To, class From, class Enable>
To matcl::convert_scalar(const From& s)
{
    using SP = typename md::promote_scalar<From>::type;
    return details::convert_scalar_helper<To,SP>::eval(SP(s));
};

};
