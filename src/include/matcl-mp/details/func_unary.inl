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

#include "matcl-mp/func_unary.h"
#include "matcl-mp/mp_float.h"
#include "matcl-mp/mp_rational.h"
#include "matcl-mp/mp_complex.h"
#include "matcl-mp/details/utils.h"
#include "matcl-core/matrix/complex_type.h"

#pragma warning(push)
#pragma warning(disable:4244)// possible loss of data

namespace matcl { namespace mp { namespace details
{

namespace md = matcl::details;

inline Integer  get_mp_scalar_int(const mp_int& s)      { return s.cast_int(); };
inline Integer  get_mp_scalar_int(const mp_float& s)    { return s.cast_int(); };
inline Integer  get_mp_scalar_int(const mp_complex& s)  { return s.real().cast_int(); };
inline Integer  get_mp_scalar_int(const mp_rational& s) { return s.cast_int(); };

inline Real     get_mp_scalar_real(const mp_int& s)     { return s.cast_float(); };
inline Real     get_mp_scalar_real(const mp_float& s)   { return s.cast_float(); };
inline Complex  get_mp_scalar_real(const mp_complex& s) { return s.cast_complex(); };
inline Real     get_mp_scalar_real(const mp_rational& s){ return s.cast_float(); };

template<class To, class From>
struct convert_between_mp;

template<>
struct convert_between_mp<mp_int, mp_int>
{
    using To    = mp_int;
    using From  = mp_int;
    static const To& eval(const From& s)   { return s; };
};
template<>
struct convert_between_mp<mp_int, mp_float>
{
    using To    = mp_int;
    using From  = mp_float;
    static To eval(const From& s)   { return s.cast_mp_int(); };
};
template<>
struct convert_between_mp<mp_int, mp_rational>
{
    using To    = mp_int;
    using From  = mp_rational;
    static To eval(const From& s)   { return s.cast_mp_int(); };
};
template<>
struct convert_between_mp<mp_int, mp_complex>
{
    using To    = mp_int;
    using From  = mp_complex;
    static To eval(const From& s)   { return s.cast_mp_int(); };
};

template<>
struct convert_between_mp<mp_float, mp_int>
{
    using To    = mp_float;
    using From  = mp_int;
    static To eval(const From& s)   { return s.cast_mp_float(); };
};
template<>
struct convert_between_mp<mp_float, mp_float>
{
    using To    = mp_float;
    using From  = mp_float;
    static const To& eval(const From& s)   { return s; };
};
template<>
struct convert_between_mp<mp_float, mp_rational>
{
    using To    = mp_float;
    using From  = mp_rational;
    static To eval(const From& s)   { return s.cast_mp_float(); };
};
template<>
struct convert_between_mp<mp_float, mp_complex>
{
    using To    = mp_float;
    using From  = mp_complex;
    static To eval(const From& s)   { return s.cast_mp_float(); };
};

template<>
struct convert_between_mp<mp_rational, mp_int>
{
    using To    = mp_rational;
    using From  = mp_int;
    static To eval(const From& s)   { return s.cast_mp_rational(); };
};
template<>
struct convert_between_mp<mp_rational, mp_float>
{
    using To    = mp_rational;
    using From  = mp_float;
    static To eval(const From& s)   { return s.cast_mp_rational(); };
};
template<>
struct convert_between_mp<mp_rational, mp_rational>
{
    using To    = mp_rational;
    using From  = mp_rational;
    static const To& eval(const From& s)   { return s; };
};
template<>
struct convert_between_mp<mp_rational, mp_complex>
{
    using To    = mp_rational;
    using From  = mp_complex;
    static To eval(const From& s)   { return s.cast_mp_rational(); };
};

template<>
struct convert_between_mp<mp_complex, mp_int>
{
    using To    = mp_complex;
    using From  = mp_int;
    static To eval(const From& s)   { return s.cast_mp_complex(); };
};
template<>
struct convert_between_mp<mp_complex, mp_float>
{
    using To    = mp_complex;
    using From  = mp_float;
    static To eval(const From& s)   { return s.cast_mp_complex(); };
};
template<>
struct convert_between_mp<mp_complex, mp_rational>
{
    using To    = mp_complex;
    using From  = mp_rational;
    static To eval(const From& s)   { return s.cast_mp_complex(); };
};
template<>
struct convert_between_mp<mp_complex, mp_complex>
{
    using To    = mp_complex;
    using From  = mp_complex;
    static const To& eval(const From& s)   { return s; };
};

template<class To, class From>
struct convert_mp_scalar_impl
{
    static To eval(const From& s)
    {
        return convert_between_mp<To,From>::eval(s);
    }
};

template<class From>
struct convert_mp_scalar_impl<Integer, From>    
{
    static Integer eval(const From& s)  { return get_mp_scalar_int(s); };
};

template<class From>
struct convert_mp_scalar_impl<Real, From>    
{
    static Real eval(const From& s)     { return convert_scalar<Real>(get_mp_scalar_real(s)); };
};

template<class From>
struct convert_mp_scalar_impl<Float, From>    
{
    static Float eval(const From& s)    
    { 
        Real tmp = convert_mp_scalar_impl<Real,From>::eval(s);
        return convert_scalar<Float>(tmp); 
    };
};

template<class From>
struct convert_mp_scalar_impl<Complex, From>
{
    static Complex eval(const From& s)
    { 
        return convert_scalar<Complex>(get_mp_scalar_real(s)); 
    };
};

template<class From>
struct convert_mp_scalar_impl<Float_complex, From>    
{
    static Float_complex eval(const From& s)
    { 
        return convert_scalar<Float_complex>(get_mp_scalar_real(s)); 
    };
};

template<class To, class From>
struct convert_scalar_impl
{
    static To eval(const From& s)
    {
        return convert_mp_scalar_impl<To, From>::eval(s);
    }
};

template<class To>
struct convert_scalar_impl<To, Integer> { static To eval(Integer s) { return To(s); } };

template<class To>
struct convert_scalar_impl<To, Float>   { static To eval(Float s)   { return To(s); } };

template<class To>
struct convert_scalar_impl<To, Real>    { static To eval(Real s)    { return To(s); } };

template<>
struct convert_scalar_impl<mp_complex, Complex> 
{ 
    static mp_complex eval(const Complex& s)    { return mp_complex(s); } 
};
template<>
struct convert_scalar_impl<mp_rational, Complex>
{ 
    static mp_rational eval(const Complex& s)   { return mp_rational(real(s)); }
};
template<>
struct convert_scalar_impl<mp_float, Complex>
{ 
    static mp_float eval(const Complex& s)      { return mp_float(real(s)); }
};
template<>
struct convert_scalar_impl<mp_int, Complex>
{ 
    static mp_int eval(const Complex& s)        { return mp_int(real(s)); }
};

template<>
struct convert_scalar_impl<mp_complex, Float_complex>
{ 
    static mp_complex eval(const Float_complex& s)  { return mp_complex(s); }
};
template<>
struct convert_scalar_impl<mp_rational, Float_complex>
{ 
    static mp_rational eval(const Float_complex& s)  { return mp_rational(real(s)); }
};
template<>
struct convert_scalar_impl<mp_float, Float_complex>
{ 
    static mp_float eval(const Float_complex& s)  { return mp_float(real(s)); }
};
template<>
struct convert_scalar_impl<mp_int, Float_complex>
{ 
    static mp_int eval(const Float_complex& s)  { return mp_int(real(s)); }
};

}};}

#pragma warning(pop)

namespace matcl
{

template<class To, class From>
To matcl::convert_scalar(const From& s, typename mp::details::enable_mp_bin<To,From,void*>::type)
{
    return matcl::mp::details::convert_scalar_impl<To, From>::eval(s);
}

};
