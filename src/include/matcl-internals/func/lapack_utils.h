/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas-lapack/blas/blas.h"

namespace matcl { namespace details
{

template<class T> struct lusol_value_type                   { using type = T;} ;
template<>        struct lusol_value_type<Complex>          { using type = std::complex<double>;} ;
template<>        struct lusol_value_type<Float_complex>    { using type = std::complex<float>;} ;

template<class T> struct lapack_value_type                  { using type = T;} ;
template<>        struct lapack_value_type<Complex>         { using type = std::complex<double>;} ;
template<>        struct lapack_value_type<Float_complex>   { using type = std::complex<float>;} ;

template<bool cond> struct check_size_impl            {};
template<>          struct check_size_impl<true>      { static void eval(){}; };

template<class T1, class T2>
struct check_size : public check_size_impl<sizeof(T1) == sizeof(T2)> {};

inline lapack::i_type* lap(Integer* ptr)
{
    check_size<lapack::i_type,Integer>::eval();
    return reinterpret_cast<lapack::i_type*>(ptr);
};

inline lapack::d_type* lap(Real* ptr)
{
    return reinterpret_cast<lapack::d_type*>(ptr);
};

inline lapack::s_type* lap(Float* ptr)
{
    return reinterpret_cast<lapack::s_type*>(ptr);
};

inline lapack::z_type* lap(Complex* ptr)
{
    return reinterpret_cast<lapack::z_type*>(ptr);
};

inline lapack::c_type* lap(Float_complex* ptr)
{
    return reinterpret_cast<lapack::c_type*>(ptr);
};

inline const lapack::i_type* lap(const Integer* ptr)
{
    check_size<lapack::i_type,Integer>::eval();
    return reinterpret_cast<const lapack::i_type*>(ptr);
};

inline const lapack::d_type* lap(const Real* ptr)
{
    return reinterpret_cast<const lapack::d_type*>(ptr);
};

inline const lapack::s_type* lap(const Float* ptr)
{
    return reinterpret_cast<const lapack::s_type*>(ptr);
};

inline const lapack::z_type* lap(const Complex* ptr)
{
    return reinterpret_cast<const lapack::z_type*>(ptr);
};

inline const lapack::c_type* lap(const Float_complex* ptr)
{
    return reinterpret_cast<const lapack::c_type*>(ptr);
};

inline const lapack::d_type& lap(const Real& ptr)
{
    return *reinterpret_cast<const lapack::d_type*>(&ptr);
};

inline const lapack::s_type& lap(const Float& ptr)
{
    return *reinterpret_cast<const lapack::s_type*>(&ptr);
};

inline const lapack::z_type& lap(const Complex& ptr)
{
    return *reinterpret_cast<const lapack::z_type*>(&ptr);
};

inline const lapack::c_type& lap(const Float_complex& ptr)
{
    return *reinterpret_cast<const lapack::c_type*>(&ptr);
};

};};