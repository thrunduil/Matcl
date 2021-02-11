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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-internals/base/vector.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/matrix/matrix.h"

#include <vector>

namespace matcl { namespace raw { namespace details
{
    template<class V> class spptr_helper;
    template<class V> struct array_helper;

};};};

namespace matcl { namespace details
{

//indices are 1-based
inline void pos2ind(Integer pos, Integer s, Integer& i, Integer& j)
{
    --pos;
    i = (pos%s);
    j = (pos/s);
};

template<class T>
struct default_value_impl
{
    static T eval(ti::ti_empty )
    {
        return T();
    };
};

template<>
struct default_value_impl<Object>
{
    static Object eval(ti::ti_object ti)
    {
        return Object(ti);
    };
};

template<class T>
struct default_value_impl<matcl::raw::details::array_helper<T>>
{
    using value_type = matcl::raw::details::array_helper<T>;
    static value_type eval(ti::ti_type<T> ti)
    {
        return value_type(ti);
    };
};

template<class V, class S>
struct default_value_impl<matcl::raw::Matrix<V,S>>
{
    static matcl::raw::Matrix<V,S> eval(ti::ti_type<V> ti)
    {
        return matcl::raw::Matrix<V,S>(ti);
    };
};

template<class V>
struct default_value_impl<matcl::raw::details::spptr_helper<V>>
{
    using value_type = matcl::raw::details::spptr_helper<V>;

    static value_type eval(ti::ti_type<V> ti)
    {
        return value_type(ti);
    };
};

template<class T> 
inline T default_value(typename ti::get_ti_type<T>::type ti)
{
    return default_value_impl<T>::eval(ti);
};

template<class V> inline V one_value(ti::ti_type<V> )          { return V(1); };
template<> inline Object   one_value<Object>(ti::ti_object ti) { return Object::make_one(ti); };


template<class val_type,class struct_type>
struct zero_matrix
{
    using matrix_type = raw::Matrix<val_type,struct_type>;
    static matrix_type eval(Integer m,Integer n)
    {};
};

template<class val_type>
struct zero_matrix<val_type,struct_dense>
{
    using matrix_type   = raw::Matrix<val_type,struct_dense>;
    using tinfo_t       = ti::ti_type<val_type>;

    static matrix_type eval(tinfo_t ret_ti, Integer m,Integer n)
    {
        val_type Z = default_value<val_type>(ret_ti);
        matrix_type out(ret_ti,Z,m,n);
        out.get_struct().set(predefined_struct_type::diag);
        return out;
    };
};

template<class val_type>
struct zero_matrix<val_type,struct_sparse>
{
    using matrix_type   = raw::Matrix<val_type,struct_sparse>;
    using tinfo_t       = ti::ti_type<val_type>;

    static matrix_type eval(tinfo_t ret_ti, Integer m,Integer n)
    {
        return matrix_type(ret_ti, m,n);
    };
};

template<class val_type>
struct zero_matrix<val_type,struct_banded>
{
    using matrix_type   = raw::Matrix<val_type,struct_banded>;
    using tinfo_t       = ti::ti_type<val_type>;

    static matrix_type eval(tinfo_t ret_ti,Integer m,Integer n)
    {
        matrix_type out(ret_ti,m,n);
        val_type Z = default_value<val_type>(ret_ti);

        Integer out_size = out.impl_size();
        val_type* ptr_out = out.rep_ptr();

        for  (Integer i = 0; i < out_size; ++i)
            ptr_out[i] = Z;

        out.get_struct().set(predefined_struct_type::diag);
        return out;
    };
};

template<class val_type,class struct_type>
struct MATCL_MATREP_EXPORT unit_matrix
{};

template<class val_type>
struct MATCL_MATREP_EXPORT unit_matrix<val_type,struct_dense>
{
    using matrix_type   = raw::Matrix<val_type,struct_dense>;
    using tinfo_t       = ti::ti_type<val_type>;

    static matrix_type eval(tinfo_t ret_ti, Integer m);
};

template<>
struct MATCL_MATREP_EXPORT unit_matrix<Object,struct_dense>
{
    using matrix_type   = raw::Matrix<Object,struct_dense>;
    using tinfo_t       = ti::ti_type<Object>;

    static matrix_type eval(tinfo_t ret_ti, Integer m);
};

template<class val_type>
struct MATCL_MATREP_EXPORT unit_matrix<val_type,struct_sparse>
{
    using matrix_type   = raw::Matrix<val_type,struct_sparse>;
    using tinfo_t       = ti::ti_type<val_type>;

    static matrix_type eval(tinfo_t ret_ti, Integer m);
};

template<>
struct MATCL_MATREP_EXPORT unit_matrix<Object,struct_sparse>
{
    using matrix_type   = raw::Matrix<Object,struct_sparse>;
    using tinfo_t       = ti::ti_type<Object>;

    static matrix_type eval(tinfo_t ret_ti, Integer m);
};

template<class val_type>
struct MATCL_MATREP_EXPORT unit_matrix<val_type,struct_banded>
{
    using matrix_type   = raw::Matrix<val_type,struct_banded>;
    using tinfo_t       = ti::ti_type<val_type>;

    static matrix_type eval(tinfo_t ret_ti, Integer m);
};

template<>
struct MATCL_MATREP_EXPORT unit_matrix<Object,struct_banded>
{
    using matrix_type   = raw::Matrix<Object,struct_banded>;
    using tinfo_t       = ti::ti_type<Object>;

    static matrix_type eval(tinfo_t ret_ti, Integer m);
};

template<class T>
struct ti_initializer 
{ };

template<>
struct ti_initializer<Object>
{
    private:
        ti::ti_object   m_ti;

    public:
        ti_initializer(ti::ti_object ti)
            :m_ti(ti)
        {};

        template<class T>
        void initialize(T* ptr) const
        {
            new(ptr) T(m_ti);
        };
};

template<class T, class T_TI = ti::ti_type<T>>
class workspace2 : public std::vector<T>
{
    private:
        using base_type = std::vector<T>;

    public:
        workspace2(ti::ti_empty)                 : base_type(){};
        workspace2(ti::ti_empty, size_t size)    : base_type(size){};
};

template<class T>
class workspace2<T,ti::ti_object> : public ::data_struct::vector<T,ti_initializer<Object>>
{
    private:
        using base_type = ::data_struct::vector<T,ti_initializer<Object>>;

    public:
        workspace2(ti::ti_object ti)  
            : base_type(ti_initializer<Object>(ti))
        {};
        workspace2(ti::ti_object ti,size_t size)
            : base_type(ti_initializer<Object>(ti),size)
        {};
};

template<>
class workspace2<Integer, ti::ti_type<Integer>> : public pod_workspace<Integer>
{
    public:
        using base_type = pod_workspace<Integer>;

    public:
        workspace2(size_t size)                 : base_type(size){};
        workspace2(ti::ti_empty, size_t size)   : base_type(size){};
        workspace2(ti::ti_empty)                : base_type(){};
};

template<>
class workspace2<Real, ti::ti_type<Real>> : public pod_workspace<Real>
{
    public:
        using base_type = pod_workspace<Real>;

    public:
        workspace2(size_t size)                 : base_type(size){};
        workspace2(ti::ti_empty, size_t size)   : base_type(size){};
        workspace2(ti::ti_empty)                : base_type(){};
};

template<>
class workspace2<Float, ti::ti_type<Float>> : public pod_workspace<Float>
{
    public:
        using base_type = pod_workspace<Float>;

    public:
        workspace2(size_t size)                 : base_type(size){};
        workspace2(ti::ti_empty, size_t size)   : base_type(size){};
        workspace2(ti::ti_empty)                : base_type(){};
};

template<>
class workspace2<Complex, ti::ti_type<Complex>> : public pod_workspace<pod_type<Complex>>
{
    private:
        using T     = Complex;
        using pod_T = pod_type<Complex>;

    public:
        using base_type = pod_workspace<pod_T>;

    public:
        workspace2(size_t size)                 : base_type(size){};
        workspace2(ti::ti_empty, size_t size)   : base_type(size){};
        workspace2(ti::ti_empty)                : base_type(){};

        T*          ptr()                                   { return (T*)base_type::ptr(); };
        const T*    ptr() const                             { return (const T*)base_type::ptr(); };

        T&          operator[](size_t pos)                  { return ptr()[pos]; };
        const T&    operator[](size_t pos) const            { return ptr()[pos]; };

        void        resize(size_t new_size)                 { base_type::resize(new_size); };
        void        resize(size_t new_size, const T& val)   { base_type::resize(new_size, (pod_T&)val); };
};

template<>
class workspace2<Float_complex, ti::ti_type<Float_complex>> : public pod_workspace<pod_type<Float_complex>>
{
    private:
        using T     = Float_complex;
        using pod_T = pod_type<Float_complex>;

    public:
        using base_type = pod_workspace<pod_T>;

    public:
        workspace2(size_t size)                 : base_type(size){};
        workspace2(ti::ti_empty, size_t size)   : base_type(size){};
        workspace2(ti::ti_empty)                : base_type(){};

        T*          ptr()                                   { return (T*)base_type::ptr(); };
        const T*    ptr() const                             { return (const T*)base_type::ptr(); };

        T&          operator[](size_t pos)                  { return ptr()[pos]; };
        const T&    operator[](size_t pos) const            { return ptr()[pos]; };

        void        resize(size_t new_size)                 { base_type::resize(new_size); };
        void        resize(size_t new_size, const T& val)   { base_type::resize(new_size, (pod_T&)val); };
};

template<class Lhs, class Rhs>
struct has_trivial_assignment
{};

template<class V, class S, class Rhs>
struct has_trivial_assignment<matcl::raw::Matrix<V,S>,Rhs>
{
    using matrix_type = matcl::raw::Matrix<V,S>;

    static bool eval(const matrix_type&, const Rhs&)
    {
        return true;
    };
};

template<class Rhs> struct has_trivial_assignment<Integer,Rhs>
{
    static bool eval(Integer, const Rhs&)           { return true; }; 
};

template<class Rhs> struct has_trivial_assignment<Real,Rhs>
{ 
    static bool eval(Real, const Rhs&)              { return true; }; 
};

template<class Rhs> struct has_trivial_assignment<Float,Rhs>
{ 
    static bool eval(Float, const Rhs&)             { return true; }; 
};

template<class Rhs> struct has_trivial_assignment<Complex,Rhs>
{ 
    static bool eval(Complex, const Rhs&)           { return true; }; 
};

template<class Rhs> struct has_trivial_assignment<Float_complex,Rhs>
{ 
    static bool eval(Float_complex, const Rhs&)     { return true; }; 
};

template<class Rhs> struct has_trivial_assignment<Object,Rhs>
{
    static bool eval(const Object& val, const Rhs& v)
    { 
        return ti::has_trivial_assignment(val.get_type(), ti::get_ti(v)); 
    }; 
};

template<class S, class Rhs>
struct has_trivial_assignment<matcl::raw::Matrix<Object,S>,Rhs>
{
    using matrix_type = matcl::raw::Matrix<Object,S>;

    static bool eval(const matrix_type& mat, const Rhs& val)
    {
        return ti::has_trivial_assignment(mat.get_type(), ti::get_ti(val)); 
    };
};

};};