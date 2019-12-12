/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

#include "mkgen/details/mkgen_fwd.h"

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      scalar_data
//-----------------------------------------------------------------------
// base class for ct_scalar arrays
template<class Data>
struct scalar_data
{
    //check arguments
    template<class Dummy>
    using check     = typename details::check_scalar_data_impl<Data, Dummy>::type;

    // ct_scalar arrays must implement:
    // template<class Subs_Context>
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);
};

//-----------------------------------------------------------------------
//                      scal_data_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating value_scalar
template<class Tag>
struct scal_data_value_tag
{
    //check arguments
    template<class Dummy>
    using check     = typename details::check_scalar_data_tag_impl<Tag, Dummy>::type;

    // Tag arrays must implement:
    // template<class Subs_Context>
    // static void print(std::ostream& os, int prior)
    // {
    //     Data::print<Subs_Context>(os,prior);
    // };
    // 
    // template<class Value_type>
    // static Value_type eval<Value_type>();
};

//-----------------------------------------------------------------------
//                      scal_data_evaled
//-----------------------------------------------------------------------
template<class Data, class Tag>
struct scal_data_evaled : mkd::scalar_data<scal_data_evaled<Data, Tag>>
{
    // check arguments
    using check1    = typename details::check_scalar_data_impl<Data, void>::type;
    using check2    = typename details::check_computation_tag<Tag>::type;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Data::print<Subs_Context>(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage& ls)
    {
        return Val(ls.get_scalar<Tag>());
    };

    /*
    TODO
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        return Data::accept<Visitor>(vis);
    };
    */
};

//------------------------------------------------------------------------------
//                      scal_data_rational
//------------------------------------------------------------------------------
template<Integer N, Integer D>
struct scal_data_rational : mkd::scalar_data<scal_data_rational<N, D>>
{
    static_assert(D != 0, " denominator cannot be zero");

    template<class Subs_Context>
    static void print(std::ostream& os,int)
    { 
        os << N; 
        if constexpr(D != 1)
            os << "/" << D;
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage&)
    { 
        return Val(N) / Val(D); 
    };

    //TODO
    /*
    template<class Visitor>
    static void accept(Visitor&)
    { 
        return; 
    };
    */
};

//------------------------------------------------------------------------------
//                      scal_data_value
//------------------------------------------------------------------------------
template<class Tag, class Value_type>
struct scal_data_value : mkd::scalar_data<scal_data_value<Tag, Value_type>>
{
    //check
    using check1    = typename mkd::check_valid_scalar_data_tag<Tag>::type;

    template<class Subs_Context>
    static void print(std::ostream& os,int prior)
    { 
        Tag::print<Subs_Context>(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage&)
    { 
        Value_type v    = Tag::eval<Value_type>();
        return Val(v); 
    };

    // needed by get_scalar_value
    template<class Val>
    inline_lev_1
    static Val eval()
    { 
        Value_type v    = Tag::eval<Value_type>();
        return Val(v); 
    };

    //TODO
    /*
    template<class Visitor>
    static void accept(Visitor&)
    { 
        return; 
    };
    */
};

//------------------------------------------------------------------------------
//                      scalar_expr_data
//------------------------------------------------------------------------------

template<class Expr>
struct scalar_expr_data : mkd::scalar_data<scalar_expr_data<Expr>>
{
    // check
    using check_expr    = typename details::check_scalar_data_impl<Expr, void>::type;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Expr::print<Subs_Context>(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage& ls)
    {
        return Expr::eval<Val, Local_Storage>(ls);
    };
};

//------------------------------------------------------------------------------
//                      matrix element as scalar
//------------------------------------------------------------------------------
template<class T, Integer Row, Integer Col>
struct scalar_mat_elem_2 
{
    static_assert(dependent_false<T>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, Integer Row, Integer Col>
struct scalar_mat_elem_2<ct_matrix<M, N, Array_t, Deps>, Row, Col>
             : mkd::scalar_data<scalar_mat_elem_2<ct_matrix<M, N, Array_t, Deps>, Row, Col>>
{
    using elem_type = typename get_array_elem<Array_t,Row, Col>::type;

    // check
    using check_elem    = typename details::check_scalar_data_impl<elem_type, void>::type;
    
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        elem::print<Subs_Context>(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        return elem::eval<Val>(ls);
    };
};

template<class T, Integer Pos>
struct scalar_mat_elem_1 
{
    static_assert(dependent_false<T>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, Integer Pos>
struct scalar_mat_elem_1<ct_matrix<M, N, Array_t, Deps>, Pos>
            : mkd::scalar_data<scalar_mat_elem_1<ct_matrix<M, N, Array_t, Deps>, Pos>>
{
    static const Integer col = (Pos-1)/M + 1;
    static const Integer row = Pos - (col-1) * M;

    using elem_type     = typename get_array_elem<Array_t, row, col>::type;

    // check
    using check_elem    = typename details::check_scalar_data_impl<elem_type, void>::type;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {        
        elem_type::print<Subs_Context>(os,prior);
    };
    
    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        return elem_type::eval<Val>(ls);
    };
    
    // TODO
    /*
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        return elem_type::accept<Visitor>(vis);
    };
    */
};

//TODO
template<class Array, class Deps>
struct scalar_ctrans_array : mkd::scalar_data<scalar_ctrans_array<Array, Deps>>
{};

template<class Tag, class Array, class Deps>
struct scalar_ufunc_array : mkd::scalar_data<scalar_ufunc_array<Tag, Array, Deps>>
{};

template<class Tag, class Array1, class Array2>
struct scalar_bfunc_array: mkd::scalar_data<scalar_bfunc_array<Tag, Array1, Array2>>
{};

}}}
