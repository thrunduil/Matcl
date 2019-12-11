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

//TODO: remove this type, not needed
template<class Data>
struct scalar_data
{    
    using check1    = typename details::check_scalar_data_impl<Data>::type;

    using base_data = Data;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Data::print<Subs_Context>(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage& ls)
    {
        return Data::eval<Val, Local_Storage>(ls);
    };
};

//-----------------------------------------------------------------------
//                      scal_data_evaled
//-----------------------------------------------------------------------
template<class Data, class Tag>
struct scal_data_evaled
{
    // check arguments
    using check1    = typename details::check_scalar_data_impl<Data>::type;
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
struct scal_data_rational
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
struct scal_data_value
{
    //TODO: add checks

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
//                      matrix element as scalar
//------------------------------------------------------------------------------
template<class T, Integer Row, Integer Col>
struct scalar_mat_elem_2 
{
    static_assert(dependent_false<T>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, Integer Row, Integer Col>
struct scalar_mat_elem_2<ct_matrix<M, N, Array_t, Deps>, Row, Col>
{
    // TODO
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        using elem = typename get_array_elem<Array_t,Row, col>::type;
        elem::print<Subs_Context>(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        using elem = typename get_array_elem<Array_t, Row, col>::type;
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
{
    //TODO
    static const Integer col = (Pos-1)/M + 1;
    static const Integer row = Pos - (col-1) * M;

    using matrix_type   = ct_matrix<M,N,Array_t, Deps>;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        using elem = typename get_array_elem<Array_t,row, col>::type;
        elem::print<Subs_Context>(os,prior);
    };
    
    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        using elem = typename get_array_elem<Array_t,row, col>::type;
        return elem::eval<Val>(ls);
    };
    
    // TODO
    /*
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using elem = typename get_array_elem<Array_t,row, col>::type;
        return elem::accept<Visitor>(vis);
    };
    */
};

//TODO
template<class Array, class Deps>
struct scalar_ctrans_array{};

template<class Tag, class Array, class Deps>
struct scalar_ufunc_array{};

template<class Tag, class Array1, class Array2>
struct scalar_bfunc_array{};

}}}
