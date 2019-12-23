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
    using check_scalar_data = typename details::check_scalar_data_impl<Data, Dummy>::type;

    // ct_scalar arrays must implement:
    // template<class Subs_Context>
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);
    //
    // template<class Visitor>
    // static void accept(Visitor& vis);
    //
    // template<class Void>
    // using simplify;
    //
    // template<class Void>
    // static const bool is_simplified;

    //TODO: remove or implement
    // append to Arr_List all arrays required by this scalar
    template<Integer Step, class Arr_List>
    using get_arrays    = Arr_List;
};

//------------------------------------------------------------------------------
//                      scal_data_rational
//------------------------------------------------------------------------------
// represent rational value N / D, where D > 0
template<Integer N, Integer D>
struct scal_data_rational : mkd::scalar_data<scal_data_rational<N, D>>
{
    static_assert(D > 0, " denominator must be positive");

    using this_type = scal_data_rational<N, D>;

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

    template<class Val>
    static constexpr Val value()
    {
        return Val(double(N) / double(D)); 
    };

    template<class Visitor>
    static void accept(Visitor&)
    { 
        return; 
    };

    template<class Void>
    using simplify  = this_type;

    static constexpr bool is_simplified()   { return true; };

    // TODO: remove?
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_lev_1
    static void eval_loop(Ret& ret, Integer off, const Local_Storage& cont)
    {
        mkd::eval_loop_scalar<Loop_Storage, this_type>::eval<Ret>(ret,off,cont);
    };
};

//------------------------------------------------------------------------------
//                      scal_data_const_value
//------------------------------------------------------------------------------
template<class Tag, class Value_type>
struct scal_data_const_value : mkd::scalar_data<scal_data_const_value<Tag, Value_type>>
{
    //check
    using check1    = typename mkd::check_valid_const_data_tag<Tag>::type;

    using this_type = scal_data_const_value<Tag, Value_type>;

    template<class Subs_Context>
    static void print(std::ostream& os,int prior)
    { 
        Tag::print(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage&)
    { 
        return Val(Tag::value<Value_type>());
    };

    template<class Val>
    static constexpr Val value()
    {
        return Tag::value<Value_type>();
    };

    template<class Visitor>
    static void accept(Visitor&)
    { 
        return; 
    };

    template<class Void>
    using simplify  = this_type;

    static constexpr bool is_simplified()   { return true; };
};

//------------------------------------------------------------------------------
//                      scal_data_value
//------------------------------------------------------------------------------
template<class Tag, class Value_type>
struct scal_data_value : mkd::scalar_data<scal_data_value<Tag, Value_type>>
{
    //check
    using check1    = typename mkd::check_valid_data_tag<Tag>::type;

    using this_type = scal_data_value<Tag, Value_type>;

    template<class Subs_Context>
    static void print(std::ostream& os,int prior)
    { 
        Tag::print(os,prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage&)
    { 
        return Val(Tag::value<Value_type>());
    };

    template<class Val>
    static Val value()
    {
        return Tag::value<Value_type>();
    };

    template<class Visitor>
    static void accept(Visitor&)
    { 
        return; 
    };

    template<class Void>
    using simplify  = this_type;

    static constexpr bool is_simplified()   { return true; };
};

//------------------------------------------------------------------------------
//                      scal_data_gen_value
//------------------------------------------------------------------------------
template<class Tag>
struct scal_data_gen_value : mkd::scalar_data<scal_data_gen_value<Tag>>
{
    //check
    using check2    = typename mkd::check_valid_gen_data_tag<Tag>::type;

    using this_type = scal_data_gen_value<Tag>;

    template<class Subs_Context>
    static void print(std::ostream& os,int prior)
    { 
        Tag::print(os, prior);
    };

    template<class Val, class Local_Storage>
    inline_lev_1
    static Val eval(const Local_Storage& ls)
    { 
        return Tag::eval<Val, Local_Storage>(ls);
    };

    template<class Visitor>
    static void accept(Visitor&)
    { 
        return; 
    };

    template<class Void>
    using simplify  = this_type;

    static constexpr bool is_simplified()   { return true; };
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

    using this_type = scal_data_evaled<Data, Tag>;

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

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        return Data::accept<Visitor>(vis);
    };

    template<class Void>
    using simplify  = this_type;

    static constexpr bool is_simplified()   { return true; };

    //TODO: remove?
    // append to Arr_List all arrays required by this scalar
    template<Integer Step, class Arr_List>
    using get_arrays    = Arr_List;

    // TODO: remove?
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_lev_1
    static void eval_loop(Ret& ret, Integer off, const Local_Storage& cont)
    {
        mkd::eval_loop_scalar<Loop_Storage, Data>::eval<Ret>(ret,off,cont);
    };
};

//------------------------------------------------------------------------------
//                      matrix element as scalar
//------------------------------------------------------------------------------
template<class T, Integer Row, Integer Col>
struct scalar_mat_elem_2 
{
    static_assert(md::dependent_false<T>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, Integer Row, Integer Col>
struct scalar_mat_elem_2<ct_matrix<M, N, Array_t, Deps>, Row, Col>
             : mkd::scalar_data<scalar_mat_elem_2<ct_matrix<M, N, Array_t, Deps>, Row, Col>>
{
    using elem_type = typename Array_t::template get_element<Row, Col>::type;
    using this_type = scalar_mat_elem_2<ct_matrix<M, N, Array_t, Deps>, Row, Col>;

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

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        return elem_type::accept<Visitor>(vis);
    };

    //TODO
    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };
};

template<class T, Integer Pos>
struct scalar_mat_elem_1 
{
    static_assert(md::dependent_false<T>::value, "class T must be ct_matrix");
};

template<Integer M, Integer N, class Array_t, class Deps, Integer Pos>
struct scalar_mat_elem_1<ct_matrix<M, N, Array_t, Deps>, Pos>
            : mkd::scalar_data<scalar_mat_elem_1<ct_matrix<M, N, Array_t, Deps>, Pos>>
{
    static const Integer col = (Pos-1)/M + 1;
    static const Integer row = Pos - (col-1) * M;

    using elem_type     = typename Array_t::template get_element<row, col>::type;
    using this_type     = scalar_mat_elem_1<ct_matrix<M, N, Array_t, Deps>, Pos>;

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
    
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        return elem_type::accept<Visitor>(vis);
    };

    //TODO
    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };
};

//TODO
template<class Tag, class Array, class Deps>
struct scalar_ufunc_array : mkd::scalar_data<scalar_ufunc_array<Tag, Array, Deps>>
{};

template<class Tag, class Array1, class Array2>
struct scalar_bfunc_array: mkd::scalar_data<scalar_bfunc_array<Tag, Array1, Array2>>
{};

//------------------------------------------------------------------------------
//                      predefined scalar_data
//------------------------------------------------------------------------------
using zero_sd   = scal_data_rational<0, 1>;
using one_sd    = scal_data_rational<1, 1>;
using mone_sd   = scal_data_rational<-1, 1>;

//------------------------------------------------------------------------------
//                      isa functions
//------------------------------------------------------------------------------
template<class T> 
struct is_value_scalar_data             {static const bool value = false; };

template<Integer N, Integer D> 
struct is_value_scalar_data<mkd::scal_data_rational<N,D>>
                                        {static const bool value = true; };

template<class Tag, class Val> 
struct is_value_scalar_data<mkd::scal_data_const_value<Tag,Val>>
                                        {static const bool value = true; };

template<class Tag, class Val> 
struct is_value_scalar_data<mkd::scal_data_value<Tag,Val>>
                                        {static const bool value = true; };

// return true if T is scalar_data that represents value 0
template<class T>
struct is_scalar_data_zero              { static const bool value = false; };

template<>
struct is_scalar_data_zero<mkd::scal_data_rational<0, 1>>
                                        { static const bool value = true; };

template<class Tag, class Val>
struct is_scalar_data_zero<mkd::scal_data_const_value<Tag, Val>>
                                        { static const bool value = (Tag::value<Val>() == Val(0)); };

// return true if T is scalar_data that represents value 1
template<class T>
struct is_scalar_data_one               { static const bool value = false; };

template<>
struct is_scalar_data_one<mkd::scal_data_rational<1, 1>>
                                        { static const bool value = true; };

template<class Tag, class Val>
struct is_scalar_data_one<mkd::scal_data_const_value<Tag, Val>>
                                        { static const bool value = (Tag::value<Val>() == Val(1)); };


// return true if T is scalar_data that represents value -1
template<class T>
struct is_scalar_data_mone              { static const bool value = false; };

template<>
struct is_scalar_data_mone<mkd::scal_data_rational<-1, 1>>
                                        { static const bool value = true; };

template<class Tag, class Val>
struct is_scalar_data_mone<mkd::scal_data_const_value<Tag, Val>>
                                        { static const bool value = (Tag::value<Val>() == Val(-1)); };

}}}

namespace matcl { namespace mkgen 
{

//-----------------------------------------------------------------------
//                      scal_data_const_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating const_value_scalar
template<class Tag>
struct scal_data_const_value_tag
{
    //check arguments
    template<class Dummy>
    using check_scal_data_const_value_tag
                    = typename details::check_const_data_tag_impl<Tag, Dummy>::type;

    // Tag arrays must implement:
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val>
    // static constexpr Val value();
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
    using check_scal_data_value_tag = typename details::check_data_tag_impl<Tag, Dummy>::type;

    // Tag arrays must implement:
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val>
    // static Val value();
};

//-----------------------------------------------------------------------
//                      scal_data_gen_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating gen_scalar
template<class Tag>
struct scal_data_gen_value_tag
{
    //check arguments
    template<class Dummy>
    using check_scal_data_gen_value_tag = typename details::check_gen_data_tag_impl<Tag, Dummy>::type;

    // Tag arrays must implement:
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);
};

}}