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
#include "mkgen/matrix/base_types.h"

namespace matcl { namespace mkgen { namespace details
{

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
template<Tag_scalar_const_value Tag, class Value_type>
struct scal_data_const_value : mkd::scalar_data<scal_data_const_value<Tag, Value_type>>
{
    using this_type = scal_data_const_value<Tag, Value_type>;

    template<class Subs_Context>
    static void print(std::ostream& os,int prior)
    { 
        (void)prior;
        os << "const(" << Tag :: value<Value_type>() << ")";
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
template<Tag_scalar_value Tag, class Value_type>
struct scal_data_value : mkd::scalar_data<scal_data_value<Tag, Value_type>>
{
    using this_type = scal_data_value<Tag, Value_type>;

    template<class Subs_Context>
    static void print(std::ostream& os,int prior)
    { 
        (void)prior;
        os << "scalar(" << Tag :: value<Value_type>() << ")";
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
template<Tag_scalar_gen_value Tag>
struct scal_data_gen_value : mkd::scalar_data<scal_data_gen_value<Tag>>
{
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
template<Scal_data Data, Tag_comp Tag>
struct scal_data_evaled : mkd::scalar_data<scal_data_evaled<Data, Tag>>
{
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
struct is_value_scalar_data<mkd::scal_data_const_value<Tag, Val>>
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
