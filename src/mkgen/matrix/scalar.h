/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2019
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
#include "mkgen/TODO/matrix/ct_matrix_details.h"
#include "mkgen/details/matrix/scalar_checks.h"
#include <iosfwd>

namespace matcl { namespace mkgen
{

namespace mkd = matcl::mkgen::details;

// compile time scalar value, which stores symbolic element in Data.
// Deps is a type representing dependencies from runtime values; must be
// specialization of dps type; Data must be specialization of scalar_data<T>
template<class Data, class Deps>
class ct_scalar
{
    private:
        // check arguments
        using check1    = typename mkd::check_scalar_data<Data>::type;
        using check2    = typename mkd::check_scalar_deps<Deps>::type;

    public:
        // Data argument
        using data_type     = Data;

        // Deps argument
        using deps_type     = Deps;

        // type of this class
        using this_type     = ct_scalar<Data, Deps>;

        // print scalar
        template<class Subs_Context>
        static void print(std::ostream& os, int prior)
        {
            Data::print<Subs_Context>(os,prior);
        };

        // perform computations, return type is a specialization of ct_scalar. Tag must be a
        // fresh unique tag. If value of this scalar is already known this has no effect. 
        // However if this scalar represents some computation (for example multiplication
        // of two scalars), then such computation will be performed explicitly (otherwise it
        // will be performed each time when given scalar is used)
        template<class Tag>
        static auto compute() -> typename make_evaled_scalar<Data, Deps, Tag>::type;

    public:
        //-----------------------------------------------------------------------
        //                      internal use only
        //-----------------------------------------------------------------------

        // append to Arr_List all arrays required by this scalar
        template<Integer Step, class Arr_List>
        using get_arrays    = typename mkd::get_arrays_scalar<Data, Deps, Step, Arr_List>::type;

        //get value of type Val associated with this element; called internally
        template<class Val, class Local_Storage>
        inline_lev_1
        static Val eval(const Local_Storage& ls)
        {
            return Data::eval<Val>(ls);
        };

        // evaluate expression using data in arrays starting at position off and store
        // result in ret variable. TODO
        template<class Loop_Storage, class Ret, class Local_Storage>
        inline_lev_1
        static void eval_loop(Ret& ret, Integer off, const Local_Storage& cont)
        {
            mkd::eval_loop_scalar<Loop_Storage, Data, Deps>::eval<Ret>(ret,off,cont);
        };

        //TODO
        /*
        template<class Visitor>
        static void accept(Visitor& vis)
        {
            return Data::accept<Visitor>(vis);
        };
        */
};

//------------------------------------------------------------------------------
//                      Value scalars
//------------------------------------------------------------------------------
// represent integer value Value
template<Integer Value>
using integer_scalar = ct_scalar<mkd::scalar_data<mkd::scal_data_rational<Value, 1>>, 
                                    empty_deps>;

// represent rational value N / D
template<Integer N, Integer D>
using rational_scalar = ct_scalar<mkd::scalar_data<mkd::scal_data_rational<N,D>>, 
                                    empty_deps>;

// scalar storing a value of type Value_type defined by the tag Tag; value
// cannot depend on external data
template<class Tag, class Value_type>
using value_scalar  = ct_scalar<mkd::scalar_data<mkd::scal_data_value<Tag, Value_type>>, 
                                    empty_deps>;

//------------------------------------------------------------------------------
//                      Predefined scalars
//------------------------------------------------------------------------------
using zero  = integer_scalar<0>;
using one   = integer_scalar<1>;
using two   = integer_scalar<2>;
using mone  = integer_scalar<-1>;
using mtwo  = integer_scalar<-2>; 
using half  = rational_scalar<1,2>;

//------------------------------------------------------------------------------
//                      isa function
//------------------------------------------------------------------------------
 
// return true if T is a scalar
template<class T> 
struct is_scalar                        {static const bool value = false; };

template<class T, class D> 
struct is_scalar<ct_scalar<T,D>>        {static const bool value = true; };

// return true if T is a scalar with compile time known value
// TODO: remove?
template<class T> 
struct is_value_scalar                  {static const bool value = false; };

template<Integer N, Integer D> 
struct is_value_scalar<ct_scalar<mkd::scalar_data<mkd::scal_data_rational<N,D>>, 
                                    empty_deps>>
                                        {static const bool value = true; };

//------------------------------------------------------------------------------
//                      get_scalar_value
//------------------------------------------------------------------------------

// get value stored in scalar Scal of type Val.
// Scal must be a value scalar, i.e. integer_scalar, rational_scalar, or value_scalar
template<class Scal>
struct get_scalar_value
{
    static_assert(details::dependent_false<Scal>::value,
                  "this type should not be instantiated");

    template<class Val>
    static Val get();
};

}}

#include "mkgen/details/matrix/scalar_impl.h"