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

#include "mkgen/matrix/scalar.h"
#include "mkgen/details/matrix/scalar_data.h"
#include "mkgen/details/utils/utils.h"

//TODO: move to details
namespace matcl { namespace mkgen
{

//-----------------------------------------------------------------------
//                          make_evaled_scalar
//-----------------------------------------------------------------------
template<class Scalar_data, class Deps, class Tag>
struct make_evaled_scalar
{
    static_assert(details::dependent_false<Scalar_data>::value, 
                "this type should not be instantiated");
};

template<class Scalar_data, class ... Deps, class Tag>
struct make_evaled_scalar<Scalar_data, dps<Deps...>, Tag>
{
    using scalar    = ct_scalar<Scalar_data, dps<Deps...>>;
    using new_dep   = scalar_dep<Tag>;

    using base_data = typename Scalar_data::base_data;
    using array_base= mkd::scal_data_evaled<base_data, Tag>;
    using array     = mkd::scalar_data<array_base>;
    
    using type      = ct_scalar<array, dps<new_dep, Deps...>>;

    template<class Val, class Local_Storage>
    inline_lev_1
    friend Val eval_scalar(Tag, Local_Storage* ls)
    {
        return scalar::eval<Val>(*ls);
    };
};

template<Integer N, Integer D, class Tag>
struct make_evaled_scalar<mkd::scalar_data<mkd::scal_data_rational<N,D>>, 
                            empty_deps, Tag>
{
    // nothing to compute
    using data      = mkd::scal_data_rational<N,D>;
    using type      = ct_scalar<mkd::scalar_data<data>, empty_deps>;
};

template<class Data, class Tag_data, class Deps, class Tag>
struct make_evaled_scalar<mkd::scalar_data<mkd::scal_data_evaled<Data, Tag_data>>, 
                            Deps, Tag>
{
    // nothing to compute
    using data      = mkd::scal_data_evaled<Data, Tag_data>;
    using type      = ct_scalar<mkd::scalar_data<data>, Deps>;
};

template<class Data, class Value_type, class Deps, class Tag>
struct make_evaled_scalar<mkd::scalar_data<mkd::scal_data_value<Data, Value_type>>, 
                            Deps, Tag>
{
    // nothing to compute
    using data      = mkd::scal_data_value<Data,Value_type>;
    using type      = ct_scalar<mkd::scalar_data<data>, Deps>;
};

//------------------------------------------------------------------------------
//                      get_scalar_value
//------------------------------------------------------------------------------
template<Integer N, class Deps>
struct get_scalar_value<ct_scalar<mkd::scalar_data<mkd::scal_data_rational<N,1>>,
                                    Deps>>
{
    template<class Val>
    static Val get()    { return Val(N); };
};

template<Integer N, Integer D, class Deps>
struct get_scalar_value<ct_scalar<mkd::scalar_data<mkd::scal_data_rational<N,D>>,
                                    Deps>>
{
    template<class Val>
    static Val get()    { return Val(N) / Val(D); };
};

template<class Tag, class VT, class Deps>
struct get_scalar_value<ct_scalar<mkd::scalar_data<mkd::scal_data_value<Tag, VT>>,
                                    Deps>>
{
    using rep   = mkd::scal_data_value<Tag, VT>;

    template<class Val>
    static Val get()    { return rep::eval<Val>(); };
};

}}

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                         get_arrays_scalar
//-----------------------------------------------------------------------

//TODO: 
template<class Data, class Deps, Integer Step, class Arr_List>
struct get_arrays_scalar
{
    using scalar    = ct_scalar<Data,Deps>;
    using item      = array_item<scalar, Step, array_item_scalar>;
    using type      = typename push_back<Arr_List,item> :: type;
};

template<Integer N, Integer D, class Deps, Integer Step, class Arr_List>
struct get_arrays_scalar<mkd::scalar_data<mkd::scal_data_rational<N,D>>,
                            Deps, Step, Arr_List>
{
    using type      = Arr_List;    
};

template<class Data, class VT, class Deps, Integer Step, class Arr_List>
struct get_arrays_scalar<scal_data_value<Data, VT>, Deps,Step,Arr_List>
{
    using type      = Arr_List;
};

//-----------------------------------------------------------------------
//                         eval_loop_scalar
//-----------------------------------------------------------------------
template<class Loop_Storage, class Data, class Deps>
struct eval_loop_scalar
{
    template<class Ret, class Local_Storage>
    inline_lev_1
    static void eval(Ret& ret, Integer off, const Local_Storage& cont)
    {
        (void)off;
        ret = Data::eval<Ret, Local_Storage>(cont);
    };
};

}}}
