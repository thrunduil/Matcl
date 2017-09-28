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

#include "details/function_traits.h"
#include "details/dynamic_function_utils.h"

namespace matcl { namespace dynamic { namespace details
{

template<typename Fun, class base_type_ = void>
class dynamic_function : public get_base_type<base_type_>::type
{
    protected:
        using function_traits   = details::function_traits<Fun>;
        using function_tag      = typename function_traits::function_type;
        using function          = typename details::proper_func_type<Fun,function_tag>::type;
        using input_type        = typename function_traits::input_type;
        using return_type       = typename function_traits::return_type;
        using base_type         = typename details::get_base_type<base_type_>::type;

        static const int arg_size = matcl::dynamic::details::size<input_type>::value;

    private:
        function            m_fun;

    public:
        dynamic_function(Fun f);

        //for nonmember functions
        template<class data_constructor>
        void                eval(data_constructor) const;

        //for member functions
        template<class class_type, class data_constructor>
        void                eval(class_type object, data_constructor) const;
};

};};};

#include "details/dynamic_function.inl"