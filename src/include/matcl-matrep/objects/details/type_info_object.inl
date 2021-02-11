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

#include "matcl-matrep/objects/details/type_info_object.h"

namespace matcl { namespace ti
{

namespace details
{
    MATCL_MATREP_EXPORT ti_object  unify_ti_assign_obj(ti_object t1, ti_object t2);
    MATCL_MATREP_EXPORT ti_object  unify_ti_obj(ti_object t1, ti_object t2);
    MATCL_MATREP_EXPORT ti_object  get_return_ti_obj(const dynamic::function_name& fn, 
                                              ti_object t1, ti_object t2);
    MATCL_MATREP_EXPORT ti_object  get_return_ti_obj(const dynamic::function_name& fn, ti_object t1);

    inline ti_object to_obj(ti_object t)    { return t; };
    inline ti_object to_obj(ti_int)         { return predefined::get_ti_int(); };
    inline ti_object to_obj(ti_float)       { return predefined::get_ti_float(); };
    inline ti_object to_obj(ti_real)        { return predefined::get_ti_real(); };
    inline ti_object to_obj(ti_float_compl) { return predefined::get_ti_float_complex(); };
    inline ti_object to_obj(ti_compl)       { return predefined::get_ti_complex(); };

    template<class TI_r, class TI1, class TI2>
    struct unify_ti_assign_impl
    {
        static TI_r eval(TI1, TI2)
        { 
            return TI_r(); 
        };
    };

    template<class TI1, class TI2>
    struct unify_ti_assign_impl<ti_object,TI1,TI2>
    {
        static ti_object eval(TI1 t1, TI2 t2)
        { 
            return unify_ti_assign_obj(to_obj(t1),to_obj(t2)); 
        };
    };

    template<class TI_r, class TI1, class TI2>
    struct unify_ti_impl
    {
        static TI_r eval(TI1, TI2)
        { 
            return TI_r(); 
        };
    };

    template<class TI1, class TI2>
    struct unify_ti_impl<ti_object,TI1,TI2>
    {
        static ti_object eval(TI1 t1, TI2 t2)
        { 
            return unify_ti_obj(to_obj(t1),to_obj(t2)); 
        };
    };

    template<class TI_r, class TI1, class TI2>
    struct get_return_ti_impl
    {
        static TI_r eval(const matcl::dynamic::function_name&, TI1, TI2)
        { 
            return TI_r(); 
        };
    };

    template<class TI1, class TI2>
    struct get_return_ti_impl<ti_object,TI1,TI2>
    {
        static ti_object eval(const matcl::dynamic::function_name& fn, TI1 t1, TI2 t2)
        { 
            return get_return_ti_obj(fn,to_obj(t1),to_obj(t2)); 
        };
    };

    template<class TI_r, class TI1>
    struct get_return_ti1_impl
    {
        static TI_r eval(const matcl::dynamic::function_name&, TI1)
        { 
            return TI_r(); 
        };
    };

    template<class TI1>
    struct get_return_ti1_impl<ti_object,TI1>
    {
        static ti_object eval(const matcl::dynamic::function_name& fn, TI1 t1)
        { 
            return get_return_ti_obj(fn,to_obj(t1)); 
        };
    };

    template<class TI2>
    struct has_trivial_assignment_impl
    {
        static bool eval(ti_object t1, TI2 t2)
        {
            return dynamic::operations::has_trivial_assignment(t1,to_obj(t2));
        };
    };

    template<>
    struct has_trivial_assignment_impl<ti_object>
    {
        static bool eval(ti_object t1, ti_object t2)
        {
            return dynamic::operations::has_trivial_assignment(t1,t2);
        };
    };
};

template<class TI2>
bool ti::has_trivial_assignment(ti_object t1, TI2 t2)
{
    return details::has_trivial_assignment_impl<TI2>::eval(t1,t2);
};

inline predefined::Type predefined::get_ti_int()
{
    return matcl::dynamic::predefined::type_int();
}

inline predefined::Type predefined::get_ti_real()
{
    return matcl::dynamic::predefined::type_real();
}

inline predefined::Type predefined::get_ti_float()
{
    return matcl::dynamic::predefined::type_float();
}

inline predefined::Type predefined::get_ti_complex()
{
    return matcl::dynamic::predefined::type_complex();
}

inline predefined::Type predefined::get_ti_float_complex()
{
    return matcl::dynamic::predefined::type_float_complex();
}

template<class TI_ret, class TI1, class TI2>
inline TI_ret ti::unify_ti_assign(TI1 t1, TI2 t2)
{
    return details::unify_ti_assign_impl<TI_ret,TI1,TI2>::eval(t1,t2);
};

template<class TI_ret, class TI1, class TI2>
inline TI_ret unify_ti(TI1 t1, TI2 t2)
{
    return details::unify_ti_impl<TI_ret,TI1,TI2>::eval(t1,t2);
};

template<class TI_ret, class TI1, class TI2>
TI_ret ti::get_return_ti(const matcl::dynamic::function_name& fn, TI1 t1, TI2 t2)
{
    return details::get_return_ti_impl<TI_ret,TI1,TI2>::eval(fn,t1,t2);
}

template<class TI_ret, class TI1>
TI_ret ti::get_return_ti(const matcl::dynamic::function_name& fn, TI1 t1)
{
    return details::get_return_ti1_impl<TI_ret,TI1>::eval(fn,t1);
};

};};
