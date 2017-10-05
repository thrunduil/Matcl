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

#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/object.h"
#include "matcl-core/general/thread.h"

#include <vector>

namespace matcl { namespace dynamic { namespace details
{

class evaler
{
    private:
        using func_vec  = std::vector<function>;
        using asize_t   = matcl::atomic<size_t>;

    public:
        mutable asize_t m_refcount;
        Type*           m_arg_ti;
        Type            m_ret_ti;
        int             m_args_size;        

    public:
        evaler();

        virtual ~evaler() {};

        void            increase_refcount() const;
        bool            decrease_refcount() const;

        virtual function    make_converter(int n_deduced, const Type deduced[],
                                Type deduced_ret, const func_vec& arg_conv) const = 0;

        // if function returns true, ret is a new object; otherwise ret is one
        // of existing argument (only fun_conv_null converter returns false) 
        virtual bool        make_eval(const object** args, object& ret) const = 0;

        // equivalent to make_eval, but return is ignored
        virtual void        make_eval(const object** args) const = 0;
};

struct evaler_traits
{
    private:
        using value_type    = evaler;

    public:
        static void     check_free(const value_type* p);
        static void     copy(const value_type* p);
        static void     assign(const value_type* to, const value_type* from);
};

inline void evaler_traits::check_free(const value_type* p)
{
    if (p && p->decrease_refcount())
        delete p;
};

inline void evaler_traits::copy(const value_type* p)
{
    if (p) 
        p->increase_refcount();
}

inline void evaler_traits::assign(const value_type* to, const value_type* from)
{
    copy(from);
    check_free(to);
};

inline evaler::evaler()
    :m_refcount(0)
{};

inline void evaler::increase_refcount() const
{
    ++m_refcount;
};

inline bool evaler::decrease_refcount() const
{
    return ((--m_refcount ) == 0);
}

template<class T>
struct get_object_type
{
    static Type eval(const T&)  
    { 
        return T::get_static_type(); 
    };
};

// no static type information
template<>
struct get_object_type<object>
{
    static Type eval(const object& obj)
    { 
        return obj.get_type(); 
    };
};


//get type of given argument; add reference if required
template<class Object>
struct get_arg_type_eval
{
    static Type eval(const Object& obj)
    {
        static_assert(is_object<Object>::value, "object type required");
        return get_object_type<Object>::eval(obj);
    };
};

template<class Object>
struct get_arg_type_eval<const Object&>
{
    static Type eval(const Object& obj)
    {
        static_assert(is_object<Object>::value, "object type required");
        return get_object_type<Object>::eval(obj);
    };
};

template<class Object>
struct get_arg_type_eval<Object&&>
{
    static Type eval(const Object& obj)
    {
        static_assert(is_object<Object>::value, "object type required");
        return get_object_type<Object>::eval(obj);
    };
};

template<class Object>
struct get_arg_type_eval<const Object&&>
{
    static Type eval(const Object& obj)
    {
        static_assert(is_object<Object>::value, "object type required");
        return get_object_type<Object>::eval(obj);
    };
};

template<class Object>
struct get_arg_type_eval<Object&>
{
    static Type eval(const Object&)
    {
        static_assert(is_object<Object>::value, "object type required");
        return mark_reference_type<Object>().get();
    };
};

template<>
struct MATCL_DYN_EXPORT get_arg_type_eval<object&>
{
    static Type eval(const object& obj);
};

};}};
