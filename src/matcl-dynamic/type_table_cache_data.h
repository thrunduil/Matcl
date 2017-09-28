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
#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/type.h"
#include "matcl-dynamic/function.h"
#include <boost/functional/hash.hpp>

#pragma warning(push)
#pragma warning(disable: 4146) // unary minus operator applied to unsigned

namespace matcl { namespace dynamic { namespace details
{

enum class converter_type : int
{
    conv_implicit,
    conv_explicit,
    conv_cast
};

inline bool allow_explicit(converter_type ctype)
{
    return ctype == converter_type::conv_explicit
            || ctype == converter_type::conv_cast;
};

inline bool allow_cast(converter_type ctype)
{
    return ctype == converter_type::conv_cast;
};

struct unifier_info     { Type t1; Type t2; Type unif;};
struct convert_info     { Type t1; Type t2; converter_type ctype; function f;};
struct assign_info      { Type t1; Type t2; function f; };
struct overload_info    { function_name f; int n_args; const Type* t; function func;};
struct toverload_info   { function_name f; int n_templ; const Type* templ; int n_args; 
                            const Type* t; function func;};

class unifier_impl
{
    private:
        Type            m_t1;
        Type            m_t2;
        Type            m_unifier;
        size_t          m_hash;

    public:
        unifier_impl(unifier_info in);

        static size_t   eval_hash(unifier_info in);
        bool            equal(unifier_info in) const;

        Type            get() const             { return m_unifier; };
        std::size_t     hash_value() const      { return m_hash; };
};

class overload_impl
{
    private:
        function_name   m_name;
        int             n_args;
        Type*           m_args;
        function        m_overload;
        size_t          m_hash;

    public:
        overload_impl(const overload_info& in);
        ~overload_impl();

        static size_t   eval_hash(const overload_info& in);
        bool            equal(const overload_info& in) const;

        function        get() const             { return m_overload; };
        std::size_t     hash_value() const      { return m_hash; };

    private:
        overload_impl(const overload_impl&) = delete;
        overload_impl& operator=(const overload_impl&) = delete;
};

class toverload_impl
{
    private:
        function_name   m_name;
        int             n_template;
        Type*           m_templates;
        int             n_args;
        Type*           m_args;
        function        m_overload;
        size_t          m_hash;

    public:
        toverload_impl(const toverload_info& in);
        ~toverload_impl();

        static size_t   eval_hash(const toverload_info& in);
        bool            equal(const toverload_info& in) const;

        function        get() const             { return m_overload; };
        std::size_t     hash_value() const      { return m_hash; };

    private:
        toverload_impl(const toverload_impl&) = delete;
        toverload_impl& operator=(const toverload_impl&) = delete;
};

class convert_impl
{
    private:
        Type            m_to;
        Type            m_from;
        converter_type  m_type;
        function        m_converter;
        size_t          m_hash;

    public:
        convert_impl(convert_info in);

        static size_t   eval_hash(convert_info in);
        bool            equal(convert_info in) const;

        function        get() const             { return m_converter; };
        std::size_t     hash_value() const      { return m_hash; };

};

class assign_impl
{
    private:
        Type            m_to;
        Type            m_from;
        function        m_assigner;
        size_t          m_hash;

    public:
        assign_impl(assign_info in);

        static size_t   eval_hash(assign_info in);
        bool            equal(assign_info in) const;

        function        get() const             { return m_assigner; };
        std::size_t     hash_value() const      { return m_hash; };
};

template<class T>
struct simple_ptr_traits
{
    static void check_free(const T*) {};
    static void copy(const T*) {};
};

template<class T>
struct simple_ptr
{
    using value_type    = T*;
    using type_traits   = simple_ptr_traits<T>;

    const T*            m_ptr;

    static simple_ptr   make(const T* p) { return simple_ptr{p}; };
};

//----------------------------------------------------------------------
//                  unifier_impl
//----------------------------------------------------------------------
inline size_t unifier_impl::eval_hash(unifier_info in)
{
    size_t seed = in.t1.hash_value();
    boost::hash_combine(seed,1);
    boost::hash_combine(seed,in.t2.hash_value());

    return seed;
};

inline bool unifier_impl::equal(unifier_info in) const
{
    return m_t1 == in.t1 && m_t2 == in.t2;
}

//----------------------------------------------------------------------
//                  unifier_impl
//----------------------------------------------------------------------
inline unifier_impl::unifier_impl(unifier_info in)
    :m_t1(in.t1), m_t2(in.t2), m_unifier(in.unif)
{
    m_hash = eval_hash(in);
};

//----------------------------------------------------------------------
//                  overload_impl
//----------------------------------------------------------------------
inline overload_impl::overload_impl(const overload_info& in)
    :m_name(in.f), n_args(in.n_args), m_overload(in.func)
{
    m_hash  = eval_hash(in);
    m_args  = new Type[n_args];

    for (Integer i = 0; i < n_args; ++i)
        m_args[i] = in.t[i];
};

inline overload_impl::~overload_impl()
{
    delete[] m_args;
};

inline size_t overload_impl::eval_hash(const overload_info& in)
{
    size_t seed = in.f.hash_value();
    boost::hash_combine(seed,in.n_args);

    for (int i = 0; i < in.n_args; ++i)
    {
        boost::hash_combine(seed,size_t(-i));
        boost::hash_combine(seed,in.t[i].hash_value());
    };

    return seed;
};

inline bool overload_impl::equal(const overload_info& in) const
{
    if (m_name != in.f)         return false;
    if (n_args != in.n_args)    return false;

    for (int i = 0; i < n_args; ++i)
    {
        if (m_args[i] != in.t[i])
            return false;
    };

    return true;
};

//----------------------------------------------------------------------
//                  toverload_impl
//----------------------------------------------------------------------
inline toverload_impl::toverload_impl(const toverload_info& in)
    :m_name(in.f), n_args(in.n_args), m_overload(in.func), n_template(in.n_templ)
{
    m_hash      = eval_hash(in);
    m_args      = new Type[n_args];
    m_templates = new Type[n_template];

    for (Integer i = 0; i < n_args; ++i)
        m_args[i] = in.t[i];

    for (Integer i = 0; i < n_template; ++i)
        m_templates[i] = in.templ[i];
};

inline toverload_impl::~toverload_impl()
{
    delete[] m_args;
    delete[] m_templates;
};

inline size_t toverload_impl::eval_hash(const toverload_info& in)
{
    size_t seed = in.f.hash_value();
    boost::hash_combine(seed,in.n_args);

    for (int i = 0; i < in.n_args; ++i)
    {
        boost::hash_combine(seed,size_t(-i));
        boost::hash_combine(seed,in.t[i].hash_value());
    };

    for (int i = 0; i < in.n_templ; ++i)
    {
        boost::hash_combine(seed,size_t(6-i));
        boost::hash_combine(seed,in.templ[i].hash_value());
    };
    return seed;
};

inline bool toverload_impl::equal(const toverload_info& in) const
{
    if (m_name != in.f)             return false;
    if (n_args != in.n_args)        return false;
    if (n_template != in.n_templ)   return false;

    for (int i = 0; i < n_args; ++i)
    {
        if (m_args[i] != in.t[i])
            return false;
    };

    for (int i = 0; i < n_template; ++i)
    {
        if (m_templates[i] != in.templ[i])
            return false;
    };

    return true;
};

//----------------------------------------------------------------------
//                  convert_impl
//----------------------------------------------------------------------
inline convert_impl::convert_impl(convert_info in)
    :m_to(in.t1), m_from(in.t2), m_type(in.ctype), m_converter(in.f)
{
    m_hash = eval_hash(in);
};

inline size_t convert_impl::eval_hash(convert_info in)
{
    size_t seed = in.t1.hash_value();    
    boost::hash_combine(seed,in.t2.hash_value());
    boost::hash_combine(seed,in.ctype);

    return seed;
};

inline bool convert_impl::equal(convert_info in) const
{
    return m_to == in.t1 && m_from == in.t2 && m_type == in.ctype;
}

//----------------------------------------------------------------------
//                  assign_impl
//----------------------------------------------------------------------
inline assign_impl::assign_impl(assign_info in)
    :m_to(in.t1), m_from(in.t2), m_assigner(in.f)
{
    m_hash = eval_hash(in);
};

inline size_t assign_impl::eval_hash(assign_info in)
{
    size_t seed = in.t1.hash_value();
    boost::hash_combine(seed,1);
    boost::hash_combine(seed,in.t2.hash_value());

    return seed;
};

inline bool assign_impl::equal(assign_info in) const
{
    return m_to == in.t1 && m_from == in.t2;
}

};};};

#pragma warning(pop)