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

#include "matcl-dynamic/object.h"
#include "matcl-core/details/mpl.h"

namespace matcl { namespace dynamic { namespace details
{

template<class T>
struct enable_const
{
    using type = const typename std::decay<T>::type&;
};

template<class T>
struct enable_const<const T&>
{
    using type = const typename std::decay<T>::type&;
};

template<class T>
struct enable_const<T&&>
{
    using type = const typename std::decay<T>::type&;
};

template<class T>
struct enable_const<const T&&>
{
    using type = const typename std::decay<T>::type&;
};

template<class T>
struct enable_const<T&>
{};

template<class T>
struct enable_nconst
{};

template<class T>
struct enable_nconst<T&>
{
    using type = T&;
};

template<class T>
struct enable_nconst<const T&>
{};

template <class FuncTraits>
struct data_constructor
{
    using return_type       = typename FuncTraits::return_type;

    static const 
    int n_inputs            = details::size<typename FuncTraits::input_type>::value;

    const object**          m_args;
    object*                 m_ret;

    template<int N, class T>
    typename enable_const<T>::type  get() const;

    template<int N, class T>
    typename enable_nconst<T>::type get() const;

    void      set_return(object&& arg)      {m_ret->reset(std::move(arg)); };
    void      set_return(const object& arg) {m_ret->reset(arg); };

    template<class Ty>
    void      set_return(const Ty& arg)     {m_ret->reset(object(arg)); };

    template<class Ty>
    void      set_return(Ty&& arg)          {m_ret->reset(object(std::move(arg))); };
};

template <class Base_constr>
struct data_constructor_no_ret : Base_constr
{
    void      set_return(object&&)      {};
    void      set_return(const object&) {};

    template<class Ty>
    void      set_return(const Ty&)     {};

    template<class Ty>
    void      set_return(Ty&)           {};
};

};};};