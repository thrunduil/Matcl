/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-dynamic/details/utils.h"
#include "matcl-core/details/mpl.h"

#include <type_traits>

namespace matcl { namespace details
{

template<bool val>
struct bool_type                                                { static const bool value = val;  };
                                                                                                                                    
template<bool cond,class T1,class T2> struct lazy_select_if     { using type = typename T1::type; };
template<class T1,class T2> struct lazy_select_if<false,T1,T2>  { using type = typename T2::type; };
                                                                  
template<bool cond,class T> struct lazy_enable_if               { using type = typename T::type;  };
template<class T>           struct lazy_enable_if<false,T>      {};


template<class T>   struct lazy_type                            { using type = T;                 };
template<int val>   struct integer_type                         { static const int value = val;   };

//disable type deduction                                          
template<class T>   struct hide_type                            { using type = T;  };

template<class T> struct hide_object_type
{                
    using type = typename select_if< is_object<T>::value, Object, T >::type;
};

//enabler for rvalue reference types
template<class T>   struct enable_rvalue                        {};
template<class T>   struct enable_rvalue<T&&>                   { using type = T&&;};

template<class T>   struct enable_lvalue                        {};
template<class T>   struct enable_lvalue<T&>                    { using type = T&; };

};};