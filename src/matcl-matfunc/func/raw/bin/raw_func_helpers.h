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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/lib_functions/manip.h"

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

template<class str> struct str_to_code{};
template<> struct str_to_code<struct_dense>		{ static const int value = 0; };
template<> struct str_to_code<struct_sparse>	{ static const int value = 1; };
template<> struct str_to_code<struct_banded>	{ static const int value = 2; };

template<int val> struct code_to_str{};
template<> struct code_to_str<0>				{ using type = struct_dense ;  };
template<> struct code_to_str<1>				{ using type = struct_sparse;  };
template<> struct code_to_str<2>				{ using type = struct_banded;  };

template <class str_1,class str_2>
struct max_struct
{
    static const int code_1 = str_to_code<str_1>::value;
    static const int code_2 = str_to_code<str_2>::value;

    using type = typename code_to_str<(code_1 > code_2)? code_1 : code_2>::type;
};

template <class str_1,class str_2>
struct min_struct
{
    static const int code_1 = str_to_code<str_1>::value;
    static const int code_2 = str_to_code<str_2>::value;

    using type = typename code_to_str<(code_1 < code_2)? code_1 : code_2>::type;
};

//ZZ: func(0,0), ZN: func(0,x), NZ: func(x,0) for any x!= 0
//value: true - always 0, false- nonzero
template<class str_1,class str_2,bool ZZ,bool ZN,bool NZ>
struct stor_struct_cons
{};
/*
ZZ			Z	Z	Z	Z	N	N	N	N
ZN			Z	Z	N	N	Z	Z	N	N
NZ			Z	N	Z	N	Z	N	Z	N
                                        
GM	GM		GM	GM	GM	GM	GM	GM	GM	GM
GM	SM		SM	GM	SM	GM	SM	GM	SM	GM
GM	B		B	GM	B	GM	SM	GM	SM	GM
                                        
SM	GM		SM	SM	GM	GM	SM	SM	GM	GM
SM	SM		SM	SM	SM	SM	GM	GM	GM	GM
SM	B		B	SM	B	SM	GM	GM	GM	GM
                                        
B	GM		B	B	GM	GM	SM	SM	GM	GM
B	SM		B	B	SM	SM	GM	GM	GM	GM
B	B		B	B	B	B	GM	GM	GM	GM										
*/

template<class str_1,class str_2,bool ZZ,bool ZN,bool NZ>
struct stor_struct_cons_ms
{
    using type  = typename matcl::details::select_if
                <
                    ZN == true,
                    str_1,
                    struct_dense
                >::type;
};
/*
ZZ			Z	Z	Z	Z	N	N	N	N
ZN			Z	Z	N	N	Z	Z	N	N
NZ			Z	N	Z	N	Z	N	Z	N
                                        
GM	S		GM	GM	GM	GM	GM	GM	GM	GM
SM	S		SM	SM	GM	GM	SM	SM	GM	GM
B	S		B	B	GM	GM	B	B	GM	GM
*/

template<class str_1,class str_2,bool ZZ,bool ZN,bool NZ>
struct stor_struct_cons_sm
{
    using type  = typename matcl::details::select_if
                <
                    NZ == true,
                    str_2,
                    struct_dense
                >::type;
};
/*
ZZ			Z	Z	Z	Z	N	N	N	N
ZN			Z	Z	N	N	Z	Z	N	N
NZ			Z	N	Z	N	Z	N	Z	N
                                        
S	GM		GM	GM	GM	GM	GM	GM	GM	GM
S	SM		SM	GM	SM	GM	SM	GM	SM	GM
S	B		B	GM	B	GM	B	GM	B	GM
*/

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,true,true,true>
{
    using type = typename max_struct<str_1,str_2>::type;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,true,true,false>
{
    using type = str_1;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,true,false,true>
{
    using type = str_2;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,true,false,false>
{
    using type = typename min_struct<str_1,str_2>::type;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,false,true,true>
{
    using type  = typename matcl::details::select_if
                <
                    std::is_same<str_1,struct_dense>::value 
                        && !std::is_same<str_2,struct_dense>::value
                    ||std::is_same<str_2,struct_dense>::value 
                        && !std::is_same<str_1,struct_dense>::value,
                    struct_sparse,
                    struct_dense
                >::type;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,false,true,false>
{
    using type  = typename matcl::details::select_if
                <
                    !std::is_same<str_1,struct_dense>::value 
                        && std::is_same<str_2,struct_dense>::value,
                    struct_sparse,
                    struct_dense
                >::type;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,false,false,true>
{
    using type  = typename matcl::details::select_if
                <
                    std::is_same<str_1,struct_dense>::value 
                        && !std::is_same<str_2,struct_dense>::value,
                    struct_sparse,
                    struct_dense
                >::type;
};

template<class str_1,class str_2>
struct stor_struct_cons<str_1,str_2,false,false,false>
{
    using type = struct_dense;
};

template<class T1, class T2, bool Real_ret>
struct unify_types_real : md::unify_types<T1,T2>
{};

template<class T1, class T2>
struct unify_types_real<T1, T2, true>
{
     using type0    = typename md::unify_types<T1,T2>::type;
     using type     = typename md::real_type<type0>::type;
};

template<class M1, class M2, bool ZZ, bool ZN, bool NZ,bool int_ret, bool real_ret>
struct ret_type_sm_cons
{
    using val_1         = M1;
    using val_2         = typename M2::value_type;
    using str_1         = struct_dense;
    using str_2         = typename M2::struct_type;
    using ret_struct    = typename stor_struct_cons_sm<str_1,str_2,ZZ,ZN,NZ>::type;

    static const bool is_obj    = md::is_object<val_1>::value 
                                || md::is_object<val_2>::value;

    using int_or_obj    = typename matcl::details::select_if
                        <
                            is_obj,
                            Object,
                            Integer
                        >::type;

    using ret_value     = typename matcl::details::select_if
                        <
                            int_ret,
                            int_or_obj,
                            typename unify_types_real<val_1,val_2,real_ret>::type
                        >::type;

    using type          = Matrix<ret_value,ret_struct>;
};

template<class M1, class M2, bool ZZ, bool ZN, bool NZ,bool int_ret, bool real_ret>
struct ret_type_ms_cons
{
    using val_1         = typename M1::value_type;
    using val_2         = M2;
    using str_1         = typename M1::struct_type;
    using str_2         = struct_dense;
    using ret_struct    = typename stor_struct_cons_ms<str_1,str_2,ZZ,ZN,NZ>::type;
    
    static const bool is_obj    = md::is_object<val_1>::value 
                                || md::is_object<val_2>::value;

    using int_or_obj    = typename matcl::details::select_if
                        <
                            is_obj,
                            Object,
                            Integer
                        >::type;

    using ret_value     = typename matcl::details::select_if
                        <
                            int_ret,
                            int_or_obj,
                            typename unify_types_real<val_1,val_2,real_ret>::type
                        >::type;

    using type          = Matrix<ret_value,ret_struct>;
};

template<class M1, class M2, bool ZZ, bool ZN, bool NZ,bool int_ret, bool real_ret>
struct ret_type_constructor
{
    using val_1         = typename M1::value_type;
    using val_2         = typename M2::value_type;
    using str_1         = typename M1::struct_type;
    using str_2         = typename M2::struct_type;
    using ret_struct    = typename stor_struct_cons<str_1,str_2,ZZ,ZN,NZ>::type;

    static const bool is_obj    = md::is_object<val_1>::value 
                                || md::is_object<val_2>::value;

    using int_or_obj    = typename matcl::details::select_if
                        <
                            is_obj,
                            Object,
                            Integer
                        >::type;

    using ret_value     = typename matcl::details::select_if
                        <
                            int_ret,
                            int_or_obj,
                            typename unify_types_real<val_1, val_2,real_ret>::type
                        >::type;

    using type          = Matrix<ret_value,ret_struct>;
};

template<class mat_type>
struct correct_int_ret
{
    using value_type_0  = typename mat_type::value_type;
    using struct_type   = typename mat_type::struct_type;

    using value_type    = typename matcl::details::select_if
                        <
                            std::is_same<value_type_0,Integer>::value,
                            Real,
                            value_type_0
                        >::type;

    using type          = Matrix<value_type,struct_type>;
};

template<class Mat, bool Is_obj>
struct correct_object_type
{
    using type = Mat;
};

template<class Mat>
struct correct_object_type<Mat, true>
{
    using struct_type   = typename Mat::struct_type;
    using type          = raw::Matrix<Object,struct_type>;
};

}}}
