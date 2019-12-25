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

#include "mkgen/matrix/matrix.h"
#include "mkgen/details/matrix/virtual_matrix.h"
#include "matcl-core/details/mpl.h"

//TODO

namespace matcl { namespace mkgen { namespace details
{

//--------------------------------------------------------
//              forward declarations
//--------------------------------------------------------
template<class Dep_List>
struct print_matrix_deps;

template<class Dep>
struct print_dependency;

template<class Matrix_List, Integer Mat_Rows>
struct print_matrix_elems_expand;

template<Integer M, Integer N, class Colon, class Array, Integer Row, 
            Integer Col, Integer Mat_Rows>
struct print_matrix_elems_impl;

//--------------------------------------------------------
//              utils
//--------------------------------------------------------

// print n spaces
inline void print_whitespace(std::ostream& os, int n)
{
    if (n > 0)
        os << std::string(n, ' ');
};

//--------------------------------------------------------
//              print priorities
//--------------------------------------------------------
static const int prior_start    = 0;
static const int prior_assign   = 10;
static const int prior_plus     = 20;
static const int prior_mult     = 30;
static const int prior_pow      = 40;

//--------------------------------------------------------
//              print_matrix_list
//--------------------------------------------------------
template<class Matrix_List>
struct print_matrix_list
{
    static_assert(md::dependent_false<Matrix_List>::value, 
                "this type should not be instantiated");
};

template<class Mat, class ... Matrix>
struct print_matrix_list<list::list<Mat, Matrix...>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        print_matrix_elems<Mat>::eval<Subs_Context>(os, details::prior_start, nspaces);

        if constexpr(sizeof...(Matrix) > 0)
            os << "\n";

        using list_type = list::list<Matrix...>;
        using printer   = print_matrix_list<list_type>;
        printer::eval<Subs_Context>(os, nspaces);
    };
};

template<>
struct print_matrix_list<list::list<>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        (void)os;
        (void)nspaces;
    };
};

//--------------------------------------------------------
//              print_matrix_deps
//--------------------------------------------------------
template<class Deps_List>
struct print_matrix_deps
{
    static_assert(md::dependent_false<Deps_List>::value, 
                "this type should not be instantiated");
};

template<>
struct print_matrix_deps<dps<>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {};
};

template<class Dep1, class... Deps>
struct print_matrix_deps<dps<Dep1, Deps...>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        print_dependency<Dep1>::eval<Subs_Context>(os, nspaces);        

        if constexpr(sizeof...(Deps) > 0)
        {
            using printer   = print_matrix_deps<dps<Deps...>>;
            printer::eval<Subs_Context>(os, nspaces);
        };
    };
};

//--------------------------------------------------------
//              print_dependency
//--------------------------------------------------------
template<class Dep>
struct print_dependency
{
    static_assert(md::dependent_false<Dep>::value, 
                "this type should not be instantiated");
};

template<class Tag, Integer Size>
struct print_dependency<dep<Tag, Size, dep_return>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {     
        (void)os;
        (void)nspaces;
    };
};

template<class Tag, Integer Size, class Type>
struct print_dependency< dep<Tag, Size, Type>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {        
        bool need_init = tag_need_initialization(Subs_Context(), Tag());

        if (need_init == true)
        {
            os << "\n";

            print_whitespace(os, nspaces);
            os << "and" << "\n";

            os << "\n";

            tag_printer<Subs_Context>(Tag(), os, nspaces);
        };
    };
};

//--------------------------------------------------------
//              print_matrix_elems_expand
//--------------------------------------------------------
template<class Matrix_List, Integer Mat_Rows>
struct print_matrix_elems_expand
{
    static_assert(md::dependent_false<Matrix_List>::value, 
                "this type should not be instantiated");
};

template<Integer Mat_Rows>
struct print_matrix_elems_expand<list::list<>, Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int nspaces)
    {
        (void)os;
        (void)prior;
        (void)nspaces;
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int nspaces)
    {
        (void)os;
        (void)nspaces;
    };
};

template<class Mat_Colon, class ... Mats, Integer Mat_Rows>
struct print_matrix_elems_expand<list::list<Mat_Colon, Mats...>, Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int nspaces)
    {
        using colon_assign          = typename list::elem_at_pos<Mat_Colon,0>::type;
        using matrix_type           = typename list::elem_at_pos<Mat_Colon,1>::type;
        static const Integer M      = matrix_type::rows;
        static const Integer N      = matrix_type::cols;
        using code_gen              = typename Subs_Context::code_gen;

        static const bool allow_mstep = code_gen::simd_allow_negative_step;
        static const bool is_cont   = ( N == 1 && (M >= 4 || M >= 2 && code_gen::simd_half_allow) )
                                    && simd_enable<Subs_Context,matrix_type>::value;

        if constexpr(is_cont == true)
        {
            print_whitespace(os, nspaces);
            os << "loop vectorized \n";
        };    

        using array_t               = typename matrix_type::array_type;

        print_matrix_elems_impl<M,N,colon_assign,array_t,0,0,Mat_Rows>
                ::eval<Subs_Context>(os,prior,nspaces);

        if constexpr(sizeof...(Mats) > 0)
            os << "\n";

        using printer   = print_matrix_elems_expand<list::list<Mats...>,Mat_Rows>;
        printer::eval<Subs_Context>(os,prior,nspaces);
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int nspaces)
    {
        using colon_assign      = typename list::elem_at_pos<Mat_Colon,0>::type;
        using matrix_type       = typename list::elem_at_pos<Mat_Colon,1>::type;

        using code_gen          = typename Subs_Context::code_gen;
        using subs              = decltype(get_substitution(Subs_Context(),Tag()));
        using colon             = typename subs::colon;
        static const Integer M  = matrix_type::rows;
        static const Integer N  = matrix_type::cols;
        using array_t           = typename matrix_type::array_type;

        static const bool allow_mstep = code_gen::simd_allow_negative_step;
        static const bool is_cont
                            = ( N == 1 && (M >= 4 || M >= 2 && code_gen::simd_half_allow) )
                            && simd_enable<Subs_Context,matrix_type>::value
                            ;

        if constexpr(is_cont == true)
        {
            print_whitespace(os, nspaces);
            os << "loop vectorized \n";
        };

        print_matrix_elems_impl<M,N,colon_assign,array_t,0,0,Mat_Rows>
                    ::eval_dep<Tag,Subs_Context>(os,nspaces);

        if constexpr(sizeof...(Mats) > 0)
            os << "\n";

        using printer = print_matrix_elems_expand<list::list<Mats...>,Mat_Rows>;
        printer::eval_dep<Tag,Subs_Context>(os,nspaces);
    };
};

//--------------------------------------------------------
//              print_matrix_elems_impl
//--------------------------------------------------------
template<Integer M, Integer N, class Colon, class Array, Integer Row, 
            Integer Col, Integer Mat_Rows>
struct print_matrix_elems_impl
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int nspaces)
    {
        (void)prior;

        print_whitespace(os, nspaces);

        using elem = typename Array::template get_element<Row + 1, Col + 1>::type;
        elem::print<Subs_Context>(os,details::prior_start);

        os << "\n";

        print_matrix_elems_impl<M,N,Colon,Array,Row+1,Col,Mat_Rows>
            ::eval<Subs_Context>(os, details::prior_start, nspaces);
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int nspaces)
    {
        using subs  = decltype(get_substitution(Subs_Context(),Tag()));

        static const Integer pos0   = Row+1 + Col*M;
        static const Integer pos    = colon_func::index<pos0,Colon>::value;

        print_whitespace(os, nspaces);
        print_subs(subs(), pos, os, Mat_Rows);
        os << " = ";

        using elem = typename Array::template get_element<Row + 1, Col + 1>::type;

        elem::print<Subs_Context>(os,details::prior_start);

        os << "\n";

        print_matrix_elems_impl<M,N,Colon,Array,Row+1,Col,Mat_Rows>
                    ::eval_dep<Tag,Subs_Context>(os, nspaces);
    };
};

template<Integer M, Integer N, class Colon, class Ret_Tag, Integer Row, 
                Integer Col, Integer Mat_Rows>
struct print_matrix_elems_impl<M, N, Colon, mkd::empty_array<Ret_Tag>, 
                                    Row, Col, Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int nspaces)
    {
        (void)prior;

        print_whitespace(os, nspaces);
        Ret_Tag::print(os,details::prior_start);
        os <<"\n";
    };
};

template<Integer M, Integer N, class Colon, class Array, Integer Col, 
                Integer Mat_Rows>
struct print_matrix_elems_impl<M, N, Colon, Array, M, Col, Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int nspaces)
    {
        print_matrix_elems_impl<M,N,Colon,Array,0,Col+1,Mat_Rows>
                    ::eval<Subs_Context>(os, prior, nspaces);
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int nspaces)
    {
        print_matrix_elems_impl<M,N,Colon,Array,0,Col+1,Mat_Rows>
            ::eval_dep<Tag,Subs_Context>(os, nspaces);
    };
};

template<Integer M, Integer N, class Colon, class Array, Integer Mat_Rows>
struct print_matrix_elems_impl<M, N, Colon, Array, 0, N, Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int, int)
    {
        (void)os;
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int)
    {
        (void)os;
    };
};

}}}

namespace matcl { namespace mkgen
{

//--------------------------------------------------------
//              print matrix
//--------------------------------------------------------

template<Integer M, Integer N, Mat_array Array_t, DPS Deps>
template<class Subs_Context>
void ct_matrix<M, N, Array_t, Deps>::print(std::ostream& os, int nspaces)
{
    using matrix_type           = ct_matrix<M, N, Array_t, Deps>;
    using expand_vm             = typename mkd::expand_virtual_matrix<matrix_type>::type;

    mkd::print_whitespace(os, nspaces);
    os << "matrix " << M << " x " << N << " with elements:" << "\n";

    mkd::print_whitespace(os, nspaces);
    os << "[" << "\n";

    mkd::print_matrix_list<expand_vm>::eval<Subs_Context>(os, nspaces);

    mkd::print_whitespace(os, nspaces);
    os << "]" << "\n";

    using deps_temp  = typename Subs_Context::deps_temp;

    if constexpr(std::is_same<deps_temp, empty_deps>::value == false)
        mkd::print_matrix_deps<deps_temp>::eval<Subs_Context>(os, nspaces);
};

//--------------------------------------------------------
//              print_matrix_elems
//--------------------------------------------------------
template<class Mat>
struct print_matrix_elems
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int nspaces)
    {
        using expand_vm = typename mkd::expand_virtual_matrix2<Mat>::type;
        Integer size    = list::size<expand_vm>::value;
        
        if (size > 1)
        {
            mkd::print_whitespace(os, nspaces);
            os << "[" << "\n";
        };

        mkd::print_matrix_elems_expand<expand_vm,Mat::rows>::eval<Subs_Context>(os,prior,nspaces);        

        if (size > 1)
        {
            mkd::print_whitespace(os, nspaces);
            os << "]" << "\n";
        };
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int nspaces)
    {
        using expand_vm = typename mkd::expand_virtual_matrix2<Mat>::type;
        Integer size    = list::size<expand_vm>::value;
        
        if (size > 1)
        {
            mkd::print_whitespace(os, nspaces);
            os << "[" << "\n";
        };

        mkd::print_matrix_elems_expand<expand_vm,Mat::rows>
                    ::eval_dep<Tag,Subs_Context>(os,nspaces);

        if (size > 1)
        {
            mkd::print_whitespace(os, nspaces);
            os << "]" << "\n";
        };
    };
};

}}