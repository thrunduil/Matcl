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

#include "mkgen/details/matrix/matrix_printer.h"
#include "mkgen/TODO/expression/mat_assign.h"

namespace matcl { namespace mkgen { namespace details
{

//--------------------------------------------------------
//              forward declarations
//--------------------------------------------------------
template<class Subject, class Assignments>
struct print_computations_assignments;

template<Integer M, Integer N, class Array, class Deps, class Assign_Type>
struct print_comp_assing_1;

template<Integer M, Integer N, class Array, class Deps, class Colon,
                    class Mat, Integer Pos, Integer Size>
struct print_comp_assing_1_impl;

//--------------------------------------------------------
//              print_computations_assignments
//--------------------------------------------------------

template<class Subject, class Assignments>
struct print_computations_assignments
{
    static_assert(md::dependent_false<Subject>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, Mat_array Array, DPS Deps, class Assign_Type,
                class... Items>
struct print_computations_assignments<ct_matrix<M, N, Array,Deps>,
                    list::list<Assign_Type, Items...>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        using matrix = ct_matrix<M,N,Array,Deps>;

        print_comp_assing_1<M,N,Array,Deps, Assign_Type>
                    ::eval<Subs_Context>(os,nspaces);

        using printer = print_computations_assignments<matrix, list::list<Items...>>;
        printer::eval<Subs_Context>(os,nspaces);
    };
};

template<class Subject>
struct print_computations_assignments<Subject, list::list<>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int)
    {
        (void)os;
    };
};

//--------------------------------------------------------
//              print_comp_assing_1
//--------------------------------------------------------

template<Integer M, Integer N, class Array, class Deps, class Assign_Type>
struct print_comp_assing_1
{
    static_assert(md::dependent_false<Array>::value, 
                "this type should not be instantiated");
};

template<Integer M, Integer N, Mat_array Array, DPS Deps, Integer Pos, class Scalar>
struct print_comp_assing_1<M,N,Array,Deps, mkgen::assign_colon_scal<Pos,Scalar>>
{
    using matrix = ct_matrix<M,N,Array,Deps>;

    static const Integer col    = (Pos-1)/M + 1;
    static const Integer row    = Pos - (col-1) * M;

    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        using elem = typename Array :: template get_element<row,col>::type;
        print_whitespace(os,nspaces);
        elem::print<Subs_Context>(os, details::prior_start);

        os << " = ";

        Scalar::print<Subs_Context>(os, details::prior_start);
        os <<   "\n";
    };
};

template<Integer M, Integer N, class Array, class Deps, class Colon, class Mat>
struct print_comp_assing_1<M,N,Array,Deps, assign_colon<Colon,Mat>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        static const Integer size   = colon_func::size<Colon,M*N>::value;
        print_comp_assing_1_impl<M,N,Array,Deps,Colon,Mat,0,size>::eval<Subs_Context>(os,nspaces);
    };
};

//--------------------------------------------------------
//              print_comp_assing_1_impl
//--------------------------------------------------------
template<Integer M, Integer N, class Array, class Deps, class Colon,
                    class Mat, Integer Pos, Integer Size>
struct print_comp_assing_1_impl
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        static const Integer pos    = colon_func::index<Pos+1,Colon>::value;
        using mat_array             = typename Mat::array_type;
        using elem                  = typename mat_array :: template get_element<Pos+1,1>::type;

        print_comp_assing_1<M,N,Array,Deps,assign_colon_scal<pos,elem>>
                    ::eval<Subs_Context>(os,nspaces);

        print_comp_assing_1_impl<M,N,Array,Deps,Colon,Mat,Pos+1,Size>
                    ::eval<Subs_Context>(os,nspaces);
    };
};

template<Integer M, Integer N, class Array, class Deps, class Colon, 
                class Mat, Integer Size>
struct print_comp_assing_1_impl<M,N,Array,Deps,Colon,Mat,Size,Size>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        (void)os;
        (void)nspaces;
    };
};

}}}

namespace matcl { namespace mkgen
{

//--------------------------------------------------------
//              struct print_computations_elems
//--------------------------------------------------------

template<class Tag, class Computation>
struct print_computations_elems
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int nspaces)
    {
        mkd::print_whitespace(os,nspaces);
        os << "computation ";
        Tag::print(os,details::prior_start);
        os << ":" << "\n";

        using subject       = typename Computation::subject;
        using assignments   = typename Computation::assignments;

        mkd::print_computations_assignments<subject,assignments>
            ::eval<Subs_Context>(os, nspaces);
    };
};

}}