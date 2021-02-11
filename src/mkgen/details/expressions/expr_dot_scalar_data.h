/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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
#include "mkgen/details/expressions/expr_plus.h"
#include "matcl-simd/simd.h"
#include "matcl-scalar/lib_functions/func_binary.h"

//TODO
namespace matcl { namespace mkgen { namespace details
{

//----------------------------------------------------------------------------------
//                              forward declarations
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class List_1, class List_2>
struct expr_dot_arrays;

template<class List_1, class List_2, Integer Simd_Size>
struct make_dot_evaler;

template<class List_1, class List_2, Integer Block_Size>
struct make_blocks;

template<class Blocks, Integer N_Blocks>
struct dot_evaler;

template<Integer Block_Size, Integer Pos, class ... Elems>
struct split_list;

template<class List1, class List2, Integer Size>
struct dot_evaler_impl;

template<class List_1, class List_2>
struct expand_dot;

//----------------------------------------------------------------------------------
//                              expr_dot_scalar_data
//----------------------------------------------------------------------------------
//TODO: add checks
template<class List_1, class List_2>
struct expr_dot_scalar_data : public mkd::scalar_data<expr_dot_scalar_data<List_1, List_2>>
{
    using this_type     = expr_dot_scalar_data<List_1, List_2>;

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_dot_arrays<Step, Arr_List, List_1, List_2> :: type;

    //TODO
    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        using expanded  = typename expand_dot<List_1, List_2>::type;
        expanded::print<Subs_Context>(os,prior);
    };

    template<class Ret, class Local_Storage>
    inline_expr
    static Ret eval(const Local_Storage& ls)
    {
        using expanded  = typename expand_dot<List_1, List_2>::type;
        return expanded::eval<Ret>(ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using expanded  = typename expand_dot<List_1, List_2>::type;
        return expanded::accept(vis);
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        using value_type                = typename Local_Storage::value_type;
        static const Integer simd_size  = sizeof(Ret) / sizeof(value_type);
        using dot_evaler                = typename make_dot_evaler<List_1, List_2, simd_size>::type;

        dot_evaler::eval<Loop_Storage>(ret,offset,cont);
    };
};

//----------------------------------------------------------------------------------
//                              expand_dot
//----------------------------------------------------------------------------------
template<class List_1, class List_2>
struct expand_dot
{
    static_assert(md::dependent_false<List_1>::value,
                  "this type should not be instantiated");
};

template<>
struct expand_dot<list::list<>,list::list<>>
{
    using type      = zero_sd;
};

template<class Elem_1, class ... Elems_1, class Elem_2, class ... Elems_2>
struct expand_dot<list::list<Elem_1,Elems_1...>,list::list<Elem_2, Elems_2...>>
{
    using sum_2     = typename expand_dot<list::list<Elems_1...>,list::list<Elems_2...>>::type;
    using mult      = typename make_mult_root<Elem_1,Elem_2>::type;
    using type      = typename make_plus_root<mult,sum_2>::type;
};

//----------------------------------------------------------------------------------
//                              expr_dot_arrays
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class List_1, class List_2>
struct expr_dot_arrays
{
    static_assert(md::dependent_false<List_1>::value,
                  "this type should not be instantiated");
};

template<Integer Step, class Arr_List, class ... Elems>
struct expr_dot_arrays_list
{
    using type = Arr_List;
};

template<Integer Step, class Arr_List, class Elem, class ... Elems>
struct expr_dot_arrays_list<Step,Arr_List, Elem, Elems...>
{
    using arr_1     = typename Elem::template get_arrays<Step,Arr_List>;
    using type      = typename expr_dot_arrays_list<Step, arr_1, Elems...>::type;
};

template<Integer Step, class Arr_List, class ... Elems_1, class ... Elems_2>
struct expr_dot_arrays<Step, Arr_List, list::list<Elems_1...>, list::list<Elems_2...>>
{
    using arr       = typename expr_dot_arrays_list<Step, Arr_List, Elems_1...>::type;
    using type      = typename expr_dot_arrays_list<0, arr, Elems_2...>::type;
};

//----------------------------------------------------------------------------------
//                              make_dot_evaler
//----------------------------------------------------------------------------------
template<class List_1, class List_2, Integer Simd_Size>
struct make_dot_evaler
{
    static const Integer size       = list::size<List_2>::value;
    static const Integer block_size = (Simd_Size > 4) ? Simd_Size : 4;
    static const Integer n_blocks   = (size - 1) / block_size + 1;

    using blocks                    = typename make_blocks<List_1, List_2, block_size>::type;

    using type                      = dot_evaler<blocks,n_blocks>;
};

//----------------------------------------------------------------------------------
//                              make_blocks
//----------------------------------------------------------------------------------
template<class List_1, class List_2, Integer Block_Size>
struct make_blocks
{
    static_assert(md::dependent_false<List_1>::value, 
                "this type should not be instantiated");
};

template<class ... List_1, class ... List_2, Integer Block_Size>
struct make_blocks<list::list<List_1...>,list::list<List_2...>, Block_Size>
{
    using split_1   = split_list<Block_Size, 0, List_1...>;
    using split_2   = split_list<Block_Size, 0, List_2...>;

    using list_1    = typename split_1::front;
    using list_2    = typename split_2::front;

    using rem_1     = typename split_1::tail;
    using rem_2     = typename split_2::tail;

    using elem      = list::list<list_1,list_2>;
    using next      = typename make_blocks<rem_1,rem_2,Block_Size>::type;

    using type      = typename list::push_front<next, elem>::type;
};

template<Integer Block_Size>
struct make_blocks<list::list<>,list::list<>, Block_Size>
{
    using type      = list::list<>;
};

//----------------------------------------------------------------------------------
//                              split_list
//----------------------------------------------------------------------------------
template<Integer Block_Size, Integer Pos, class ... Elems>
struct split_list
{
    using front = list::list<>;
    using tail  = list::list<>;
};

template<Integer Block_Size, Integer Pos, class Elem, class ... Elems>
struct split_list<Block_Size,Pos,Elem,Elems...>
{
    using next  = split_list<Block_Size,Pos+1,Elems...>;
    using tail  = typename next::tail;
    using front = typename list::push_front<typename next::front,Elem>::type;
};
template<Integer Block_Size, class ... Elems>
struct split_list<Block_Size,Block_Size,Elems...>
{
    using tail  = list::list<Elems...>;
    using front = list::list<>;
};

//----------------------------------------------------------------------------------
//                              dot_evaler
//----------------------------------------------------------------------------------
//TODO
template<class Blocks, Integer N_Blocks>
struct dot_evaler
{
    static_assert(md::dependent_false<Blocks>::value, 
                "not implemented");
};

template<class Blocks>
struct dot_evaler<Blocks, 1>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        using block                 = typename list::elem_at_pos<Blocks,0>::type;
        using list_1                = typename list::elem_at_pos<block,0>::type;
        using list_2                = typename list::elem_at_pos<block,1>::type;
        using value_type            = typename Local_Storage::value_type;

        static const Integer size   = list::size<list_1>::value;

        ret = dot_evaler_impl<list_1,list_2,size>::eval<Loop_Storage,Ret>(offset,cont);

        //using expanded  = typename expand_dot<list_1, list_2>::type;
        //return expanded::eval_loop<Loop_Storage,Ret>(ret,offset,cont);
    };
};

//----------------------------------------------------------------------------------
//                              dot_evaler_impl
//----------------------------------------------------------------------------------
template<class List1, class List2, Integer Size>
struct dot_evaler_impl
{
    static_assert(md::dependent_false<List1>::value, 
                "not implemented");
};

template<class List1, class List2>
struct dot_evaler_impl<List1, List2, 1>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static Ret eval(Integer offset, const Local_Storage& cont)
    {
        using elem_11   = typename list::elem_at_pos<List1,0>::type;
        using elem_21   = typename list::elem_at_pos<List2,0>::type;

        Ret scal1;
        Ret val1;

        elem_21::eval_loop<Loop_Storage,Ret>(scal1,offset,cont);
        elem_11::eval_loop<Loop_Storage,Ret>(val1,offset,cont);

        Ret val = val1 * scal1 + val2*scal2;
        return val;
    };
};

template<class List1, class List2>
struct dot_evaler_impl<List1, List2, 2>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static Ret eval(Integer offset, const Local_Storage& cont)
    {
        using elem_11   = typename list::elem_at_pos<List1,0>::type;
        using elem_12   = typename list::elem_at_pos<List1,1>::type;

        using elem_21   = typename list::elem_at_pos<List2,0>::type;
        using elem_22   = typename list::elem_at_pos<List2,1>::type;

        Ret scal1, scal2;
        Ret val1, val2;

        elem_21::eval_loop<Loop_Storage,Ret>(scal1,offset,cont);
        elem_22::eval_loop<Loop_Storage,Ret>(scal2,offset,cont);

        elem_11::eval_loop<Loop_Storage,Ret>(val1,offset,cont);
        elem_12::eval_loop<Loop_Storage,Ret>(val2,offset,cont);

        /*
        Ret tmp1    = val1 * scal1;        
        tmp1        = fma_f(val2, scal2, tmp1);
        return tmp1;
        */

        Ret val = val1 * scal1 + val2*scal2;
        return val;
    };
};

template<class List1, class List2>
struct dot_evaler_impl<List1, List2, 3>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static Ret eval(Integer offset, const Local_Storage& cont)
    {
        using elem_11   = typename list::elem_at_pos<List1,0>::type;
        using elem_12   = typename list::elem_at_pos<List1,1>::type;
        using elem_13   = typename list::elem_at_pos<List1,2>::type;

        using elem_21   = typename list::elem_at_pos<List2,0>::type;
        using elem_22   = typename list::elem_at_pos<List2,1>::type;
        using elem_23   = typename list::elem_at_pos<List2,2>::type;

        Ret scal1, scal2, scal3;
        Ret val1, val2, val3;

        elem_21::eval_loop<Loop_Storage,Ret>(scal1,offset,cont);
        elem_22::eval_loop<Loop_Storage,Ret>(scal2,offset,cont);
        elem_23::eval_loop<Loop_Storage,Ret>(scal3,offset,cont);

        elem_11::eval_loop<Loop_Storage,Ret>(val1,offset,cont);
        elem_12::eval_loop<Loop_Storage,Ret>(val2,offset,cont);
        elem_13::eval_loop<Loop_Storage,Ret>(val3,offset,cont);

        Ret tmp1    = val1 * scal1;        
        tmp1        = fma_f(val2, scal2, tmp1);
        Ret tmp2    = val3 * scal3;
        Ret val     = tmp1 + tmp2;
        return val;

        //Ret val = (val1 * scal1 + val2*scal2) + (val3*scal3);
        //return val;
    };
};

template<class List1, class List2>
struct dot_evaler_impl<List1, List2, 4>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static Ret eval(Integer offset, const Local_Storage& cont)
    {
        using elem_11   = typename list::elem_at_pos<List1,0>::type;
        using elem_12   = typename list::elem_at_pos<List1,1>::type;
        using elem_13   = typename list::elem_at_pos<List1,2>::type;
        using elem_14   = typename list::elem_at_pos<List1,3>::type;

        using elem_21   = typename list::elem_at_pos<List2,0>::type;
        using elem_22   = typename list::elem_at_pos<List2,1>::type;
        using elem_23   = typename list::elem_at_pos<List2,2>::type;
        using elem_24   = typename list::elem_at_pos<List2,3>::type;

        Ret scal1, scal2, scal3, scal4;
        Ret val1, val2, val3, val4;

        elem_21::eval_loop<Loop_Storage,Ret>(scal1,offset,cont);
        elem_22::eval_loop<Loop_Storage,Ret>(scal2,offset,cont);
        elem_23::eval_loop<Loop_Storage,Ret>(scal3,offset,cont);
        elem_24::eval_loop<Loop_Storage,Ret>(scal4,offset,cont);

        elem_11::eval_loop<Loop_Storage,Ret>(val1,offset,cont);
        elem_12::eval_loop<Loop_Storage,Ret>(val2,offset,cont);
        elem_13::eval_loop<Loop_Storage,Ret>(val3,offset,cont);
        elem_14::eval_loop<Loop_Storage,Ret>(val4,offset,cont);

        Ret tmp1    = val1 * scal1;        
        tmp1        = fma_f(val2, scal2, tmp1);
        Ret tmp2    = val3 * scal3;
        tmp2        = fma_f(val4, scal4, tmp2);
        Ret val     = tmp1 + tmp2;
        return val;

        //Ret val     = (val1*scal1 + val2*scal2) + (val3*scal3 + val4*scal4);
        //return val;        
    };
};


}}}