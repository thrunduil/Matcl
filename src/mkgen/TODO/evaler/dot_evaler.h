#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"

#include "matcl-scalar/lib_functions/func_binary.h"

namespace matcl { namespace mkgen
{

template<Integer Block_Size, Integer Pos, class ... Elems>
struct split_list
{
    using front = list<>;
    using tail  = list<>;
};
template<Integer Block_Size, Integer Pos, class Elem, class ... Elems>
struct split_list<Block_Size,Pos,Elem,Elems...>
{
    using next  = split_list<Block_Size,Pos+1,Elems...>;
    using tail  = typename next::tail;
    using front = typename push_front<typename next::front,Elem>::type;
};
template<Integer Block_Size, class ... Elems>
struct split_list<Block_Size,Block_Size,Elems...>
{
    using tail  = list<Elems...>;
    using front = list<>;
};

template<class List_1, class List_2, Integer Block_Size>
struct make_blocks
{
    static_assert(details::dependent_false<List_1>::value, 
                "this type should not be instantiated");
};
template<class ... List_1, class ... List_2, Integer Block_Size>
struct make_blocks<list<List_1...>,list<List_2...>, Block_Size>
{
    using split_1   = split_list<Block_Size, 0, List_1...>;
    using split_2   = split_list<Block_Size, 0, List_2...>;

    using list_1    = typename split_1::front;
    using list_2    = typename split_2::front;

    using rem_1     = typename split_1::tail;
    using rem_2     = typename split_2::tail;

    using elem      = list<list_1,list_2>;
    using next      = typename make_blocks<rem_1,rem_2,Block_Size>::type;

    using type      = typename push_front<next, elem>::type;
};
template<Integer Block_Size>
struct make_blocks<list<>,list<>, Block_Size>
{
    using type      = list<>;
};

template<class List1, class List2, Integer Size>
struct dot_evaler_impl
{
    static_assert(details::dependent_false<List1>::value, 
                "not implemented");
};

template<class List1, class List2>
struct dot_evaler_impl<List1, List2, 1>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static Ret eval(Integer offset, const Local_Storage& cont)
    {
        using elem_11   = typename get_elem_at_pos<List1,0>::type;
        using elem_21   = typename get_elem_at_pos<List2,0>::type;

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
        using elem_11   = typename get_elem_at_pos<List1,0>::type;
        using elem_12   = typename get_elem_at_pos<List1,1>::type;

        using elem_21   = typename get_elem_at_pos<List2,0>::type;
        using elem_22   = typename get_elem_at_pos<List2,1>::type;

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
        using elem_11   = typename get_elem_at_pos<List1,0>::type;
        using elem_12   = typename get_elem_at_pos<List1,1>::type;
        using elem_13   = typename get_elem_at_pos<List1,2>::type;

        using elem_21   = typename get_elem_at_pos<List2,0>::type;
        using elem_22   = typename get_elem_at_pos<List2,1>::type;
        using elem_23   = typename get_elem_at_pos<List2,2>::type;

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
        using elem_11   = typename get_elem_at_pos<List1,0>::type;
        using elem_12   = typename get_elem_at_pos<List1,1>::type;
        using elem_13   = typename get_elem_at_pos<List1,2>::type;
        using elem_14   = typename get_elem_at_pos<List1,3>::type;

        using elem_21   = typename get_elem_at_pos<List2,0>::type;
        using elem_22   = typename get_elem_at_pos<List2,1>::type;
        using elem_23   = typename get_elem_at_pos<List2,2>::type;
        using elem_24   = typename get_elem_at_pos<List2,3>::type;

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

template<class Blocks, Integer N_Blocks>
struct dot_evaler
{
    static_assert(details::dependent_false<Blocks>::value, 
                "not implemented");
};
template<class Blocks>
struct dot_evaler<Blocks,1>
{
    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_expr
    static void eval(Ret& ret, Integer offset, const Local_Storage& cont)
    {
        using block                 = typename get_elem_at_pos<Blocks,0>::type;
        using list_1                = typename get_elem_at_pos<block,0>::type;
        using list_2                = typename get_elem_at_pos<block,1>::type;
        using value_type            = typename Local_Storage::value_type;

        static const Integer size   = list_size<list_1>::value;

        ret = dot_evaler_impl<list_1,list_2,size>::eval<Loop_Storage,Ret>(offset,cont);

        //using expanded  = typename expand_dot<list_1, list_2>::type;
        //return expanded::eval_loop<Loop_Storage,Ret>(ret,offset,cont);
    };
};

template<class List_1, class List_2, Integer Simd_Size>
struct make_dot_evaler
{
    static const Integer size       = list_size<List_2>::value;
    static const Integer block_size = (Simd_Size > 4) ? Simd_Size : 4;
    static const Integer n_blocks   = (size - 1) / block_size + 1;

    using blocks                    = typename make_blocks<List_1, List_2, block_size>::type;

    using type                      = dot_evaler<blocks,n_blocks>;
};

}}