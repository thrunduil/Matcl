#pragma once

#include "mkgen/mkgen_fwd.h"
#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/utils/simd_utils.h"
#include "mkgen/matrix/dependency.h"
#include "mkgen/TODO/evaler/loop_promotion.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              is_continuous
//----------------------------------------------------------------------------------
template<class Code_Gen, class Modif>
struct is_continuous
{};

template<class Tag>
struct is_tag_continuous
{
    static const bool value = Tag::is_continuous;
};

template<class Code_Gen, class Tag, class Colon, Integer R, Integer C>
struct is_continuous<Code_Gen, modif2<Tag,Colon, R, C>>
{
    static const bool allow_mone_step  = Code_Gen::simd_allow_negative_step;
    static const bool value = ((colon_func::step<Colon>::value == 1 
                                    || allow_mone_step && colon_func::step<Colon>::value == -1) 
                               && is_tag_continuous<Tag>::value)
                            || colon_func::size<Colon,R*C>::value == 1;
};

//----------------------------------------------------------------------------------
//                              simd_enable
//----------------------------------------------------------------------------------
//now ony operations of type func(x,y) can be vectorized
template<class Subs_Context, class Stored_Matrix>
struct simd_enable
{
    static_assert(md::dependent_false<Stored_Matrix>::value,
                  "this type should not be instantiated");
};

template<class Subs_Context, class Array_T>
struct enable_vectorization_array
{
    static_assert(md::dependent_false<Array_T>::value,
                  "this type should not be instantiated");
};

template<class Subs_Context, Integer M, Integer N, Mat_array Array_T, class Deps>
struct simd_enable<Subs_Context,ct_matrix<M,N,Array_T,Deps>>
{
    using codegen           = typename Subs_Context::code_gen;
    static const bool value = codegen::simd_enable
                            && enable_vectorization_array<Subs_Context,Array_T>::value;
};

template<class Subs_Context, class Tag, class... Assign_List>
struct enable_vectorization_array<Subs_Context,mkd::virtual_array<Tag,Assign_List...>> 
{
    static const bool value = false;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_assign_array<M, N, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array2>::value;
};
template<class Subs_Context, Integer M,Integer N,class Array1,class Colon, 
    Integer M2, Integer N2, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_assign_array_colon<M, N, Array1, Colon, M2, N2, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_scal_assign_array<M, N, Array1, Array2>> 
{
    static const bool value = false;
};

template<class Subs_Context, class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer Step2>
struct enable_vectorization_array<Subs_Context, mkd::sub_array_2<Array_t, Offset1, Offset2, Step1, Step2>> 
{
    static const bool value = false;
};

template<class Subs_Context, class Ret_Tag>
struct enable_vectorization_array<Subs_Context, mkd::empty_array<Ret_Tag>> 
{
    static const bool value = false;
};

template<class Subs_Context>
struct enable_vectorization_array<Subs_Context, call_array_type> 
{
    static const bool value = false;
};

template<class Subs_Context, class Tag, Integer M, Integer N, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_bfunc_array<Tag, M, N, Array1, Array2>>
{
    static const bool value = false;
};

template<class Subs_Context, class Tag, Integer M, Integer N, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_scal_bfunc_array<Tag, M, N, Array1, Array2>>
{
    static const bool value = false;
};

template<class Subs_Context, class Tag, Integer M, Integer N, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::scal_mat_bfunc_array<Tag, M, N, Array1, Array2>> 
{
    static const bool value = false;
};

template<class Subs_Context, class Array_t, Integer Offset, Integer Step>
struct enable_vectorization_array<Subs_Context, mkd::sub_array_1<Array_t, Offset, Step>> 
{
    using code_gen  = typename Subs_Context::code_gen;
    static const bool allow_mone_step  = code_gen::simd_allow_negative_step;

    static const bool value = (Step == 1 || allow_mone_step && Step == -1) 
                            && enable_vectorization_array<Subs_Context,Array_t>::value;
};

template<class Subs_Context, Integer K, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_mult_array<K, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_plus_array<M, N, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value 
                                && enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mult_array<Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value 
                                && enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::div_array<Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value 
                                && enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_minus_array<M, N, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value 
                                && enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mult_rows_array<Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value 
                                && enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, class Array1, class Scalar2>
struct enable_vectorization_array<Subs_Context, mkd::mat_scal_mult_array<Array1, Scalar2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, class Array1, class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mult_cols_array<Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, class Array1, class Scal2>
struct enable_vectorization_array<Subs_Context, mkd::div_array_mat_scal<Array1, Scal2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, class Array2, class Scal1>
struct enable_vectorization_array<Subs_Context, mkd::div_array_scal_mat<Array2, Scal1>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array2>::value;
};

template<class Subs_Context, Integer M, Integer N, class Array>
struct enable_vectorization_array<Subs_Context, mkd::mat_uminus_array<M, N, Array>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array>::value;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_scal_plus_array<M, N, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::mat_scal_minus_array<M, N, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, Integer M,Integer N,class Array1,class Array2>
struct enable_vectorization_array<Subs_Context, mkd::scal_mat_minus_array<M, N, Array1, Array2>> 
{
    static const bool value = enable_vectorization_array<Subs_Context,Array1>::value;
};

template<class Subs_Context, class Tag, Integer Rows, Integer Cols>
struct enable_vectorization_array<Subs_Context, mkd::temp_output_array<Tag, Rows, Cols>> 
{
    using code_gen          = typename Subs_Context::code_gen;
    using ret_subs          = decltype(get_substitution(Subs_Context(), Tag()));
    static const bool value = is_continuous<code_gen,ret_subs>::value;
};

template<class Subs_Context, class Tag>
struct enable_vectorization_array<Subs_Context, mkd::gen_array<Tag>>
{
    static const bool value = Tag::is_continuous;
};

template<class Subs_Context, class Tag, class Array, class Deps>
struct enable_vectorization_array<Subs_Context, details::scalar_ufunc_array<Tag, Array, Deps>> 
{
    static const bool value = true;
};

template<class Subs_Context, class Tag, Integer Mat_Rows, Integer Mat_Cols, bool Force>
struct enable_vectorization_array<Subs_Context, mkd::mat_temp_array<Tag, Mat_Rows, Mat_Cols, Force>> 
{
    using code_gen          = typename Subs_Context::code_gen;
    using ret_subs          = decltype(get_substitution(Subs_Context(), Tag()));
    static const bool value = is_continuous<code_gen,ret_subs>::value;
};

template<class Subs_Context, class Tag>
struct enable_vectorization_array<Subs_Context, mkd::output_array<Tag>> 
{
    using code_gen          = typename Subs_Context::code_gen;
    using ret_subs          = decltype(get_substitution(Subs_Context(), Tag()));
    static const bool value = is_continuous<code_gen,ret_subs>::value;
};

template<class Subs_Context, class Tag, Integer M, Integer N, class Array>
struct enable_vectorization_array<Subs_Context, mkd::mat_ufunc_array<Tag, M, N, Array>> 
{
    static const bool value = Tag::is_continuous;
};

template<class Subs_Context, Tag_matrix_const_data Tag, class Val>
struct enable_vectorization_array<Subs_Context, mkd::matrix_array_const_value<Tag, Val>> 
{
    static const bool value = false;
};

template<class Subs_Context, Tag_matrix_data Tag, class Val>
struct enable_vectorization_array<Subs_Context, mkd::matrix_array_value<Tag, Val>> 
{
    static const bool value = false;
};

//----------------------------------------------------------------------------------
//                              loop_context
//----------------------------------------------------------------------------------

template<class Array_Tags_List>
struct loop_context
{
    using array_tags_list   = Array_Tags_List;
};

template<class Elem>
struct make_loop_context
{
    using array_collector   = typename Elem::template get_arrays<1,list::list<>>;
    using type              = loop_context<array_collector>;
};

//----------------------------------------------------------------------------------
//                              value_getter, value_setter
//----------------------------------------------------------------------------------
template<class Val, class Ret, bool Is_Aligned, Integer Step, Integer Offset>
struct value_getter
{
    static_assert(md::dependent_false<Val>::value,
                  "this type should not be instantiated");
};

template<class Val, class Ret, bool Is_Aligned, Integer Offset>
struct value_getter<Val,Ret,Is_Aligned,0,Offset>
{
    static_assert(sizeof(Val) != sizeof(Ret), "invalid getter, check const modifiers");

    static const Integer vec_size = Ret::vector_size;
    
    inline_lev_1
    static void eval(Ret& ret, const Val* arr, Integer offset)
    {
        (void)offset;
        ret = Ret(arr[Offset]);
    };
};

template<class Val, class Ret, bool Is_Aligned, Integer Offset>
struct value_getter<Val,Ret,Is_Aligned,1,Offset>
{
    static_assert(sizeof(Val) != sizeof(Ret), "invalid getter, check const modifiers");

    static const Integer vec_size = Ret::vector_size;
    
    inline_lev_1
    static void eval(Ret& ret, const Val* arr, Integer offset)
    {
        ret = Ret::load(arr + Offset + offset*vec_size, std::integral_constant<bool,Is_Aligned>());
    };
};

template<class Val, class Ret, bool Is_Aligned, Integer Offset>
struct value_getter<Val,Ret,Is_Aligned,-1, Offset>
{
    static_assert(sizeof(Val) != sizeof(Ret), "invalid getter, check const modifiers");

    static const Integer vec_size = Ret::vector_size;
    
    inline_lev_1
    static void eval(Ret& ret, const Val* arr, Integer offset)
    {
        using aligned_type  = std::integral_constant<bool,Is_Aligned>;

        ret = details::load_reverse<Ret, aligned_type>
                    ::eval(arr + Offset - offset*vec_size);
        //ret = Ret::load_reverse(arr + Offset - offset*vec_size, aligned_type());
    };
};

template<class Val, bool Is_Aligned, Integer Offset>
struct value_getter<Val,Val,Is_Aligned,1, Offset>
{
    inline_lev_1
    static void eval(Val& ret, const Val* arr, Integer off)
    {
        ret = arr[Offset + off];
    };
};

template<class Val, bool Is_Aligned, Integer Offset>
struct value_getter<Val,Val,Is_Aligned,0, Offset>
{
    inline_lev_1
    static void eval(Val& ret, const Val* arr, Integer off)
    {
        (void)off;
        ret = arr[Offset];
    };
};

template<class Val, bool Is_Aligned, Integer Offset>
struct value_getter<Val,Val,Is_Aligned,-1, Offset>
{
    inline_lev_1
    static void eval(Val& ret, const Val* arr, Integer off)
    {
        ret = arr[Offset-off];
    };
};

template<class Val, class Ret, bool Is_Aligned, Integer Step, Integer Offset>
struct value_setter
{
    static_assert(md::dependent_false<Val>::value,
                  "this type should not be instantiated");
};

template<class Val, class Ret, bool Is_Aligned, Integer Offset>
struct value_setter<Val,Ret,Is_Aligned,1,Offset>
{
    static_assert(sizeof(Val) != sizeof(Ret), "invalid getter, check const modifiers");

    using val_nc                    = typename std::remove_cv<Val>::type;
    static const Integer vec_size   = Ret::vector_size;
    
    inline_lev_1
    static void eval(val_nc* arr, Integer offset, const Ret& val)
    {
        val.store(arr + Offset + offset*vec_size,std::integral_constant<bool,Is_Aligned>());
    };
};

template<class Val, class Ret, bool Is_Aligned, Integer Offset>
struct value_setter<Val,Ret,Is_Aligned,-1, Offset>
{
    static_assert(sizeof(Val) != sizeof(Ret), "invalid getter, check const modifiers");

    static const Integer vec_size = Ret::vector_size;
    
    inline_lev_1
    static void eval(Val* arr, Integer offset, const Ret& val)
    {
        using aligned_type  = std::integral_constant<bool,Is_Aligned>;

        //val.store_reverse(arr + Offset - offset*vec_size,aligned_type());
        details::store_reverse<Ret, aligned_type>
                ::eval(arr + Offset - offset*vec_size, val);
        
    };
};

template<class Val, bool Is_Aligned, class Ret, Integer Step, Integer Offset>
struct value_setter<Val,Ret,Is_Aligned,Step,Offset>
{
    static_assert(sizeof(Val) != sizeof(Ret), "invalid getter, check const modifiers");

    static const Integer vec_size = Ret::vector_size;
    
    inline_lev_1
    static void eval(Val* arr, Integer offset, const Ret& val)
    {
        val.scatter<Step>(arr + Offset + (offset*vec_size) * Step);
    };
};

template<class Val, bool Is_Aligned, Integer Offset>
struct value_setter<Val,Val,Is_Aligned,1,Offset>
{
    inline_lev_1
    static void eval(Val* arr, Integer offset, const Val& val)
    {
        arr[Offset + offset] = val;
    };
};
template<class Val, bool Is_Aligned, Integer Offset>
struct value_setter<Val,Val,Is_Aligned,-1, Offset>
{
    inline_lev_1
    static void eval(Val* arr, Integer offset, const Val& val)
    {
        arr[Offset-offset] = val;
    };
};
template<class Val, bool Is_Aligned, Integer Step, Integer Offset>
struct value_setter<Val,Val,Is_Aligned,Step, Offset>
{
    inline_lev_1
    static void eval(Val* arr, Integer offset, const Val& val)
    {
        arr[Offset+Step*offset] = val;
    };
};

//----------------------------------------------------------------------------------
//                             loop_context_data
//----------------------------------------------------------------------------------
template<class Val, class Data_Provider, class Subs_Context, class Array_Tags_List>
struct make_loop_context_info
{
    static_assert(md::dependent_false<Array_Tags_List>::value,
                  "this type should not be instantiated");
};

template<class Val, class Aligned_Root_Type, Integer Step, Integer Offset, class Dep>
struct loop_context_data
{
    using dep   = Dep;

    using align_off             = typename get_offset_alignment<Val,Offset,Step>::type;
    using align_type            = typename link_alignment<Aligned_Root_Type, align_off>::type;
    static const Integer step   = Step;
    static const Integer offset = Offset;
};
template<class Val, class Elem>
struct loop_context_data_scalar
{
    static const Integer step   = 0;
    using align_type            = align_full;
    using dep                   = loop_context_data_scalar;
};

template<class Val, class Aligned, Integer Step, Integer Offset>
struct loop_data_ret
{
    template<class Ret>
    inline_lev_1
    static void set_ret(Val* arr, Ret& val, int offset)
    {
        static const bool is_align = is_aligned<Aligned, sizeof(Ret)/sizeof(Val)>::value;
        return value_setter<Val,Ret,is_align,Step,Offset>::eval(arr, offset, val); 
    };
};

template<class Val, Integer Step, Integer Offset>
struct loop_data_ret_unaligned
{
    template<class Ret>
    inline_lev_1
    static void set_ret(Val* arr, Ret& val, int offset)
    {
        return value_setter<Val,Ret,false,Step,Offset>::eval(arr, offset, val); 
    };
};

template<class Data_Info>
struct check_has_negative_step
{
    static_assert(md::dependent_false<Data_Info>::value,
                  "this type should not be instantiated");
};
template<class Data_Info>
struct get_align_type
{
    static_assert(md::dependent_false<Data_Info>::value,
                  "this type should not be instantiated");
};

template<class Item, class ... Loop_Data>
struct check_has_negative_step<list::list<Item, Loop_Data...>>
{
    static const bool value = (Item::step != 1 && Item::step != 0) 
                            || check_has_negative_step<list::list<Loop_Data...>>::value;
};
template<>
struct check_has_negative_step<list::list<>>
{
    static const bool value = false;
};

template<class Item, class ... Loop_Data>
struct get_align_type<list::list<Item, Loop_Data...>>
{
    using align_tail    = typename get_align_type<list::list<Loop_Data...>>::type;
    using align_elem    = typename Item::align_type;
    using type          = typename link_alignment<align_tail,align_elem>::type;
};

template<>
struct get_align_type<list::list<>>
{
    using type          = align_full;
};

template<class Data_Provider, class Colon, Integer Offset, class Tag, class Ret_Tag>
struct get_ret_offset
{
    static const Integer pos1   = colon_func::index<Offset + 1, Colon>::value;
    static const Integer pos2   = Ret_Tag::get_offset(1,1);
    static const Integer value  = pos2+pos1-1;
};
template<class Data_Provider, class Colon, Integer Offset, class Tag>
struct get_ret_offset<Data_Provider,Colon,Offset,Tag,Tag>
{
    using offset_info           = decltype(get_temp_offset_info(Tag()));
    static const Integer value  = Offset + offset_info::offset;
};

template<class Temp_Elem, class Modif, class Tag, class Ret_Tag>
struct make_dep_subs
{
    using type = dep<Ret_Tag,0,dep_extern>;
};
template<class Temp_Elem, class Modif, class Tag>
struct make_dep_subs<Temp_Elem,Modif,Tag,Tag>
{
    using type = dep<Tag,Temp_Elem::rows * Temp_Elem::cols, dep_temp>;
};

template<class Elem>
struct make_dep
{
    using type  = dep<Elem,0,dep_extern>;
};
template<class Tag, Integer Row, Integer Col>
struct make_dep<mkd::element<Tag,Row,Col>>
{
    using type = dep<Tag,0,dep_extern>;
};
template<class Tag, Integer Rows, Integer Cols, Integer Row, Integer Col>
struct make_dep<get_temporary<Tag,Rows,Cols,Row,Col>>
{
    using type = dep<Tag,Rows*Cols,dep_temp>;
};

template<class Val, class Data_Provider, class Subs_Context, class Elem>
struct make_loop_context_data
{
    static_assert(md::dependent_false<Val>::value,
                  "this type should not be instantiated");
};

template<class Val, class Data_Provider, class Subs_Context, class Elem, Integer Step>
struct make_loop_context_data<Val,Data_Provider,Subs_Context,
    details::array_item<Elem,Step, details::array_item_extern>>
{
    static const Integer offset = Elem::get_offset;
    using aligned_r             = typename Elem::root_align_type;
 
    using dep           = typename make_dep<Elem>::type;
    using type          = loop_context_data<Val,aligned_r,Step,offset,dep>;
};

template<class Val, class Data_Provider, class Subs_Context, class Elem, Integer Step>
struct make_loop_context_data<Val,Data_Provider,Subs_Context, details::array_item<Elem,Step, details::array_item_temp>>
{
    using tag           = typename Elem::tag;
    using ret_subs      = decltype(get_substitution(Subs_Context(), tag()));
    using colon         = typename ret_subs::colon;
    using ret_tag       = typename ret_subs::tag;

    static const Integer offset1    = Elem::get_offset;
    static const Integer ret_step   = colon_func::step<colon>::value;
    static const Integer offset     = get_ret_offset<Data_Provider,colon,offset1, tag, ret_tag>::value;
    using align_r                   = typename ret_tag::root_align_type;

    using dep           = typename make_dep_subs<Elem, ret_subs, tag, ret_tag>::type;
    using type          = loop_context_data<Val,align_r,Step * ret_step, offset, dep>;
};
template<class Val, class Data_Provider, class Subs_Context, class Elem, Integer Step>
struct make_loop_context_data<Val,Data_Provider,Subs_Context, details::array_item<Elem,Step,details::array_item_scalar>>
{    
    using type          = loop_context_data_scalar<Val,Elem>;
};

template<class Info, class List_Arr, class List_Scal>
struct collect_tags
{};
template<class Val, class Aligned_Root, Integer Step, Integer Offset, class Dep, 
        class ... Elems, class ... Arr, class List_Scal>
struct collect_tags<list::list<loop_context_data<Val,Aligned_Root,Step,Offset,Dep>, Elems...>, list::list<Arr...>, List_Scal>
{
    using collector     = collect_tags<list::list<Elems...>,list::list<Dep,Arr...>, List_Scal>;
    using array_list    = typename collector::array_list;
    using scalar_list   = typename collector::scalar_list;
};
template<class Val, class Elem,  class ... Elems, class List_Arr, class ... Scal>
struct collect_tags<list::list<loop_context_data_scalar<Val,Elem>, Elems...>, List_Arr, list::list<Scal...> >
{
    using collector     = collect_tags<list::list<Elems...>, List_Arr, list::list<Elem,Scal...> >;
    using array_list    = typename collector::array_list;
    using scalar_list   = typename collector::scalar_list;
};

template<class List_Arr, class List_Scal>
struct collect_tags<list::list<>, List_Arr, List_Scal>
{
    using array_list    = List_Arr;
    using scalar_list   = List_Scal;
};

template<class Val, class Data_Provider, class Subs_Context, class ...Elems>
struct make_loop_context_info<Val, Data_Provider, Subs_Context, list::list<Elems...>>
{
    using info          = list::list<typename make_loop_context_data<Val, Data_Provider, 
                                Subs_Context, Elems>::type...>;
    using collector     = collect_tags<info, list::list<>, list::list<>>;
    using tags_array    = typename collector::array_list;
    using tags_scal     = typename collector::scalar_list;
    using unique_arrays = typename list::unique_list<tags_array>::type;
    using unique_scal   = typename list::unique_list<tags_scal>::type;
};

template<class Array_List>
struct remove_step
{
    static_assert(md::dependent_false<Array_List>::value,
                  "this type should not be instantiated");
};
template<class ... Elems>
struct remove_step<list::list<Elems...>>
{
    using type = list::list<typename Elems::elem...>;
};

template<class Val, class Ret_Align, Integer Ret_Step, Integer Ret_Offset>
struct make_ret_storage
{
    using type = loop_data_ret_unaligned<Val,Ret_Step,Ret_Offset>;
};

template<class Val, class Ret_Align, Integer Ret_Offset>
struct make_ret_storage<Val,Ret_Align,1, Ret_Offset>
{
    using type = loop_data_ret<Val,Ret_Align,1,Ret_Offset>;
};
template<class Val, class Ret_Align, Integer Ret_Offset>
struct make_ret_storage<Val,Ret_Align,-1, Ret_Offset>
{
    using type = loop_data_ret<Val,Ret_Align,-1,Ret_Offset>;
};

//----------------------------------------------------------------------------------
//                             loop_storage
//----------------------------------------------------------------------------------
template<class Loop_Context, class Val, class Local_Storage, class Ret_Dep, class Ret_Align, 
        Integer Ret_Step, Integer Ret_Offset>
struct loop_storage
{
    static_assert(std::is_same<Val, typename std::remove_cv<Val>::type>::value == true, "invalid Val");

    using value_type            = Val;
    using array_tags_list       = typename Loop_Context::array_tags_list;
    using elems_list            = typename remove_step<array_tags_list>::type;
    using data_provider         = typename Local_Storage::data_provider;
    using subs_context          = typename Local_Storage::subs_context;
    using context_info          = make_loop_context_info<Val,data_provider, subs_context, array_tags_list>;
    using data_info             = typename context_info::info;
    using unique_arrays         = typename context_info::unique_arrays;
    using unique_scalar         = typename context_info::unique_scal;

    using ret_storage           = typename make_ret_storage<Val, Ret_Align, Ret_Step, Ret_Offset>::type;

    static const Integer size   = list::size<array_tags_list>::value;
    static const bool has_negative_step = check_has_negative_step<data_info>::value;
    using aligned_type          = typename get_align_type<data_info>::type;

    template<class Ret, class Elem, class Local_Storage>
    inline_lev_1
    static void get_value(Local_Storage& ls, Ret& ret, int off)
    {
        static const Integer pos    = list::elem_pos<elems_list,Elem>::value;
        using info                  = typename list::elem_at_pos<data_info,pos>::type;
        using dep                   = typename info::dep;
        static const Integer pos_d  = list::elem_pos<unique_arrays,dep>::value;

        using align_type            = typename info::align_type;
        static const Integer step   = info::step;
        static const Integer offset = info::offset;
        static const bool is_align  = is_aligned<align_type, sizeof(Ret)/sizeof(Val)>::value;

        const Val* arr              = ls.get_array<dep>();//m_data.elem<pos_d>()
        return value_getter<Val,Ret,is_align,step,offset>::eval(ret, arr, off);
    };
    
    template<class Ret, class Elem, class Local_Storage>
    inline_lev_1
    static void get_value_array(Local_Storage& ls, Ret& ret, int off, const Val* arr)
    {
        static const Integer pos    = list::elem_pos<elems_list,Elem>::value;
        using info                  = typename list::elem_at_pos<data_info,pos>::type;
        using dep                   = typename info::dep;
        static const Integer pos_d  = list::elem_pos<unique_arrays,dep>::value;

        using align_type            = typename info::align_type;
        static const Integer step   = info::step;
        static const Integer offset = info::offset;
        static const bool is_align  = is_aligned<align_type, sizeof(Ret)/sizeof(Val)>::value;

        return value_getter<Val,Ret,is_align,step,offset>::eval(ret, arr, off);
    };
    
    template<class Ret, class Local_Storage>
    inline_lev_1
    static void set_ret(Local_Storage& ls, Ret& val, int offset)
    {
        static_assert(std::is_same<Ret, typename std::remove_cv<Ret>::type>::value == true, "invalid Ret");

        const Val* arr  = ls.get_array<Ret_Dep>();
        ret_storage::set_ret<Ret>(const_cast<Val*>(arr), val, offset);
    };
};

//----------------------------------------------------------------------------------
//                              loop_unrolling
//----------------------------------------------------------------------------------
template<Integer Rows>
struct do_unrolling
{
    static const bool value     = (Rows <= 8);
};

template<Integer N, Integer Pos, Integer Start, class Val, class Elem, class Loop_Storage>
struct loop_unrolling_impl
{
    template<class Local_Storage>
    inline_force
    static void eval(Local_Storage& ls)
    {
        {
            Val ret;
            Elem::eval_loop<Loop_Storage, Val>(ret, Start + Pos, ls);
            Loop_Storage::set_ret<Val>(ls, ret, Start + Pos);
        };

        loop_unrolling_impl<N, Pos+1, Start, Val, Elem, Loop_Storage>
            ::eval(ls);
    };
};

template<Integer N, Integer Start, class Val, class Elem, class Loop_Storage>
struct loop_unrolling_impl<N,N,Start,Val,Elem,Loop_Storage>
{
    template<class Local_Storage>
    static void eval(Local_Storage& ls)
    {
        (void)ls;
    };
};

template<Integer Rows_Start, Integer Rows_End, class Val, class Elem, class Loop_Storage, 
    bool Do_unrolling = do_unrolling<Rows_End-Rows_Start+1>::value>
struct loop_unrolling
{
    template<class Local_Storage>
    inline_loop
    static void eval(Local_Storage& ls)
    {   
        static const Integer vec_size   = vector_size<Val>::value;
        static const Integer step       = vec_size;
        static const Integer start      = (Rows_Start - 1)/vec_size;

        for (Integer i = Rows_Start, j = start; i <= Rows_End ; i += step, ++j)
        {
            Val ret;
            Elem::eval_loop<Loop_Storage, Val>(ret, j, ls);
            ls.set_ret<Val>(ret,j);
        };

    };
};
template<Integer Row_Start, Integer Row_End, class Val, class Elem, class Loop_Storage>
struct loop_unrolling<Row_Start, Row_End,Val,Elem,Loop_Storage,true>    
{
    template<class Local_Storage>
    inline_lev_1
    static void eval(Local_Storage& ls)
    {   
        static const Integer vec_size   = details::simd_vector_size<Val>::value;
        static const Integer step       = vec_size;
        static const Integer start      = (Row_Start - 1)/vec_size;
        static const Integer N          = (Row_End - Row_Start + 1) / step;

        loop_unrolling_impl<N, 0, start, Val, Elem, Loop_Storage>
            ::eval(ls);
    };
};

template<class Val, class Simd_Half, Integer N_Rem, bool Select_Half = (N_Rem <= Simd_Half::vector_size)>
struct make_simd_rem2
{
    using type = Simd_Half;
};

template<class Val, class Simd_Half, Integer N_Rem>
struct make_simd_rem2<Val,Simd_Half,N_Rem,false>
{
    //using type = simd<Val, type_avx_limit<N_Rem>>;
    using type = Simd_Half;
};

template<class Val, class Simd, class Simd_Half, Integer Rows>
struct make_simd_rem
{
    static const Integer n_rem  = Rows % Simd::vector_size;
    using type                  = typename make_simd_rem2<Val, Simd_Half, n_rem>::type;
};

//----------------------------------------------------------------------------------
//                              loop_evaler
//----------------------------------------------------------------------------------
template<class Loop_Context, Integer Rows, class Elem, class Ret_Dep, class Ret_Align, 
    Integer Ret_Step, Integer Ret_Offset>
struct loop_evaler
{
    template<class Val, class Local_Storage>
    inline_loop
    static void eval(Local_Storage& st)
    {
        using subs_context  = typename Local_Storage::subs_context;
        using code_gen      = typename subs_context::code_gen;
        using simd_t        = typename code_gen::simd_type;
        static const int Bits
                            = code_gen::simd_bits;
        using simd_type     = matcl::simd::simd<Val, Bits, simd_t>;
        using simd_half_t   = typename simd_type::simd_half;
        using simd2_type    = typename make_simd_rem<Val,simd_type,simd_half_t,Rows>::type;        

        using storage       = loop_storage<Loop_Context, Val, Local_Storage, Ret_Dep, 
                                Ret_Align, Ret_Step, Ret_Offset>;

        static const bool has_negative_step = storage::has_negative_step;
        using aligned_type                  = typename storage::aligned_type;
        
        static const bool allow_full    = true
                                        && (has_negative_step == false || code_gen::simd_allow_negative_step )
                                        && (is_align_full<aligned_type>::value || code_gen::simd_allow_unaligned );

        static const bool allow_half    = code_gen::simd_half_allow
                                        && (has_negative_step == false || code_gen::simd_half_allow_negative_step )
                                        && (is_align_half<aligned_type>::value || code_gen::simd_half_allow_unaligned );

        static const Integer Rows_1_s   = 1;
        static const Integer Rows_1_e   = allow_full ? (Rows / simd_type::vector_size) * simd_type::vector_size : 0;
        static const Integer Rows_2_s   = Rows_1_e + 1;
        static const Integer Rows_2_e   = allow_half? (Rows / simd2_type::vector_size) 
                                            * simd2_type::vector_size : Rows_1_e;
        static const Integer Rows_3_s   = Rows_2_e + 1;
        static const Integer Rows_3_e   = Rows;

        loop_unrolling<Rows_1_s, Rows_1_e, simd_type, Elem, storage>::eval(st);
        loop_unrolling<Rows_2_s, Rows_2_e, simd2_type, Elem, storage>::eval(st);
        loop_unrolling<Rows_3_s, Rows_3_e, Val, Elem, storage>::eval(st);
    };
};

}}

//#pragma warning(pop)