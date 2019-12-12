#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/utils/alignment.h"
#include "mkgen/TODO/evaler/dependency.h"

namespace matcl { namespace mkgen
{

//base tag for temp variables
struct temp_tag_base 
{
    static const bool is_continuous     = true;
    static const Integer step           = 1;
    using root_align_type               = align_full;

    template<Integer Row, Integer Col>
    using   get_offset                  = std::integral_constant<Integer,Row - 1>;
};

//----------------------------------------------------------------------------------
//                              mat_temporary
//----------------------------------------------------------------------------------
template<class Tag, Integer Mat_Rows, Integer Mat_Cols, bool Force>
struct mat_temp_array{};

template<class Tag, Integer Mat_Rows, Integer Mat_Cols, Integer Row, Integer Col, bool Force>
struct get_array_elem<mat_temp_array<Tag, Mat_Rows,Mat_Cols,Force>, Row, Col>
{
    using type = get_temporary<Tag,Mat_Rows,Mat_Cols, Row,Col>;
};

template<class Subs_Context, class Dep>
dps<> get_child_deps(Dep)
{
    return dps<>();
};

template<class Dep>
struct is_extern_dep
{
    static const bool value = false;
};
template<class Tag>
struct is_extern_dep<dep<Tag,0,dep_extern>>
{
    static const bool value = true;
};

template<class Subs_Context, class New_Deps, class Collected>
struct collect_deps
{
    static_assert(details::dependent_false<New_Deps>::value,
                  "this type should not be instantiated");
};
template<class Deps_All, class Collected>
struct collect_deps_temp
{
    static_assert(details::dependent_false<Deps_All>::value,
                  "this type should not be instantiated");
};

template<class Subs_Context, bool Is_1, class Dep, class Collected>
struct collect_deps2
{
    using type      = Collected;
};
template<class Subs_Context, class Dep, class Collected>
struct collect_deps2<Subs_Context, false, Dep, Collected>
{
    using deps1     = decltype(get_child_deps<Subs_Context>(Dep()));
    using collect2  = typename collect_deps<Subs_Context, deps1, Collected>::type;
    using type      = typename merge_deps<dps<Dep>, collect2>::type;
};

template<class Subs_Context, class Dep, class ... Deps, class Collected>
struct collect_deps<Subs_Context,dps<Dep, Deps...>, Collected>
{    
    using collected2        = typename collect_deps<Subs_Context,dps<Deps...>, Collected>::type;
    static const bool is_1  = list::is_member<Dep, collected2>::value;
    using type              = typename collect_deps2<Subs_Context, is_1, Dep, collected2>::type;
};

template<class Subs_Context, class Collected>
struct collect_deps<Subs_Context,dps<>, Collected>
{
    using type = Collected;
};

template<class Dep, class ... Deps, class ...Collected>
struct collect_deps_temp<dps<Dep,Deps...>, dps<Collected...>>
{
    using type =  typename collect_deps_temp<dps<Deps...>,dps<Collected...,Dep>>::type;
};
template<class Tag, Integer Size, class ... Deps, class ...Collected>
struct collect_deps_temp<dps<dep<Tag,Size,dep_extern>,Deps...>, dps<Collected...>>
{
    using type =  typename collect_deps_temp<dps<Deps...>,dps<Collected...>>::type;
};
template<class Tag, Integer Size, class ... Deps, class ...Collected>
struct collect_deps_temp<dps<dep<Tag,Size,dep_scalar>,Deps...>, dps<Collected...>>
{
    using type =  typename collect_deps_temp<dps<Deps...>,dps<Collected...>>::type;
};
template<class Collected>
struct collect_deps_temp<dps<>,Collected>
{
    using type = Collected;
};

template<class Mat, class Tag, bool Is_Temp, bool Force>
struct mat_temporary2
{
    using type = Mat;
};
template<Integer M1, Integer N1, class Deps, class Tag1, Integer Mat_Rows, Integer Mat_Cols, class Tag>
struct mat_temporary2<ct_matrix<M1,N1,mat_temp_array<Tag1, Mat_Rows, Mat_Cols,false>,Deps>,Tag,true,true>
{    
    using array = mat_temp_array<Tag1, Mat_Rows, Mat_Cols,true>;
    using type  = ct_matrix<M1, N1, array, Deps>;
};
template<class Mat, class Tag, bool Force>
struct mat_temporary2<Mat, Tag, false, Force>
{
    static const Integer M1     = Mat::rows;
    static const Integer N1     = Mat::cols;
    using deps1                 = typename Mat::dps_type; 
    using new_dep               = dep<Tag, M1 * N1, dep_temp>;
    using deps                  = dps<new_dep>;

    using array_type            = mat_temp_array<Tag, M1, N1, Force>;
    using type                  = ct_matrix<M1, N1, array_type,deps>;

    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    inline_initializer
    friend void tag_initializer(Tag, Local_Storage& ls, Data_Provider& dp, Temp_Storage* ts)
    {
        (void)ts;
        (void)dp;

        using subs_context      = typename Temp_Storage::subs_context;
        using val_type          = typename Temp_Storage::val_type;
        init_temporary_impl<Mat>(subs_context(), new_dep(), ls);
    };    
    template<class Visitor, class Temp_Storage>
    friend void tag_initializer_accept(Tag, Visitor& vis, Temp_Storage* ts)
    {
        using subs_context      = typename Temp_Storage::subs_context;
        init_temporary_impl_accept<Mat>(subs_context(), new_dep(), vis);
    };

    template<class Subs_Context>
    friend void tag_printer(Tag, std::ostream& os, int tabs)
    {
        print_matrix_elems<Mat>::eval_dep<Tag, Subs_Context>(os, tabs);
    };

    template<class Subs_Context>
    friend Mat get_stored_matrix(Tag)
    {
        return Mat();
    };

    template<class Subs_Context>
    friend typename deps1 get_child_deps(new_dep)
    {
        return deps1();
    };
};

template<class Mat, bool With_Forced>
struct is_temporary_mat
{
    static const bool value = false;
};
template<class Array_t, bool With_Forced>
struct is_temporary_mat_array
{
    static const bool value = false;
};
template<class Tag, Integer Mat_Rows, Integer Mat_Cols,bool Force>
struct is_temporary_mat_array<mat_temp_array<Tag, Mat_Rows, Mat_Cols,Force>, true>
{
    static const bool value = true;
};
template<class Tag, Integer Mat_Rows, Integer Mat_Cols,bool Force>
struct is_temporary_mat_array<mat_temp_array<Tag, Mat_Rows, Mat_Cols,Force>, false>
{
    static const bool value = (Force == false);
};

template<bool With_Forced, class ... Args>
struct is_temporary_virtual_array
{
    static_assert(details::dependent_false_var<Args...>::value, 
                "this type should not be instantiated");
};
template<class Item, bool With_Forced>
struct is_temporary_virtual_array_elem
{
    static_assert(details::dependent_false<Item>::value, 
                "this type should not be instantiated");
};
template<Integer M1, Integer N1, Integer M2, Integer N2, class Array2, class Colon_1, bool With_Forced>
struct is_temporary_virtual_array_elem<assign_item<M1, N1, M2, N2, Array2, Colon_1>, With_Forced >
{
    static const bool value = is_temporary_mat_array<Array2,With_Forced>::value; 
};
template<Integer M1, Integer N1, class scalar1, class Colon_1, bool With_Forced>
struct is_temporary_virtual_array_elem<assign_item_scalar<M1, N1, scalar1, Colon_1>, With_Forced>
{
    static const bool value = false; 
};

template<bool With_Forced>
struct is_temporary_virtual_array<With_Forced>
{
    static const bool value = true;
};
template<bool With_Forced, class Item, class... Args>
struct is_temporary_virtual_array<With_Forced, Item, Args...>
{
    static const bool value = is_temporary_virtual_array_elem<Item,With_Forced>::value
                            && is_temporary_virtual_array<With_Forced, Args...>::value;
};

template<class Tag, class ... Args, bool With_Forced>
struct is_temporary_mat_array<mkd::virtual_array<Tag, Args...>, With_Forced >
{
    //virtual array of temporary arrays is also temporary
    static const bool value = is_temporary_virtual_array<With_Forced, Args...>::value;
};

template<Integer M1, Integer N1, class Array1, class Deps1, bool With_Forced>
struct is_temporary_mat<ct_matrix<M1, N1, Array1, Deps1>,With_Forced>
{
    static const bool value = is_temporary_mat_array<Array1,With_Forced>::value;
};

template<class Mat, class Tag, bool Force>
struct mat_temporary
{
    static const bool is_temp   = is_temporary_mat<Mat,true>::value;
    using type                  = typename mat_temporary2<Mat,Tag,is_temp,Force>::type;
};

//----------------------------------------------------------------------------------
//                              get_temporary
//----------------------------------------------------------------------------------
template<class Tag, class Subs_Context>
struct get_temporary_elem
{
    using type = decltype(get_substitution(Subs_Context(),Tag()));
};

template<class Colon>
struct get_pos_colon2
{
    static_assert(details::dependent_false<Colon>::value,
                  "this type should not be instantiated");
};
template<>
struct get_pos_colon2<colon_all>
{
    inline_lev_1
    static void eval(Integer pos0, Integer& row, Integer& col, Integer mat_rows)
    {
        Integer pos     = pos0;
        col             = (pos-1) / mat_rows + 1;
        row             = pos - (col-1)*mat_rows;        
    };
};
template<Integer Start, Integer End>
struct get_pos_colon2<colon2<Start, End>>
{
    inline_lev_1
    static void eval(Integer pos0, Integer& row, Integer& col, Integer mat_rows)
    {
        Integer pos     = Start + (pos0 - 1);
        col             = pos / mat_rows + 1;
        row             = pos - (col - 1)*mat_rows;
    };
};
template<Integer Start, Integer Step, Integer End>
struct get_pos_colon2<colon3<Start, Step, End>>
{
    inline_lev_1
    static void eval(Integer pos0, Integer& row, Integer& col, Integer mat_rows)
    {
        Integer pos     = Start + Step * (pos0 - 1);
        col             = (pos-1) / mat_rows + 1;
        row             = pos - (col - 1)*mat_rows;
    };
};

template<Integer Mat_Rows>
struct get_mat_rows_val
{
    inline_lev_1
    static Integer eval(Integer mat_rows2)
    {
        (void)mat_rows2;
        return Mat_Rows;
    };
};
template<>
struct get_mat_rows_val<-1>
{
    inline_lev_1
    static Integer eval(Integer mat_rows2)
    {
        return mat_rows2;
    };
};

template<class Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols>
void print_subs(modif2<Tag,Colon,Mat_Rows,Mat_Cols>, Integer pos0, std::ostream& os, Integer mat_rows2)
{
    Integer row, col;
    get_pos_colon2<Colon>::eval(pos0, row, col, get_mat_rows_val<Mat_Rows>::eval(mat_rows2));

    Tag::print(os, details::prior_start);
    os << "[";
    os << row << "," << col;
    os << "]";
};

template<Integer Step, class Arr_List, class Elem>
struct get_temporary_array
{
    //substitution not taken into account
    using item      = details::array_item<Elem,Step, details::array_item_temp>;
    using type      = typename list::push_back<Arr_List,item> :: type;
};

template <class Tag, Integer Mat_Rows, Integer Mat_Cols, Integer Row, Integer Col>
struct get_temporary
{
    using tag                   = Tag;
    static const Integer rows   = Mat_Rows;
    static const Integer cols   = Mat_Cols;
    static const Integer pos    = (Col - 1) * Mat_Rows + Row;

    static const Integer        get_offset  = pos - 1;

    template<Integer Step, class Arr_List>
    using get_arrays    = typename get_temporary_array<Step, Arr_List, get_temporary> :: type;

    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        (void)prior;

        using subs      = typename get_temporary_elem<Tag,Subs_Context>::type;
        Integer pos2    = Row + (Col-1)*Mat_Rows;

        print_subs(subs(), pos2, os, Mat_Rows);
    };

    template<class Ret, class Local_Storage>
    inline_lev_1
    static Ret eval(const Local_Storage& ls)
    { 
        using dep = dep<Tag, Mat_Rows * Mat_Cols, dep_temp>;
        return ls.get_temp<dep,pos>();
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    { 
        vis.visit_load();
    };
    template<class Val, class Local_Storage>
    inline_lev_1
    static void assign(const Val& v, const Local_Storage& ls)
    { 
        using dep   = dep<Tag, Mat_Rows * Mat_Cols, dep_temp>;
        Val& elem   = const_cast<Val&>(ls.get_temp<dep,pos>());
        elem        = v;        
    };

    template<class Loop_Storage, class Ret, class Local_Storage>
    inline_lev_1
    static void eval_loop(Ret& ret, Integer offset, const Local_Storage& cont)
    { 
        return Loop_Storage::get_value<Ret, get_temporary>(cont, ret, offset);
    };

    template<class Visitor>
    static void accept_assign(Visitor& vis)
    {
        vis.visit_store();
    };
};

}}