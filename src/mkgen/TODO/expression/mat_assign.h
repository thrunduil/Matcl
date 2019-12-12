#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/evaler/dependency.h"
#include "mkgen/details/matrix/colon_func.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              expr_assign
//----------------------------------------------------------------------------------
template<Integer Step, class Arr_List, class T1, class T2>
struct expr_assign_arrays
{
    using arr_1 = typename T1::template get_arrays<Step, Arr_List>;
    using type  = typename T2::template get_arrays<Step, arr_1>;
};

template<class T1, class T2>
struct expr_assign
{
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        if (prior >= details::prior_assign)
            os << "(";

        T1::print<Subs_Context>(os, details::prior_assign);

        os << "=";

        T2::print<Subs_Context>(os, details::prior_start);

        if (prior >= details::prior_assign)
            os << ")";
    };
    template<class Ret, class Local_Storage>
    inline_lev_1
    static void eval(const Local_Storage& ls)
    {
        Ret val = T2::eval<Ret>(ls);
        return T1::assign<Ret>(val, ls);
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        T2::accept<Visitor>(vis);
        T1::accept_assign<Visitor>(vis);
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_assign_arrays<Step, Arr_List, T1, T2> :: type;
};

//----------------------------------------------------------------------------------
//                              mat_assign_array
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array1, class Array2>
struct mat_scal_assign_array{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct get_array_elem<mkd::mat_assign_array<M,N,Array1,Array2>, Row, Col>
{
    static_assert(Row <= M && Col <= N, "invalid elem");

    using elem_1    = typename get_array_elem<Array1, Row, Col>::type;
    using elem_2    = typename get_array_elem<Array2, Row, Col>::type;
    using new_item  = expr_assign<elem_1,elem_2>;
    using type      = new_item;
};

template<Integer M, Integer N, class Array1, class Array2, class Deps2,
        Integer Row, Integer Col>
struct get_array_elem<mat_scal_assign_array<M,N,Array1, ct_scalar<Array2,Deps2>>, Row, Col>
{
    static_assert(Row <= M && Col <= N, "invalid elem");

    using elem_1    = typename get_array_elem<Array1, Row, Col>::type;
    using elem_2    = typename ct_scalar<Array2,Deps2>;
    using new_item  = expr_assign<elem_1,elem_2>;
    using type      = new_item;
};

//----------------------------------------------------------------------------------
//                              mat_assign_array_colon
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Array1, class Colon, Integer M2, Integer N2, class Array2, 
        Integer Row, Integer Col>
struct get_array_elem<mkd::mat_assign_array_colon<M,N,Array1,Colon, M2, N2, Array2>, Row, Col>
{
    static_assert(Row <= M2 && Col <= N2, "invalid elem");

    static const Integer Pos    = (Col-1)*M2 + Row;
    static const Integer Pos_1  = colon_func::index<Pos,Colon>::value;
    static const Integer Col_1  = (Pos_1 - 1) / M + 1;
    static const Integer Row_1  = (Pos_1 - 1) % M + 1;

    using elem_1    = typename get_array_elem<Array1, Row_1, Col_1>::type;
    using elem_2    = typename get_array_elem<Array2, Row, Col>::type;
    using new_item  = expr_assign<elem_1,elem_2>;
    using type      = new_item;
};

//----------------------------------------------------------------------------------
//                              mat_assign
//----------------------------------------------------------------------------------
template<Integer M1_M2, Integer N1_N2, class Array1, class Deps1, class Array2, class Deps2>
struct mat_assign<ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_assign_array<M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1, class Array2, class Deps2>
struct mat_assign<ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2,Deps2>>
{
    using array_type    = mat_scal_assign_array<M1, N1, Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, Integer M2, Integer N2, class Array1, class Deps1, class Array2, class Deps2>
struct mat_assign<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid matrix operation, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              mat_virtual_assign_1
//----------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Colon_1>
struct is_assignment_valid_scalar{};
template<Integer M1, Integer N1, class Colon_1, Integer M2, Integer N2>
struct is_assignment_valid{};

template<Integer M1, Integer N1>
struct is_assignment_valid_scalar<M1,N1,colon_all>
{
    static const bool value = true;
};
template<Integer M1, Integer N1, Integer M2, Integer N2>
struct is_assignment_valid<M1,N1,colon_all, M2, N2>
{
    static const bool value = (M2 == 1 || N2 == 1) && M2 * N2 == M1 * N1;
};

template<Integer M1, Integer N1, Integer Pos>
struct is_assignment_valid_scalar<M1,N1,colon<Pos>>
{
    static const bool value = (Pos > 0 && Pos <= M1*N1);
};
template<Integer M1, Integer N1, Integer Pos, Integer M2, Integer N2>
struct is_assignment_valid<M1,N1,colon<Pos>, M2, N2>
{
    static const bool value = (M2 == 1 && N2 == 1) && Pos > 0 && Pos <= M1 * N1;
};

template<Integer M1, Integer N1, Integer Start, Integer End>
struct is_assignment_valid_scalar<M1,N1,colon2<Start,End>>
{
    static const bool value = (Start > 0 && End >= Start && End <= M1*N1);
};
template<Integer M1, Integer N1, Integer Start, Integer End, Integer M2, Integer N2>
struct is_assignment_valid<M1,N1,colon2<Start,End>, M2, N2>
{
    static const Integer size     = End - Start + 1;
    static const bool value = (Start > 0 && End >= Start && End <= M1*N1) && (M2 == 1 || N2 == 1) 
                                && M2 * N2 == size;
};

template<Integer M1, Integer N1, Integer Start, Integer Step, Integer End>
struct is_assignment_valid_scalar<M1,N1,colon3<Start,Step,End>>
{
    static const Integer size = M1 * N1;
    static const bool value = (Start > 0 && Start <= size && End > 0 && End <= size
                                && Step != 0);
};
template<Integer M1, Integer N1, Integer Start, Integer Step, Integer End, Integer M2, Integer N2>
struct is_assignment_valid<M1,N1,colon3<Start,Step,End>,M2,N2>
{
    static const Integer size       = M1 * N1;
    static const Integer colon_size = (End - Start) / Step + 1;
    static const bool value         = (Start > 0 && Start <= size && End > 0 && End <= size
                                        && Step != 0) && (M2 == 1 || N2 == 1) && (M2 * N2 == colon_size);
};

template<Integer M1, Integer N1, Integer M2, Integer N2, class Array2, class Colon_1>
struct assign_item
{
    static_assert(is_assignment_valid<M1,N1,Colon_1,M2,N2>::value, "invalid assignment");
};

template<Integer M1, Integer N1, class Scalar, class Colon_1>
struct assign_item_scalar
{
    static_assert(is_assignment_valid_scalar<M1,N1,Colon_1>::value, "invalid assignment");
};

template<Integer Row, Integer Col, class... Items>
struct get_assignment
{};

template<Integer Row, Integer Col, class Colon_1, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon
{};

template<Integer Row, Integer Col, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row,Col,colon_all,Mat_Row,Mat_Col>
{
    static_assert(Col <= Mat_Col && Col > 0 && Row <= Mat_Row * Mat_Col 
                  && Row > 0, "position out of matrix range");
    static const bool value = true;
};
template<Integer Row, Integer Col, Integer Pos, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row,Col,colon<Pos>,Mat_Row,Mat_Col>
{
    static_assert(Col <= Mat_Col && Col > 0 && Row <= Mat_Row * Mat_Col 
                  && Row > 0, "position out of matrix range");

    static const Integer pos_col = (Pos - 1) / Mat_Row + 1;
    static const Integer pos_row = Pos - (pos_col - 1) * Mat_Row;

    static const bool value = (pos_col == Col && pos_row == Row);
};
template<Integer Row, Integer Col, Integer Start, Integer End, Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row,Col,colon2<Start,End>,Mat_Row,Mat_Col>
{
    static_assert(Col <= Mat_Col && Col > 0 && Row <= Mat_Row * Mat_Col 
                  && Row > 0, "position out of matrix range");

    static const Integer pos = (Col - 1) * Mat_Row + Row;

    static const bool value = (pos >= Start && pos <= End);
};
template<Integer Row, Integer Col, Integer Start, Integer Step, Integer End, 
                    Integer Mat_Row, Integer Mat_Col>
struct is_in_colon<Row,Col,colon3<Start,Step,End>,Mat_Row,Mat_Col>
{
    static_assert(Col <= Mat_Col && Col > 0 && Row <= Mat_Row * Mat_Col 
                  && Row > 0, "position out of matrix range");

    static const Integer pos    = (Col - 1) * Mat_Row + Row;
    static const Integer dif    = pos - Start;


    static const bool value     = dif == 0 || (dif > 0 && Step > 0 && dif % Step == 0 && pos >= Start && pos <= End)
                                || (dif < 0 && Step < 0 && ((-dif) % (-Step) == 0) && pos >= End && pos <= Start);
};

template<Integer Row, Integer Col, class Item>
struct is_member_assign
{
    static_assert(details::dependent_false<Item>::value, "invalid item type, this should not happened");
};

template<Integer Row, Integer Col, Integer M1, Integer N1, Integer M2, Integer N2, class Array2,
        class Colon_1>
struct is_member_assign<Row,Col,assign_item<M1,N1,M2,N2,Array2,Colon_1>>
{
    static const bool value = is_in_colon<Row, Col, Colon_1, M1, N1>::value;
};
template<Integer Row, Integer Col, Integer M1, Integer N1, class Scalar, class Colon_1>
struct is_member_assign<Row,Col,assign_item_scalar<M1,N1,Scalar,Colon_1>>
{
    static const bool value = is_in_colon<Row,Col, Colon_1, M1, N1>::value;
};

template<Integer Pos, class Colon>
struct get_relative_pos{};

template<Integer Pos>
struct get_relative_pos<Pos,colon_all>
{
    static const Integer value = Pos;
};
template<Integer Pos, Integer Sel>
struct get_relative_pos<Pos,colon<Sel>>
{
    static_assert(Pos == Sel, "invalid element");
    static const Integer value = 1;
};
template<Integer Pos, Integer Start, Integer End>
struct get_relative_pos<Pos,colon2<Start,End>>
{
    static_assert(Pos >= Start && Pos <= End, "invalid element");
    static const Integer value = Pos - Start + 1;
};
template<Integer Pos, Integer Start, Integer Step, Integer End>
struct get_relative_pos<Pos, colon3<Start, Step, End>>
{
    static_assert(Step > 0 && (Pos >= Start && Pos <= End)
                  || Step < 0 && (Pos >= End && Pos <= Start) , "invalid element");

    static const Integer dif    = Pos - Start;

    static_assert(dif == 0 || dif > 0 && Step > 0 && dif % Step == 0 
                  || dif < 0 && Step < 0 && (-dif) % (-Step) == 0, "invalid element");
   
    static const Integer steps  = (dif == 0 ? 0 : dif / Step);
    static const Integer value  = steps + 1;
};

template<Integer Row, Integer Col, class Item>
struct get_assignment_item_impl
{};
template<Integer Row, Integer Col, Integer M1, Integer N1, Integer M2, Integer N2, class Array2,
        class Colon_1>
struct get_assignment_item_impl<Row,Col, assign_item<M1,N1,M2,N2,Array2,Colon_1>>
{
    static const Integer pos    = (Col - 1) * M1 + Row;
    static const Integer rel_pos= get_relative_pos<pos,Colon_1>::value;

    static_assert(M2 == 1 || N2 == 1, "rhs of assignment should be a vector");
    static const Integer size   = M2 * N2;
    
    static_assert(rel_pos > 0 && rel_pos <= size, "invalid element");

    static const Integer col2   = (rel_pos-1)/M2;
    static const Integer row2   = rel_pos - col2 * M2;

    using type = typename get_array_elem<Array2, row2, col2+1>::type;
};
template<Integer Row, Integer Col, Integer M1, Integer N1, class Scalar, class Colon_1>
struct get_assignment_item_impl<Row,Col, assign_item_scalar<M1,N1,Scalar,Colon_1>>
{
    using type = Scalar;
};

template<bool Found, Integer Row, Integer Col, class Item, class ... Items>
struct get_assignment_item
{
    using type = typename get_assignment<Row,Col,Items...>::type;
};
template<Integer Row, Integer Col, class Item, class ... Items>
struct get_assignment_item<true, Row, Col, Item, Items...>
{
    using type = typename get_assignment_item_impl<Row,Col, Item>::type;
};

template<Integer Row, Integer Col, class Item, class ... Items>
struct get_assignment<Row,Col,Item,Items...>
{
    static const bool found = is_member_assign<Row,Col,Item>::value;
    using type = typename get_assignment_item<found, Row, Col, Item, Items...>::type;
};
template<Integer Row, Integer Col>
struct get_assignment<Row,Col>
{
    static_assert(details::dependent_value_false<Integer,Row>::value,
                  "element is not assigned to given virtual matrix");
};

template<Integer Row, Integer Col, class Tag, class ... Args>
struct get_array_elem<mkd::virtual_array<Tag, Args...> , Row, Col>
{
    using type = typename get_assignment<Row, Col, Args...>::type;
};

template<Integer M1, Integer N1, class Tag, class Deps1,
            Integer M2, Integer N2, class Array2, class Deps2, class Colon_1, class ... Args>
struct mat_virtual_assign_1<ct_matrix<M1,N1,mkd::virtual_array<Tag, Args...>,Deps1>, 
        ct_matrix<M2,N2,Array2,Deps2>, Colon_1>
{
    using item      = assign_item<M1, N1, M2, N2, Array2, Colon_1>;
    using array_type= mkd::virtual_array<Tag, item, Args...>;
    using deps      = typename link_deps<Deps1, Deps2>::type;
    using type      = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, class Deps1, class Tag, class Array2, class Deps2,
    class Colon_1, class ... Args>
struct mat_virtual_assign_1<ct_matrix<M1,N1,mkd::virtual_array<Tag, Args...>,Deps1>, 
                            ct_scalar<Array2,Deps2>, Colon_1>
{
    using item          = assign_item_scalar<M1, N1, ct_scalar<Array2,Deps2>, Colon_1>;
    using array_type    = mkd::virtual_array<Tag, item, Args...>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, class Array1, class Deps1,
            Integer M2, Integer N2, class Array2, class Deps2, class Colon_1>
struct mat_virtual_assign_1<ct_matrix<M1,N1,Array1,Deps1>, 
        ct_matrix<M2,N2,Array2,Deps2>, Colon_1>
{
    static_assert(details::dependent_false<Array1>::value, "virtual matrix is required");
};

//----------------------------------------------------------------------------------
//                              compl_assign_1
//----------------------------------------------------------------------------------
template<Integer Pos, class Scalar>
struct assign_colon_scal{};

template<class Colon, class Matrix>
struct assign_colon{};

template<class Tag, Integer M, Integer N, class Array, class Deps,
        class Assignments, class RHS_Tag, class Scal_Deps, Integer Pos>
struct comp_assign_1<computation<Tag, ct_matrix<M,N,Array,Deps>,Assignments>, 
                    ct_scalar<RHS_Tag,Scal_Deps>, colon<Pos> >
{
    static_assert(Pos >= 1 && Pos <= M * N, "invalid position");

    using scalar            = ct_scalar<RHS_Tag,Scal_Deps>;
    using old_matrix        = ct_matrix<M,N,Array,Deps>;

    using new_assign        = assign_colon_scal<Pos,scalar>;
    using new_assignments   = typename list::push_back<Assignments,new_assign>::type;
    using new_deps          = typename link_deps<Deps, Scal_Deps>::type;
    using new_matrix        = ct_matrix<M,N,Array,new_deps>;

    using type              = computation<Tag,new_matrix,new_assignments>;
};

template<class Tag, Integer M, Integer N, class Array, class Deps,
        class Assignments, Integer M2, Integer N2, class Array2, class Deps2, 
        class Colon>
struct comp_assign_1<computation<Tag, ct_matrix<M,N,Array,Deps>,Assignments>, 
                    ct_matrix<M2,N2,Array2,Deps2>, Colon >
{
    static_assert(N2 == 1 && colon_func::size<Colon, M * N>::value == M2
                  && colon_func::first<Colon,M*N>::value >= 1
                  && colon_func::last<Colon,M*N>::value <= M*N, "invalid assignment");

    using old_matrix        = ct_matrix<M,N,Array,Deps>;
    using rhs_matrix        = ct_matrix<M2,N2,Array2,Deps2>;
    using colon_type        = Colon;

    using new_assign        = assign_colon<colon_type,rhs_matrix>;
    using new_assignments   = typename list::push_back<Assignments,new_assign>::type;
    using new_deps          = typename link_deps<Deps, Deps2>::type;
    using new_matrix        = ct_matrix<M,N,Array,new_deps>;

    using type              = computation<Tag,new_matrix,new_assignments>;
};

template<class Comp>
struct make_comp_result
{
    static_assert(details::dependent_false<Comp>::value, "this type should not be instantiated");
};
template<class Tag, Integer M, Integer N, class Array, class Deps, class Assignments_List>
struct make_comp_result<computation<Tag, ct_matrix<M,N,Array,Deps>, Assignments_List>>
{
    using new_dep           = dep<Tag,0,dep_computation>;
    using new_deps          = typename link_deps<dps<new_dep>, Deps>::type;
    using computation_type  = computation<Tag, ct_matrix<M, N, Array, Deps>, Assignments_List>;

    using type              = ct_matrix<M,N,Array,new_deps>;

    static_assert(list::is_member<new_dep,Deps>::value == false, "computation tag in use");    

    template<class Local_Storage, class Data_Provider, class Temp_Storage>
    inline_initializer
    friend void tag_initializer(Tag, Local_Storage& ls, Data_Provider& dp, Temp_Storage* ts)
    {
        (void)dp;
        ts->init_computation<Tag, computation_type, Local_Storage>(ls);
    };    
    template<class Visitor, class Temp_Storage>
    friend void tag_initializer_accept(Tag, Visitor& vis, Temp_Storage* ts)
    {
        ts->init_computation_accept<Tag, computation_type, Visitor>(vis);
    };

    template<class Subs_Context>
    friend void tag_printer(Tag, std::ostream& os, int tabs)
    {        
        print_computations_elems<Tag, computation_type> 
            ::eval<Subs_Context>(os, tabs);
    };
};

//----------------------------------------------------------------------------------
//                              assign_elem
//----------------------------------------------------------------------------------
template<class Matrix, Integer Row, Integer Col>
struct assign_elem
{
    static_assert(details::dependent_false<Matrix>::value, 
                "this type should not be instantiated");
};
template<Integer M, Integer N, class Array, class Deps, Integer Row, Integer Col>
struct assign_elem<ct_matrix<M,N,Array,Deps>,Row,Col>
{
    template<class Val, class Local_Storage>
    inline_expr
    static void eval(const Val& elem, const Local_Storage& ls)
    {
        using elem_t = typename get_array_elem<Array,Row,Col>::type;
        elem_t::assign<Val>(elem,ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using elem_t = typename get_array_elem<Array,Row,Col>::type;
        elem_t::accept_assign<Visitor>(vis);
    };
};

//----------------------------------------------------------------------------------
//                              make_assign_info
//----------------------------------------------------------------------------------
template<class Mat>
struct make_return_dep
{
    static_assert(details::dependent_false<Mat>::value, 
                "this type should not be instantiated");
};
template<Integer M, Integer N, class Tag, class Dep>
struct make_return_dep<ct_matrix<M,N, mkd::output_array<Tag>,Dep>>
{
    using type          = extern_dep<Tag>;
};
template<Integer M, Integer N, class Output_Tag, Integer MR, Integer MC, class Dep>
struct make_return_dep<ct_matrix<M,N, mkd::temp_output_array<Output_Tag, MR, MC>,Dep>>
{
    using type          = temp_dep<Output_Tag, MR, MC>;
};


template<class LHS, class RHS, class Colon>
struct assign_info
{
    using rhs                           = RHS;
    using colon_type                    = Colon;
    using return_dep                    = typename make_return_dep<LHS>::type;
    using return_tag                    = typename return_dep::tag;

    //TODO: subs colon
    static const Integer step           = return_tag::step;

    static const Integer return_step    = colon_func::step<colon_type>::value * step;
    static const Integer return_offset  = colon_func::offset<colon_type>::value *step 
                                        + typename return_tag::template get_offset<1,1>::value;

    using ret_root_align                = typename return_tag::root_align_type;
};

template<class Mat>
struct make_assign_info
{
    static_assert(details::dependent_false<Mat>::value, 
                "this type should not be instantiated");
};
template<Integer M, Integer N, class Array1, class Array2, class Deps>
struct make_assign_info<ct_matrix<M,N,mkd::mat_assign_array<M,N,Array1,Array2>,Deps>>
{
    using lhs   = ct_matrix<M,N,Array1,Deps>;
    using rhs   = ct_matrix<M,N,Array2,Deps>;
    using type  = assign_info<lhs,rhs,colon_all>;
};

template<Integer M, Integer N, class Array1, Integer M2, Integer N2, class Colon, class Array2, class Deps>
struct make_assign_info<ct_matrix<M2,N2, mkd::mat_assign_array_colon<M,N,Array1,Colon,M2,N2,Array2>,Deps>>
{
    using lhs   = ct_matrix<M,N,Array1,Deps>;
    using rhs   = ct_matrix<M2,N2,Array2,Deps>;
    using type  = assign_info<lhs,rhs,Colon>;
};

}}