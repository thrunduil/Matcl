#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/dependency.h"
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
struct expr_assign : public mkd::scalar_data<expr_assign<T1, T2>>
{
    using this_type     = expr_assign<T1, T2>;

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
    static Ret eval(const Local_Storage& ls)
    {
        Ret val = T2::eval<Ret>(ls);
        T1::assign<Ret>(val, ls);

        return val;
    };

    template<class Visitor>
    static void accept(Visitor& vis)
    {
        T2::accept<Visitor>(vis);
        T1::accept_assign<Visitor>(vis);
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename expr_assign_arrays<Step, Arr_List, T1, T2> :: type;

    //TODO
    template<class Void>
    using simplify      = this_type;

    static constexpr bool is_simplified()   { return true; };
};

template<class Array, Integer Row, Integer Col>
struct mat_assign_array_colon_get_elem
{};

template<Integer M, Integer N, class Array1, class Colon, Integer M2, Integer N2, class Array2, 
        Integer Row, Integer Col>
struct mat_assign_array_colon_get_elem<mkd::mat_assign_array_colon<M,N,Array1,Colon, M2, N2, Array2>, Row, Col>
{
    static_assert(Row <= M2 && Col <= N2, "invalid elem");

    static const Integer Pos    = (Col-1)*M2 + Row;
    static const Integer Pos_1  = colon_func::index<Pos,Colon>::value;
    static const Integer Col_1  = (Pos_1 - 1) / M + 1;
    static const Integer Row_1  = (Pos_1 - 1) % M + 1;

    using elem_1    = typename Array1 :: template get_element<Row_1, Col_1>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;
    using new_item  = expr_assign<elem_1,elem_2>;
    using type      = new_item;
};

template<class Array, Integer Row, Integer Col>
struct mat_assign_array_get_elem
{};

template<Integer M, Integer N, class Array1, class Array2, Integer Row, Integer Col>
struct mat_assign_array_get_elem<mkd::mat_assign_array<M,N,Array1,Array2>, Row, Col>
{
    static_assert(Row <= M && Col <= N, "invalid elem");

    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename Array2 :: template get_element<Row, Col>::type;
    using new_item  = expr_assign<elem_1,elem_2>;
    using type      = new_item;
};

//----------------------------------------------------------------------------------
//                              mat_assign_array
//----------------------------------------------------------------------------------

template<class Array, Integer Row, Integer Col>
struct mat_scal_assign_array_get_elem
{};

template<Integer M, Integer N, class Array1, Scal_data Array2, DPS Deps2,
        Integer Row, Integer Col>
struct mat_scal_assign_array_get_elem<mkd::mat_scal_assign_array<M,N, Array1, 
                            ct_scalar<Array2, Deps2>>, Row, Col>
{
    static_assert(Row <= M && Col <= N, "invalid elem");

    using elem_1    = typename Array1 :: template get_element<Row, Col>::type;
    using elem_2    = typename ct_scalar<Array2, Deps2>;
    using new_item  = expr_assign<elem_1,elem_2>;
    using type      = new_item;
};

//----------------------------------------------------------------------------------
//                              mat_assign
//----------------------------------------------------------------------------------
template<Integer M1_M2, Integer N1_N2, Mat_array Array1, DPS Deps1, Mat_array Array2, DPS Deps2>
struct mat_assign<ct_matrix<M1_M2,N1_N2,Array1,Deps1>, ct_matrix<M1_M2,N1_N2,Array2,Deps2>>
{
    using array_type    = mkd::mat_assign_array<M1_M2, N1_N2, Array1, Array2>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1_M2, N1_N2, array_type,deps>;
};

template<Integer M1, Integer N1, Mat_array Array1, DPS Deps1, 
            Scal_data Array2, DPS Deps2>
struct mat_assign<ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Array2, Deps2>>
{
    using array_type    = mkd::mat_scal_assign_array<M1, N1, Array1, ct_scalar<Array2,Deps2>>;
    using deps          = typename link_deps<Deps1, Deps2>::type;
    using type          = ct_matrix<M1, N1, array_type,deps>;
};

template<Integer M1, Integer N1, Integer M2, Integer N2, Mat_array Array1, DPS Deps1, 
            Mat_array Array2, DPS Deps2>
struct mat_assign<ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>>
{
    static_assert(M1 == M2 && N1 == N2, "invalid matrix operation, check size of matrices");
};

//----------------------------------------------------------------------------------
//                              compl_assign_1
//----------------------------------------------------------------------------------
template<Integer Pos, class Scalar>
struct assign_colon_scal{};

template<class Colon, class Matrix>
struct assign_colon{};

template<class Tag, Integer M, Integer N, Mat_array Array, DPS Deps,
        class Assignments, Scal_data RHS_Tag, DPS Scal_Deps, Integer Pos>
struct comp_assign_1<computation<Tag, ct_matrix<M,N,Array, Deps>,Assignments>, 
                    ct_scalar<RHS_Tag, Scal_Deps>, colon<Pos> >
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

template<class Tag, Integer M, Integer N, Mat_array Array, DPS Deps,
        class Assignments, Integer M2, Integer N2, Mat_array Array2, DPS Deps2, 
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
    static_assert(md::dependent_false<Comp>::value, "this type should not be instantiated");
};
template<class Tag, Integer M, Integer N, Mat_array Array, DPS Deps, class Assignments_List>
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
    static_assert(md::dependent_false<Matrix>::value, 
                "this type should not be instantiated");
};
template<Integer M, Integer N, Mat_array Array, DPS Deps, Integer Row, Integer Col>
struct assign_elem<ct_matrix<M,N,Array,Deps>,Row,Col>
{
    template<class Val, class Local_Storage>
    inline_expr
    static void eval(const Val& elem, const Local_Storage& ls)
    {
        using elem_t = typename Array :: template get_element<Row,Col>::type;
        elem_t::assign<Val>(elem,ls);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        using elem_t = typename Array :: template get_element<Row,Col>::type;
        elem_t::accept_assign<Visitor>(vis);
    };
};

//----------------------------------------------------------------------------------
//                              make_assign_info
//----------------------------------------------------------------------------------
template<class Mat>
struct make_return_dep
{
    static_assert(md::dependent_false<Mat>::value, 
                "this type should not be instantiated");
};
template<Integer M, Integer N, class Tag, DPS Dep>
struct make_return_dep<ct_matrix<M,N, mkd::output_array<Tag>, Dep>>
{
    using type          = extern_dep<Tag>;
};
template<Integer M, Integer N, class Output_Tag, Integer MR, Integer MC, DPS Dep>
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
                                        + return_tag::get_offset(1,1);

    using ret_root_align                = typename return_tag::root_align_type;
};

template<class Mat>
struct make_assign_info
{
    static_assert(md::dependent_false<Mat>::value, 
                "this type should not be instantiated");
};
template<Integer M, Integer N, class Array1, class Array2, DPS Deps>
struct make_assign_info<ct_matrix<M,N,mkd::mat_assign_array<M,N,Array1,Array2>,Deps>>
{
    using lhs   = ct_matrix<M,N,Array1,Deps>;
    using rhs   = ct_matrix<M,N,Array2,Deps>;
    using type  = assign_info<lhs,rhs,colon_all>;
};

template<Integer M, Integer N, class Array1, Integer M2, Integer N2, class Colon, class Array2, DPS Deps>
struct make_assign_info<ct_matrix<M2,N2, mkd::mat_assign_array_colon<M,N,Array1,Colon,M2,N2,Array2>,Deps>>
{
    using lhs   = ct_matrix<M,N,Array1,Deps>;
    using rhs   = ct_matrix<M2,N2,Array2,Deps>;
    using type  = assign_info<lhs,rhs,Colon>;
};

}}