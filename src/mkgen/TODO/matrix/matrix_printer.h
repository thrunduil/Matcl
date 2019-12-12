#pragma once

#include "mkgen/TODO/matrix/ct_matrix.h"
#include "mkgen/TODO/expression/ct_matrix_expr.inl"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/evaler/dependency.h"
#include "mkgen/TODO/matrix/ct_matrix_details.h"

namespace matcl { namespace mkgen
{

//----------------------------------------------------------------------------------
//                              matrix printer
//----------------------------------------------------------------------------------
template<Integer M, Integer N, class Colon, class Array, Integer Row, Integer Col, Integer Mat_Rows>
struct print_matrix_elems_impl
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int tabs)
    {
        (void)prior;

        print_tabs(os, tabs);

        using elem = typename get_array_elem<Array,Row+1,Col+1>::type;
        elem::print<Subs_Context>(os,details::prior_start);

        os << "\n";

        print_matrix_elems_impl<M,N,Colon,Array,Row+1,Col,Mat_Rows>
            ::eval<Subs_Context>(os, details::prior_start, tabs);
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int tabs)
    {
        using subs  = decltype(get_substitution(Subs_Context(),Tag()));

        static const Integer pos0   = Row+1 + Col*M;
        static const Integer pos    = get_pos_colon<pos0,Colon>::value;

        print_tabs(os, tabs);
        print_subs(subs(), pos, os, Mat_Rows);
        os << " = ";

        using elem = typename get_array_elem<Array,Row+1,Col+1>::type;
        elem::print<Subs_Context>(os,details::prior_start);

        os << "\n";

        print_matrix_elems_impl<M,N,Colon,Array,Row+1,Col,Mat_Rows>::eval_dep<Tag,Subs_Context>(os, tabs);
    };
};

template<class Matrix_List>
struct print_matrix_list
{
    static_assert(details::dependent_false<Matrix_List>::value, 
                "this type should not be instantiated");
};

template<class Mat, class ... Matrix>
struct print_matrix_list<list::list<Mat, Matrix...>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        print_matrix_elems<Mat>::eval<Subs_Context>(os, details::prior_start, tabs);

        if constexpr(sizeof...(Matrix) > 0)
        {
            os << "\n";
        };

        using list_type = list::list<Matrix...>;
        using printer   = print_matrix_list<list_type>;
        printer::eval<Subs_Context>(os, tabs);
    };
};
template<>
struct print_matrix_list<list::list<>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        (void)os;
        (void)tabs;
    };
};

template<class Matrix_List, Integer Mat_Rows>
struct print_matrix_elems_expand
{
    static_assert(details::dependent_false<Matrix_List>::value, 
                "this type should not be instantiated");
};
template<Integer Mat_Rows>
struct print_matrix_elems_expand<list::list<>, Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int tabs)
    {
        (void)os;
        (void)prior;
        (void)tabs;
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int tabs)
    {
        (void)os;
        (void)tabs;
    };
};
template<class Mat_Colon, class ... Mats, Integer Mat_Rows>
struct print_matrix_elems_expand<list::list<Mat_Colon, Mats...>,Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int tabs)
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
            print_tabs(os, tabs);
            os << "loop vectorized \n";
        };    

        using array_t               = typename matrix_type::array_type;
        print_matrix_elems_impl<M,N,colon_assign,array_t,0,0,Mat_Rows>::eval<Subs_Context>(os,prior,tabs);

        if constexpr(sizeof...(Mats) > 0)
        {
            os << "\n";
        };

        using printer   = print_matrix_elems_expand<list::list<Mats...>,Mat_Rows>;
        printer::eval<Subs_Context>(os,prior,tabs);
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int tabs)
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
            print_tabs(os, tabs);
            os << "loop vectorized \n";
        };

        print_matrix_elems_impl<M,N,colon_assign,array_t,0,0,Mat_Rows>::eval_dep<Tag,Subs_Context>(os,tabs);

        if constexpr(sizeof...(Mats) > 0)
        {
            os << "\n";
        };

        using printer = print_matrix_elems_expand<list::list<Mats...>,Mat_Rows>;
        printer::eval_dep<Tag,Subs_Context>(os,tabs);
    };
};

template<class Mat>
struct print_matrix_elems
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int tabs)
    {
        using expand_vm = typename expand_virtual_matrix2<Mat>::type;
        Integer size    = list::size<expand_vm>::value;
        
        if (size > 1)
        {
            print_tabs(os, tabs);
            os << "[" << "\n";
        };

        print_matrix_elems_expand<expand_vm,Mat::rows>::eval<Subs_Context>(os,prior,tabs);        

        if (size > 1)
        {
            print_tabs(os, tabs);
            os << "]" << "\n";
        };
    };

    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int tabs)
    {
        using expand_vm = typename expand_virtual_matrix2<Mat>::type;
        Integer size    = list::size<expand_vm>::value;
        
        if (size > 1)
        {
            print_tabs(os, tabs);
            os << "[" << "\n";
        };

        print_matrix_elems_expand<expand_vm,Mat::rows>::eval_dep<Tag,Subs_Context>(os,tabs);

        if (size > 1)
        {
            print_tabs(os, tabs);
            os << "]" << "\n";
        };
    };
};

template<Integer M, Integer N, class Colon, class Ret_Tag, Integer Row, Integer Col, Integer Mat_Rows>
struct print_matrix_elems_impl<M,N,Colon,empty_array<Ret_Tag>,Row,Col,Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int tabs)
    {
        (void)prior;

        print_tabs(os, tabs);
        Ret_Tag::print(os,details::prior_start);
        os <<"\n";
    };
};

template<Integer M, Integer N, class Colon, class Array, Integer Col, Integer Mat_Rows>
struct print_matrix_elems_impl<M,N,Colon, Array,M,Col,Mat_Rows>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int prior, int tabs)
    {
        print_matrix_elems_impl<M,N,Colon,Array,0,Col+1,Mat_Rows>::eval<Subs_Context>(os, prior, tabs);
    };
    template<class Tag, class Subs_Context>
    static void eval_dep(std::ostream& os, int tabs)
    {
        print_matrix_elems_impl<M,N,Colon,Array,0,Col+1,Mat_Rows>::eval_dep<Tag,Subs_Context>(os, tabs);
    };
};
template<Integer M, Integer N, class Colon, class Array, Integer Mat_Rows>
struct print_matrix_elems_impl<M,N,Colon,Array,0,N,Mat_Rows>
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

template<class Dep>
struct print_dependency
{
    static_assert(details::dependent_false<Dep>::value, 
                "this type should not be instantiated");
};

template<class Dep_List>
struct print_matrix_deps;

template<class Tag, Integer Size>
struct print_dependency< dep<Tag,Size,dep_return>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {        
    };
};
template<class Tag, Integer Size, class Type>
struct print_dependency< dep<Tag,Size,Type>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {        
        bool need_init = tag_need_initialization(Subs_Context(), Tag());

        if (need_init == true)
        {
            os << "\n";

            print_tabs(os, tabs);
            os << "and" << "\n";

            os << "\n";

            tag_printer<Subs_Context>(Tag(), os, tabs);
        };
    };
};

template<class Deps_List>
struct print_matrix_deps
{
    static_assert(details::dependent_false<Deps_List>::value, 
                "this type should not be instantiated");
};

template<>
struct print_matrix_deps<dps<>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {};
};
template<class Dep1, class... Deps>
struct print_matrix_deps<dps<Dep1, Deps...>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        print_dependency<Dep1>::eval<Subs_Context>(os, tabs);        

        if constexpr(sizeof...(Deps) > 0)
        {
            using printer   = print_matrix_deps<dps<Deps...>>;
            printer::eval<Subs_Context>(os, tabs);
        };
    };
};

template<class Subject, class Assignments>
struct print_computations_assignments
{};

template<Integer M, Integer N, class Array, class Deps, class Assign_Type>
struct print_comp_assing_1
{};
template<Integer M, Integer N, class Array, class Deps, Integer Pos, class Scalar>
struct print_comp_assing_1<M,N,Array,Deps,assign_colon_scal<Pos,Scalar>>
{
    using matrix = ct_matrix<M,N,Array,Deps>;

    static const Integer col            = (Pos-1)/M + 1;
    static const Integer row            = Pos - (col-1) * M;

    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        using elem = typename get_array_elem<Array,row,col>::type;
        print_tabs(os,tabs);
        elem::print<Subs_Context>(os, details::prior_start);

        os << " = ";

        Scalar::print<Subs_Context>(os, details::prior_start);
        os <<   "\n";
    };
};

template<Integer M, Integer N, class Array, class Deps, class Colon, class Mat, Integer Pos, Integer Size>
struct print_comp_assing_1_impl
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        static const Integer pos    = get_pos_colon<Pos+1,Colon>::value;
        using mat_array             = typename Mat::array_type;
        using elem                  = typename get_array_elem<mat_array,Pos+1,1>::type;

        print_comp_assing_1<M,N,Array,Deps,assign_colon_scal<pos,elem>>::eval<Subs_Context>(os,tabs);

        print_comp_assing_1_impl<M,N,Array,Deps,Colon,Mat,Pos+1,Size>::eval<Subs_Context>(os,tabs);
    };
};
template<Integer M, Integer N, class Array, class Deps, class Colon, class Mat, Integer Size>
struct print_comp_assing_1_impl<M,N,Array,Deps,Colon,Mat,Size,Size>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        (void)os;
        (void)tabs;
    };
};

template<Integer M, Integer N, class Array, class Deps, class Colon, class Mat>
struct print_comp_assing_1<M,N,Array,Deps,assign_colon<Colon,Mat>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        static const Integer size   = get_size_colon<Colon,M*N>::value;
        print_comp_assing_1_impl<M,N,Array,Deps,Colon,Mat,0,size>::eval<Subs_Context>(os,tabs);
    };
};
template<Integer M, Integer N, class Array, class Deps, class Assign_Type, class... Items>
struct print_computations_assignments<ct_matrix<M,N,Array,Deps>,list::list<Assign_Type,Items...>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        using matrix = ct_matrix<M,N,Array,Deps>;

        print_comp_assing_1<M,N,Array,Deps, Assign_Type>::eval<Subs_Context>(os,tabs);

        using printer = print_computations_assignments<matrix, list::list<Items...>>;
        printer::eval<Subs_Context>(os,tabs);
    };
};
template<class Subject>
struct print_computations_assignments<Subject,list::list<>>
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int)
    {
        (void)os;
    };
};

template<class Tag, class Computation>
struct print_computations_elems
{
    template<class Subs_Context>
    static void eval(std::ostream& os, int tabs)
    {
        print_tabs(os,tabs);
        os << "computation ";
        Tag::print(os,details::prior_start);
        os << ":" << "\n";

        using subject       = typename Computation::subject;
        using assignments   = typename Computation::assignments;

        print_computations_assignments<subject,assignments>
            ::eval<Subs_Context>(os, tabs);
    };
};

//--------------------------------------------------------
//              print matrix
//--------------------------------------------------------
inline void print_tabs(std::ostream& os, int tabs)
{
    if (tabs > 0)
        os << std::string(tabs, ' ');
};
template<Integer M, Integer N, class Array_t, class Deps>
template<class Subs_Context>
void ct_matrix<M,N,Array_t,Deps>::print(std::ostream& os, int tabs)
{
    using matrix_type           = ct_matrix<M,N,Array_t,Deps>;
    using expand_vm             = typename expand_virtual_matrix<matrix_type>::type;

    print_tabs(os, tabs);
    os << "matrix " << M << " x " << N << " with elements:" << "\n";

    print_tabs(os, tabs);
    os << "[" << "\n";

    print_matrix_list<expand_vm>::eval<Subs_Context>(os, tabs);

    print_tabs(os, tabs);
    os << "]" << "\n";

    using deps_temp  = typename Subs_Context::deps_temp;

    if constexpr(std::is_same<deps_temp, empty_deps>::value == false)
        print_matrix_deps<deps_temp>::eval<Subs_Context>(os, tabs);
};

}}