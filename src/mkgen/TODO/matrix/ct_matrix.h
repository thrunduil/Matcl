#pragma once

#include <iosfwd>

#include "mkgen/matrix/scalar.h"
#include "mkgen/details/mkgen_fwd.h"
#include "mkgen/mkgen_fwd.h"

namespace matcl { namespace mkgen
{
//------------------------------------------------------------------------------
//                      Overwiev
//------------------------------------------------------------------------------
// This header implements compile time matrix, that stores symbolic elements as
// types. all numerical operations are represented as types, no data is stored
// or created. After defining computation one can create optimizer by defining
// appropriate operations on types or creare evaluators. types defined in this 
// header are intended for creating highly optimized small kernels with complete
// loop unrolling.

//------------------------------------------------------------------------------
//                      EXAMPLE
//------------------------------------------------------------------------------
// form GEMM product: C = A*B + a*C
//
// matrix dimensions:
//     static const Integer M = 4;
//     static const Integer K = 4;
//     static const Integer N = 4;
// 
// tags for matrices and scalars:
//     struct tag_A{ static void print(std::ostream& os) { os << "A"; }};
//     struct tag_B{ static void print(std::ostream& os) { os << "B"; }};
//     struct tag_C{ static void print(std::ostream& os) { os << "C"; }};
//     struct tag_a{ static void print(std::ostream& os) { os << "a"; }};
// 
// Define expression GEMM_expr:
//
//     using GEMM_expr = 
//         mat_assign
//         <
//             ct_matrix<M,N,Array<tag_C>>,
// 
//             mat_plus
//             <
//                 mat_mult< ct_matrix<M,K,Array<tag_A>>, 
//                             ct_matrix<K,N,Array<tag_B>> > :: type,
//                 mat_mult< ct_scalar<elem<tag_a>>, 
//                             ct_matrix<M,N,Array<tag_C>> > :: type
//             > :: type
//
//         >::type;
// 
// GEMM_expr::print(std::cout);

//------------------------------------------------------------------------------
//                      EXAMPLE EVAL
//------------------------------------------------------------------------------
// The same problem, but tags must be modified:
// 
// struct dp_GEMM;
//
// struct tag_A
// { 
//     static void     print(std::ostream& os) { os << "A"; }        
//     template<Integer pos>
//     static double&  get_write(dp_GEMM* dp)  { return dp->ptr_A[pos]; }
//     template<Integer pos>
//     static double   get_read(dp_GEMM* dp)   { return dp->ptr_A[pos]; }
// };
// struct tag_B
// { 
//     static void     print(std::ostream& os) { os << "B"; }
//     template<Integer pos>
//     static double&  get_write(dp_GEMM* dp)  { return dp->ptr_B[pos]; }
//     template<Integer pos>
//     static double   get_read(dp_GEMM* dp)   { return dp->ptr_B[pos]; }
// };
// struct tag_C
// { 
//     static void     print(std::ostream& os) { os << "B"; }
//     template<Integer pos>
//     static double&  get_write(dp_GEMM* dp)  { return dp->ptr_C[pos]; }
//     template<Integer pos>
//     static double   get_read(dp_GEMM* dp)   { return dp->ptr_C[pos]; }
// };
// 
// struct tag_a
// { 
//     static void     print(std::ostream& os) { os << "a"; }
//     static double   get(dp_GEMM* dp)        { return dp->scal_a; };   
// };
//
// data provider:
//
// struct dp_GEMM
// {
//     double*         ptr_A;
//     double*         ptr_B;
//     double*         ptr_C;
//     double          scal_a;
// 
//     template<class Tag, class Val, Integer row, Integer col, Integer pos>
//     Val&    get_elem_write()    { return Tag::get_write<pos>(this); }
// 
//     template<class Tag, class Val, Integer row, Integer col, Integer pos>
//     Val     get_elem_read()     { return Tag::get_read<pos>(this); }
// 
//     template<class Tag, class Val>
//     Val     get_scalar()        { return Tag::get(this); }
// };
//
// evaluation:
//     dp_GEMM dp; <... initialization ...>
//     expr_evaler<GEMM_expr, dp_GEMM, double>::eval(dp);

//------------------------------------------------------------------------------
//                      Basic Types
//------------------------------------------------------------------------------
// This class represents a generic element of a matrix tagged with Tag at row
// row, column col.
template<class Tag, Integer Row, Integer Col>
struct element : public mkd::scalar_data<element<Tag, Row, Col>>
{
    //TODO: add checks
    static const Integer    get_offset      = Tag::get_offset(Row,Col);
    using                   root_align_type = typename Tag::root_align_type;

    using tag       = Tag;

    // print element
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        (void)prior;

        Tag::print(os, details::prior_start);
        os << "[" << Row << "," << Col << "]";
    };

    //store result of type Val of computation in Array supplied by data_provided
    template<class Val, class Local_Storage>
    inline_lev_1
    static void assign(const Val& v, const Local_Storage& ls)
    {
        const_cast<Val&>(ls.get_extern<Tag,Row,Col>()) = v;
    };    

    template<class Val, class Data_Provider, class Temp_Storage>
    inline_lev_1
    static const Val* get_array(const Data_Provider& dp, const Temp_Storage* ev)
    {        
        return Tag::get_data_ptr<Val,0>(dp);
    };

    //get data of type Val associated with this element supplied by data_provided
    template<class Val, class Local_Storage>
    inline_lev_1
    static const Val& eval(const Local_Storage& ls)
    {
        return ls.get_extern<Tag,Row,Col>();
    };

    template<class Loop_Storage, class Val, class Local_Storage>
    inline_lev_1
    static void eval_loop(Val& ret, Integer off, const Local_Storage& cont)
    {
        return Loop_Storage::get_value<Val,element>(cont, ret, off);
    };

    template<Integer Step, class Arr_List>
    using get_arrays    = typename list::push_back<Arr_List, details::array_item<element,Step * Tag::step,
                            details::array_item_extern>> :: type;

    template<class Visitor>
    static void accept_assign(Visitor& vis)
    {
        vis.visit_store();
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        vis.visit_load();
    };
};

template<class Elem>
struct get_offset_elem_step
{};

template<class Tag, Integer R, Integer C>
struct get_offset_elem_step<element<Tag,R,C>>
{
    static const Integer value = element<Tag,R,C>::get_offset;
};
template<class Tag, Integer MR, Integer MC, Integer R, Integer C>
struct get_offset_elem_step<get_temporary<Tag,MR,MC,R,C>>
{
    using elem = get_temporary<Tag,MR,MC,R,C>;
    static const Integer value = elem::get_offset;
};
template<class Elem, Integer Step>
struct get_offset_elem_step<element_step<Elem,Step>>
{
    static const Integer value = Elem::get_offset;
};
template<class T1, class T2>
struct get_offset_elem_step<expr_plus<T1,T2>>
{
    //TODO: remove this
    static const Integer value = 0;
};

template<class Elem, Integer Step>
struct element_step : public mkd::scalar_data<element_step<Elem, Step>>
{
    // print element
    template<class Subs_Context>
    static void print(std::ostream& os, int prior)
    {
        Elem::print<Subs_Context>(os,prior);
    };

    //store result of type Val of computation in Array supplied by data_provided
    template<class Val, class Local_Storage>
    inline_lev_1
    static void assign(const Val& v, const Local_Storage& ls)
    {
        Elem::assign<Val>(v,ls);
    };

    static const Integer    get_offset = get_offset_elem_step<Elem>::value;

    template<class Val, class Data_Provider, class Temp_Storage>
    inline_lev_1
    static const Val* get_array(const Data_Provider& dp, const Temp_Storage* ev)
    {     
        return Elem::get_array<Val>(dp,ev);
    };

    //get data of type Val associated with this element supplied by data_provided
    template<class Val, class Local_Storage>
    inline_expr
    static Val eval(const Local_Storage& ls)
    {
        return Elem::eval<Val>(ls);
    };

    template<class Loop_Storage, class Val, class Local_Storage>
    inline_expr
    static void eval_loop(Val& ret, Integer off, const Local_Storage& cont)
    {
        return Elem::eval_loop<Loop_Storage,Val>(ret,off,cont);
    };

    template<Integer Step_loc, class Arr_List>
    using get_arrays    = typename Elem::template get_arrays<Step*Step_loc,Arr_List>;

    template<class Visitor>
    static void accept_assign(Visitor& vis)
    {
        Elem::accept_assign<Visitor>(vis);
    };
    template<class Visitor>
    static void accept(Visitor& vis)
    {
        Elem::accept<Visitor>(vis);
    };
};

//------------------------------------------------------------------------------
//                      Numerical Operations
//------------------------------------------------------------------------------
// matrix or scalar multiplication
template<class Mat, class Mat2>
struct mat_mult {};

// matrix or scalar element by element division
template<class Mat, class Mat2>
struct mat_div {};

// Mat .* (D * J), J = ones(1, cols)
template<class Mat, class D>
struct make_mult_rows {};

// Mat .* (J * D), J = ones(rows, 1)
template<class Mat, class D>
struct make_mult_cols {};

// Mat1 .* mat2
template<class Mat1, class Mat2>
struct make_mult_mat {};

template<class Tag, template<class Arg> class Func, class Mat1>
struct make_call_inline {};

template<class Mat>
struct is_value_matrix;

template<class Tag, class Func, class Mat1, bool Is_Value = is_value_matrix<Mat1>::value>
struct make_call_external {};

// matrix or scalar addition
template<class Mat1, class Mat2>
struct mat_plus {};

// matrix or scalar substraction
template<class Mat1, class Mat2>
struct mat_minus {};

// matrix transposition
template<class Mat1>
struct mat_trans {};

// matrix conjugate transposition
template<class Mat1>
struct mat_ctrans {};

// assignment operator
template<class Mat1, class Mat2>
struct mat_assign {};

// assignment in computation
template<class Subject, class Mat, class Colon_1>
struct comp_assign_1
{
    static_assert(md::dependent_false<Subject>::value, "this type should not be instantiated");
};

// Create temporary buffer 
template<class Mat, class Tag, bool Force>
struct mat_temporary;

// unary minus
template<class Mat>
struct unary_minus {};

// unary function with Tag tag applied to matrix or scalar Mat
template<class Tag, class Mat>
struct func_unary {};

// binary element-wise function with tag tag applied to matrices or scalars
// Mat1 and mat2
template<class Tag, class Mat1, class Mat2>
struct func_bin {};

//------------------------------------------------------------------------------
//                      Working with matrices
//------------------------------------------------------------------------------
// get element of matrix matrix_type at row row and column col
template<class Matrix_Type, Integer Row, Integer Col>
struct get_elem{};

// Create evaluator for expression expr. it makes sense to create evaluators only
// for assing expression. Needed data are provided by data_provider. all computations
// are performed on values of type Val.
template<class Tag, class Code_Gen, class Ret_Matrix, class Matrix, class Val>
struct expr_evaler;

//TODO
#if 1
//------------------------------------------------------------------------------
//                      operator +
//------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        operator+(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_plus<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1, class Data2, class Deps2>
auto        operator+(ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_plus<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

template<class Data1,class Deps1, Integer M2, Integer N2, class Array2, class Deps2>
auto        operator+(ct_scalar<Data1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_plus<ct_scalar<Data1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<class Data1, class Deps1, class Data2, class Deps2>
auto        operator+(ct_scalar<Data1,Deps1>,ct_scalar<Data2,Deps2>)
                                        -> typename mat_plus<ct_scalar<Data1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

//------------------------------------------------------------------------------
//                      operator -
//------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        operator-(ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_minus<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1, class Data2, class Deps2>
auto        operator-(ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_minus<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

template<class Data1,class Deps1,Integer M2, Integer N2, class Array2, class Deps2>
auto        operator-(ct_scalar<Data1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_minus<ct_scalar<Data1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<class Data1, class Deps1, class Data2, class Deps2>
auto        operator-(ct_scalar<Data1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_minus<ct_scalar<Data1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

//------------------------------------------------------------------------------
//                      operator *
//------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        operator*(ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_mult<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1, class Data2, class Deps2>
auto        operator*(ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_mult<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

template<class Data1,class Deps1, Integer M2, Integer N2, class Array2, class Deps2>
auto        operator*(ct_scalar<Data1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_mult<ct_scalar<Data1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<class Data1, class Deps1, class Data2, class Deps2>
auto        operator*(ct_scalar<Data1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_mult<ct_scalar<Data1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

//------------------------------------------------------------------------------
//                      operator /
//------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        operator/(ct_matrix<M1,N1,Array1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_div<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1, class Data2, class Deps2>
auto        operator/(ct_matrix<M1,N1,Array1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_div<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

template<class Data1,class Deps1, Integer M2, Integer N2, class Array2, class Deps2>
auto        operator/(ct_scalar<Data1,Deps1>, ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename mat_div<ct_scalar<Data1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<class Data1, class Deps1, class Data2, class Deps2>
auto        operator/(ct_scalar<Data1,Deps1>, ct_scalar<Data2,Deps2>)
                                        -> typename mat_div<ct_scalar<Data1,Deps1>,
                                                    ct_scalar<Data2,Deps2>>::type;

//------------------------------------------------------------------------------
//                      mult elem by elem
//------------------------------------------------------------------------------
template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        mult_rows(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename make_mult_rows<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        mult_cols(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename make_mult_cols<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;

template<Integer M1, Integer N1, class Array1, class Deps1,
         Integer M2, Integer N2, class Array2, class Deps2>
auto        mult(ct_matrix<M1,N1,Array1,Deps1>,ct_matrix<M2,N2,Array2,Deps2>)
                                        -> typename make_mult_mat<ct_matrix<M1,N1,Array1,Deps1>,
                                                    ct_matrix<M2,N2,Array2,Deps2>>::type;
#endif

//------------------------------------------------------------------------------
//                      call
//------------------------------------------------------------------------------
template<class Tag, template<class Arg> class Func,
         Integer M1, Integer N1, class Array1, class Deps1>
auto        call_inline(ct_matrix<M1, N1, Array1, Deps1>)
                                        -> typename make_call_inline<Tag, Func, 
                                                ct_matrix<M1, N1, Array1, Deps1>> :: type;

template<class Tag, class Func,
         Integer M1, Integer N1, class Array1, class Deps1>
auto        call_external(ct_matrix<M1, N1, Array1, Deps1>)
                                        -> typename make_call_external<Tag, Func, 
                                                ct_matrix<M1, N1, Array1, Deps1>> :: type;
}}