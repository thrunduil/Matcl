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
//                      Numerical Operations
//------------------------------------------------------------------------------

template<class Tag, template<class Arg> class Func, class Mat1>
struct make_call_inline {};

template<class Mat>
struct is_value_matrix;

template<class Tag, class Func, class Mat1, bool Is_Value = is_value_matrix<Mat1>::value>
struct make_call_external {};

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

//------------------------------------------------------------------------------
//                      call
//------------------------------------------------------------------------------
template<class Tag, template<class Arg> class Func,
         Integer M1, Integer N1, Mat_array Array1, DPS Deps1>
auto        call_inline(ct_matrix<M1, N1, Array1, Deps1>)
                                        -> typename make_call_inline<Tag, Func, 
                                                ct_matrix<M1, N1, Array1, Deps1>> :: type;

template<class Tag, class Func,
         Integer M1, Integer N1, Mat_array Array1, DPS Deps1>
auto        call_external(ct_matrix<M1, N1, Array1, Deps1>)
                                        -> typename make_call_external<Tag, Func, 
                                                ct_matrix<M1, N1, Array1, Deps1>> :: type;
}}