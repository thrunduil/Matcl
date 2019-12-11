#pragma once

#include <iosfwd>

#include "mkgen/TODO/utils/utils.h"
#include "mkgen/matrix/scalar.h"
#include "mkgen/TODO/matrix/ct_matrix_details.h"

namespace matcl { namespace mkgen
{

// access matrix element at position Pos (1-based)
template<Integer Pos>
struct colon{};

// access all elements
struct colon_all{};

// access matrix elements from index Start to index End (1-based),
// i.e. in Matlab's notation Start : End
template<Integer Start, Integer End>
struct colon2{};

// access matrix elements from index Start to index End with step Step
// (1-based), i.e. in Matlab's notation  Start : Step : End
template<Integer Start, Integer Step, Integer End>
struct colon3{};

// compile time matrix with size M x N storing symbolic elements in array Array_t.
//
// matrix 1x1 is not a scalar: for example multiplying 1x1 matrix and 2x2 matrix will
// produce an error. For a matrix with generic elements Array_t type is basically array<Tag>
// with unique Tag for each matrix. For matrices containing symbolic expressions Array_t is
// any type for which a template get_array_elem described below is specialized.
// Deps is a type representing dependencies from runtime values; must be specialization
// of dps type
template<Integer M, Integer N, class Array_t, class Deps>
struct ct_matrix
{
    public:
        // number of rows
        static const Integer    rows        = M;

        // number of columns
        static const Integer    cols        = N;

        // number of elements
        static const Integer    size        = M * N;

        // type storing elements of this matrix
        using array_type                    = Array_t;

        // type representing dependencies from runtime values
        using dps_type                      = Deps;

    public:
        // get element at position Pos as a scalar
        template<Integer Pos>
        static auto elem(colon<Pos>)        -> ct_scalar<details::scalar_data<details::scalar_mat_elem_1<ct_matrix, Pos>>, Deps>;

        // get element at row Row and column Col
        template<Integer Row, Integer Col>
        static auto elem(colon<Row>, colon<Col>) 
                                            -> ct_scalar<details::scalar_data<details::scalar_mat_elem_2<ct_matrix, Row, Col>>, 
                                                    Deps>;

        // get submatrix Mat(Colon_1, Colon_2)
        template<class Colon_1, class Colon_2>
        static auto sub(Colon_1, Colon_2)   -> typename details::submatrix_maker_2<ct_matrix, 
                                                        Colon_1, Colon_2>::type;

        // get submatrix Mat(Colon_1)
        template<class Colon_1>
        static auto sub(Colon_1)            -> typename details::submatrix_maker_1<ct_matrix, 
                                                        Colon_1>::type;

        // build virtual matrix, this type must be a virtual_matrix
        //TODO
        template<class Colon_1, class Mat>
        static auto assign_1(Colon_1, Mat)  -> typename mat_virtual_assign_1<ct_matrix, Mat, Colon_1>::type;

        // store current results in temporary matrix (placed on the stack) if cond == true,
        // otherwise return this matrix.
        template<class Tag, bool Force = false>
        static auto make_temp()             -> typename mat_temporary<ct_matrix, Tag, Force>::type;

        // print matrix
        template<class Subs_Context>
        static void print(std::ostream& os, int tabs);

        //TODO: add compute function
};

}}

#include "mkgen/TODO/matrix/matrix.h"
