/*
* Quern: a sparse QR library.
* This code is in the public domain.
* - Robert Bridson
*/

#pragma once

#include "matcl-matrep/details/utils.h"

namespace quern
{

namespace md = matcl::details;

// Function return values: nonzero indicates an error
#define QUERN_OK 0
#define QUERN_INPUT_ERROR 1
#define QUERN_OUT_OF_MEMORY 2

template<class Val>
struct make_val_givens
{
    using VR    = typename md::real_type<Val>::type;
    using type  = std::pair<VR,Val>;
    using ret   = type;

    static bool is_swap(VR c)
    {
        return c == VR(2.0);
    }
    static ret eval_one()
    {
        return ret(VR(2.0),Val());
    };
    static ret eval(const Val& c, const Val& s)
    {
        return ret(real(c),s);
    };

    static void encode(Val* Q_value, int j, const ret& qvalue)
    {
        Q_value[j*2 + 0] = qvalue.first;
        Q_value[j*2 + 1] = qvalue.second;
    };
    static void decode(const Val* Q_value, int j, VR& c, Val& s)
    {
        c = matcl::real(Q_value[j*2 + 0]);
        s = Q_value[j*2 + 1];
    };
};

// Takes a compressed sparse row format m*n matrix (m>=n) and outputs the upper
// triangular R factor from QR, also in compressed sparse row format. If
// row_order is non-null, it should contain the order in which rows of A should be
// taken; if it is null, the natural ordering is used. Memory for R is
// allocated internally; use QUERN_free_result to free.
template<class Val>
int QUERN_compute_qr_without_q(int m, int n, const int* A_row_start, 
        const int* A_column_index, const Val* A_value, const int* row_order,
        int** ptr_R_row_start, int** ptr_R_column_index, Val** ptr_R_value);

// Takes a compressed sparse row format m*n matrix (m>=n) and outputs the QR
// factors. The storage for Q encodes a sequence of Givens rotations and row
// swaps; it is not directly usable as a standard sparse matrix. R, however,
// is stored in standard compressed sparse row format. If row_order is non-null
// it should contain the order in which rows of A should be taken; if it is
// null, the natural ordering is used. Memory for Q and R is allocated
// internally; use QUERN_free_result to free.
template<class Val>
int QUERN_compute_qr(int m, int n, const int* A_row_start, const int* A_column_index,
        const Val* A_value, const int* row_order, int** ptr_Q_row_start,
        int** ptr_Q_column_index, Val** ptr_Q_value, int** ptr_R_row_start,
        int** ptr_R_column_index, Val** ptr_R_value);

// Free the memory allocated during QR factorization (for either Q or R).
// After calling this, do not try to access the factor again!
template<class Val>
void QUERN_free_result(int* row_start, int* column_index, Val* value);

};