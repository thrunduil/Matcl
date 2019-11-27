/*
* Quern: a sparse QR library.
* This code is in the public domain.
* - Robert Bridson
*
*   Modifications:
*       2016, Pawel Kowal:  template version    
*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>

#include "matcl-linalg/decompositions/quern/include/quern.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-linalg/decompositions/quern/include/quern_list.h"
#include "matcl-linalg/decompositions/givens.h"

namespace quern
{

using namespace matcl;

template<class Val>
struct SparseEntry
{
   int      index;
   Val      value;
   SparseEntry() {}

   SparseEntry(int index_, const Val& value_)
       : index(index_), value(value_)
   {}
};

template<class Val>
static bool copy_row (int nnz, const int* index, const Val* value, quern::List<SparseEntry<Val>>& x)
{
    for (int i = nnz-1; i >= 0; --i) 
    {
        if( value[i] != Val() )
        {
            if(!x.push_front(SparseEntry<Val>(index[i], value[i])))
                return false;
        }
    }
    return true;
}

template<class Val>
static bool apply_givens(quern::List<SparseEntry<Val>>& x,
                         quern::List<SparseEntry<Val>>& y, int diagonal, 
                         typename md::real_type<Val>::type& c, Val& s)
{
    (void)diagonal;
    assert(!x.empty() && x.front().index == diagonal && x.front().value);
    assert(!y.empty() && y.front().index == diagonal && y.front().value);

    // find the rotation we need
    Val a   = x.front().value;
    Val b   = y.front().value;

    Val r;
    givens(a, b, c, s, r);
    Val s2  = -conj(s);

    // rotate the start of each list
    x.front().value = r;
    y.pop_front();
   
    // then update the rest (x_new = c*x+s*y and y_new = s2*x+c*y)
    quern::ListIterator<SparseEntry<Val>> p = x.begin(), q = y.begin();

    ++p; // skip the first value we already took care of

    while(p.still_going() && q.still_going())
    {
        if(p->index==q->index)
        {
            Val xnew    = c*p->value  + s*q->value;
            Val ynew    = s2*p->value + c*q->value;

            if (xnew != Val())
            {
                p->value    = xnew;
                ++p;
            }
            else
                x.erase(p);

            if (ynew != Val())
            {
                q->value    = ynew;
                ++q;
            }
            else
                y.erase(q);

        }
        else if (p->index < q->index)
        {
            int k       = p->index;
            Val xnew    = c * p->value;
            Val ynew    = s2 * p->value;

            p->value    = xnew;
            ++p;

            if ( !y.insert(q, SparseEntry<Val>(k, ynew)))
                return false;

            ++q;
        }
        else
        {
            int k       = q->index;
            Val xnew    = s*q->value;
            Val ynew    = c*q->value;

            if( !x.insert(p, SparseEntry<Val>(k, xnew))) 
                return false;
         
            ++p;
            q->value    = ynew;
            ++q;
        }
    };

    if(p.still_going())
    {
        do
        {
            int k       = p->index;
            Val xnew    = c*p->value;
            Val ynew    = s2*p->value;
            p->value    = xnew;
            ++p;

            if ( !y.insert(q, SparseEntry<Val>(k, ynew)))
                return false;

            ++q;
        }
        while(p.still_going());
    }
    else if (q.still_going())
    {
        do
        {
            int k       = q->index;
            Val xnew    = s*q->value;
            Val ynew    = c*q->value;

            if ( !x.insert(p, SparseEntry<Val>(k, xnew)))
                return false;

            ++p;
            q->value    = ynew;
            ++q;
        }
        while(q.still_going());
    }
    return true;
}

template<class Val>
static void givens(const Val& a, const Val& b, typename md::real_type<Val>::type& c, Val& s, Val& r)
{
    return matcl::construct_givens(a,b,c,s, r);
}

template<class Val>
int quern::QUERN_compute_qr(int m, int n, const int* A_row_start, const int* A_column_index,
    const Val* A_value, const int* row_order, int** ptr_Q_row_start, int** ptr_Q_column_index,
    Val** ptr_Q_value, int** ptr_R_row_start, int** ptr_R_column_index, Val** ptr_R_value)
{
    if ( m<=0 || n<=0 || !A_row_start || !A_column_index || !A_value
        || !ptr_Q_row_start || !ptr_Q_column_index || !ptr_Q_value
        || !ptr_R_row_start || !ptr_R_column_index || !ptr_R_value )
    {
        return QUERN_INPUT_ERROR;
    };

    using VR                = typename md::real_type<Val>::type;
    using Val_givens        = typename make_val_givens<Val>::type;
    using SparseVector_R    = quern::List<SparseEntry<Val>>;
    using SparseVector_Q    = quern::List<SparseEntry<Val_givens>>;

    n            = std::min(n,m);

    // set up lists for dynamically building Q and R
    quern::Pool<SparseEntry<Val>>           pool_R;
    quern::Pool<SparseEntry<Val_givens>>    pool_Q;

    SparseVector_Q* Q  = (SparseVector_Q*)std::malloc( m * sizeof(SparseVector_Q) );

    if(!Q) 
       return QUERN_OUT_OF_MEMORY;

    SparseVector_R* R = (SparseVector_R*)std::malloc(m*sizeof(SparseVector_R));
    if(!R)
    {
        std::free(Q);
        return QUERN_OUT_OF_MEMORY;
    }

    for(int i=0; i<m; ++i)
    {
        Q[i].init(&pool_Q);
        R[i].init(&pool_R);
    }

    // do the Givens QR
    SparseVector_R row;
    row.init(&pool_R);

    VR c;
    Val s;

    for(int a = 0; a < m; ++a)
    {
        int i   = (row_order ? row_order[a] : a);

        if(!copy_row(A_row_start[i+1]-A_row_start[i], A_column_index+A_row_start[i],
                   A_value+A_row_start[i], row))
        {
            std::free(Q);
            std::free(R);
            return QUERN_OUT_OF_MEMORY;
        }

        quern::ListIterator<SparseEntry<Val_givens>> q = Q[a].begin();

        while( !row.empty() && row.front().index < a && row.front().index < n)
        {
            int j = row.front().index;

            if ( R[j].empty() || R[j].front().index > j)
            { 
                // swap?
                R[j].swap(row);
                Q[a].insert(q, SparseEntry<Val_givens>(j, make_val_givens<Val>::eval_one()));
                ++q;
            }
            else
            { 
                // use Givens
                if( !apply_givens<Val>(R[j], row, j, c, s))
                {
                    std::free(Q);
                    std::free(R);
                    return QUERN_OUT_OF_MEMORY;
                }

                Q[a].insert(q, SparseEntry<Val_givens>(j, make_val_givens<Val>::eval(c,s)));
                ++q;
            }
        }

        if(a<n)
        {
            R[a].swap(row);
            assert(R[a].empty() || R[a].front().index>=a);
        }
    };

    // transfer Q's lists to static CSR form
    int* Q_row_start = (int*)std::malloc((m+1)*sizeof(int));

    if(!Q_row_start)
    {
        std::free(Q);
        std::free(R);
        Q = nullptr;
        R = nullptr;
        return QUERN_OUT_OF_MEMORY;
    }

    Q_row_start[0] = 0;

    for(int i=0; i<m; ++i)
        Q_row_start[i+1] = Q_row_start[i] + Q[i].size();

    int Qnnz            = Q_row_start[m];
    int* Q_column_index = (int*)std::malloc(Qnnz*sizeof(int));

    if(!Q_column_index)
    {
        std::free(Q);
        std::free(R);
        std::free(Q_row_start);

        Q           = nullptr;
        R           = nullptr;
        Q_row_start = nullptr;
        return QUERN_OUT_OF_MEMORY;
    }
   
    static const int num_val_Q  = sizeof(Val_givens)/sizeof(VR);
    Val* Q_value                = (Val*)std::malloc(Qnnz*sizeof(Val)*2);

    if (!Q_value)
    {
        std::free(Q);
        std::free(R);
        std::free(Q_row_start);
        std::free(Q_column_index);

        Q               = nullptr;
        R               = nullptr;
        Q_row_start     = nullptr;
        Q_column_index  = nullptr;
        return QUERN_OUT_OF_MEMORY;
    }

    int j = 0;
    for(int i = 0; i < m; ++i)
    {
        quern::ListIterator<SparseEntry<Val_givens>> q;

        for (q = Q[i].begin(); q.still_going(); ++q)
        {
            Q_column_index[j]   = q->index;

            make_val_givens<Val>::encode(Q_value, j, q->value);
            ++j;
        }
    }

    std::free(Q);
    Q = nullptr;

    // transfer R's lists to static CSR form
    int* R_row_start = (int*)std::malloc((n+1)*sizeof(int));
    if(!R_row_start)
    {
        std::free(Q);
        std::free(R);
        std::free(Q_row_start);
        std::free(Q_column_index);
        std::free(Q_value);

        Q               = nullptr;
        R               = nullptr;
        Q_row_start     = nullptr;
        Q_column_index  = nullptr;
        Q_value         = nullptr;
        return QUERN_OUT_OF_MEMORY;
    }

    R_row_start[0] = 0;

    for(int i=0; i<n; ++i)
        R_row_start[i+1] = R_row_start[i] + R[i].size();

    int Rnnz = R_row_start[n];

    int* R_column_index = (int*)std::malloc(Rnnz*sizeof(int));

    if(!R_column_index)
    {
        std::free(Q);
        std::free(R);
        std::free(Q_row_start);
        std::free(Q_column_index);
        std::free(Q_value);
        std::free(R_row_start);

        Q               = nullptr;
        R               = nullptr;
        Q_row_start     = nullptr;
        Q_column_index  = nullptr;
        Q_value         = nullptr;
        R_row_start     = nullptr;
        return QUERN_OUT_OF_MEMORY;
    }

    Val* R_value = (Val*)std::malloc(Rnnz*sizeof(Val));
    
    if(!R_value)
    {
        std::free(Q);
        std::free(R);
        std::free(Q_row_start);
        std::free(Q_column_index);
        std::free(Q_value);
        std::free(R_row_start);
        std::free(R_column_index);

        Q               = nullptr;
        R               = nullptr;
        Q_row_start     = nullptr;
        Q_column_index  = nullptr;
        Q_value         = nullptr;
        R_row_start     = nullptr;
        R_column_index  = nullptr;
        return QUERN_OUT_OF_MEMORY;
    }

    j = 0;

    for(int i = 0; i < n; ++i)
    {
        quern::ListIterator<SparseEntry<Val>> p;

        for( p = R[i].begin(); p.still_going(); ++p)
        {
            R_column_index[j]   = p->index;
            R_value[j]          = p->value;
            ++j;
        }
    }

    std::free(R);
    R = nullptr;

    *ptr_Q_row_start    = Q_row_start;
    *ptr_Q_column_index = Q_column_index;
    *ptr_Q_value        = Q_value;
    *ptr_R_row_start    = R_row_start;
    *ptr_R_column_index = R_column_index;
    *ptr_R_value        = R_value;
    return QUERN_OK;
}

template<class Val>
int quern::QUERN_compute_qr_without_q(int m, int n, const int* A_row_start, const int* A_column_index,
        const Val* A_value, const int* row_order, int** ptr_R_row_start, int** ptr_R_column_index,
        Val** ptr_R_value)
{
    if( m <= 0 || n <= 0 || !A_row_start || !A_column_index || !A_value
        || !ptr_R_row_start || !ptr_R_column_index || !ptr_R_value)
    {
        return QUERN_INPUT_ERROR;
    };

    using VR                = typename md::real_type<Val>::type;
    using Val_givens        = typename make_val_givens<Val>::type;
    using SparseVector_R    = quern::List<SparseEntry<Val>>;
    using SparseVector_Q    = quern::List<SparseEntry<Val_givens>>;

    n            = std::min(n,m);

    // set up lists for dynamically building R
    quern::Pool<SparseEntry<Val>>           pool_R;

    SparseVector_R* R = (SparseVector_R*)std::malloc(m*sizeof(SparseVector_R));

    if(!R)
        return QUERN_OUT_OF_MEMORY;

    for(int i = 0; i < m; ++i) 
       R[i].init(&pool_R);

    // do the Givens QR
    SparseVector_R row;
    row.init(&pool_R);

    VR c;
    Val s;

    for (int a = 0; a < m; ++a)
    {
        int i   = (row_order ? row_order[a] : a);

        if( !copy_row(A_row_start[i+1]-A_row_start[i], A_column_index+A_row_start[i],
            A_value+A_row_start[i], row))
        {
             std::free(R);
            return QUERN_OUT_OF_MEMORY;
        }

        while( !row.empty() && row.front().index < a && row.front().index < n)
        {
            int j = row.front().index;

            if ( R[j].empty() || R[j].front().index > j)
            { 
                // swap?
                R[j].swap(row);
            }
            else
            { 
                // use Givens
                if(!apply_givens<Val>(R[j], row, j, c, s))
                {
                    std::free(R);
                    return QUERN_OUT_OF_MEMORY;
                }
            }
        }

        if (a < n)
        {
            R[a].swap(row);
            assert(R[a].empty() || R[a].front().index>=a);
        }
    }

    // transfer R's lists to static CSR form
    int* R_row_start = (int*)std::malloc((n+1)*sizeof(int));
    if (!R_row_start)
    {
        std::free(R);
        return QUERN_OUT_OF_MEMORY;
    }

    R_row_start[0] = 0;

    for (int i = 0; i < n; ++i)
        R_row_start[i+1] = R_row_start[i] + R[i].size();

    int Rnnz = R_row_start[n];

    int* R_column_index = (int*)std::malloc(Rnnz*sizeof(int));
    
    if(!R_column_index)
    {
        std::free(R);
        std::free(R_row_start);
        return QUERN_OUT_OF_MEMORY;
    }

    Val* R_value = (Val*)std::malloc(Rnnz*sizeof(Val));
    
    if(!R_value)
    {
        std::free(R);
        std::free(R_row_start);
        std::free(R_column_index);
        return QUERN_OUT_OF_MEMORY;
    }

    int j = 0;
    for(int i = 0; i < n; ++i)
    {
        quern::ListIterator<SparseEntry<Val>> p;

        for (p = R[i].begin(); p.still_going(); ++p)
        {
            R_column_index[j]   = p->index;
            R_value[j]          = p->value;
            ++j;
        }
    }

    std::free(R);

    *ptr_R_row_start    = R_row_start;
    *ptr_R_column_index = R_column_index;
    *ptr_R_value        = R_value;

    return QUERN_OK;
}

template <class Val>
void quern::QUERN_free_result(int* row_start, int* column_index, Val* value)
{
   std::free(row_start);
   std::free(column_index);
   std::free(value);
}

template int quern::QUERN_compute_qr(int m, int n, const int* A_row_start, const int* A_column_index,
        const Real* A_value, const int* row_order, int** ptr_Q_row_start, int** ptr_Q_column_index,
        Real** ptr_Q_value, int** ptr_R_row_start, int** ptr_R_column_index, Real** ptr_R_value);

template int quern::QUERN_compute_qr(int m, int n, const int* A_row_start, const int* A_column_index,
        const Float* A_value, const int* row_order, int** ptr_Q_row_start, int** ptr_Q_column_index,
        Float** ptr_Q_value, int** ptr_R_row_start, int** ptr_R_column_index, Float** ptr_R_value);

template int quern::QUERN_compute_qr(int m, int n, const int* A_row_start, const int* A_column_index,
        const Complex* A_value, const int* row_order, int** ptr_Q_row_start, int** ptr_Q_column_index,
        Complex** ptr_Q_value, int** ptr_R_row_start, int** ptr_R_column_index, Complex** ptr_R_value);

template int quern::QUERN_compute_qr(int m, int n, const int* A_row_start, const int* A_column_index,
        const Float_complex* A_value, const int* row_order, int** ptr_Q_row_start, int** ptr_Q_column_index,
        Float_complex** ptr_Q_value, int** ptr_R_row_start, int** ptr_R_column_index, Float_complex** ptr_R_value);

template int quern::QUERN_compute_qr_without_q(int m, int n, const int* A_row_start, const int* A_column_index,
        const Real* A_value, const int* row_order, int** ptr_R_row_start, int** ptr_R_column_index,
        Real** ptr_R_value);
template int quern::QUERN_compute_qr_without_q(int m, int n, const int* A_row_start, const int* A_column_index,
        const Float* A_value, const int* row_order, int** ptr_R_row_start, int** ptr_R_column_index,
        Float** ptr_R_value);
template int quern::QUERN_compute_qr_without_q(int m, int n, const int* A_row_start, const int* A_column_index,
        const Complex* A_value, const int* row_order, int** ptr_R_row_start, int** ptr_R_column_index,
        Complex** ptr_R_value);
template int quern::QUERN_compute_qr_without_q(int m, int n, const int* A_row_start, const int* A_column_index,
        const Float_complex* A_value, const int* row_order, int** ptr_R_row_start, int** ptr_R_column_index,
        Float_complex** ptr_R_value);

template void quern::QUERN_free_result(int* row_start, int* column_index, Real* value);
template void quern::QUERN_free_result(int* row_start, int* column_index, Float* value);
template void quern::QUERN_free_result(int* row_start, int* column_index, Complex* value);
template void quern::QUERN_free_result(int* row_start, int* column_index, Float_complex* value);


};
