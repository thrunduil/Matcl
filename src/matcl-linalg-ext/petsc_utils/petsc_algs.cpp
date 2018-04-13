/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#if 0
TODO

#include "petsc_algs.h"
#include "matcl-linalg/petsc_utils/petsc_option.h"
#include "matcl-linalg/petsc_utils/petsc_objects.h"
#include "matcl-linalg/petsc_utils/petsc_matrix.h"

#pragma warning (push)
#pragma warning (disable:4127)
#pragma warning (disable:4101)
#pragma warning (disable:4100)

#include "petsc/private/matimpl.h"

#pragma warning (pop)

#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/func_unary.h"

namespace matcl { namespace petsc
{

namespace mr = matcl::raw;

static void reorder_impl(const Matrix& mat, petsc_ordering ord_type, permvec& pr,
                         permvec& pc, bool sym)
{
    std::string ord_str     = petsc::converter_order().to_string((Integer)ord_type);
    ::MatOrderingType type  = ord_str.c_str();

    using petsc_mat         = std::shared_ptr<petsc::petsc_sparse_struct>;

    IS rperm;
    IS cperm;
    petsc_mat pmat;

    petsc::mmlib_to_sparse_struct(mat, true, pmat);
    PetscErrorCode err = ::MatGetOrdering(*pmat, type, &rperm, &cperm);

    petsc::IS_to_perm_vector(rperm, pr);

    if (sym == false)
        petsc::IS_to_perm_vector(cperm, pc);

    ::ISDestroy(&rperm);
    ::ISDestroy(&cperm);
    (void)err;
};

permvec petsc::reorder_sym(const Matrix& mat, petsc_ordering ord_type)
{
    permvec ret, pc;
    reorder_impl(mat, ord_type, ret, pc, true);
    return ret;
};

tuple<permvec,permvec>
petsc::reorder_nsym(const Matrix& mat, petsc_ordering ord_type)
{
    permvec pr, pc;
    reorder_impl(mat, ord_type, pr, pc, false);
    return tuple<permvec,permvec>(pc,pr);
};

template<class Val>
static Matrix make_coarsening(PetscCoarsenData* agg_lists)
{
    Integer N           = agg_lists->size;
    Integer n_coarse    = 0;

    for (Integer i = 0; i < N; ++i)
    {
        PetscCDIntNd* n = agg_lists->array[i];
        if (n) 
            ++n_coarse;
    };

    using Mat_S     = raw::Matrix<Val, struct_sparse>;
    using ccs_rep   = raw::details::sparse_ccs<Val>;
    Mat_S coarse(ti::ti_empty(), N, n_coarse, N);

    ccs_rep rep     = coarse.rep();
    Integer nz      = 0;

    Integer* d_c    = rep.ptr_c();
    Integer* d_r    = rep.ptr_r();
    Val* d_x        = rep.ptr_x();

    Integer head    = 0;

    for (Integer j = 0; j < n_coarse; ++j)
    {
        d_c[j]      = nz;

        //find new head
        while(agg_lists->array[head] == nullptr)
            ++head;

        PetscCDIntNd* el    = agg_lists->array[head];
        ++head;

        while(el)
        {
            Integer r       = el->gid;
            d_r[nz]         = r;
            d_x[nz]         = Val(1.0);
            el              = el->next;
            ++nz;            
        };
    };

    d_c[n_coarse]   = nz;

    return Matrix(coarse, false);
};

Matrix petsc::coarsen(const Matrix& A0, const permvec& perm, coarsen_type ct)
{
    Matrix A(A0);

    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    MatCoarsen crs;
    MatCoarsenType type;
    
    switch(ct)
    {
        case coarsen_type::mis:
            type = MATCOARSENMIS;
            break;
        case coarsen_type::hem:
        default:
            type = MATCOARSENHEM;
            break;
    }

    MatCoarsenCreate(PETSC_COMM_SELF,&crs);
    MatCoarsenSetType(crs,type);
    
    using petsc_matrix_ptr  = std::shared_ptr<details::petsc_matrix>;

    value_code vc           = A.get_value_code();
    bool is_compl           = matrix_traits::is_float_complex(vc);
    value_code vcr          = matrix_traits::real_value_type(vc);

    if (vc == value_code::v_object)
        throw error::object_value_type_not_allowed("coarsen");

    // petsc allows for real type only
    if (is_compl)
        A                   = abs(A);

    petsc_matrix_ptr pm_A   = details::petsc_matrix::create(A);

    std::shared_ptr<petsc::petsc_IS> is_perm;
    Matrix iperm = perm.to_matrix();
    petsc::mmlib_to_IS(iperm.get_impl<mr::integer_dense>(), A.rows(), is_perm); 
    MatCoarsenSetAdjacency(crs, pm_A->get_Aksp());
    MatCoarsenSetGreedyOrdering(crs, *is_perm);

    //nonoverlapping coarsening; effectively only this case is implemented
    MatCoarsenSetStrictAggs(crs, PETSC_TRUE);

    MatCoarsenApply(crs);

    PetscCoarsenData *agg_lists;
    MatCoarsenGetData(crs, &agg_lists);

    Matrix coarse;
    
    switch (vcr)
    {
        case value_code::v_float:
            coarse  = make_coarsening<Float>(agg_lists);
            break;
        default:
            coarse  = make_coarsening<Real>(agg_lists);
            break;
    };

    PetscCDDestroy(agg_lists);
    MatCoarsenDestroy(&crs);

    return coarse;
};

}};

#endif
