/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "linsolve_objects_decomp.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/decompositions/chol.h"
#include "matcl-linalg/decompositions/lu.h"
#include "matcl-linalg/decompositions/svd.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/decompositions/ldl.h"
#include "matcl-linalg/utils/optim_params.h"

namespace matcl { namespace details
{

struct make_linsolve_impl
{
    static linsolve_obj eval(const Matrix& A, const matcl::options& opts)
    {
        linsolve_obj ret = build(A, opts);
        ret = matcl::linsolve_extended_sol(ret, opts);
        return ret;
    };

    static linsolve_obj build(const Matrix& A, const matcl::options& opts)
    {
        struct_code sc  = A.get_struct_code();
        struct_flag sf  = A.get_struct();

        using data_ptr = linsolve_obj::linsolve_data_ptr;

        if (is_unitary(sf))
            return linsolve_unitary(A);

        Integer ld      = get_ld(A, -1);
        Integer ud      = get_ud(A, -1);

        if (ld == 0 && ud == 0)
            return linsolve_diag(A, opts);

        bool rank_revealing = opts.get_option<bool>(opt::linsolve::use_rr());

        if ((ld == 0 || ud == 0) && rank_revealing == false)
            return linsolve_triang(A, opts);

        if (rank_revealing == true)
        {
            if (sc == struct_code::struct_dense)
                return make_rr_dense(A, sf, opts);
            else
                return make_rr_nondense(A, sf, opts);
        }
        else
        {
            if (sc == struct_code::struct_dense)
                return make_standard_dense(A, ld, ud, sf, opts);
            else if (sc == struct_code::struct_banded)
                return make_standard_band(A, ld, ud, sf, opts);
            else
                return make_standard_sparse(A, ld, ud, sf, opts);
        };
    };

    static void change_pivot_partial(matcl::options& opts, const option& default_rr_pivot)
    {
        using pivot_type    = opt::linsolve::pivot_type;

        bool has_opt_piv    = opts.has_option(opt::linsolve::pivot());                    

        if (has_opt_piv == false)
        {
            opts.set(default_rr_pivot);
        }
        else
        {        
            Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
            pivot_type piv      = (pivot_type)piv_int;

            if (piv == pivot_type::partial)
                opts.set(default_rr_pivot);
        };
    };

    static linsolve_obj make_rr_nondense(const Matrix& A, struct_flag sf, const matcl::options& opts)
    {
        bool is_posdef  = matcl::is_posdef(sf) || matcl::is_semi_posdef(sf);

        using pivot_type    = opt::linsolve::pivot_type;

        matcl::option default_rr_pivot  = opt::linsolve::pivot((Integer)pivot_type::rook);

        if (is_posdef == false)
        {
            //possible factorizations: LU
                    
            options opts2       = opts;

            //change partial pivoting to default_rr_pivot
            change_pivot_partial(opts2, default_rr_pivot);
            return linsolve_lu(A, opts2);
        }
        else
        {
            return linsolve_chol(A, true, opts);
        }
    };

    static linsolve_obj make_rr_dense(const Matrix& A, struct_flag sf, const matcl::options& opts)
    {
        bool is_posdef  = matcl::is_posdef(sf) || matcl::is_semi_posdef(sf);

        if (is_posdef == true)
            return linsolve_chol(A, true, opts);

        using pivot_type    = opt::linsolve::pivot_type;

        options opts2       = opts;

        matcl::option default_rr_pivot  = opt::linsolve::pivot((Integer)pivot_type::rook);

        change_pivot_partial(opts2, default_rr_pivot);
        return linsolve_lu(A, opts2);
    };

    static linsolve_obj make_standard_sparse(const Matrix& A, Integer ld, Integer ud, struct_flag sf, 
                                                  const matcl::options& opts)
    {
        //possible factorizations: LU, CHOL

        bool is_posdef  = matcl::is_posdef(sf);

        if (is_posdef == true)
            return linsolve_chol(A, false, opts);

        Integer N           = A.rows();
        Integer max_lu      = Integer(Real(N) * sparse_qr_lu_threashold());
        bool small_ludiags  = ld < max_lu || ud < max_lu;
        bool is_tridiag     = ld == 1 && ud == 1;

        if (small_ludiags && is_tridiag == false)
            return linsolve_qr(A, false, opts);

        return linsolve_lu(A, opts);
    };

    static linsolve_obj make_standard_band(const Matrix& A, Integer ld, Integer ud, struct_flag sf, 
                                                  const matcl::options& opts)
    {
        //possible factorizations: LU, CHOL 

        bool is_posdef  = matcl::is_posdef(sf);

        if (is_posdef == true)
            return linsolve_chol(A, false, opts);

        Integer N           = A.rows();
        Integer max_lu      = Integer(Real(N) * band_qr_lu_threashold());
        bool small_ludiags  = ld < max_lu || ud < max_lu;
        bool is_tridiag     = ld == 1 && ud == 1;

        if (small_ludiags && is_tridiag == false)
            return linsolve_qr(A, false, opts);

        return linsolve_lu(A, opts);
    };

    static linsolve_obj make_standard_dense(const Matrix&  A, Integer ld, Integer ud,
                                               struct_flag sf, const matcl::options& opts)
    {
        //possible factorizations: LU, CHOL, QR, LDL

        bool is_posdef  = matcl::is_posdef(sf);
        bool is_sym     = sf.has_sym_flag();
        bool is_her     = sf.has_her_flag();

        if (is_posdef == true)
            return linsolve_chol(A, false, opts);

        if (is_sym == true || is_her == true)
        {
            bool can_use_ldl    = opts.get_option<bool>(opt::linsolve::can_use_ldl());

            if (can_use_ldl == true)
            {
                bool is_real    = matrix_traits::is_float_real(A.get_value_code());
                if (is_her == true || is_real)
                    return linsolve_ldl(A, false, opts);
                else
                    return linsolve_ldl(A, true,opts);
            }
        };

        Integer N           = A.rows();
        Integer max_lu      = Integer(Real(N) * dense_qr_lu_threashold());
        bool small_ludiags  = ld < max_lu || ud < max_lu;
        bool is_tridiag     = ld == 1 && ud == 1;

        if (small_ludiags && is_tridiag == false)
            return linsolve_qr(A, false, opts);

        return linsolve_lu(A, opts);
    }

    static Real dense_qr_lu_threashold()    { return linalg_optim_params::dense_qr_lu_threashold(); };
    static Real band_qr_lu_threashold()     { return linalg_optim_params::band_qr_lu_threashold(); };
    static Real sparse_qr_lu_threashold()   { return linalg_optim_params::sparse_qr_lu_threashold(); };
};

};}

namespace matcl
{

linsolve_obj matcl::make_linsolve_obj(const Matrix& A, const matcl::options& opts)
{
    return details::make_linsolve_impl::eval(A, opts);
};
linsolve_obj matcl::make_linsolve_obj(Matrix&& A, const matcl::options& opts)
{
    return details::make_linsolve_impl::eval(std::move(A), opts);
};

};
