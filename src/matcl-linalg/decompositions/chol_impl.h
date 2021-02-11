/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#pragma once

#include "extern/cholmod/cholmod.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/inplace.h"

namespace matcl { namespace details
{

namespace mr = matcl::raw;

template<class V, class S>
struct chol_diag
{
	using Mat           = raw::Matrix<V,S>;
	using DenseMatrix   = raw::Matrix<V,struct_dense>;

	static void eval(chol2_return_type& ret, const Mat& A)
	{
		DenseMatrix Ub = raw::converter<DenseMatrix,decltype(A.get_diag())>
                            ::eval(A.get_diag()).make_explicit();

		Integer K   = Ub.size();
		V* ptr      = Ub.ptr();

        using VR    = md::real_type<V>::type;

        permvec p   = permvec::identity(A.rows());

		for(Integer i = 0; i < K; ++i)
		{
            VR val      = matcl::real(ptr[i]);

            if (val <= VR(0.0))
            {
                ret = chol2_return_type(matcl::bdiag(Matrix(Ub, false), 0), p, i);
                return;
            };

			ptr[i]      = raw::details::sqrt_helper<V>::eval(val);
        }
        
		ret = chol2_return_type(matcl::bdiag(Matrix(Ub, false), 0), p, K);
        return;
	}
};

template<class V, class S> 
struct chol_impl{};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V> 
struct chol_impl<V,struct_dense> 
{
    using Mat = raw::Matrix<V,struct_dense>;

    static void eval(chol2_return_type& ret, const Mat& A, bool upper, const options&)
    {
		Integer N = A.rows();

		if (A.get_struct().is_id())
        {
            ret = chol2_return_type(Matrix(A,false), permvec::identity(N), N);
            return;
        }

		if (N == 0)
        {
            ret = chol2_return_type(Matrix(A,false),details::pv_constructor::make(izeros(1,0)), N);
            return;
        }

		if (is_diag(Matrix(A,false)))
			return chol_diag<V,struct_dense>::eval(ret, A);

		char uplo = (upper ? 'U' : 'L');
		
		matcl::lapack::i_type info;

        Mat Ac          = A.make_unique();
        Ac.get_struct().reset();

        Integer N_suc   = N;

		lapack::potrf(&uplo, N, lap(Ac.ptr()), Ac.ld(), &info);
        
		if (info < 0)
        {
            throw error::error_general("invalid argument passed to potrf");
			//throw error::error_nonposdef();
        };

        if (info > 0)
            N_suc       = info - 1;

        if (upper)
        {
            Matrix U;
            matcl::raw::inplace::make_triu<Mat>::eval(U,Ac,0);
            ret = chol2_return_type(U, permvec::identity(N), N_suc); 
        }
        else
        {
            Matrix U;
            matcl::raw::inplace::make_tril<Mat>::eval(U,Ac,0);
            ret = chol2_return_type(U, permvec::identity(N), N_suc);
        }

        return;
    };
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V> 
struct chol_impl<V,struct_banded> 
{
    using VR        = typename md::real_type<V>::type;
	using Mat       = raw::Matrix<V,struct_banded>;
    using Mat_D     = raw::Matrix<V,struct_dense>;
    using Mat_DR    = raw::Matrix<VR,struct_dense>;

    static void eval(chol2_return_type& ret, const Mat& A, bool upper, const options&)
	{
		Integer N = A.rows();

		if (A.get_struct().is_id())
        {
            ret = chol2_return_type(Matrix(A,false), permvec::identity(N), N);
            return;
        }

		if (N == 0)
        {
            ret = chol2_return_type(Matrix(A,false), details::pv_constructor::make(izeros(1,0)), N);
            return;
        }

		if (is_diag(Matrix(A,false)))
			return chol_diag<V,struct_banded>::eval(ret, A);

        if (A.has_diag(0) == false)
        {
            value_code vc   = matrix_traits::value_code<VR>::value;
            Matrix U        = spzeros(A.rows(), A.rows(), 0, vc);

            ret = chol2_return_type(U, permvec::identity(N), 0);
            return;
        }

        Integer udiags  = A.number_superdiagonals();
		char uplo       = (upper ? 'U' : 'L');		        
        Mat Ac          = A.make_unique();

        Ac.get_struct().reset_value();

        Integer N_suc   = N;
        Integer offset  = upper? 0 : Ac.first_elem_diag(0);
        V* A_ptr        = Ac.rep_ptr() + offset;

		matcl::lapack::i_type info;        

		lapack::pbtrf(&uplo, N, udiags, lap(A_ptr), Ac.ld(), &info);

		if (info < 0)
			throw error::error_general("invalid argument passed to pbtrf");

        if (info > 0)
            N_suc       = info - 1;

        if (upper)
        {
            Matrix U;
            matcl::raw::inplace::make_triu<Mat>::eval(U,Ac,0);
            ret = chol2_return_type(U, permvec::identity(N), N_suc);
        }
        else
        {
            Matrix U;
            matcl::raw::inplace::make_tril<Mat>::eval(U,Ac,0);
            ret = chol2_return_type(U, permvec::identity(N), N_suc);
        }

        return;
	};
};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
namespace cholmod_helpers
{

std::enable_if<sizeof(Integer)==sizeof(SuiteSparse_long), void*>::type cholmod_type_guard;

//Convert matcl raw sparse real matrix to `cholmod_sparse' format.
//NOTE #1: remember to clean up the returned cholmod_sparse object! (by call to `cholmod_free_sparse').
//     It is not refcounted! 
//NOTE #2: it will work only on symmetric matrices.
static cholmod_sparse* matcl_sparse_symreal_to_cholmod(const raw::Matrix<Real,struct_sparse> &M, 
                                                        cholmod_common& cc)
{
    int stype = -1;

    auto sprep      = M.rep();
    Integer nrow    = sprep.rows();
    Integer ncol    = sprep.cols();
    Integer nnz     = sprep.nnz();

    cholmod_sparse* ret = cholmod_l_allocate_sparse(nrow, ncol, nnz,
        1, // TRUE if columns of A sorted, FALSE otherwise
        1, // TRUE if A will be packed, FALSE otherwise
       stype, CHOLMOD_REAL, &cc);

    const Integer *ptr_c = sprep.ptr_c(),
                  *ptr_r = sprep.ptr_r();
    const Real    *ptr_x = sprep.ptr_x();

    Integer *target_c = static_cast<Integer*>(ret->p),
            *target_r = static_cast<Integer*>(ret->i);
    Real    *target_x = static_cast<Real*   >(ret->x);

    if(ret != NULL)
    {
        // Matrix is assumed to be symmetric, copy to cholmod just the tril part.
        Integer pos = 0;
        target_c[0] = 0;

        for (Integer c = 0; c < ncol; ++c)
        {
            Integer k       = ptr_c[c];
            Integer last_k  = ptr_c[c+1];

            //ignore upper-triangular part
            for (; k < last_k; ++k)
            {
                Integer r   = ptr_r[k];

                if (r < c)
                    continue;

                break;
            };

            for (; k < last_k; ++k)
            {
                target_x[pos] = ptr_x[k];
                target_r[pos] = ptr_r[k];
                ++pos;
            };

            target_c[c+1]   = pos;
        };
    }

    return ret;
}

static matcl::Matrix cholmod_to_matcl_sparse_real(const cholmod_sparse* A, Integer N_suc)
{
    matcl::Matrix ret;
    if(A != NULL)
    {
        Integer nrow  = (int)A->nrow;
        Integer nzmax = (int)A->nzmax;

        ret     = matcl::make_real_sparse_noinit(nrow, nrow, nzmax);

        using rep_type  = sparse_matrix<Real,false>;
        rep_type sprepu(ret);

        Integer *target_c = sprepu.ptr_c(),
                *target_r = sprepu.ptr_r();
        Real    *target_x = sprepu.ptr_x();
        nzmax             = sprepu.nzmax();

        if (N_suc == 0)
        {
            for (Integer i = 0; i <= nrow; ++i)
                target_c[i] = 0;

            return ret;
        };
        
        const Integer *ptr_c = static_cast<Integer*>(A->p),
                      *ptr_r = static_cast<Integer*>(A->i);
        const Real    *ptr_x = static_cast<Real*   >(A->x);

        assert(A->sorted || "CHOLMOD returned CCS matrix with unsorted row indices. matcl requires that.");

        for(size_t ind = 0; ind < (size_t)nzmax; ++ind)
        {
            target_x[ind] = ptr_x[ind];
            target_r[ind] = ptr_r[ind];
        }

        for(size_t ind = 0; ind <= (size_t)N_suc; ++ind)
            target_c[ind] = ptr_c[ind];
        for(Integer ind = N_suc + 1; ind <= nrow; ++ind)
            target_c[ind] = target_c[N_suc];

        //returned matrix should be unsymmetric (stype = 0)
        if(     A->stype > 0) ret = trans(tril(ret,1)) + ret;
        else if(A->stype < 0) ret = trans(triu(ret,1)) + ret;
    }

    return ret;
}

// Convert matcl raw sparse complex matrix to `cholmod_sparse' format.
//   NOTE #1: remember to clean up the returned cholmod_sparse object! (by call to `cholmod_free_sparse').
//      It is not refcounted! 
// NOTE #2: it will work only on hermitian matrices.
static cholmod_sparse* matcl_sparse_hercomplex_to_cholmod(const raw::Matrix<Complex,struct_sparse> &M,
                                                          cholmod_common& cc)
{
    int stype = -1;

    auto sprep      = M.rep();
    Integer nrow    = sprep.rows();
    Integer ncol    = sprep.cols();
    Integer nnz     = sprep.nnz();

    cholmod_sparse* ret = cholmod_l_allocate_sparse(nrow, ncol, nnz,
                            1,  // TRUE if columns of A sorted, FALSE otherwise
                            1,   // TRUE if A will be packed
                           stype, CHOLMOD_COMPLEX, &cc);

    const Integer *ptr_c = sprep.ptr_c(),
                  *ptr_r = sprep.ptr_r();
    const Complex *ptr_x = sprep.ptr_x();

    Integer *target_c = static_cast<Integer*>(ret->p),
            *target_r = static_cast<Integer*>(ret->i);
    Complex *target_x = static_cast<Complex*>(ret->x);

    if(ret != NULL)
    {
        // Matrix is assumed to be hermitian, copy to cholmod just the tril part.
        Integer pos = 0;
        target_c[0] = 0;

        for (Integer c = 0; c < ncol; ++c)
        {
            Integer k       = ptr_c[c];
            Integer last_k  = ptr_c[c+1];

            //ignore upper-triangular part
            for (; k < last_k; ++k)
            {
                Integer r   = ptr_r[k];

                if (r < c)
                    continue;

                break;
            };

            for (; k < last_k; ++k)
            {
                target_x[pos] = ptr_x[k];
                target_r[pos] = ptr_r[k];
                ++pos;
            };

            target_c[c+1]   = pos;
        };
    }

    return ret;
}

static matcl::Matrix cholmod_to_matcl_sparse_complex(const cholmod_sparse* A, Integer N_suc)
{
    matcl::Matrix ret;

    if(A != NULL)
    {
        Integer nrow  = (int)A->nrow;        
        Integer nzmax = (int)A->nzmax;

        ret     = matcl::make_complex_sparse_noinit(nrow, nrow, nzmax);

        using rep_type  = sparse_matrix<Complex,false>;
        rep_type sprepu(ret);

        Integer *target_c = sprepu.ptr_c(),
                *target_r = sprepu.ptr_r();
        Complex *target_x = sprepu.ptr_x();
        nzmax             = sprepu.nzmax();

        if (N_suc == 0)
        {
            for (Integer i = 0; i <= nrow; ++i)
                target_c[i] = 0;

            return ret;
        }

        const Integer *ptr_c = static_cast<Integer*>(A->p),
                      *ptr_r = static_cast<Integer*>(A->i);
        const Real    *ptr_x = static_cast<Real*   >(A->x);

        assert(A->sorted || "CHOLMOD returned CCS matrix with unsorted row indices. matcl requires that.");

        for(size_t ind = 0; ind < (size_t)nzmax; ++ind)
        {
            std::complex<Real> tempc(ptr_x[2*ind], ptr_x[2*ind+1]);
            target_x[ind] = Complex(tempc);
            target_r[ind] = ptr_r[ind];
        }
        for(size_t ind = 0; ind <= (size_t)N_suc; ++ind)
            target_c[ind] = ptr_c[ind];
        for(Integer ind = N_suc + 1; ind <= nrow; ++ind)
            target_c[ind] = target_c[N_suc];

        //returned matrix should be unsymmetric (stype = 0)
        if(     A->stype > 0) ret = ctrans(tril(ret,1)) + ret;
        else if(A->stype < 0) ret = ctrans(triu(ret,1)) + ret;
    }

    return ret;
}

static permvec read_permutation(const cholmod_factor* F)
{
    Integer n = (int)F->n;    

    if(F != nullptr)
    {
        Matrix ret      = izeros(1, n);
        using rep_type  = dense_matrix<Integer,false>;

        rep_type dru(ret);

        Integer* target = dru.ptr();
        for(size_t ind = 0; ind < (size_t)n; ++ind)
        {
            // matcl indexing starts with 1.
            target[ind] = 1 + static_cast<const Integer*>(F->Perm)[ind];
        }

        return details::pv_constructor::make(ret);
    }
    else
    {
        return permvec::identity(n);
    };
}

static void cholmod_error_handler(int status, const char *file, int line, const char *message)
{
    (void)line;
    (void)file;

    if(status == 1) // nonposdef status code
        return;
    //if(status == 1) // nonposdef status code
    //    throw error::error_nonposdef();

    if(status > 0)
        error::get_global_messanger_linalg()->warning_cholmod(message);
    else if(status < 0)
        throw error::error_cholmod(message);
}

}

template<> 
struct chol_impl<Real,struct_sparse> 
{
    using Mat   = raw::Matrix<Real,struct_sparse>;

    static void eval(chol2_return_type& ret, const Mat& A, bool upper, const options& opts)
    {
        if(A.rows() != A.cols())
            throw error::square_matrix_required(A.rows(), A.cols());

        Integer n = A.rows();

		if (A.get_struct().is_id())
        {
            ret = chol2_return_type(Matrix(A,false), permvec::identity(n), n);
            return;
        };

		if (is_diag(Matrix(A,false)))
			return chol_diag<Real,struct_sparse>::eval(ret, A);

        cholmod_common c;
        cholmod_l_start(&c);

        // Verbosity level of cholmod. Uncomment for tests.
        //c.print = 5; 
        
        c.itype         = CHOLMOD_LONG;
        c.final_ll      = true; // Leave the result in LL' form
        c.error_handler = &cholmod_helpers::cholmod_error_handler;

        namespace opt   = opt::linsolve;
        opt::chol_ordering_type ord  = (opt::chol_ordering_type)opts.get_option<Integer>(opt::chol_ordering());

        switch (ord)
        {
            case opt::chol_ordering_type::default_val:
            default:
                break;
            case opt::chol_ordering_type::natural:
                c.method[0].ordering    = CHOLMOD_NATURAL;
                c.nmethods              = 1;
                break;
            /*
            metis is broken in cholmod
            case opt::ordering_type::metis:
                c.method[0].ordering    = CHOLMOD_METIS;
                c.postorder   = true;
                c.nmethods              = 1;
                break;
            */
            case opt::chol_ordering_type::amd:
                c.method[0].ordering    = CHOLMOD_AMD;
                c.postorder             = true;
                c.nmethods              = 1;
                break;                
            case opt::chol_ordering_type::cholmod_nested:
                c.method[0].ordering    = CHOLMOD_NESDIS;
                c.postorder             = true;
                c.nmethods              = 1;
                break;        
        };

        cholmod_sparse *cholmA  = cholmod_helpers::matcl_sparse_symreal_to_cholmod(A, c);        
        cholmod_factor *fac     = cholmod_l_allocate_factor(n, &c);       
        cholmod_sparse *L;

        if (cholmA == nullptr)
        {
            cholmod_l_free_sparse(&cholmA, &c);
            cholmod_l_finish(&c);
            throw error::error_general("CHOLMOD: unable to allocate matrix.");
        }

        // Analyze and permute.
        fac = cholmod_l_analyze(cholmA, &c);

        // Factorize.
        cholmod_l_factorize(cholmA, fac, &c);

        Integer N_suc       = (int)fac->minor;
        L                   = cholmod_l_factor_to_sparse(fac, &c);

        Matrix ret_factor   = cholmod_helpers::cholmod_to_matcl_sparse_real(L, N_suc);
        permvec p           = cholmod_helpers::read_permutation(fac);

        cholmod_l_free_factor(&fac, &c);
        cholmod_l_free_sparse(&cholmA, &c);
        cholmod_l_free_sparse(&L, &c);
        cholmod_l_finish(&c);

        // matcl expects U'U factorization rather than LL'. We need to transpose ret_factor.
        if (upper)
        {
            ret_factor = trans(ret_factor);
            ret_factor.add_struct(predefined_struct_type::triu);
        }
        else
        {
            ret_factor.add_struct(predefined_struct_type::tril);
        };

        ret = chol2_return_type(ret_factor, p, N_suc);
        return;
    };
};

template<> 
struct chol_impl<Float,struct_sparse> 
{
    using Mat   = raw::Matrix<Float,struct_sparse>;
    using Mat_R = raw::Matrix<Real,struct_sparse>;

    static void eval(chol2_return_type& ret, const Mat& A, bool upper, const options& opts)
    {
        //cholmod is not implemented for single precision
        Mat_R tmp   = mr::converter<Mat_R,Mat>::eval(A);
        chol_impl<Real,struct_sparse>::eval(ret, tmp, upper, opts);

        // convert back to single precision
        auto scope  = error::enable_warnings(false);
        Matrix res  = matcl::convert_value(ret.get<1>(), value_code::v_float);
        ret         = chol2_return_type(res, ret.get<2>(), ret.get<3>());
    };
};

template<> 
struct chol_impl<Complex,struct_sparse> 
{
    using Mat = raw::Matrix<Complex,struct_sparse>;

    static void eval(chol2_return_type& ret, const Mat& A, bool upper, const options& opts)
    {
        if(A.rows() != A.cols())
            throw error::square_matrix_required(A.rows(), A.cols());

        Integer n = A.rows();

		if (A.get_struct().is_id())
        {
            ret = chol2_return_type(Matrix(A,false), permvec::identity(n),n);
            return;
        }

		if (is_diag(Matrix(A,false)))
			return chol_diag<Complex,struct_sparse>::eval(ret, A);

        cholmod_common c;
        cholmod_l_start(&c);

        // Verbosity level of cholmod. Uncomment me for tests.
        //c.print = 5;

        c.itype         = CHOLMOD_LONG;
        c.final_ll      = true; // Leave the result in LL' form
        c.error_handler = &cholmod_helpers::cholmod_error_handler;

        namespace opt   = opt::linsolve;
        opt::chol_ordering_type ord  = (opt::chol_ordering_type)opts.get_option<Integer>(opt::chol_ordering());

        switch (ord)
        {
            case opt::chol_ordering_type::default_val:
                break;
            case opt::chol_ordering_type::natural:
                c.method[0].ordering    = CHOLMOD_NATURAL;
                c.nmethods              = 1;
                break;
            /*
            ///metis is broken in cholmod
            case opt::ordering_type::metis:
                c.method[0].ordering    = CHOLMOD_METIS;
                c.postorder             = true;
                c.nmethods              = 1;
                break;
            */
            case opt::chol_ordering_type::amd:
                c.method[0].ordering    = CHOLMOD_AMD;
                c.postorder             = true;
                c.nmethods              = 1;
                break;                
            case opt::chol_ordering_type::cholmod_nested:
                c.method[0].ordering    = CHOLMOD_NESDIS;
                c.postorder             = true;
                c.nmethods              = 1;
                break;    
            default:
                break;
        };

        cholmod_sparse *cholmA  = cholmod_helpers::matcl_sparse_hercomplex_to_cholmod(A, c);        
        cholmod_factor *fac     = cholmod_l_allocate_factor(n, &c);
        cholmod_sparse *L;

        if (cholmA == NULL)
        {
            cholmod_l_free_sparse(&cholmA, &c);
            cholmod_l_finish(&c);
            throw error::error_general("CHOLMOD: unable to allocate matrix.");
        }

        // Analyze and permute.
        fac = cholmod_l_analyze(cholmA, &c);

        // Factorize.
        cholmod_l_factorize(cholmA, fac, &c);

        Integer N_suc       = (int)fac->minor;
        L = cholmod_l_factor_to_sparse(fac, &c);

        Matrix ret_factor   = cholmod_helpers::cholmod_to_matcl_sparse_complex(L,N_suc);
        permvec p           = cholmod_helpers::read_permutation(fac);

        cholmod_l_free_factor(&fac, &c);
        cholmod_l_free_sparse(&cholmA, &c);
        cholmod_l_free_sparse(&L, &c);
        cholmod_l_finish(&c);

        // matcl expects U'U factorization rather than LL'. We need to hermite conjugate ret_factor.
        if (upper)
        {
            ret_factor = ctrans(ret_factor);
            ret_factor.add_struct(predefined_struct_type::triu);
        }
        else
        {
            ret_factor.add_struct(predefined_struct_type::tril);
        }

        ret = chol2_return_type(ret_factor, p, N_suc);
        return;
    };
};

template<> 
struct chol_impl<Float_complex,struct_sparse> 
{
    using Mat       = raw::Matrix<Float_complex,struct_sparse>;
    using Mat_R     = raw::Matrix<Complex,struct_sparse>;

    static void eval(chol2_return_type& ret, const Mat& A, bool upper, const options& opts)
    {
        //cholmod is not implemented for single precision
        Mat_R tmp   = mr::converter<Mat_R,Mat>::eval(A);
        chol_impl<Complex,struct_sparse>::eval(ret, tmp, upper, opts);

        // convert back to single precision
        auto scope  = error::enable_warnings(false);
        Matrix res  = matcl::convert_value(ret.get<1>(), value_code::v_float_complex);
        ret         = chol2_return_type(res, ret.get<2>(), ret.get<3>());
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------

template<class V, class S>
struct chol
{
    using M = raw::Matrix<V,S>;
    static void eval(chol2_return_type& ret, const M& A, bool upper, const options& opts)
    {
        return chol_impl<V,S>::eval(ret, A, upper, opts);
    };
};

template<class S>
struct chol<Integer,S>
{
    using M = raw::Matrix<Integer,S>;
    static void eval(chol2_return_type& ret, const M& A, bool upper, const options& opts)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return chol<Real,S>::eval(ret, AC, upper, opts);
    };
};

};};
