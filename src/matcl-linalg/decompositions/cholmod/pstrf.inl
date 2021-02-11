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
#include "pstrf.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matrep/lib_functions/func_binary.h"

namespace matcl { namespace details
{

//===========================================================================================
//                      LOCAL FUNCTIONS
//===========================================================================================

template<class T> inline bool is_gt_c(T a,T b)    
{ 
    using TR = typename details::real_type<T>::type;
	TR tmp = abs2(a) - abs2(b);

    if (tmp > 0)
		return true;

	if (tmp < 0)
		return false;

	return imag(a) > imag(b); 
};

template<class T> inline bool is_gt(T a,T b)
{ 
    return a > b; 
};
template<> inline bool is_gt(Float_complex a,Float_complex b)    
{ 
    return is_gt_c(a,b); 
};
template<> inline bool is_gt(Complex a,Complex b)    
{ 
    return is_gt_c(a,b); 
};

template<class T>	
Integer maxloc(const T* mat, Integer size)
{
	Integer pos_max				    = 0;

    if (size == 0)
		return pos_max;

	T val_max						= mat[0];	
	
    for (Integer i = 1;i<size;i++)
	{
		if (is_gt(mat[i],val_max))
		{
			val_max					= mat[i];
			pos_max					= i;
		};
	};
	return pos_max;
};

template<class T>	
Integer minloc(const T* mat, Integer size)
{
	Integer pos_min					= 0;

	if (size == 0)
		return pos_min;

	T val_min						= mat[0];	

    for (Integer i = 1;i<size;i++)
	{
		if (mat[i] < val_min)
		{
			val_min					= mat[i];
			pos_min					= i;
		};
	};
	
    return pos_min;
};

template<class T>	
Integer maxabsloc(const T* mat, Integer size)
{
	Integer pos_max				    = 0;

    if (size == 0)
		return pos_max;

	T val_max						= abs(mat[0]);	
	
    for (Integer i = 1;i<size;i++)
	{
		if (abs(mat[i])>val_max)
		{
			val_max					= abs(mat[i]);
			pos_max					= i;
		};
	};
	
    return pos_max;
};

// i-th gershgorin circle of a symmetric matrix
template<class T> 
typename details::real_type<T>::type gershgorin_circle(const cholmod_matrix<T>& A,Integer j)
{
	if (A.rows() != A.cols())
		throw std::runtime_error("symmetrix matrix is required");

	Integer N			= A.rows();

	if (N == 0)
		return 0;

	using real_type = typename details::real_type<T>::type;

	real_type val	= 0;
    if (A.triang_type() == triang_part::UPPER)
	{
		for (Integer i = j+1; i < N; i++)
			val			= val + abs(A(j,i));

		for (Integer i = 0; i < j; i++)
			val			= val + abs(A(i,j));
	}
	else
	{
		for (Integer i = j+1; i < N; i++)
			val			= val + abs(A(i,j));

		for (Integer i = 0; i < j; i++)
			val			= val + abs(A(j,i));
	};

	return val;
};

inline void conj(Integer , Real* , Integer ) {};
inline void conj(Integer , Float* , Integer ) {};
inline void conj(Integer n_elem, Float_complex* A, Integer lda)
{
	for(Integer i=0;i<n_elem;++i)
	{
        A[0]    = conj(A[0]);
		A	    +=lda;
	};
};
inline void conj(Integer n_elem, Complex* A, Integer lda)
{
	for(Integer i=0;i<n_elem;++i)
	{
		A[0]    = conj(A[0]);
		A		+=lda;
	};
};

inline void swap_conj(Integer n_elem, Real* A, Integer lda, Real* B, Integer ldb)
{
    lapack::swap(n_elem, A, lda, B, ldb);
};
inline void swap_conj(Integer n_elem, Float* A, Integer lda, Float* B, Integer ldb)
{
	lapack::swap(n_elem, A, lda, B, ldb);
};
inline void swap_conj(Integer n_elem, Float_complex* A, Integer lda, Float_complex* B, Integer ldb)
{
    using VL = lapack_value_type<Float_complex>::type;

	lapack::swap<VL>(n_elem, lap(A), lda, lap(B), ldb);
	conj(n_elem, A,lda);
	conj(n_elem, B,ldb);
};
inline void swap_conj(Integer n_elem, Complex* A, Integer lda, Complex* B, Integer ldb)
{
    using VL = lapack_value_type<Complex>::type;

	lapack::swap<VL>(n_elem, lap(A), lda, lap(B), ldb);
	conj(n_elem,A,lda);
	conj(n_elem,B,ldb);
};

//===========================================================================================
//                      MAIN FUNCTIONS
//===========================================================================================
template<class T>
inline pstrf<T>::pstrf(cholmod_matrix<T>& A, void* workspace, const cholmod_options& opt, pstrf_output<T>& out)
	:A(A),workspace(workspace),opt(opt),out(out)
{
	init_local_vars();
	init_output();
	check_arg();
};

template<class T>
inline void pstrf<T>::init_output()
{
	out.step_phase_I				= 0;
	out.step_phase_II				= 0;
	out.exit_block_size				= n;
	out.tol_used					= (real_type)opt.tol();
	out.norm_E						= 0;
};
template<class T>
inline void pstrf<T>::init_local_vars()
{
    upper							= (A.triang_type() == triang_part::UPPER);
	n								= A.rows();
	lda								= A.ld();
};

template<class T>
inline void pstrf<T>::check_arg() const
{
	opt.check();
};
template<class T>
void pstrf<T>::make()
{
    using VL = typename lapack_value_type<T>::type;

    Integer* piv                     = out.m_piv.get();

    for (Integer i = 0; i < n; ++i) 
		piv[i]					    = i; 

	make(true,piv);

	Integer n_step_I                = out.step_phase_I;
	Integer M                       = A.rows() - n_step_I;

	if (!(opt.corr_alg() == opt.E_GMW || opt.corr_alg() == opt.E_SE) || n_step_I == A.rows())
		return;

	cholmod_matrix<T> A_sav			= A;
	A								= cholmod_matrix<T>(&A(n_step_I,n_step_I),M,M,A.ld(),A.triang_type());
	n								= M;
	piv							    = piv + n_step_I;

	make(false,piv);

	// permute remaining part of the matrix
	for ( Integer i =0; i<M; i++)
	{
		Integer tmp				    = piv[i];
		if (tmp == i)
			continue;

		if (upper)
			lapack::swap<VL>(n_step_I, lap(&A_sav(0,n_step_I+i)), 1, lap(&A_sav(0,n_step_I+tmp)),1);
		else
			lapack::swap<VL>(n_step_I, lap(&A_sav(n_step_I+i,0)), lda, lap(&A_sav(n_step_I+tmp,0)), lda);
	};

	for ( Integer i =0; i<M; i++)
		piv[i]					    = piv[i] + n_step_I;
};
template<class T>
void pstrf<T>::make(bool phase_I, Integer* piv)
{
    using VL = typename lapack_value_type<T>::type;

	if (n == 0)
		return;

    VL ONE                          = 1.;
	work							= static_cast<real_type*>(workspace);
	void* workspace_tmp				= work + 2*n;
	work_sav						= static_cast<T*>(workspace_tmp);
	bool b_exit						= false;

	real_type mu_used				= (real_type)opt.mu();
	if (mu_used < 0)
	{
		if (opt.corr_alg() == opt.E_SE )
			mu_used					= (real_type)0.1;
		else
			mu_used					= (real_type)0.5;
	};
	
	Integer nb						= get_block_size(upper == 1);

	if (nb < 1 || nb >= n )
		nb							= n;

	if (phase_I == true)
	{
		compute_stopping_value();
	}
	else
	{
		if (opt.corr_alg() == opt.E_GMW)
			init_correction_alg_GMW();
		else
			init_correction_alg_SE();

		if (opt.piv_type() == opt.PIV_GERSHGORIN)
			init_gershgorin_bounds();
	};		

	if (upper) 
	{
		//  Compute the Cholesky factorization P' * A * P = U' * U 
	    for (Integer k = 0; k < n; k += nb) 
		{
			// Account for last block not being NB wide 
			Integer jb				= lapack::minimum(nb,n - k);
			real_type diag_saved	= 0;

			for (Integer i = k; i < n; ++i) 
			    work[i]				= 0;

			for (j = k; j < k + jb; ++j) 
			{
				for (Integer i = j; i < n; ++i) 
				{
					if (j > k) 
					{
						T val		= A(j-1,i);
						work[i]		+= abs2(val);
					}
					work[n + i]		= real(A(i,i)) - work[i];
				};

				pvt					= find_pivot(phase_I);
				max_diag			= work[n + pvt];

				if (phase_I == true)
				{
					Integer tmp_min	= minloc(&work[n + j], n - j);
					min_diag		= work[n + j + tmp_min];
				};

				// first and second test on the pivot
				if (phase_I == true)
				{
					if (max_diag <= stop_value || min_diag < -mu_used*max_diag || isnan(max_diag)) 
					{
						b_exit		= true;
						break;
					};
				};		

				if (j != pvt ) 
				{
					//  Pivot OK, so can now swap pivot rows and columns
					diag_saved		= real(A(pvt,pvt));

					A(pvt,pvt)		= A(j,j);
					lapack::swap<VL>(j, lap(&A(0,j)), 1, lap(&A(0,pvt)), 1);
					if (pvt+1 < n) 
					{
						Integer  n_elements = n - pvt-1;
						lapack::swap<VL>(n_elements, lap(&A(j,pvt+1)), lda, lap(&A(pvt,pvt+1)), lda);
					};
					swap_conj(pvt - j-1, &A(j,j+1), lda, &A(j+1,pvt), 1);
					conj(1, &A(j,pvt), lda);

					//  Swap dot products
		
					real_type val	= work[j];
					work[j]			= work[pvt];
					work[pvt]		= val;

					// Swap gershgorin bounds
					if (phase_I == false && opt.piv_type() == opt.PIV_GERSHGORIN )
					{
						T val2			= work_sav[j];
						work_sav[j]		= work_sav[pvt];
						work_sav[pvt]	= val2;
					};
				};
				piv[j]			        = pvt;

				// update column J
				if (j+1 < n) 
				{			
					Integer tmp_1	= j - k;
					Integer tmp_2	= n - j - 1;

					conj(j,&A(0,j),1);

					if (phase_I == true)
						lapack::lacpy<VL>("All",1,tmp_2,lap(&A(j,j+1)),lda,lap(&work_sav[j+1]),1);
					
					lapack::gemv<VL>("Trans", tmp_1, tmp_2, -ONE, lap(&A(k,j+1)), lda, lap(&A(k,j)),
                                     1, ONE, lap(&A(j,j+1)), lda);

					conj(j,&A(0,j),1);
				};

				if (phase_I == false)
				{
					if (opt.corr_alg() == opt.E_GMW)
						correct_max_diag_GMW();
					else
						correct_max_diag_SE();

					if ( opt.piv_type() == opt.PIV_GERSHGORIN )
						update_gerschgorin_bounds();
				};

				// scale column J
				if (j+1 < n) 
				{
                    if (max_diag != real_type(0.))
                    {
					    real_type val	= real_type(1.) / sqrt(max_diag);
					    lapack::scal<VL>(n - j - 1, val, lap(&A(j,j+1)), lda);
                    }
				};

				if (phase_I == true)
				{
					// compute mimimal element on diagonal in next step
					min_diag				= max_diag_abs;
					for (Integer i = j+1; i < n; ++i) 
					{
						T& tmp				= A(j,i);
						real_type tmp2		= work[i] + abs2(tmp);
						tmp2				= real(A(i,i)) - tmp2;
						if (tmp2 < min_diag)
						{
							min_diag		= tmp2;
						};
					};

					// third test of the pivot
					if (min_diag < - mu_used * max_diag_abs)
					{
						if (j+1 < n) 
							lapack::lacpy<VL>("All",1,n-j-1,lap(&work_sav[j+1]),1,lap(&A(j,j+1)),lda);

						// undo the last permutation

						if (j == pvt) 
						{
							b_exit	= true;
							break;
						};
			
						A(pvt,pvt)			= diag_saved;
						lapack::swap<VL>(j, lap(&A(0,j)), 1, lap(&A(0,pvt)),1);

                        if (pvt+1 < n) 
							lapack::swap<VL>(n-pvt-1, lap(&A(j,pvt+1)), lda, lap(&A(pvt, pvt+1)), lda);

						swap_conj(pvt - j-1, &A(j,j+1), lda, &A(j+1,pvt), 1);
						conj(1, &A(j,pvt), lda);

						piv[j]	    = j;

						b_exit		= true;
						break;
					};
				};			

				max_diag			= sqrt(max_diag);
				A(j,j)				= max_diag;
			};

			if (j-k > 0 && j < n) 
				lapack::herk<VL>("Upper", "Conj Trans", n-j, j-k, -1.0, lap(&A(k,j)), 
                                    lda, 1.0, lap(&A(j,j)), lda);

			if (b_exit == true)
			    break;
		};
	}
	else
	{
		//  Compute the Cholesky factorization P' * A * P = U' * U 
	    for (Integer k = 0; k < n; k += nb) 
		{
			// Account for last block not being NB wide 
			Integer jb				= lapack::minimum(nb,n - k);
			real_type diag_saved	= 0;

			for (Integer i = k; i < n; ++i) 
			    work[i]				= 0;

			for (j = k; j < k + jb; ++j) 
			{
				for (Integer i = j; i < n; ++i) 
				{
					if (j > k) 
					{
						T val		= A(i,j-1);
						work[i]		+= abs2(val);
					}
					work[n + i]		= real(A(i,i)) - work[i];
				};

				pvt					= find_pivot(phase_I);
				max_diag			= work[n + pvt];

				if (phase_I == true)
				{
					Integer tmp_min	= minloc(&work[n + j], n - j);
					min_diag		= work[n + j + tmp_min];
				};

				// first and second test on the pivot
				if (phase_I == true)
				{
					if (max_diag <= stop_value || min_diag < -mu_used*max_diag || isnan(max_diag)) 
					{
						b_exit		= true;
						break;
					};
				};							

				if (j != pvt ) 
				{
					//  Pivot OK, so can now swap pivot rows and columns
					diag_saved		= real(A(pvt,pvt));

					A(pvt,pvt)		= A(j,j);
					lapack::swap<VL>(j, lap(&A(j,0)), lda, lap(&A(pvt,0)), lda);

                    if (pvt+1 < n) 
						lapack::swap<VL>(n-pvt-1, lap(&A(pvt+1,j)), 1, lap(&A(pvt+1,pvt)), 1);

					swap_conj(pvt - j-1, &A(j+1,j), 1, &A(pvt,j+1), lda);
					conj(1, &A(pvt,j), 1);

					//  Swap dot products
		
					real_type val	= work[j];
					work[j]			= work[pvt];
					work[pvt]		= val;

					// Swap gershgorin bounds
					if (phase_I == false && opt.piv_type() == opt.PIV_GERSHGORIN )
					{
						T val2			= work_sav[j];
						work_sav[j]		= work_sav[pvt];
						work_sav[pvt]	= val2;
					};
				};
				piv[j]			        = pvt;

				// update column J
				if (j+1 < n) 
				{			
					conj(j,&A(j,0),lda);

					if (phase_I == true)
						lapack::lacpy<VL>("All",n-j-1, 1, lap(&A(j+1,j)),1, lap(&work_sav[j+1]),1);
					
					lapack::gemv<VL>("No Trans", n-j-1, j - k, -ONE, lap(&A(j+1,k)), lda, lap(&A(j,k)),
                                     lda, ONE, lap(&A(j+1,j)), 1);

					conj(j,&A(j,0),lda);
				};

				if (phase_I == false)
				{
					if (opt.corr_alg() == opt.E_GMW)
						correct_max_diag_GMW();
					else
						correct_max_diag_SE();

					if ( opt.piv_type() == opt.PIV_GERSHGORIN )
						update_gerschgorin_bounds();
				};

				// scale column J
				if (j+1 < n) 
				{			
                    if (max_diag != real_type(0.0))
                    {
					    real_type val	= real_type(1.) / sqrt(max_diag);
					    lapack::scal<VL>(n-j-1, val, lap(&A(j+1,j)), 1);
                    }
				};

				if (phase_I == true)
				{
					// compute mimimal element on diagonal in next step
					min_diag		= max_diag_abs;
					for (Integer i = j+1; i < n; ++i) 
					{
						T& tmp				= A(i,j);
						real_type tmp2		= work[i] + abs2(tmp);
						tmp2				= real(A(i,i)) - tmp2;

						if (tmp2 < min_diag)
							min_diag		= tmp2;
					};

					// third test of the pivot
					if (min_diag < - mu_used * max_diag_abs)
					{
						if (j+1 < n) 
							lapack::lacpy<VL>("All",n-j-1,1,lap(&work_sav[j+1]),1,lap(&A(j+1,j)),1);

						// undo the last permutation

						if (j == pvt) 
						{
							b_exit	= true;
							break;
						};
			
						A(pvt,pvt)			= diag_saved;
						lapack::swap<VL>(j, lap(&A(j,0)), lda, lap(&A(pvt,0)),lda);

						if (pvt+1 < n) 
							lapack::swap<VL>(n-pvt-1, lap(&A(pvt+1,j)), 1, lap(&A(pvt+1, pvt)), 1);

						swap_conj(pvt-j-1, &A(j+1,j), 1, &A(pvt,j+1), lda);
						conj(1, &A(pvt,j), 1);

						piv[j]	    = j;

						b_exit		= true;
						break;
					};
				};

				max_diag			= sqrt(max_diag);
				A(j,j)				= max_diag;
			};

			if (j-k > 0 && j < n) 
				lapack::herk<VL>("Lower", "No Trans", n-j, j-k, -1.0, lap(&A(j,k)), 
                                 lda, 1.0, lap(&A(j,j)), lda);

			if (b_exit == true)
				break;
		};
	};

	if (phase_I == true)
	{
		out.step_phase_I			= j;
	}
	else		
	{
		out.step_phase_II			= j;
		out.exit_block_size			= n-j;
	};

    return;
};
template<class T>
void pstrf<T>::compute_stopping_value()
{
    using VL    = typename lapack_value_type<T>::type;
    using VR    = typename md::real_type<T>::type;

    for (Integer i = 0; i < n; ++i) 
		work[i]						= real(A(i,i));

	max_diag_abs					= abs(work[0]);

    for (Integer i = 1; i < n; ++i) 
	{
		real_type tmp				= work[i];

		if (abs(tmp) > max_diag_abs)
			max_diag_abs			= abs(tmp);
    };

    VR eps                          = lapack::lamch<VL>("Epsilon");

	m_corr_alg_data.tau				= pow(eps,(real_type)(1./3.0));
	m_corr_alg_data.tau_bar			= m_corr_alg_data.tau*m_corr_alg_data.tau;;

    if (opt.tol() > 0) 
	{
		stop_value					= (real_type)opt.tol();
		out.tol_used				= stop_value;
		return;
	};

	if (opt.corr_alg() == opt.E_SE)
		stop_value					= m_corr_alg_data.tau_bar * max_diag_abs;
	else if (opt.corr_alg() == opt.E_GMW && opt.corr_type() == opt.TYPE_I)
		stop_value					= n * eps * max_diag_abs;
	else if (opt.corr_alg() == opt.E_GMW && opt.corr_type() == opt.TYPE_II)
		stop_value					= m_corr_alg_data.tau_bar * max_diag_abs;
	else
		stop_value					= n * eps * max_diag_abs;

	out.tol_used					= stop_value;
};

template<class T>
Integer pstrf<T>::get_block_size(bool upper2)
{
    using VL = typename lapack_value_type<T>::type;

	const char * uplo				= (upper2)? "U" : "L";
	lapack::details::lapack_code c  = lapack::details::lapack_type<VL>::value;
	Integer N			            = n;
	Integer int_1		            = 1;
	Integer int_m1		            = -1;

	switch (c)
	{
        case lapack::details::LAPACK_C:
			return matcl::lapack::ilaenv(int_1, "CPOTRF", uplo, N, int_m1, int_m1, int_m1);
		case lapack::details::LAPACK_D:
			return matcl::lapack::ilaenv(int_1, "DPOTRF", uplo, N, int_m1, int_m1, int_m1);
		case lapack::details::LAPACK_S:
			return matcl::lapack::ilaenv(int_1, "SPOTRF", uplo, N, int_m1, int_m1, int_m1);
		case lapack::details::LAPACK_Z:
			return matcl::lapack::ilaenv(int_1, "ZPOTRF", uplo, N, int_m1, int_m1, int_m1);
		default:
		{
			//matcl_assert(0,"unknown lapack code");
			throw error::error_general("invalid case");
		}
	};	
};
template<class T>
void pstrf<T>::init_correction_alg_GMW()
{
    using VL = typename lapack_value_type<T>::type;

	// maximum magnitute of off diagonals
	real_type xi_hat		= 0;

	for (Integer i = 0; i<n; i++)
	{
		if (upper)
		{
			for (Integer k = i + 1; k<n; k++)
			{
				real_type tmp	= abs(A(i,k));

				if (tmp > xi_hat)
					xi_hat		= tmp;
			};
		}
		else
		{
			for (Integer k = i + 1; k<n; k++)
			{
				real_type tmp	= abs(A(k,i));
				if (tmp > xi_hat)
					xi_hat		= tmp;
			};
		};
	};

	if (n > 1)
	{
		if (opt.corr_type() == opt.TYPE_I)
			m_corr_alg_data.beta	= xi_hat/sqrt(static_cast<real_type>(n*n-1));
		else
			m_corr_alg_data.beta	= xi_hat/sqrt(static_cast<real_type>(n*n-n));
	}

	if (m_corr_alg_data.beta == 0)
		m_corr_alg_data.beta		= lapack::lamch<VL>("Epsilon");
};

template<class T>
void pstrf<T>::init_correction_alg_SE()
{};

template<class T>
void pstrf<T>::correct_max_diag_GMW()
{
	real_type new_diag				= 0;	
	real_type c_norm				= get_c_norm_inf();

	if (opt.corr_type() == opt.TYPE_I)
	{
		new_diag					= lapack::maximum<real_type>(stop_value,max_diag + m_corr_alg_data.delta);
		new_diag					= lapack::maximum<real_type>(new_diag,-max_diag);
		new_diag					= lapack::maximum<real_type>(new_diag,c_norm*c_norm/m_corr_alg_data.beta);
	}
	else
	{
		new_diag					= lapack::maximum<real_type>(stop_value,max_diag + m_corr_alg_data.delta);
		new_diag					= lapack::maximum<real_type>(new_diag,c_norm*c_norm/m_corr_alg_data.beta);					
	};

	real_type d_diag				= new_diag - max_diag;

	if (opt.nondecreasing_strategy() == true)
		m_corr_alg_data.delta		= d_diag;

	if (d_diag > out.norm_E)
		out.norm_E					= d_diag;

	if (out.norm_E == 0)
	{
		out.step_phase_I++;
		out.step_phase_II--;
	};
	
	max_diag						= new_diag;
};

template<class T>
typename pstrf<T>::real_type pstrf<T>::get_c_norm_inf()
{
    using VL = typename lapack_value_type<T>::type;

	if (upper)
		return lapack::lange<VL>("1", 1, n-j-1, lap(&A(j,j+1)), lda, lap(&work[n]));
	else
		return lapack::lange<VL>("I", n-j-1, 1, lap(&A(j+1,j)), 1, lap(&work[n]));
};

template<class T>
typename pstrf<T>::real_type pstrf<T>::get_c_norm_1()
{
    using VL = typename lapack_value_type<T>::type;

	if (upper)
		return lapack::lange<VL>("I", 1, n-j-1, lap(&A(j,j+1)), lda, lap(&work[n]));
	else
		return lapack::lange<VL>("1", n-j-1, 1, lap(&A(j+1,j)), 1, lap(&work[n]));
};

template<class T>
void pstrf<T>::correct_max_diag_SE()
{
	real_type d_diag				= get_delta_SE();	
	real_type new_diag				= max_diag + d_diag;

	if (opt.nondecreasing_strategy() == true)
		m_corr_alg_data.delta		= d_diag;

	if (d_diag > out.norm_E)
		out.norm_E					= d_diag;

	if (out.norm_E == 0)
	{
		out.step_phase_I++;
		out.step_phase_II--;
	};
	
	max_diag						= new_diag;
};

template<class T>
typename pstrf<T>::real_type pstrf<T>::get_delta_SE()
{
	if (n-j == 2)
		return get_delta_SE_block22();

	if (n-j == 1 && n>=2)
		return m_corr_alg_data.delta;

	real_type d_diag;	
	real_type c_norm;

	if (n-j > 1)
		c_norm						= get_c_norm_1();
	else
		c_norm						= -m_corr_alg_data.tau/(1-m_corr_alg_data.tau)* max_diag;

	if (opt.corr_type() == opt.TYPE_II)
	{
		d_diag						= -max_diag + lapack::maximum(c_norm,m_corr_alg_data.tau_bar*max_diag_abs);
	}
	else
	{
		d_diag						= lapack::maximum(c_norm,m_corr_alg_data.tau_bar*max_diag_abs);
		d_diag						= lapack::maximum<real_type>(-max_diag*2,-max_diag + d_diag);		
	};

    d_diag							= lapack::maximum(m_corr_alg_data.delta,d_diag);
	return d_diag;
};

template<class T>
typename pstrf<T>::real_type pstrf<T>::get_delta_SE_block22()
{
    using VL = typename lapack_value_type<T>::type;

	real_type d_diag				= m_corr_alg_data.delta;

	real_type	lam_1;
	real_type	lam_2;
	T			A_12;

    if (upper)
		A_12						= A(j,j+1);
	else
		A_12						= A(j+1,j);

    lapack::lae2<VL>(work[n + j], *lap(&A_12), work[n + j+1], lam_1, lam_2);

	if (lam_1 > lam_2)
	{
		real_type tmp				= lam_1;
		lam_1						= lam_2;
		lam_2						= tmp;
	};    

	real_type tmp					= (lam_2-lam_1)*m_corr_alg_data.tau/(1-m_corr_alg_data.tau);
	tmp								= -lam_1+lapack::maximum(m_corr_alg_data.tau_bar*max_diag_abs,tmp);

	if (opt.corr_type() == opt.TYPE_II)
	{			
		d_diag						= lapack::maximum(d_diag,tmp);
	}
	else
	{
		d_diag						= lapack::maximum(d_diag,tmp);
		d_diag						= lapack::maximum<real_type>(d_diag,-2*lam_1);
	};
	
	m_corr_alg_data.delta			= d_diag;
	return d_diag;
};
template<class T>
Integer pstrf<T>::find_pivot(bool phase_I)
{
	if (phase_I == true )
		return j + maxloc(&work[n + j], n - j);				

	switch (opt.piv_type())
	{
		case opt.PIV_NONE:
			return j;
		case opt.PIV_DIAG:
			return j + maxloc(&work[n + j], n - j);	
		case opt.PIV_ABS_DIAG:
			return j + maxabsloc(&work[n + j], n - j);	
		case opt.PIV_GERSHGORIN:
			return j + maxloc(&work_sav[j], n - j);	
		default:
		{
			//matcl_assert(0,"unknown pivoting strategy");
			throw error::error_general("invalid case");
		}
	};
};
template<class T>
void pstrf<T>::init_gershgorin_bounds()
{
	for (Integer i = 0; i < n; i++)
	{
		real_type circle			= gershgorin_circle(A,i);
		work_sav[i]					= A(i,i) - circle;
	};
};
template<class T>
void pstrf<T>::update_gerschgorin_bounds()
{
	real_type c_norm                = 1.0;
    
    if (max_diag != real_type(0.0))
        c_norm                      = get_c_norm_1()/max_diag;

	if (upper == true)
	{		
		for (Integer i = j+1; i < n; i++)
		{
			real_type tmp			= abs(A(j,i))*(1-c_norm);
			T circle				= work_sav[i] + tmp;
			work_sav[i]				= circle;
		};
	}
	else
	{
		for (Integer i = j+1; i < n; i++)
		{
			real_type tmp			= abs(A(i,j))*(1-c_norm);
			T circle				= work_sav[i] + tmp;
			work_sav[i]				= circle;
		};
	};
};

};};