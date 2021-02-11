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

#include "matcl-linalg/decompositions/cholmod/cholmod_options.h"
#include "matcl-linalg/decompositions/cholmod/cholmod_matrix.h"
#include <boost/shared_array.hpp>
#include "matcl-matrep/details/utils.h"

namespace matcl { namespace details
{

template<class T>
struct pstrf_output
{
	using real_type = typename details::real_type<T>::type;
    using int_array = boost::shared_array<Integer>;

	int_array						    m_piv;
	real_type							tol_used;
	Integer							    step_phase_I;
	Integer							    step_phase_II;
	Integer							    exit_block_size;
	real_type							norm_E;
};
template<class T>
struct corr_alg_data
{
	using real_type = typename details::real_type<T>::type;

	real_type							beta;
	real_type							delta;
	real_type							tau;
	real_type							tau_bar;

	corr_alg_data()
	{
		beta							= 0;
		delta							= 0;
		tau								= 0;
		tau_bar							= 0;
	};
};

template<class T>
class pstrf
{
	private:
		using real_type = typename details::real_type<T>::type;

		cholmod_matrix<T>				A;
		void*							workspace;
		const cholmod_options&			opt;
		pstrf_output<T>&				out;

		bool    						upper;
		Integer  						n;
		Integer						    lda;

		real_type *						work;
		T *								work_sav;

		real_type						max_diag;		
		real_type						min_diag;
		real_type						max_diag_abs;
		Integer						    pvt;
		real_type						stop_value;
		corr_alg_data<T>				m_corr_alg_data;
		Integer						    j;

		void							init_local_vars();
		void							init_output();
		void							check_arg() const;
		void							compute_stopping_value();
		Integer						    get_block_size(bool upper);
		void							correct_max_diag_GMW();
		void							correct_max_diag_SE();
		real_type						get_delta_SE();
		real_type						get_delta_SE_block22();
		void							init_correction_alg_GMW();
		void							init_correction_alg_SE();
		Integer						    find_pivot(bool phase_I);
		void							init_gershgorin_bounds();
		void							update_gerschgorin_bounds();
		real_type						get_c_norm_inf();
		real_type						get_c_norm_1();

		void							make(bool phase_I, Integer* piv);

	public:
		pstrf(cholmod_matrix<T>& A, void* workspace, const cholmod_options& opt, pstrf_output<T>& out);

		void make();

    private:
        pstrf& operator=(const pstrf&) = delete;
};

};};