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
#pragma once

#include "matcl-linalg/decompositions/cholmod/cholmod_options.h"
#include "matcl-linalg/decompositions/cholmod/cholmod_matrix.h"
#include "matcl-blas-lapack/blas/details/blas_utils.h"

#include <boost/shared_array.hpp>
#include <memory>

#pragma warning(push)
#pragma warning(disable:4251)       //needs to have dll-interface to be used by clients of class

namespace matcl { namespace details
{

namespace md = matcl::details;

template<class T>
class cholmod
{
    public:
        using int_array = boost::shared_array<Integer>;

	private:
		using real_type = typename md::real_type<T>::type;
		using self_type = cholmod;        

	private:
        void			make_factorization(cholmod_matrix<T>& mat,const cholmod_options& opt);
		void			clear();
		void			init();	

		int_array		m_piv;
		Integer		    m_rank;
		real_type		m_norm_E;

	public:
        // perform decomposition of the matrix mat using options opt;
        // on return mat represent the matrix U
        cholmod(cholmod_matrix<T>& mat, const cholmod_options& opt);

        // standard destructor
		~cholmod();

        // estimated rank of the matrix
		Integer         rank() const				{ return m_rank; };

        // return permutation vector
		int_array		permutation() const			{ return m_piv;	};

        // return upper bound of norm 2 of correction matrix E
		real_type		norm_correction() const		{ return m_norm_E; };
};

};};

#pragma warning(pop)