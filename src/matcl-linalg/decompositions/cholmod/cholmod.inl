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
#include <limits>
#include <exception>
#include "pstrf.inl"
#include "matcl-linalg/decompositions/cholmod/cholmod.h"

namespace matcl { namespace details
{

template<class T>
cholmod<T>::cholmod(cholmod_matrix<T>& mat,const cholmod_options& opt)
{
	init();

	try
	{
		make_factorization(mat,opt);
	}
	catch(...)
	{
		clear();
		throw;
	};
};
template<class T>
inline cholmod<T>::~cholmod()
{};

template<class T>
inline void cholmod<T>::init()
{	
	m_piv					= int_array();
	m_norm_E				= 0;
	m_rank					= 0;
};

template<class T>
inline void cholmod<T>::clear()
{
	init();
};

template<class T>
void cholmod<T>::make_factorization(cholmod_matrix<T>& mat,const cholmod_options& opt)
{
    Integer n				= mat.rows();

	if (n != mat.cols())
		throw std::runtime_error("symmetrix matrix is required");

    if (n == 0)
    {
        m_rank      = 0;
        m_norm_E    = 0;
        return;
    };

    using real_array = boost::shared_array<real_type>;

	real_array m_work		= real_array(new real_type[(2+sizeof(T)/sizeof(real_type))*n+1]);
	m_piv					= int_array(new Integer[2*n+1]);
	
    pstrf_output<T> out;
    out.m_piv				= m_piv;

    pstrf<T> fact(mat, m_work.get(),opt,out);	
    fact.make();

    m_rank					= out.step_phase_I;
    m_norm_E				= out.norm_E;	

    matcl::lapack::int2perm(n,out.m_piv.get(), out.m_piv.get()+n, true);
    return;
};

};};