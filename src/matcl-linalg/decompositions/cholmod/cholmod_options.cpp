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
#include <stdexcept>
#include "matcl-linalg/decompositions/cholmod/cholmod_options.h"

namespace matcl { namespace details
{

cholmod_options::cholmod_options()
{
	m_tol						= -1;		
	m_mu						= -1;
	m_corr_alg					= E_SE;
	m_corr_type					= TYPE_II;
	m_nondecreasing_strategy	= true;	
	m_piv_type					= PIV_GERSHGORIN;
};

void cholmod_options::mu(Real val)
{
	m_mu						= val;
};
void cholmod_options::check() const
{
	if (m_piv_type==PIV_GERSHGORIN && m_corr_alg !=E_SE)
		throw std::runtime_error("Gerschgoring pivoting strategy is allowed only for SE correction algorithm");
};

};};