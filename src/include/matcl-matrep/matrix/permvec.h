/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/matrix/matrix.h"

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl
{

// class representing permutations based on shared pointer
class MATCL_MATREP_EXPORT permvec
{
    private:
        using impl_type     = std::shared_ptr<class pv_impl>;

    public:
        // create uninitialized permutation
        permvec();

        // destructor; release memory if reference count drops to zero
        ~permvec();

        // length of permutation
        Integer             length() const;

        // fast check if permutation is identity; only informative if is_id() == true
        bool                is_id() const;
        
        // return inverse permutation
        permvec             invperm() const;

        // concatenate two permutations
        permvec             operator()(const permvec& p) const;
		
        // convert permutation to matrix; i.e. return matrix [1:length](p), where
        // p represent this permutation
        Matrix              to_matrix() const;

        // convert permutation to matrix representing elements interchanges as in 
        // Lapack, i.e. element k and q(k) are interchanged, where q is returned
        // matrix of interchanges
        Matrix              to_interchanges_matrix() const;

        // equivalent to to_matrix().get_array<Integer>();
        const Integer*      to_array() const;

        // create identity permutation of given length
        static permvec      identity(Integer length);

        // create permutation object from matrix p; then to_matrix() == p
        static permvec      from_matrix(const Matrix& p);
		
    private:
        permvec(const impl_type& impl);

        static permvec      from_matrix_nocheck(const Matrix& p);

    private:		
        impl_type			m_data;

        friend details::pv_constructor;      
        friend class pv_impl;
};

// create inverse permutation
inline permvec invperm(const permvec& p)
{
    return p.invperm();
};

};

#pragma warning(pop)