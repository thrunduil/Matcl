/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

namespace matcl { namespace details
{

template<class M1,class val_type>
struct MATCL_MATREP_EXPORT assign_functor
{
	static Matrix& eval(Matrix& A, M1& mat, const val_type& val,Integer i,Integer j);
	static Matrix& eval(Matrix& A,M1& mat, const val_type& val,Integer i);
	static Matrix& eval(Matrix& A, M1& mat,const val_type& val,const colon& c1);
	static Matrix& eval(Matrix& A, M1& mat,const val_type& val,const colon& c1,const colon& c2);
	static Matrix& eval_diag(Matrix& A,M1& mat, const val_type& val,Integer d);
};

template<class M1,class val_type>
struct MATCL_MATREP_EXPORT assign_functor_scal
{
	static Matrix& eval(Matrix& A,const val_type& val,Integer i,Integer j);
	static Matrix& eval(Matrix& A,const val_type& val,Integer i);
	static Matrix& eval(Matrix& A,const val_type& val,const colon& c1);
	static Matrix& eval(Matrix& A,const val_type& val,const colon& c1,const colon& c2);
};

};};