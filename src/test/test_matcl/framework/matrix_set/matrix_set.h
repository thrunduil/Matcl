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

#include "test/test_matcl/framework/test_functions/unary_function.h"
#include "test/test_matcl/framework/test_functions/bin_function.h"
#include "test/test_matcl/framework/data_struct/options.h"
#include "test/test_matcl/framework/matrix_set/dynamic_mat_set.h"
#include "test/test_matcl/framework/matrix_set/matrix_rand.h"

#include <vector>
#include <map>

namespace matcl { namespace test
{

class matrix_set
{
    public:
        using int2          = std::pair<Integer,Integer>;
        using sparse_info   = std::vector<Real>;
        using band_info     = std::vector<int2>;

    protected:
        std::vector<Matrix>		            m_matrices;
        std::vector<Scalar>		            m_scalars;
        rand_matrix_ptr                     m_rand;
        sparse_info                         m_sparse_info;
        band_info                           m_band_info;

    private:
        matrix_set(const matrix_set&) = delete;
        matrix_set& operator=(const matrix_set&) = delete;

    protected:
        void		add_matrices(std::vector<Matrix>& mat, Integer m, Integer n);
        void		add_matrices_band(std::vector<Matrix>& mat, Integer m, Integer n, Integer ld, Integer ud);
        void		add_matrices_dense(std::vector<Matrix>& mat, Integer m, Integer n);
        void		add_matrices_sparse(std::vector<Matrix>& mat, Integer m, Integer n, Real d);

    public:
        matrix_set(rand_matrix_ptr rand, const sparse_info& si, const band_info& bi)    
            :m_rand(rand),m_sparse_info(si), m_band_info(bi)
        {} ;
        virtual ~matrix_set() {};

        virtual Real		make(unary_function* func, options opts) const;
        Integer             get_matrix_vector_size() const;
        Matrix		        get_matrix(int code) const;

    protected:		
        std::string		make_mat_name(const Matrix& mat, int code) const;
        std::string		make_scal_name(const Scalar& mat, int code) const;
        Real			make(unary_function* func, options opts, const std::vector<Matrix>&, const std::vector<Scalar>&) const;
        Real			make(const Matrix& mat,unary_function* func,bool show_partial_res, int code) const;
        Real			make(const Scalar& sc,unary_function* func,bool show_partial_res, int code) const;
};

class matrix_set_bin
{
    public:
        using matrix_pair   = std::pair<Matrix,Matrix>;
        using scalar_pair   = std::pair<Scalar,Scalar>;
        using int2          = std::pair<Integer,Integer>;
        using sparse_info   = std::vector<Real>;
        using band_info     = std::vector<int2>;

    protected:
        std::vector<matrix_pair>	m_matrices;
        std::vector<scalar_pair>	m_scalars;
        rand_matrix_ptr             m_rand;
        sparse_info                 m_sparse_info;
        band_info                   m_band_info;

    protected:
        Real			make(bin_function* func, options opts, const std::vector<matrix_pair>&,
                        const std::vector<scalar_pair>& ) const;
        Real			make(const matrix_pair& mat,bin_function* func,bool show_partial_res, int code) const;
        Real			make(const scalar_pair& mat,bin_function* func,bool show_partial_res, int code) const;
        std::string		make_mat_name(const matrix_pair& mat, int code) const;
        std::string		make_scal_name(const scalar_pair& mat, int code) const;

    private:
        matrix_set_bin(const matrix_set_bin&) = delete;
        matrix_set_bin& operator=(const matrix_set_bin&) = delete;

    protected:
        void			add_matrices(std::vector<matrix_pair>& mat, Integer m1, Integer n1,
                                    Integer m2, Integer n2);
        void			add_matrices_dense(std::vector<matrix_pair>& mat, Integer m1, Integer n1,
                                    Integer m2, Integer n2);
        void			add_matrices_sparse(std::vector<matrix_pair>& mat, Integer m1, Integer n1, 
                                    Integer m2, Integer n2, Real d);
        void			add_matrices_band(std::vector<matrix_pair>& mat, Integer m1, Integer n1, 
                                    Integer m2, Integer n2, Integer ld, Integer ud);

        void			add_matrices(std::vector<matrix_pair>& mat, Integer m, Integer n, Matrix nm);
        void			add_matrices_dense(std::vector<matrix_pair>& mat, Matrix nm, 
                                        Integer m, Integer n);
        void			add_matrices_sparse(std::vector<matrix_pair>& mat, Matrix nm, 
                                        Integer m, Integer n, Real d);
        void			add_matrices_band(std::vector<matrix_pair>& mat, Matrix nm, 
                                        Integer m, Integer n, Integer ld, Integer ud);
        void            add(std::vector<matrix_pair>& mat,const matrix_pair& mp);

    public:
        matrix_set_bin(rand_matrix_ptr rand, const sparse_info& si, const band_info& bi)    
            :m_rand(rand),m_sparse_info(si), m_band_info(bi)
        {} ;
        virtual ~matrix_set_bin() {};

        virtual Real			make(bin_function* func,options opts) const;
        virtual matrix_pair		get_matrix(int code) const;
        virtual scalar_pair		get_scalar(int code) const;
};

};};