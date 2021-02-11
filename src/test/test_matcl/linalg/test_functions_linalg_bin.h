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

#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "matcl-linalg/matcl_linalg.h"

#include <vector>
#include <string>

namespace matcl { namespace test
{

class linalg_bin_functions_list
{
    public:
        using matrix_pair   = matrix_set_bin::matrix_pair;
        using scalar_pair   = matrix_set_bin::scalar_pair;

    private:
        const matrix_set_bin&	m_tests;
        options					m_options;
        bool					matrices_with_nan;
        dynamic_mat_set&	    ms_dyn;

    public:
        linalg_bin_functions_list(const matrix_set_bin& t,dynamic_mat_set& dyn,bool matrices_with_nan) 
                                : m_tests(t),ms_dyn(dyn),matrices_with_nan(matrices_with_nan){};

        void			make(options opts, Integer thread_id);

        matrix_pair		get_matrix(int code) const;
        scalar_pair		get_scalar(int code) const;

    private:

        void        test_gschur();
        void        test_gen_sym_eigen();
        void        test_geigs();
        
    private:
        linalg_bin_functions_list(const linalg_bin_functions_list&) = delete;
        linalg_bin_functions_list& operator=(const linalg_bin_functions_list&) = delete;
};

class test_function_gschur : public bin_function
{
    public:

        virtual Real eval_mat(const Matrix& mat1, const Matrix& mat2, int code);
        virtual Real eval_scalar(const Scalar& s1, const Scalar& s2, int code);
};

class test_function_geigs : public bin_function
{
    public:

        virtual Real eval_mat(const Matrix& mat1, const Matrix& mat2, int code);
        virtual Real eval_scalar(const Scalar& s1, const Scalar& s2, int code);
};

class test_function_gen_sym_eigen : public bin_function
{
    public:

        virtual Real eval_mat(const Matrix& mat1, const Matrix& mat2, int code);
        virtual Real eval_scalar(const Scalar& s1, const Scalar& s2, int code);
};

};};
