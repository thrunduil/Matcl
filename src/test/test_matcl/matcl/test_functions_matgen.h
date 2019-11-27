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

#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

class matgen_functions_list
{
    public:
        void	make();

        void operator()()
        {
            make();
        };

    private:
        void    test_seq();
        void    test_zeros();
        void    test_ones();
        void    test_eye();
        void    test_rand();
        void    test_constructors();

    public:
        matgen_functions_list();
};

class test_function_construct
{
    public:
        test_function_construct()	{};

        Real make();

    private:
        Real    test_dense_ld(Integer r, Integer c, int vc);
        Real    test_other_dense_cons(Integer r, Integer c, int vc);
        Real    test_other_sparse_cons(Integer r, Integer c, Integer nz, int vc);
        Real    test_other_band_cons(Integer r, Integer c, Integer fd, Integer ld, int vc);
};

class test_function_seq
{
    public:
        test_function_seq()	{};

        Real make();

    private:
        Real    test_linspace(Real s, Real e, Integer n, value_code vc, Matrix& A);
        Real    test_logspace(Real s, Real e, Integer n, value_code vc, Matrix& A);
};

class test_function_zeros
{
    public:
        test_function_zeros(){};

        Real make();

    private:
        Real    test_zeros(Integer r, Integer c, value_code vc, Matrix& A);
        Real    test_spzeros(Integer r, Integer c, Integer nz, value_code vc, Matrix& A);
        Real    test_bzeros(Integer r, Integer c, Integer fd, Integer ld, value_code vc, Matrix& A);
};

class test_function_ones
{
    public:
        test_function_ones(){};

        Real make();

    private:
        Real    test_ones(Integer r, Integer c, value_code vc, Matrix& A);
        Real    test_spones(Integer r, Integer c, value_code vc, Matrix& A);
        Real    test_bones(Integer r, Integer c, value_code vc, Matrix& A);
};

class test_function_eye
{
    public:
        test_function_eye(){};

        Real make();

    private:
        Real    test_eye(Integer r, Integer c, value_code vc, Matrix& A);
        Real    test_speye(Integer r, Integer c, value_code vc, Matrix& A);
        Real    test_beye(Integer r, Integer c, Integer fd, Integer ld, value_code vc, Matrix& A);
};

class test_function_rand
{
    public:
        test_function_rand(){};

        Real make();

    private:
        Real    test_dense(Integer r, Integer c, value_code vc, int type, Matrix& A);
        Real    test_sparse(Integer r, Integer c, value_code vc, int type, Matrix& A);
        Real    test_band(Integer r, Integer c, Integer fd, Integer ld, value_code vc, int type, Matrix& A);
        Real    test_randperm(Integer n, Matrix& A);
};

};};
