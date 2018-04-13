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

class io_functions_list
{
	public:
		using matrix_pair   = matrix_set_bin::matrix_pair;
		using scalar_pair   = matrix_set_bin::scalar_pair;

	private:
		const matrix_set_bin&	m_tests;
		options					m_options;

	public:
		io_functions_list(const matrix_set_bin& t) : m_tests(t){};

		void			make(options opts, Integer thread_id);

		matrix_pair		get_matrix(int code) const;
		scalar_pair		get_scalar(int code) const;

	public:
		void		test_io();
        void		test_io2();
		void		test_serialize();
        void        test_io_mm();

	private:
		io_functions_list(const io_functions_list&) = delete;
		io_functions_list& operator=(const io_functions_list&) = delete;
};

class test_function_serialize : public bin_function
{
    public:
        template<class T1,class T2>
        Real eval_scal_func(T1 a, T2 b)
        {
            Real dif = 0;
            std::ostringstream ss;
            oarchive ia(ss);

            save(ia,a);
            save(ia,b);

            T1 tmp1;
            T2 tmp2;
            std::istringstream ss2(ss.str());
            iarchive ia2(ss2);

            load(ia2,tmp1);
            load(ia2,tmp2);
            
            dif += norm_1(a - tmp1);
            dif += norm_1(b - tmp2);

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_io_formatted : public bin_function
{
    public:
        template<class T1,class T2>
        Real eval_scal_func(T1 a, T2 b)
        {
            Real dif = 0;
            T1 tmp1;
            T2 tmp2;
            
            {
                std::ostringstream ss;
                ss.precision(20);
                matcl::operator<<(ss, a);
                matcl::operator<<(ss, b);

                std::istringstream ss2(ss.str());
                matcl::operator>>(ss2, tmp1);
                matcl::operator>>(ss2, tmp2);

                check_struct(tmp1);
                check_struct(tmp2);
                
                dif += norm_1(a - tmp1);
                dif += norm_1(b - tmp2);
            };

            return dif;
        };

    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_io_formatted2 : public bin_function
{
    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

class test_function_io_mm : public bin_function
{
    public:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2,int code);
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2,int code);
};

};};
