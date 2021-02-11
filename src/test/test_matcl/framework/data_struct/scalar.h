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

#include <stdexcept>

namespace matcl { namespace test
{

class Scalar
{
	private:
		Integer				m_int;
		Real				m_real;
        Float				m_float;
		Complex				m_complex;
        Float_complex       m_fcomplex;
		matcl::value_code	m_type;

	public:
		Scalar()                : m_int(0),		m_type(value_code::v_integer)	{};
		Scalar(Integer i)       : m_int(i),		m_type(value_code::v_integer)	{};
		Scalar(Real r)          : m_real(r),	m_type(value_code::v_real)		{};
        Scalar(Float r)         : m_float(r),	m_type(value_code::v_float)		{};
		Scalar(Complex c)       : m_complex(c), m_type(value_code::v_complex)	{};
        Scalar(Float_complex c) : m_fcomplex(c), m_type(value_code::v_float_complex){};

		matcl::value_code	get_value_code() const	{ return m_type; };
		Integer				get_int() const		    { check(value_code::v_integer); return m_int; };
		Real				get_real() const	    { check(value_code::v_real); return m_real; };
        Float				get_float() const	    { check(value_code::v_float); return m_float; };
		Complex				get_complex() const	    { check(value_code::v_complex); return m_complex; };
        Float_complex		get_fcomplex() const    { check(value_code::v_float_complex); return m_fcomplex; };

	private:

		void check(matcl::value_code t) const
		{
			if (t != m_type)
				throw std::runtime_error("invalid scalar type");
		};
};

};};