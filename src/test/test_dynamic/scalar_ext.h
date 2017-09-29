/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-mp-obj/mp_object.h"

namespace matcl { namespace test
{

enum class value_ext_code
{
    v_int = 0,  v_float,    v_real,
    v_fcomplex, v_complex,  v_mpint,
    v_mpfloat,  v_mpcompl,  v_mprat = 8
};

class Scalar_ext
{
	private:
        using value_code    = value_ext_code;

		Integer				m_int;
		Real				m_real;
        Float				m_float;
		Complex				m_complex;
        Float_complex       m_fcomplex;
        mp_int              m_mpint;
        mp_float            m_mpfloat;
        mp_complex          m_mpcompl;
        mp_rational         m_mprat;

		value_code	        m_type;

	public:
		Scalar_ext()                    : m_int(0),		m_type(value_code::v_int)   	{};
		Scalar_ext(Integer i)           : m_int(i),		m_type(value_code::v_int)	    {};
		Scalar_ext(Real r)              : m_real(r),	m_type(value_code::v_real)		{};
        Scalar_ext(Float r)             : m_float(r),	m_type(value_code::v_float)		{};
		Scalar_ext(Complex c)           : m_complex(c), m_type(value_code::v_complex)	{};
        Scalar_ext(Float_complex c)     : m_fcomplex(c),m_type(value_code::v_fcomplex)  {};
        Scalar_ext(const mp_int& c)     : m_mpint(c),   m_type(value_code::v_mpint)     {};
        Scalar_ext(const mp_float& c)   : m_mpfloat(c), m_type(value_code::v_mpfloat)   {};
        Scalar_ext(const mp_complex& c) : m_mpcompl(c), m_type(value_code::v_mpcompl)   {};
        Scalar_ext(const mp_rational& c): m_mprat(c),   m_type(value_code::v_mprat)     {};

		value_code	        get_value_code() const	{ return m_type; };

		Integer				get_int() const		    { check(value_code::v_int); return m_int; };
		Real				get_real() const	    { check(value_code::v_real); return m_real; };
        Float				get_float() const	    { check(value_code::v_float); return m_float; };
		Complex				get_complex() const	    { check(value_code::v_complex); return m_complex; };
        Float_complex		get_fcomplex() const    { check(value_code::v_fcomplex); return m_fcomplex; };
        const mp_int&       get_mpint() const       { check(value_code::v_mpint); return m_mpint; };
        const mp_float&     get_mpfloat() const     { check(value_code::v_mpfloat); return m_mpfloat; };
        const mp_complex&   get_mpcomplex() const   { check(value_code::v_mpcompl); return m_mpcompl; };
        const mp_rational&  get_mprational() const  { check(value_code::v_mprat); return m_mprat; };

	private:

		void check(value_code t) const
		{
			if (t != m_type)
				throw std::runtime_error("invalid scalar type");
		};
};

};};