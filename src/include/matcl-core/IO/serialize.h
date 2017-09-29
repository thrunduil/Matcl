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

#include "matcl-core/IO/archive.h"
#include "matcl-core/matrix/complex_type.h"

namespace matcl { namespace details
{

template<class Archive>
void serialize_save(Archive& ar, const struct_flag& sf, const unsigned int ver)
{
    sf.save(ar, ver);
};

template<class Archive>
void serialize_load(Archive& ar, struct_flag& sf, const unsigned int ver)
{
    sf.load(ar, ver);
};

template<class Archive>
void serialize_save(Archive& ar, Integer v, const unsigned int)
{
    ar << v;
};

template<class Archive>
void serialize_save(Archive& ar, Real v, const unsigned int)
{
    ar << v;
};

template<class Archive>
void serialize_save(Archive& ar, const Complex& v, const unsigned int)
{
    Real re = v.value.real(), im = v.value.imag();
    ar << re;
    ar << im;
};

template<class Archive>
void serialize_save(Archive& ar, const std::string& v, const unsigned int)
{
    ar << v;
}
template<class Archive>
void serialize_save(Archive& ar, const Object& v, const unsigned int)
{
    ar << v;
};

template<class Archive>
void serialize_load(Archive& ar, Integer& v, const unsigned int)
{
    ar >> v;
};

template<class Archive>
void serialize_load(Archive& ar, Real& v, const unsigned int)
{
    ar >> v;
};

template<class Archive>
void serialize_load(Archive& ar, Complex& v, const unsigned int)
{
    Real re, im;
    ar >> re;
    ar >> im;
    v = Complex(std::complex<Real>(re,im));
};

template<class Archive>
void serialize_load(Archive& ar, std::string& v, const unsigned int)
{
    ar >> v;
};

template<class Archive>
void serialize_load(Archive& ar, Object& v, const unsigned int)
{
    ar >> v;
};

template<class Archive>
void serialize_save(Archive& ar, Float v, const unsigned int)
{
    ar << v;
};

template<class Archive>
void serialize_save(Archive& ar, const Float_complex& v, const unsigned int)
{
    Float re = v.value.real(), im = v.value.imag();
    ar << re;
    ar << im;
};

template<class Archive>
void serialize_load(Archive& ar, Float& v, const unsigned int)
{
    ar >> v;
};

template<class Archive>
void serialize_load(Archive& ar, Float_complex& v, const unsigned int)
{
    Float re, im;
    ar >> re;
    ar >> im;
    v = Float_complex(std::complex<Float>(re,im));
};

};};

