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

#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/IO/archive.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl
{

template<class T>
void complex<T>::serialize(oarchive_impl & ar, const unsigned int version)
{
    (void)version;
    T re = real(value);
    T im = imag(value);
    ar << re;
    ar << im;
};

template<class T>
void complex<T>::serialize(iarchive_impl & ar, const unsigned int version)
{
    (void)version;

    T re;
    T im;
    ar >> re;
    ar >> im;

    value = std::complex<T>(re,im);
};

template struct complex<Float>;
template struct complex<Real>;

template bool operator==(const complex<Float>& arg, const complex<Float>& val);
template bool operator==(const complex<Real>& arg, const complex<Real>& val);

template bool operator!=(const complex<Float>& arg, const complex<Float>& val);
template bool operator!=(const complex<Real>& arg, const complex<Real>& val);

template bool operator>=(const complex<Float>& arg, const complex<Float>& val);
template bool operator>=(const complex<Real>& arg, const complex<Real>& val);

template bool operator<=(const complex<Float>& arg, const complex<Float>& val);
template bool operator<=(const complex<Real>& arg, const complex<Real>& val);

template bool operator>(const complex<Float>& arg, const complex<Float>& val);
template bool operator>(const complex<Real>& arg, const complex<Real>& val);

template bool operator<(const complex<Float>& arg, const complex<Float>& val);
template bool operator<(const complex<Real>& arg, const complex<Real>& val);

};
