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

#include "matcl-dynamic/predefined_type_traits.h"
#include "matcl-dynamic/details/object.inl"
#include "matcl-core/IO/archive.h"
#include "matcl-core/details/complex_details.h"
#include "matcl-core/details/IO/printer.h"
#include "matcl-core/IO/scalar_io.h"

namespace matcl { namespace dynamic
{

void object_type_traits<bool>::disp(bool t, matcl::details::printer& pr, Integer elem_width, 
                                    align_type at, Integer)
{
	pr.disp_elem(elem_width, t? "true" : "false", at, 0);
};

void object_type_traits<bool>::load_data(iarchive_impl& ar, T& val, unsigned int)
{
    ar >> val;
}

void object_type_traits<bool>::save_data(oarchive_impl& ar, const T& val, unsigned int)
{
    ar << val;
}

void object_type_traits<Integer>::disp(const Integer& t, matcl::details::printer& pr, 
                                        Integer elem_width, align_type at, Integer)
{
	pr.disp_elem(elem_width, t, at, 0);
};

void object_type_traits<Integer>::load_data(iarchive_impl& ar, T& val, unsigned int)
{
    ar >> val;
}

void object_type_traits<Integer>::save_data(oarchive_impl& ar, const T& val, unsigned int)
{
    ar << val;
}

void object_type_traits<Real>::disp(const Real& t, matcl::details::printer& pr, Integer elem_width, 
                                           align_type at, Integer)
{
	pr.disp_elem(elem_width, t, at, 0);
};

void object_type_traits<Real>::load_data(iarchive_impl& ar, T& val, unsigned int)
{
    ar >> val;
}

void object_type_traits<Real>::save_data(oarchive_impl& ar, const T& val, unsigned int)
{
    ar << val;
}

void object_type_traits<Float>::disp(const Float& t, matcl::details::printer& pr, Integer elem_width, 
                                           align_type at, Integer)
{
	pr.disp_elem(elem_width, t, at, 0);
};

void object_type_traits<Float>::load_data(iarchive_impl& ar, T& val, unsigned int)
{
    ar >> val;
}

void object_type_traits<Float>::save_data(oarchive_impl& ar, const T& val, unsigned int)
{
    ar << val;
}

void object_type_traits<Complex>::disp(const Complex& t, matcl::details::printer& pr, Integer elem_width, 
                                              align_type at, Integer)
{
	pr.disp_elem(elem_width, t, at, 0);
};

void object_type_traits<Complex>::load_data(iarchive_impl& ar, T& val, unsigned int)
{
    Real re, im;
    ar >> re;
    ar >> im;

    val = Complex(re,im);
}

void object_type_traits<Complex>::save_data(oarchive_impl& ar, const T& val, unsigned int)
{
    Real re = real(val);
    Real im = imag(val);
    ar << re;
    ar << im;
}

bool object_type_traits<Complex>::read(std::istream& is, T& val)
{
    Real re, im;
    is >> re;
    is >> im;

    val = Complex(re,im);

    if (is.fail() || is.bad())
        return false;
    else
        return true;
}

void object_type_traits<Complex>::write(std::ostream& os, const T& t)
{
    os << real(t);
    os << ' ';
    os << imag(t);
}

void object_type_traits<Float_complex>::disp(const Float_complex& t, matcl::details::printer& pr, 
                                        Integer elem_width, align_type at, Integer)
{
	pr.disp_elem(elem_width, t, at, 0);
};

void object_type_traits<Float_complex>::load_data(iarchive_impl& ar, T& val, unsigned int)
{
    Float re, im;
    ar >> re;
    ar >> im;

    val = Float_complex(re,im);
}

void object_type_traits<Float_complex>::save_data(oarchive_impl& ar, const T& val, 
                                                  unsigned int)
{
    Float re = real(val);
    Float im = imag(val);
    ar << re;
    ar << im;
}

bool object_type_traits<Float_complex>::read(std::istream& is, T& val)
{
    Float re, im;
    is >> re;
    is >> im;

    val = Float_complex(re,im);

    if (is.fail() || is.bad())
        return false;
    else
        return true;
}

void object_type_traits<Float_complex>::write(std::ostream& os, const T& t)
{
    os << real(t);
    os << ' ';
    os << imag(t);
}

void object_type_traits<std::string>::disp(const T& t, matcl::details::printer& pr, 
                                Integer elem_width, align_type at, Integer)
{
    pr.disp_elem(elem_width, t, at, 0);
};

void object_type_traits<std::string>::load_data(iarchive_impl& ar, T& val, 
                                                unsigned int)
{
    ar >> val;
};

void object_type_traits<std::string>::save_data(oarchive_impl& ar, const T& val, 
                                                unsigned int)
{
    ar << val;
};

//--------------------------------------------------------------------------
//                      Unit
//--------------------------------------------------------------------------
void object_type_traits<unit_type>::disp(const T&, matcl::details::printer& pr, Integer elem_width, 
                                    align_type at, Integer)
{
    pr.disp_elem(elem_width, "unit", at, 0);
};

void object_type_traits<unit_type>::load_data(iarchive_impl&, T&, unsigned int)
{
    //nothing to load;
};

void object_type_traits<unit_type>::save_data(oarchive_impl&, const T&, unsigned int)
{
    //nothing to save
};

//--------------------------------------------------------------------------
//                      Null
//--------------------------------------------------------------------------
void object_type_traits<null_type>::disp(const T&, matcl::details::printer& pr, Integer elem_width, 
                                    align_type at, Integer)
{
    pr.disp_elem(elem_width, "unit", at, 0);
};

void object_type_traits<null_type>::load_data(iarchive_impl&, T&, unsigned int)
{
    //nothing to load;
};

void object_type_traits<null_type>::save_data(oarchive_impl&, const T&, unsigned int)
{
    //nothing to save
};
//--------------------------------------------------------------------------
//                      Any
//--------------------------------------------------------------------------
bool object_type_traits<any_type>::is_zero(const T& val)
{
    return val.is_zero();
};

//--------------------------------------------------------------------------
//                      Type
//--------------------------------------------------------------------------
void object_type_traits<Type>::disp(const T& t, matcl::details::printer& pr, Integer elem_width,
                                    align_type at, Integer)
{
    pr.disp_elem(elem_width, t.to_string(), at, 0);
};

void object_type_traits<Type>::load_data(iarchive_impl& ar, T& val, 
                                         unsigned int version)
{
    val.serialize(ar, version);
};
void object_type_traits<Type>::save_data(oarchive_impl& ar, const T& val, 
                                         unsigned int version)
{
    val.serialize(ar, version);
};

};};
