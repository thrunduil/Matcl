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

#include "matcl-scalar/object.h"
#include "matcl-core/general/exception.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/object_type.h"
//#include "mmlib_basic/details/printer.h"
//#include "mmlib_basic/base/disp_stream_impl.h"
//#include "mmlib_basic/options/options_disp.h"
//#include "mmlib_basic/IO/matrix_io.h"

namespace matcl 
{

namespace md = matcl::details;

Integer matcl::cast_integer(const Object& v)
{
    dynamic::Type t = OInteger::get_static_type();
    return OInteger(cast(t, v), dynamic::from_object()).get();
};

Real matcl::cast_real(const Object& v)
{
    dynamic::Type t = OReal::get_static_type();
    return OReal(cast(t, v), dynamic::from_object()).get();
};

Float matcl::cast_float(const Object& v)
{
    dynamic::Type t = OFloat::get_static_type();
    return OFloat(cast(t, v), dynamic::from_object()).get();
};

Complex matcl::cast_complex(const Object& v)
{
    dynamic::Type t = OComplex::get_static_type();
    return OComplex(cast(t, v), dynamic::from_object()).get();
};

Float_complex matcl::cast_float_complex(const Object& v)
{
    dynamic::Type t = OFloat_complex::get_static_type();
    return OFloat_complex(cast(t, v), dynamic::from_object()).get();
};

Object matcl::convert_to_object(const dynamic::Type& ty, Integer v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, Float v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, Real v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, const Complex& v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, const Float_complex& v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, const Object& v)
{
    return Object(matcl::dynamic::convert(ty, v));
}

/*
struct printer_impl : public dynamic::printer
{
    details::printer* pr;

    printer_impl(details::printer* pr_) : pr(pr_){};

    virtual void disp_elem(Integer elem_width, const std::string& s, align_type at, 
                            Integer value_pos) override
    {
        if (pr)
            pr->disp_elem(elem_width, s, at, value_pos);
    }
    virtual void disp_elem(Integer elem_width, Integer r, align_type at, Integer value_pos)  override
    {
        if (pr)
            pr->disp_elem(elem_width, r, at, value_pos);
    }
    virtual void disp_elem(Integer elem_width, Real r, align_type at, Integer value_pos)  override
    {
        if (pr)
            pr->disp_elem(elem_width, r, at, value_pos);
    }
    virtual void disp_elem(Integer elem_width, Float r, align_type at,Integer value_pos)  override
    {
        if (pr)
            pr->disp_elem(elem_width, r, at, value_pos);
    }
    virtual void disp_elem(Integer elem_width, const Complex& r, align_type at,
                            Integer value_pos)  override
    {
        if (pr)
            pr->disp_elem(elem_width, r, at, value_pos);
    }
    virtual void disp_elem(Integer elem_width, const Float_complex& r, align_type at, 
                            Integer value_pos)  override
    {
        if (pr)
            pr->disp_elem(elem_width, r, at, value_pos);
    }

    virtual Integer get_precision() const override
    {
        if (pr)
            return pr->get_precision();
        else
            return matcl::options().get_option<Integer>(opt::disp::precision());
    };

    virtual std::string to_string(Integer v) override
    {
        return matcl::to_string(v);
    }
    virtual std::string to_string(Real v) override
    {
        return matcl::to_string(v);
    }
    virtual std::string to_string(Float v) override
    {
        return matcl::to_string(v);
    }
    virtual std::string to_string(const Complex& v) override
    {
        return matcl::to_string(v);
    }
    virtual std::string to_string(const Float_complex& v) override
    {
        return matcl::to_string(v);
    }
};

void matcl::disp_object(const Object& v, details::printer& pr, Integer elem_width, 
                 align_type at, Integer value_pos)
{
    printer_impl pi(&pr);
    v.disp(pi,elem_width,at,value_pos);
}; 

std::string matcl::to_string(const Object& v)
{
    printer_impl pi(nullptr);
    return v.to_string(pi);
};
*/

};
