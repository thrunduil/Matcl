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

#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/details/utils.h"

namespace matcl
{

namespace md = matcl :: details;

matcl::value_code matrix_traits::get_value_type(matcl::mat_code mt)
{
    switch (mt)
    {
        case mat_code::integer_dense:
        case mat_code::integer_sparse:
        case mat_code::integer_band:
        case mat_code::integer_scalar:
            return matcl::value_code::v_integer;

        case mat_code::float_dense:
        case mat_code::float_sparse:
        case mat_code::float_band:
        case mat_code::float_scalar:
            return matcl::value_code::v_float;

        case mat_code::real_dense:
        case mat_code::real_sparse:
        case mat_code::real_band:
        case mat_code::real_scalar:
            return matcl::value_code::v_real;

        case mat_code::float_complex_dense:						
        case mat_code::float_complex_sparse:						
        case mat_code::float_complex_band:
        case mat_code::float_complex_scalar:
            return matcl::value_code::v_float_complex;

        case mat_code::complex_dense:						
        case mat_code::complex_sparse:						
        case mat_code::complex_band:
        case mat_code::complex_scalar:
            return matcl::value_code::v_complex;

        case mat_code::object_dense:						
        case mat_code::object_sparse:						
        case mat_code::object_band:
        case mat_code::object_scalar:
            return matcl::value_code::v_object;

        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};
matcl::struct_code matrix_traits::get_struct_type(matcl::mat_code mt)
{
    switch (mt)
    {
        case mat_code::integer_dense:
        case mat_code::float_dense:
        case mat_code::real_dense:
        case mat_code::float_complex_dense:
        case mat_code::complex_dense:
        case mat_code::object_dense:
        {
            return matcl::struct_code::struct_dense;
        }			
        case mat_code::integer_scalar:
        case mat_code::real_scalar:
        case mat_code::float_scalar:
        case mat_code::float_complex_scalar:
        case mat_code::complex_scalar:
        case mat_code::object_scalar:
        {
            return matcl::struct_code::struct_scalar;
        }			
        case mat_code::integer_sparse:
        case mat_code::float_sparse:
        case mat_code::real_sparse:
        case mat_code::float_complex_sparse:
        case mat_code::complex_sparse:
        case mat_code::object_sparse:
        {
            return matcl::struct_code::struct_sparse;
        }			
        case mat_code::integer_band:
        case mat_code::float_band:
        case mat_code::real_band:
        case mat_code::float_complex_band:
        case mat_code::complex_band:
        case mat_code::object_band:
        {
            return matcl::struct_code::struct_banded;
        }
        default:			
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};
matcl::mat_code matrix_traits::get_matrix_type(matcl::value_code vt, 
                                             matcl::struct_code st)
{
    switch (vt)
    {
        case matcl::value_code::v_integer:
        {
            switch (st)
            {
                case matcl::struct_code::struct_dense: 
                    return matcl::mat_code::integer_dense;
                case matcl::struct_code::struct_sparse:
                    return matcl::mat_code::integer_sparse;
                case matcl::struct_code::struct_banded:
                    return matcl::mat_code::integer_band;
                case matcl::struct_code::struct_scalar:
                    return matcl::mat_code::integer_scalar;
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw;
                }
            };
        }
        case matcl::value_code::v_float:
        {
            switch (st)
            {
                case matcl::struct_code::struct_dense:
                    return matcl::mat_code::float_dense;
                case matcl::struct_code::struct_sparse:
                    return matcl::mat_code::float_sparse;
                case matcl::struct_code::struct_banded:
                    return matcl::mat_code::float_band;
                case matcl::struct_code::struct_scalar:
                    return matcl::mat_code::float_scalar;
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw;
                }
            };
        }

        case matcl::value_code::v_real:
        {
            switch (st)
            {
                case matcl::struct_code::struct_dense:
                    return matcl::mat_code::real_dense;
                case matcl::struct_code::struct_sparse:
                    return matcl::mat_code::real_sparse;
                case matcl::struct_code::struct_banded:
                    return matcl::mat_code::real_band;
                case matcl::struct_code::struct_scalar:
                    return matcl::mat_code::real_scalar;
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw;
                }
            };
        }
        case matcl::value_code::v_float_complex:
        {
            switch (st)
            {
                case matcl::struct_code::struct_dense:
                    return matcl::mat_code::float_complex_dense;
                case matcl::struct_code::struct_sparse:
                    return matcl::mat_code::float_complex_sparse;
                case matcl::struct_code::struct_banded:
                    return matcl::mat_code::float_complex_band;
                case matcl::struct_code::struct_scalar:
                    return matcl::mat_code::float_complex_scalar;
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw;
                }
            };
        }
        case matcl::value_code::v_complex:
        {
            switch (st)
            {
                case matcl::struct_code::struct_dense: 
                    return matcl::mat_code::complex_dense;
                case matcl::struct_code::struct_sparse:
                    return matcl::mat_code::complex_sparse;
                case matcl::struct_code::struct_banded:
                    return matcl::mat_code::complex_band;
                case matcl::struct_code::struct_scalar:
                    return matcl::mat_code::complex_scalar;
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw;
                }
            };
        }
        case matcl::value_code::v_object:
        {
            switch (st)
            {
                case matcl::struct_code::struct_dense:
                    return matcl::mat_code::object_dense;
                case matcl::struct_code::struct_sparse:
                    return matcl::mat_code::object_sparse;
                case matcl::struct_code::struct_banded:
                    return matcl::mat_code::object_band;
                case matcl::struct_code::struct_scalar:
                    return matcl::mat_code::object_scalar;
                default:
                {
                    matcl_assert(0,"unknown case");
                    throw;
                }
            };
        }
        default:
        {
            matcl_assert(0,"unknown case");
            throw;
        }
    };
};

matcl::value_code matrix_traits::real_value_type(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:  
            return matcl::value_code::v_integer;
        case matcl::value_code::v_float:
            return matcl::value_code::v_float;
        case matcl::value_code::v_real:
            return matcl::value_code::v_real;
        case matcl::value_code::v_float_complex:
            return matcl::value_code::v_float;
        case matcl::value_code::v_complex:
            return matcl::value_code::v_real;
        case matcl::value_code::v_object:
            return matcl::value_code::v_object;
        default:
            return matcl::value_code::v_real;
    };
};
matcl::value_code matrix_traits::unify_value_types(matcl::value_code v1,
                                                   matcl::value_code v2)
{
    switch (v1)
    {
        case matcl::value_code::v_integer:
        {
            switch (v2)
            {
                case matcl::value_code::v_integer:  
                {
                    using type = typename md::unify_types<Integer,Integer>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float:
                {
                    using type = typename md::unify_types<Integer,Float>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_real:
                {
                    using type = typename md::unify_types<Integer,Real>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float_complex:
                {
                    using type = typename md::unify_types<Integer,Float_complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_complex:
                {
                    using type = typename md::unify_types<Integer,Complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_object:
                {
                    using type = typename md::unify_types<Integer,Object>::type;
                    return md::value_to_code<type>::value;
                }
                default:
                {
                    return matcl::value_code::v_real;
                }
            }
        }
        case matcl::value_code::v_float:
        {
            switch (v2)
            {
                case matcl::value_code::v_integer:  
                {
                    using type = typename md::unify_types<Float,Integer>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float:
                {
                    using type = typename md::unify_types<Float,Float>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_real:
                {
                    using type = typename md::unify_types<Float,Real>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float_complex:
                {
                    using type = typename md::unify_types<Float,Float_complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_complex:
                {
                    using type = typename md::unify_types<Float,Complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_object:
                {
                    using type = typename md::unify_types<Float,Object>::type;
                    return md::value_to_code<type>::value;
                }
                default:
                {
                    return matcl::value_code::v_real;
                }
            }
        }
        case matcl::value_code::v_real:
        {
            switch (v2)
            {
                case matcl::value_code::v_integer:  
                {
                    using type = typename md::unify_types<Real,Integer>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float:
                {
                    using type = typename md::unify_types<Real,Float>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_real:
                {
                    using type = typename md::unify_types<Real,Real>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float_complex:
                {
                    using type = typename md::unify_types<Real,Float_complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_complex:
                {
                    using type = typename md::unify_types<Real,Complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_object:
                {
                    using type = typename md::unify_types<Real,Object>::type;
                    return md::value_to_code<type>::value;
                }
                default:
                {
                    return matcl::value_code::v_real;
                }
            }
        }
        case matcl::value_code::v_float_complex:
        {
            switch (v2)
            {
                case matcl::value_code::v_integer:  
                {
                    using type = typename md::unify_types<Float_complex,Integer>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float:
                {
                    using type = typename md::unify_types<Float_complex,Float>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_real:
                {
                    using type = typename md::unify_types<Float_complex,Real>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float_complex:
                {
                    using type = typename md::unify_types<Float_complex,Float_complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_complex:
                {
                    using type = typename md::unify_types<Float_complex,Complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_object:
                {
                    using type = typename md::unify_types<Float_complex,Object>::type;
                    return md::value_to_code<type>::value;
                }
                default:
                {
                    return matcl::value_code::v_real;
                }
            }
        }
        case matcl::value_code::v_complex:
        {
            switch (v2)
            {
                case matcl::value_code::v_integer:
                {
                    using type = typename md::unify_types<Complex,Integer>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float:
                {
                    using type = typename md::unify_types<Complex,Float>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_real:
                {
                    using type = typename md::unify_types<Complex,Real>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float_complex:
                {
                    using type = typename md::unify_types<Complex,Float_complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_complex:
                {
                    using type = typename md::unify_types<Complex,Complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_object:
                {
                    using type = typename md::unify_types<Complex,Object>::type;
                    return md::value_to_code<type>::value;
                }
                default:
                {
                    return matcl::value_code::v_real;
                }
            }
        }
        case matcl::value_code::v_object:
        {
            switch (v2)
            {
                case matcl::value_code::v_integer:  
                {
                    using type = typename md::unify_types<Object,Integer>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float:
                {
                    using type = typename md::unify_types<Object,Float>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_real:
                {
                    using type = typename md::unify_types<Object,Real>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_float_complex:
                {
                    using type = typename md::unify_types<Object,Float_complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_complex:
                {
                    using type = typename md::unify_types<Object,Complex>::type;
                    return md::value_to_code<type>::value;
                }
                case matcl::value_code::v_object:
                {
                    using type = typename md::unify_types<Object,Object>::type;
                    return md::value_to_code<type>::value;
                }
                default:
                {
                    return matcl::value_code::v_real;
                }
            }
        }
        default:
            return matcl::value_code::v_real;
    };
};

bool matrix_traits::is_float(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:          return false;
        case matcl::value_code::v_float:            return true;
        case matcl::value_code::v_real:             return true;
        case matcl::value_code::v_float_complex:    return true;
        case matcl::value_code::v_complex:          return true;
        case matcl::value_code::v_object:           return false;
        default:                                    return false;
    };
};
bool matrix_traits::is_real(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:          return true;
        case matcl::value_code::v_float:            return true;
        case matcl::value_code::v_real:             return true;
        case matcl::value_code::v_float_complex:    return false;
        case matcl::value_code::v_complex:          return false;
        case matcl::value_code::v_object:           return false;
        default:                                    return false;
    };
};

bool matrix_traits::is_float_real(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:          return false;
        case matcl::value_code::v_float:            return true;
        case matcl::value_code::v_real:             return true;
        case matcl::value_code::v_float_complex:    return false;
        case matcl::value_code::v_complex:          return false;
        case matcl::value_code::v_object:           return false;
        default:                                    return false;
    };
};
bool matrix_traits::is_float_complex(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:          return false;
        case matcl::value_code::v_float:            return false;
        case matcl::value_code::v_real:             return false;
        case matcl::value_code::v_float_complex:    return true;
        case matcl::value_code::v_complex:          return true;
        case matcl::value_code::v_object:           return false;
        default:                                    return false;
    };
};

value_code matrix_traits::complex_value_type(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:
            return matcl::value_code::v_complex;

        case matcl::value_code::v_float:
            return matcl::value_code::v_float_complex;

        case matcl::value_code::v_real:
            return matcl::value_code::v_complex;

        case matcl::value_code::v_float_complex:
            return matcl::value_code::v_float_complex;

        case matcl::value_code::v_complex:
            return matcl::value_code::v_complex;

        case matcl::value_code::v_object:
            return matcl::value_code::v_object;

        default:
            return matcl::value_code::v_complex;
    };
};

bool matrix_traits::is_single_precision(matcl::value_code vc)
{
    if (vc == matcl::value_code::v_float 
        || vc == matcl::value_code::v_float_complex)
    {
        return true;
    }
    else
    {
        return false;
    }
}

value_code matrix_traits::double_precision(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer:
            return matcl::value_code::v_real;

        case matcl::value_code::v_float:
            return matcl::value_code::v_real;

        case matcl::value_code::v_real:
            return matcl::value_code::v_real;

        case matcl::value_code::v_float_complex:
            return matcl::value_code::v_complex;

        case matcl::value_code::v_complex:
            return matcl::value_code::v_complex;

        case matcl::value_code::v_object:
            return matcl::value_code::v_object;

        default:
            return matcl::value_code::v_complex;
    };
}
value_code matrix_traits::single_precision(matcl::value_code v)
{
    switch (v)
    {
        case matcl::value_code::v_integer: 
            return matcl::value_code::v_float;

        case matcl::value_code::v_float:
            return matcl::value_code::v_float;

        case matcl::value_code::v_real:
            return matcl::value_code::v_float;

        case matcl::value_code::v_float_complex:
            return matcl::value_code::v_float_complex;

        case matcl::value_code::v_complex:
            return matcl::value_code::v_float_complex;

        case matcl::value_code::v_object:
            return matcl::value_code::v_object;

        default:
            return matcl::value_code::v_complex;
    };
}

};