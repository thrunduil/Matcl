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

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"

#include <memory>
#include <string>

namespace eos
{
    class portable_iarchive;
    class portable_oarchive;
};

namespace matcl { namespace dynamic
{

template<class T>
class object_type;

class object;
class Type;

}};


namespace boost { namespace serialization
{
    class access;
};};

namespace matcl
{
    using iarchive_impl     = eos::portable_iarchive;
    using oarchive_impl     = eos::portable_oarchive;

    class Matrix;
    class struct_flag;
    class disp_stream;
    class options;
    class options_notifier;
    class disp_data_provider;

    class iarchive;
    class oarchive;

    struct struct_dense;
    struct struct_sparse;
    struct struct_banded;

    using Object            = dynamic::object;
    using disp_stream_ptr   = std::shared_ptr<disp_stream>;
}

namespace matcl { namespace details
{

class disp_stream_impl;

}};

namespace matcl { namespace error
{

class exception_message;

using exception_message_ptr = std::shared_ptr<exception_message>;

}};

namespace matcl { namespace raw
{

namespace details
{};

template<class val_type,class struct_type> 
class MATCL_CORE_EXTERN Matrix;

using integer_dense         = Matrix<Integer,struct_dense>;
using real_dense            = Matrix<Real,struct_dense>;
using float_dense           = Matrix<Float,struct_dense>;
using complex_dense         = Matrix<Complex,struct_dense>;
using float_complex_dense   = Matrix<Float_complex,struct_dense>;
using object_dense          = Matrix<Object,struct_dense>;

using integer_sparse        = Matrix<Integer,struct_sparse>;
using real_sparse           = Matrix<Real,struct_sparse>;
using float_sparse          = Matrix<Float,struct_sparse>;
using complex_sparse        = Matrix<Complex,struct_sparse>;
using float_complex_sparse  = Matrix<Float_complex,struct_sparse>;
using object_sparse         = Matrix<Object,struct_sparse>;

using integer_band          = Matrix<Integer,struct_banded>;
using real_band             = Matrix<Real,struct_banded>;
using float_band            = Matrix<Float,struct_banded>;
using complex_band          = Matrix<Complex,struct_banded>;
using float_complex_band    = Matrix<Float_complex,struct_banded>;
using object_band           = Matrix<Object,struct_banded>;

};};