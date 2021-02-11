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

#pragma once

#include "matcl-matrep/general/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/general/fwd_decls.h"

#include <memory>
#include <string>

namespace matcl { namespace dynamic
{
    template<class Val> class object_type;
    
    class object;
    class Type;

};};

namespace matcl
{
    //template<class T>
    //using shared_ptr = std::shared_ptr<T>;

    class Matrix;
    class unique_matrix;
    class colon;
    class mat_row;
    class mat_col;
    class sub_matrix;
    class sub_matrix_1;
    class sub_matrix_2;

    using Object            = dynamic::object;
    using String            = dynamic::object_type<std::string>;
    using OBool             = dynamic::object_type<bool>;
    using OInteger          = dynamic::object_type<Integer>;
    using OReal             = dynamic::object_type<Real>;
    using OFloat            = dynamic::object_type<Float>;
    using OComplex          = dynamic::object_type<Complex>;
    using OFloat_complex    = dynamic::object_type<Float_complex>;
   
    class test_function;
    class disp_stream;
    class options;
    class options_notifier;
    class rand_state;
    class struct_flag;
    class user_flag;
    class user_flag_config;
    
    using disp_stream_ptr   = std::shared_ptr<disp_stream>;

    class iarchive;
    class oarchive;
    
    template<class derived> class test_function_templ;    

    template<class Val, bool Is_safe = true>
    class dense_matrix;

    template<class Val, bool Is_safe = true>
    class sparse_matrix;

    template<class Val, bool Is_safe = true>
    class band_matrix;    
    
    template<class Val> class sub_dense_matrix;

    template<class Val> class sub_sparse_matrix;
    template<class Val> class sub_sparse_matrix_1;
    template<class Val> class sub_sparse_matrix_2;

    template<class Val> class sub_band_matrix;
    template<class Val> class sub_band_matrix_1;
    template<class Val> class sub_band_matrix_2;

    template<class Val> class dense_row;

    template<class Val> class dense_col;
    template<class Val> class sparse_row;
    template<class Val> class sparse_col;

    template<class T> class object_matrix;
    template<class T> class object_row;
    template<class T> class object_col;

    template<class T> class sub_object_matrix;
    template<class T> class sub_object_matrix_1;
    template<class T> class sub_object_matrix_2;
};

namespace matcl { namespace details
{
    template<class V, class S>
    class matrix_container;

    class matrix_container_base;
    class matrix_base;

	class mat_cons_data;
    class printer;
	struct matrix_data_accesser;	    
    class registered_user_flags;

	template<class M> struct sparse_matrix_constructor_row;
	template<class M> struct dense_matrix_constructor_row;
	template<class M> struct sparse_matrix_constructor_col;
	template<class M> struct dense_matrix_constructor_col;

    class pv_constructor;
    class local_rand_state_raii;
    class enable_warnings_raii;

    class disp_stream_impl;
    class printer;
    class rand_state;
};};

namespace matcl { namespace raw 
{

template<class value_type, class struct_type> 
class MATCL_MATREP_EXPORT Matrix;

template<class value_type>
class MATCL_MATREP_EXPORT dense_matrix_base;

template<class value_type>
class MATCL_MATREP_EXPORT sparse_matrix_base;

};};

namespace matcl { namespace raw { namespace details
{

template<class value_type>	class sparse_ccs;

};};};

namespace matcl { namespace error
{

class exception_message;
using exception_message_ptr         = std::shared_ptr<exception_message>;

}};
