/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"

#include "matcl-dynamic/type.h"

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/general/thread.h"
#include "function_table.h"
#include "type_table_cache.h"
#include "matcl-dynamic/details/register_function_impl.h"

#include <map>

#pragma warning( push )
#pragma warning( disable:4251 ) //needs to have dll-interface to be used by clients of class

namespace matcl { namespace dynamic { namespace details
{

class type_table
{
    private:
        matcl::default_mutex            m_mutex_reg_functions;
        matcl::default_spinlock_mutex   m_mutex_build_type;

    public:
        enum predefined_type
        {
            ti_bool, ti_int, ti_float, ti_real, ti_fcompl, ti_compl, ti_string,
            ti_unit, ti_any
        };

    private:
        struct type_handle
        {
            Type ti;
            Type (*constructor)();
            type_handle()       : constructor(nullptr) {};
        };

        using type_map          = std::map<std::string,type_handle>;
        using type_map_ref      = std::map<Type, Type>;
        using type_vec          = std::vector<Type>;

    private:
        type_map                m_type_map;
        type_map_ref            m_type_map_ref;
        function_table			m_fun_tab;        

    private:
        type_table();

    public:
        static type_table*      get();

        Type                    get_predefined(predefined_type t);
        Type                    make_reference_type(Type t);

        Type                    unify_types(Type t1, Type t2);
        function                get_overload(const function_name& func, int n_args, const Type t[]);
        function                get_template_overload(const function_name& func, int n_templ, 
                                        const Type templates[], int n_args, const Type arg_types[]);
        function                get_converter(Type to, Type from, bool implicit);
        function                get_cast_function(Type to, Type from);
        function                get_assigner(Type to, Type from);

        //raw class_name must be given by typeid(T).name()
        Type                    get_type(const char* raw_class_name);

        //class_name must be returned by get_class_name function
        Type                    get_type_from_class(const std::string& class_name);
        static std::string		get_class_name(const char* name);                

        static void				initial_register_function(delayed_function_register* );
        static void				initial_register_function_template(delayed_function_template_register* );
        void                    register_object(const char* obj, Type (*creator)());

        //free cache if additional memory is required
        void                    free_cache();

    private:
        type_table_cache*       get_cache();
        void				    initialize();
        void                    finish_initialization();
        bool				    process_functions(error_handler& eh);
        void				    register_function(const function_name& func, const function& fun_evl, 
                                    error_handler& eh);
        void				    register_function_template(const function_name& func, 
                                    const function& fun_evl, const type_vec& types, 
                                    make_return_fptr ret, error_handler& eh);
        void                    add_predefined_functions(const std::string& type_name);
};

};};};

#pragma warning( pop )