/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Paweł Kowal 2017 - 2018
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

#include "matcl-blas-lapack/blas_loader/blas_loader.h"

namespace matcl { namespace lapack
{

class plugin_manager
{
    private:
        static ::blas_plugin*   m_plugin;
        static void*            m_lib_handle;

    public:
        static void             init_plugin() {};
        static ::blas_plugin*   plugin();
        static bool             load_plugin(const std::string& path);
        static void             initialize_plugin();

    private:
        static void             init_env();


        template<class plugins_list>
        static bool             get_plugin(const plugins_list& plugins, ::blas_plugin*& plugin,
                                    void*& lib_handle);

        template<class plugins_list>
        static bool             get_plugin_mix(const plugins_list& plugins, ::blas_plugin*& plugin, 
                                    const ::blas_plugin* plugin_cpu, const ::blas_plugin* plugin_gpu,
                                    void*& lib_handle);

        static ::blas_plugin*   get_plugin();
        static void             change_plugin(::blas_plugin* plug, void* mod);

        static void             free_library(void* dll_handle);
        static void*            load_library(const char* path);
        static void*            get_proc_address(void* mod, const char* func);  
        static bool             set_environment_variable(const char* name, const char* value);
};

inline ::blas_plugin* plugin_manager::plugin()
{
    return m_plugin;
}

}}
