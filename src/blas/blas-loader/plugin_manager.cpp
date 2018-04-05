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

#include "plugin_manager.h"
#include "blas_config.h"

#include <thread>
#include <sstream>
#include <iostream>

#ifdef _WIN64
    #define BLAS_LIBRARY_OPENBLAS           "matcl-openblas-plugin-x64-Release.dll"
    #define BLAS_LIBRARY_CLAPACK            "matcl-clapack-plugin-x64-Release.dll"
#elif defined _WIN32
    #define BLAS_LIBRARY_OPENBLAS           "matcl-openblas-plugin-Win32-Release.dll"
    #define BLAS_LIBRARY_CLAPACK            "matcl-clapack-plugin-Win32-Release.dll"
#else
    #error unknown platform
#endif

namespace matcl { namespace lapack
{

::blas_plugin* plugin_manager::m_plugin = get_plugin();
void* plugin_manager::m_lib_handle      = nullptr;

void plugin_manager::init_env()
{
    //it seems to be impossible to set MKL threads using omp_set_num_threads, we need
    //to set environmental variable 

    int n_threads = std::thread::hardware_concurrency();

    std::ostringstream os;
    os << n_threads;

    std::string n_str   = os.str();

    const char lpName[] = "OMP_NUM_THREADS";
    const char* lpValue = n_str.c_str();

    set_environment_variable(lpName, lpValue);
};

void plugin_manager::change_plugin(::blas_plugin* plug, void* mod)
{
    if (m_lib_handle)
        free_library(m_lib_handle);

    m_plugin        = plug;
    m_lib_handle    = mod;
}

::blas_plugin* plugin_manager::get_plugin()
{
    init_env();

    blas_config config;

    config.add_plugin_cpu(BLAS_LIBRARY_OPENBLAS);
    config.add_plugin_cpu(BLAS_LIBRARY_CLAPACK);

    ::blas_plugin* plugin_cpu   = nullptr;    
    void* lib_handle_cpu        = nullptr;

    bool ok_cpu = get_plugin(config.get_plugins_cpu(), plugin_cpu, lib_handle_cpu);

    if (ok_cpu == true && plugin_cpu == nullptr)
    {
        std::cerr << "Unable to load blas library plugin" << "\n";
        return nullptr;
    }

    m_lib_handle    = lib_handle_cpu;
    return plugin_cpu;

    /*
    //TODO
    ::blas_plugin* plugin_gpu   = nullptr;
    void* lib_handle_gpu        = nullptr;
    void* lib_handle_mix        = nullptr;

    bool ok_gpu = get_plugin(config.get_plugins_gpu(), plugin_gpu, lib_handle_gpu);

    if (ok_gpu == true && plugin_gpu != nullptr)
    {
        ::blas_plugin* plugin_cpu_gpu = nullptr;

        bool ok = get_plugin_mix(config.get_plugins_cpu_gpu(), plugin_cpu_gpu, 
                                    plugin_cpu, plugin_gpu, lib_handle_mix);

        if (ok == true && plugin_cpu_gpu != nullptr)
            return plugin_cpu_gpu;
        else
            return plugin_cpu;
    }
    else
    {
        return plugin_cpu;
    };
    */
}    

bool plugin_manager::load_plugin(const std::string& path_in)
{
    std::string path    = blas_config::normalize_path(path_in);

    void* mod           = nullptr;
    mod                 = load_library(path.c_str());

    if (mod == nullptr)
        return false;

    using plugin_type = blas_plugin::get_blas_plugin_type;

    plugin_type* f = (plugin_type*)get_proc_address(mod, "get_blas_plugin");    

    if (f == nullptr)
    {
        free_library(mod);
        return false;
    };

    ::blas_plugin* ret = f();

    if (ret == nullptr)
    {
        free_library(mod);
        return false;
    };

    ret->initialize_fptr();
    change_plugin(ret, mod);
    return true;
}

void plugin_manager::initialize_plugin()
{
    if (!m_plugin)
        return;

    m_plugin->force_initialization_fptr();
}

template<class plugins_list>
bool plugin_manager::get_plugin(const plugins_list& plugins, ::blas_plugin*& plugin,
                                void*& lib_handle)
{
    plugin          = nullptr;
    void* mod       = nullptr;
    lib_handle      = nullptr;

    for (auto pos = plugins.begin(); pos != plugins.end(); ++pos)
    {
        const auto& loc_plugin = *pos;
        mod = load_library(loc_plugin.c_str());

        if (mod != nullptr)
            break;
    };

    if (!mod)
        return true;

    using plugin_type = blas_plugin::get_blas_plugin_type;

    plugin_type* f = (plugin_type*)get_proc_address(mod, "get_blas_plugin");
    ::blas_plugin* ret = nullptr;

    if (f)
        ret = f();

    if (!ret)
    {
        free_library(mod);
        std::cerr << "unable to init blas plugin from blas library" << "\n";
        return false;
    }

    ret->initialize_fptr();

    plugin      = ret;
    lib_handle  = mod;
    return true;
}

template<class plugins_list>
bool plugin_manager::get_plugin_mix(const plugins_list& plugins, ::blas_plugin*& plugin, 
                            const ::blas_plugin* plugin_cpu, const ::blas_plugin* plugin_gpu,
                            void*& lib_handle)
{
    plugin          = nullptr;
    void* mod       = nullptr;
    lib_handle      = nullptr;

    for (auto pos = plugins.begin(); pos != plugins.end(); ++pos)
    {
        const auto& loc_plugin = *pos;
        mod = load_library(loc_plugin.c_str());

        if (mod != nullptr)
            break;
    };

    if (!mod)
        return true;

    using plugin_type = blas_plugin::get_blas_plugin_gpu_cpu_type;

    plugin_type* f = (plugin_type*)GetProcAddress(mod, "get_blas_plugin_gpu_cpu");
    ::blas_plugin* ret = nullptr;

    if (f)
        ret = f(plugin_cpu, plugin_gpu);

    if (!ret)
    {
        free_library(mod);
        std::cerr << "unable to init blas plugin from blas library" << "\n";
        return false;
    }

    plugin      = ret;
    lib_handle  = mod;
    return true;
}

#if defined(WIN32) || defined(WIN64)
    #include <Windows.h>

    void plugin_manager::free_library(void* dll_handle)
    {
        FreeLibrary((HINSTANCE)dll_handle);
    }

    void* plugin_manager::load_library(const char* path)
    {
        return LoadLibrary(path);
    };

    void* plugin_manager::get_proc_address(void* mod, const char* func)
    {
        return GetProcAddress((HINSTANCE)mod, func);
    };

    bool plugin_manager::set_environment_variable(const char* name, const char* value)
    {
        return SetEnvironmentVariable(name, value);
    };

#elif defined(__unix__)

    //TODO: UNIX, default blas plugins

    #include <sys/types.h>
    #include <dlfcn.h>

    void* plugin_manager::load_library(const char* path)
    {
        return dlopen(path, RTLD_LAZY);
    };

    void* plugin_manager::get_proc_address(void* mod, const char* func)
    {
        return dlsym(mod, func);
    };

    void plugin_manager::free_library(void* dll_handle)
    {
        dlclose(dll_handle);
    }

    bool plugin_manager::set_environment_variable(const char* name, const char* value)
    {
        //TODO: UNIX
        return SetEnvironmentVariable(name, value);
    };
#else
    #error unknown system
#endif

}}
