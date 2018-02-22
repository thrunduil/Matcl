/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "blas_config.h"
#include <fstream>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/filesystem/operations.hpp>

namespace matcl { namespace lapack
{

blas_config::blas_config()
{
    init_actions();
    load_config();
}

void blas_config::init_actions()
{
    {
        std::string var = "CPU";
        action func     = [=](const std::string& value)->bool 
                        { 
                            this->m_plugins_cpu.push_back(value);
                            return true;
                        };
        m_actions.insert(action_table::value_type(var, func));
    };
    {
        std::string var = "GPU";
        action func = [=](const std::string& value)->bool
        {
            this->m_plugins_gpu.push_back(value);
            return true;
        };
        m_actions.insert(action_table::value_type(var, func));
    };
    {
        std::string var = "CPU_GPU";
        action func = [=](const std::string& value)->bool
        {
            this->m_plugins_cpu_gpu.push_back(value);
            return true;
        };
        m_actions.insert(action_table::value_type(var, func));
    };
}

void blas_config::load_config()
{    
    std::wstring config_path = get_config_path();

    std::ifstream config_file(config_path);

    std::string line;

    while (std::getline(config_file, line))
        parse_line(line);
}

void blas_config::parse_line(const std::string& line)
{
    if (line.size() == 0)
        return;

    //check if comment
    if (line[0] == '%')
        return;

    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    boost::char_separator<char> sep(" \t");

    tokenizer tokens(line, sep);
    
    string_vec items;

    for (auto tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter)
        items.push_back(*tok_iter);

    //structure: VAR VALUE
    if (items.size() == 0)
        return;

    if (items.size() != 2)
    {
        std::cerr << "invalid command: " << line << "\n";
        return;
    }

    const std::string& var      = items[0];
    const std::string& value    = items[1];

    bool success                = insert_value(var, value);

    if (success == false)
    {
        std::cerr << "invalid command: " << line << "\n";
        return;
    };

    return;
};

bool blas_config::insert_value(const std::string& var, const std::string& value)
{
    auto pos = m_actions.find(var);

    if (pos == m_actions.end())
        return false;

    return pos->second(value);
};

blas_config::string_vec blas_config::get_plugins_cpu() const
{
    return m_plugins_cpu;
};

blas_config::string_vec blas_config::get_plugins_cpu_gpu() const
{
    return m_plugins_cpu_gpu;
};

blas_config::string_vec blas_config::get_plugins_gpu() const
{
    return m_plugins_gpu;
};

void blas_config::add_plugin_cpu(const std::string& path)
{
    m_plugins_cpu.push_back(path);
}

void blas_config::add_plugin_cpu_gpu(const std::string& path)
{
    m_plugins_cpu_gpu.push_back(path);
};

void blas_config::add_plugin_gpu(const std::string& path)
{
    m_plugins_gpu.push_back(path);
}

#ifndef __unix__

    #include <Windows.h>

    extern "C" IMAGE_DOS_HEADER __ImageBase;

    std::wstring blas_config::get_config_path() const
    {
        WCHAR   DllPath[MAX_PATH] = { 0 };
        GetModuleFileNameW((HINSTANCE)&__ImageBase, DllPath, _countof(DllPath));

        std::wstring path(DllPath);

        boost::filesystem::wpath p(path);
        p = p.parent_path() / "blas_config.txt";

        return p.wstring();
    }

#else
    std::wstring blas_config::get_config_path() const
    {
        //TODO: implement
        return "blas_config.txt";
    }
#endif

}}