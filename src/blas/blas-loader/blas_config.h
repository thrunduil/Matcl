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

#pragma once

#include <vector>
#include <string>
#include <map>
#include <functional>

namespace matcl { namespace lapack
{

class blas_config
{
    public:
        using string_vec    = std::vector<std::string>;
        using action        = std::function<bool(const std::string&)>;
        using action_table  = std::map<std::string, action>;

    private:
        string_vec      m_plugins_cpu;
        string_vec      m_plugins_gpu;
        string_vec      m_plugins_cpu_gpu;
        action_table    m_actions;

    public:
        blas_config();

        string_vec      get_plugins_cpu() const;
        string_vec      get_plugins_gpu() const;
        string_vec      get_plugins_cpu_gpu() const;
        void            add_plugin_cpu(const std::string& path);
        void            add_plugin_gpu(const std::string& path);
        void            add_plugin_cpu_gpu(const std::string& path);

    private:
        std::wstring    get_config_path() const;
        void            load_config();
        void            init_actions();
        void            parse_line(const std::string& line);
        bool            insert_value(const std::string& var, const std::string& value);
};

}}