/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-file/config_matcl_file.h"
#include "matcl-core/error/exception_classes.h"
#include "exception_message_file.h"

#include <exception>
#include <string>

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 
#pragma warning(disable:4275)	//non dll-interface class used as base for dll-interface class

namespace matcl { namespace error
{

class MATCL_FILE_EXPORT matcl_file_exception : public matcl_exception
{
    public:
        virtual const char* what() const throw() override;

    private:
        virtual const char* what(exception_message& em) const override;
        virtual const char* what(exception_message_file& em) const = 0;
};

class MATCL_FILE_EXPORT error_matclfile_locked : public matcl_file_exception
{
    public:
        error_matclfile_locked(){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_create_matclfile : public matcl_file_exception
{
    public:
        std::string msg;
    public:
        error_create_matclfile(const std::string& msg) : msg(msg){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_read_matclfile : public matcl_file_exception
{
    public:
        std::string msg;
    public:
        error_read_matclfile(const std::string& msg) : msg(msg){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_read_matclfile_mat_not_exist : public matcl_file_exception
{
    public:
        std::string mat_name;
    public:
        error_read_matclfile_mat_not_exist(const std::string& msg) : mat_name(msg){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_write_matclfile_mat_already_exist : public matcl_file_exception
{
    public:
        std::string mat_name;
    public:
        error_write_matclfile_mat_already_exist(const std::string& msg) : mat_name(msg){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_open_matclfile : public matcl_file_exception
{
    public:
        std::string file_name;
        std::string msg;

    public:
        error_open_matclfile(const std::string& file, const std::string& msg = "")
            : file_name(file), msg(msg)	{};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_open_matfile : public matcl_file_exception
{
    public:
        std::string file_name;
        std::string msg;

    public:
        error_open_matfile(const std::string& file, const std::string& msg = "")
            : file_name(file), msg(msg)	{};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_matfile_not_opened : public matcl_file_exception
{
    public:
        error_matfile_not_opened(){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_read_matfile : public matcl_file_exception
{
    public:
        error_read_matfile(){};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_read_matfile_var : public matcl_file_exception
{
    public:
        std::string var_name;
    public:
        error_read_matfile_var(const std::string& var_name)
            :var_name(var_name)
        {};

        virtual const char* what(exception_message_file& em) const override;
};

class MATCL_FILE_EXPORT error_write_matfile : public matcl_file_exception
{
    public:
        error_write_matfile(){};

        virtual const char* what(exception_message_file& em) const override;
};

};};

#pragma warning(pop)