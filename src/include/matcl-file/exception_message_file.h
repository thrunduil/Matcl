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
#include "matcl-core/error/safe_string_message.h"
#include <memory>
#include <string>

#pragma warning(push)
#pragma warning(disable:4251)   //'matcl::error::exception_message::m_linalg' : class needs to have dll-interface

namespace matcl { namespace error
{

class exception_message_file;
using  exception_message_file_ptr               = std::shared_ptr<exception_message_file>;

MATCL_FILE_EXPORT void                          set_global_messanger_file(exception_message_file_ptr msg);
MATCL_FILE_EXPORT exception_message_file_ptr    get_global_messanger_file();

class MATCL_FILE_EXPORT exception_message_file
{
    public:
        virtual ~exception_message_file(){};
        
        virtual const char* error_open_matfile(const std::string& file, const std::string& msg)= 0;
        virtual const char* error_matfile_not_opened() = 0;
        virtual const char* error_read_matfile() = 0;
        virtual const char* error_read_matfile_var(const std::string& var_name) = 0;
        virtual const char* error_write_matfile() = 0;
        virtual const char* error_open_mmlibfile(const std::string& file, const std::string& msg)= 0;
        virtual const char* error_mmlibfile_locked()= 0;
        virtual const char* error_create_mmlibfile(const std::string& msg)= 0;
        virtual const char* error_read_mmlibfile(const std::string& msg)= 0;
        virtual const char* error_read_mmlibfile_mat_not_exist(const std::string& mat)= 0;
        virtual const char* error_write_mmlibfile_mat_already_exist(const std::string& mat)= 0;
};

class MATCL_FILE_EXPORT default_exception_message_file : public exception_message_file
{
    private:
        safe_string_message current_message;

    public:
        virtual const char* error_open_matfile(const std::string& file, const std::string& msg) override;
        virtual const char* error_matfile_not_opened() override;
        virtual const char* error_read_matfile() override;
        virtual const char* error_read_matfile_var(const std::string& var_name) override;
        virtual const char* error_write_matfile() override;
        virtual const char* error_open_mmlibfile(const std::string& file, const std::string& msg) override;
        virtual const char* error_mmlibfile_locked() override;
        virtual const char* error_create_mmlibfile(const std::string& msg) override;
        virtual const char* error_read_mmlibfile(const std::string& msg) override;
        virtual const char* error_read_mmlibfile_mat_not_exist(const std::string& mat) override;
        virtual const char* error_write_mmlibfile_mat_already_exist(const std::string& mat) override;
};

};};

#pragma warning(pop)