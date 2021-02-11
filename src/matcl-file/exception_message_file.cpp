/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-file/exception_message_file.h"
#include <sstream>

namespace matcl { namespace error
{

static exception_message_file_ptr get_global_message_file()
{
    exception_message_file_ptr em(new default_exception_message_file());
    return em;
};

static exception_message_file_ptr messanger = get_global_message_file();

void error::set_global_messanger_file(exception_message_file_ptr msg)
{
    if (msg)
        messanger = msg;
};

exception_message_file_ptr error::get_global_messanger_file()
{
    return messanger;
}

const char* default_exception_message_file::error_open_matclfile(const std::string& file, const std::string& msg)
{
    std::ostringstream buf;

    buf << "unable to open file " << file;
    if (msg.size() > 0)
        buf <<", reason: " << msg;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_open_matfile(const std::string& file, const std::string& msg)
{
    std::ostringstream buf;

    buf << "unable to open file " << file;
    if (msg.size() > 0)
        buf <<", reason: " << msg;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_matfile_not_opened()
{
    std::ostringstream buf;

    buf << "matfile not opened";
    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_read_matfile()
{
    std::ostringstream buf;

    buf << "error while reading matfile";
    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_read_matfile_var(const std::string& var_name)
{
    std::ostringstream buf;

    buf << "error while reading matfile, matrix " << var_name << " not found, or memory error occured";
    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_write_matfile()
{
    std::ostringstream buf;

    buf << "error while writing to matfile";
    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_matclfile_locked()
{
    current_message = "matcl file is locked";
    return current_message.c_str();
};

const char* default_exception_message_file::error_create_matclfile(const std::string& msg)
{
    std::ostringstream buf;

    buf << "unable to create matcl file, reason: " << msg;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_read_matclfile(const std::string& msg)
{
    std::ostringstream buf;

    buf << "unable to read matcl file, reason: " << msg;

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_read_matclfile_mat_not_exist(const std::string& msg)
{
    std::ostringstream buf;

    buf << "matrix " << msg <<" does not exist in matcl file";

    current_message = buf.str();
    return current_message.c_str();
};

const char* default_exception_message_file::error_write_matclfile_mat_already_exist(const std::string& msg)
{
    std::ostringstream buf;

    buf << "unable to save " << msg <<" matrix, matrix already exists";

    current_message = buf.str();
    return current_message.c_str();
};

};};
