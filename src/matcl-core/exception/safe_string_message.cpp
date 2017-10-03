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

#include "matcl-core/error/safe_string_message.h"
#include "matcl-core/general/thread.h"
#include <mutex>

namespace matcl { namespace error
{
class safe_string_message::safe_string_message_impl
{
    private:
        matcl::default_mutex    msg_mutex;
        std::string             msg;

    public:
        safe_string_message_impl()
        {}

        safe_string_message_impl(const safe_string_message_impl& other)
        {
                std::unique_lock<matcl::default_mutex> lock(msg_mutex);
                msg = other.msg;
        }
        
        safe_string_message_impl& 
        operator=(const safe_string_message_impl& other)
        {
            if(this != &other)
            {
                std::unique_lock<matcl::default_mutex> lock(msg_mutex); 
                msg = other.msg;
            }
        
            return (*this);
        }
        
        safe_string_message_impl& operator=(const std::string& rhs)
        {
            std::unique_lock<matcl::default_mutex> lock(msg_mutex); 
            msg = rhs;
            return (*this);
        }
        
        safe_string_message_impl& operator+=(const std::string& rhs)
        {
            std::unique_lock<matcl::default_mutex> lock(msg_mutex); 
            msg += rhs;
            return (*this);
        }
        
        const char* c_str()
        {
            std::unique_lock<matcl::default_mutex> lock(msg_mutex); 
            return msg.c_str();
        }
};

safe_string_message::safe_string_message()
    : impl(new safe_string_message_impl()) 
{}
 
safe_string_message::safe_string_message(const safe_string_message& other)
    : impl(new safe_string_message_impl(*(other.impl))) 
{}
 
safe_string_message& 
safe_string_message::operator=(const safe_string_message& other) 
{
    *impl = *(other.impl);
    return *this;
}
safe_string_message& safe_string_message::operator=(const std::string& rhs) 
{
    *impl = rhs;
    return *this;
}
safe_string_message& safe_string_message::operator+=(const std::string& rhs) 
{
    *impl += rhs;
    return *this;
} 
const char* safe_string_message::c_str()
{
    return impl->c_str();
}

}}
