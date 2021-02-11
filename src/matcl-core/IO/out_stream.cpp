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

#include "matcl-core/IO/output_stream.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/IO/out_stream_initializer.h"

#include <iostream>

namespace matcl { namespace details
{

class out_stream_buf: public std::stringbuf
{
    private:
        using base_type = std::stringbuf;

    private:
        std::mutex	    m_mutex;

    public:
        virtual int sync() override
        {
            global_output_stream()->disp(str());
            str("");
            return 0;
        }

        virtual std::streamsize xsputn(const char_type* s, std::streamsize n) override
        {
            using scoped_lock   = std::unique_lock<std::mutex>;

            scoped_lock lock(m_mutex);
            std::streamsize out = base_type::xsputn(s, n);

            for(std::streamsize i = 0; i < n; ++i)
            {
                if (s[i] == '\n')
                {
                    sync();
                    return out;
                }
            }

            return out;
        };
};

// memory for the stream object
static typename std::aligned_storage<sizeof(std::ostream), alignof(std::ostream)>::type
out_data;

static typename std::aligned_storage<sizeof(out_stream_buf), alignof(out_stream_buf)>::type
out_stream_buf_data;

}}

namespace matcl
{

std::ostream& matcl::out_stream = reinterpret_cast<std::ostream&>(details::out_data);

};

namespace matcl { namespace details
{

void out_stream_initializer::open_global()
{
    new(&out_stream_buf_data) out_stream_buf();

    out_stream_buf* buf = reinterpret_cast<out_stream_buf*>(&out_stream_buf_data);

    new (&out_data) std::ostream(buf);
};

void out_stream_initializer::close_global()
{
    std::ostream* out_ptr   = &out_stream;
    out_stream_buf* buf_ptr = reinterpret_cast<out_stream_buf*>(&out_stream_buf_data);

    out_ptr->~basic_ostream();    
    buf_ptr->~out_stream_buf();
};

}}
