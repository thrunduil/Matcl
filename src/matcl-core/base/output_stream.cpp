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

#include "matcl-core/IO/output_stream.h"
#include <iostream>
#include <sstream>

namespace matcl 
{

output_stream_from_ostream::output_stream_from_ostream(std::ostream& os)
    :m_stream(os)
{};

void output_stream_from_ostream::disp(std::stringstream& of)
{
    of.seekg(0, std::ios::beg);    

    if (of.rdbuf())
        m_stream << of.rdbuf();

    of.seekg(0, std::ios::end);
    m_stream.flush();
};

output_stream_from_ostream_synchronized
    ::output_stream_from_ostream_synchronized(std::ostream& os)
    :m_stream(os)
{
    m_mutex = new matcl::default_mutex();
}

output_stream_from_ostream_synchronized
    ::~output_stream_from_ostream_synchronized() 
{
    delete m_mutex;
}

void output_stream_from_ostream_synchronized::disp(std::stringstream& of)
{
    std::unique_lock<matcl::default_mutex> lock(*m_mutex);
    of.seekg(0, std::ios::beg);

    if (of.rdbuf())
        m_stream << of.rdbuf();

    of.seekg(0, std::ios::end);
    m_stream.flush();
};

class global_output_stream_impl : public output_stream
{
    public:
        global_output_stream_impl()
        {
            using stream_type = output_stream_from_ostream_synchronized;
            m_stream_impl = output_stream_ptr(new stream_type(std::cout));
        };

        virtual void disp(std::stringstream& of) override
        {
            return get_stream_impl()->disp(of);
        };

        output_stream_ptr get_current()
        {
            return get_stream_impl();
        };
        void set_current(const output_stream_ptr& ptr)
        {
            if (!ptr)
                return;

            if (dynamic_cast<const global_output_stream_impl*>(ptr.get()) != nullptr)
            {
                //global stream cannot be current stream, 
                //this avoids strange circularity
                return;
            };

            m_stream_impl = ptr;
        };

    private:
        output_stream_ptr   get_stream_impl() { return m_stream_impl; }
        output_stream_ptr   m_stream_impl;
};

output_stream_ptr matcl::global_output_stream()
{
    static output_stream_ptr impl(new global_output_stream_impl());

    return impl;
};

output_stream_ptr matcl::get_current_output_stream()
{
    global_output_stream_impl* impl 
        = dynamic_cast<global_output_stream_impl*>(global_output_stream().get());

    return impl->get_current();
};

void set_current_output_stream(const output_stream_ptr& os)
{
    global_output_stream_impl* impl 
        = dynamic_cast<global_output_stream_impl*>(global_output_stream().get());

    return impl->set_current(os);
};

#include <windows.h>

Integer output_stream::get_terminal_width() const
{
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    Integer columns, rows;

    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    columns = csbi.srWindow.Right - csbi.srWindow.Left;
    rows = csbi.srWindow.Bottom - csbi.srWindow.Top;

    return columns;
};

}