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

#pragma once

#include "matcl-core/config.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/memory/alloc.h"

#include <string>
#include <iosfwd>
#include <memory>

namespace matcl 
{

// this class is responsible for printing on proper device
class MATCL_CORE_EXPORT output_stream : public matcl_new_delete
{
    public:
        virtual ~output_stream() {};

        // return width of current terminal, for example console, default
        // implementation return width of the standard output console
        virtual Integer     get_terminal_width() const;

        // display string stream
        virtual void        disp(std::stringstream& of) = 0;
};

// convert std::ostream to output_stream
// this stream is not protected by mutex and should be used in
// one thread only
class MATCL_CORE_EXPORT output_stream_from_ostream : public output_stream
{
    private:
        std::ostream&   m_stream;

    public:
        explicit output_stream_from_ostream(std::ostream& os);
        virtual ~output_stream_from_ostream() {};

        virtual void    disp(std::stringstream& of) override;

    private:
        output_stream_from_ostream(const output_stream_from_ostream&) = delete;
        output_stream_from_ostream& operator=(const output_stream_from_ostream&) = delete;
};

// convert std::ostream to output_stream
// this stream is protected by mutex
class MATCL_CORE_EXPORT output_stream_from_ostream_synchronized : public output_stream
{
    private:
        using this_type         = output_stream_from_ostream_synchronized;

    private:
        std::ostream&           m_stream;
        matcl::default_mutex*   m_mutex;

    public:
        explicit output_stream_from_ostream_synchronized(std::ostream& os);
        virtual ~output_stream_from_ostream_synchronized();

        virtual void    disp(std::stringstream& of) override;

    private:
        this_type(const this_type&) = delete;
        this_type& operator=(const this_type&) = delete;
};

using output_stream_ptr = std::shared_ptr<output_stream>;

// global output stream, which points to current output stream
// if set_current_output_stream() is called later, then returned
// output stream will point to the new output stream
MATCL_CORE_EXPORT 
output_stream_ptr   global_output_stream();

// change current output stream, this will change global output stream
// as well. Setting global output stream as current output stream will
// not have any impact.
MATCL_CORE_EXPORT 
void                set_current_output_stream(const output_stream_ptr&);

// return current output stream. if set_current_output_stream() is called,
// then returned output stream will not change
MATCL_CORE_EXPORT 
output_stream_ptr   get_current_output_stream();

}