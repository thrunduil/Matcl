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

#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/matrix/scalar_types.h"

#pragma warning(push)
#pragma warning(disable: 4127)  //conditional expression is constant
#pragma warning(disable: 4512)  //assignment operator could not be generated
#pragma warning(disable: 4244)  //'initializing' : conversion from 'const int' 
                                //to 'const boost::uint_least16_t', possible loss of data
#pragma warning(disable: 4996)  //Function call with parameters that may be unsafe - 
                                //this call relies on the caller to check that the passed values are correct.
#pragma warning(disable: 4702)  //unreachable code
#pragma warning(disable: 4100)  //unreferenced formal parameter
#pragma warning(disable: 4310)  //cast truncates constant value

#include "matcl-core/details/portable_iarchive.hpp"
#include "matcl-core/details/portable_oarchive.hpp"

#pragma warning(pop)

#include <iosfwd>

namespace matcl 
{

using iarchive_impl = eos::portable_iarchive;
using oarchive_impl = eos::portable_oarchive;

// input archive
class MATCL_CORE_EXPORT iarchive
{
    private:
        using archive_type  = eos::portable_iarchive;

        archive_type*   m_impl;
        bool            m_need_delete;

        iarchive(const iarchive&) = delete;
        iarchive& operator=(const iarchive&) = delete;

    public:
        // create input archive from input stream
        iarchive(std::istream& is);

        // construct iarchive from an existing archive; this
        // object cannot be used, when ar is destroyed
        explicit iarchive(archive_type& ar);

        // standard destructor
        ~iarchive();

        // get implementation type
        archive_type&   get()       { return *m_impl; };

        // load integer
        iarchive&       operator>>(Integer& v);

        // load double float
        iarchive&       operator>>(Real& v);

        // load single float
        iarchive&       operator>>(Float& v);

        // load double complex
        iarchive&       operator>>(Complex& v);

        // load single complex
        iarchive&       operator>>(Float_complex& v);
};

// output archive
class MATCL_CORE_EXPORT oarchive
{
    private:
        using archive_type  = eos::portable_oarchive;

        archive_type*   m_impl;
        bool            m_need_delete;

        oarchive(const oarchive&) = delete;
        oarchive& operator=(const oarchive&) = delete;

    public:
        // create output archive from output stream
        oarchive(std::ostream& os);

        // construct iarchive from an existing archive; this
        // object cannot be used, when ar is destroyed
        explicit oarchive(archive_type& ar);

        // standard destructor
        ~oarchive();

        // get implementation type
        archive_type&   get()       { return *m_impl; };

        // save integer
        oarchive&       operator<<(Integer v);

        // save double float
        oarchive&       operator<<(Real v);

        // save single float
        oarchive&       operator<<(Float v);

        // save double complex
        oarchive&       operator<<(const Complex& v);

        // save single complex
        oarchive&       operator<<(const Float_complex& v);
};

};