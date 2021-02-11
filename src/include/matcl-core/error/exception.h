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

#pragma once

#include "matcl-core/config.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/details/exception_details.h"

#include <exception>
#include <string>

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 
#pragma warning(disable:4275)	//non dll-interface class used as base for dll-interface class

namespace matcl { namespace error
{

// base class for exceptions throws in matcl
class MATCL_CORE_EXPORT matcl_exception : public std::exception
{
    protected:
        mutable std::string m_message;

    public:
        virtual const char* what() const throw() override;

    private:
        virtual const char* what(exception_message& em) const = 0;
};

// asserts fauilures are converted to exceptions not derived from std::exception
// if this exception is caught; then program should not continue further
class MATCL_CORE_EXPORT assert_exception
{
    private:
        std::string file;
        int         line;
        std::string message;

    public:
        // construct assert message
        assert_exception(const std::string& file, int line, const std::string& message);

        // construct assert message and throw exception
        static void make(const char* txt, const char* description, const char* file, int line);

        virtual std::string what() const throw();
};

// enable or disable all warnings in current scope, restore previous settings
// after exiting from current scope; if global is true, then warnings are
// enabled/disabled in all threads, otherwise only in current thread; global
// setting has precedence over local setting, i.e. if warnings are disabled 
// globally but enabled locally, then no warnings are displayed;
MATCL_CORE_EXPORT details::enable_warnings_raii
                        enable_warnings(bool enable, bool global = false);

// check if warnings are enabled
MATCL_CORE_EXPORT bool check_warnings_enabled();

// print warning (if warnings are enabled)
MATCL_CORE_EXPORT void warning(const std::string& msg);

};};

// matcl assert
#define matcl_assert(cond, description)  \
    ((cond) ? (void)0 : matcl::error::assert_exception::make(#cond, description, __FILE__, __LINE__))

#pragma warning(pop)