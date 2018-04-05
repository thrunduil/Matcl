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

#include <string>
#include <vector>
#include <memory>

#include "matcl-sqlite-cpp/config.h"
#include "matcl-sqlite-cpp/builder/enums.h"
#include "matcl-sqlite-cpp/builder/details/fwd_decls.h"
#include "matcl-sqlite-cpp/builder/column.h"
#include "matcl-sqlite-cpp/builder/bind.h"
#include "matcl-sqlite-cpp/builder/comma_list.h"

namespace matcl { namespace sql
{

struct null_t {}; // SQL null value
extern null_t null;

class MATCL_SQLITE_EXPORT_MACRO expression
{
    public:
        expression();
        expression(const column& c);
        expression(const bind_param& b);
        expression(double d);
        expression(const std::string& s);
        expression(const char* c);
        expression(const null_t);

        std::string to_str() const;

        // takes ownership
        explicit expression(details::expression_base_impl* impl);
    
        bool exists() const;

    private:
        using impl_ptr = std::shared_ptr<details::expression_base_impl>;
        impl_ptr m_impl;
};

// unary

MATCL_SQLITE_EXPORT_MACRO expression operator-(const expression& e);
MATCL_SQLITE_EXPORT_MACRO expression operator+(const expression& e);
MATCL_SQLITE_EXPORT_MACRO expression operator~(const expression& e);
MATCL_SQLITE_EXPORT_MACRO expression operator!(const expression& e);

// binary

MATCL_SQLITE_EXPORT_MACRO expression concat    (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator* (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator/ (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator% (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator+ (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator- (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator<<(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator>>(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator& (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator| (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator< (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator<=(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator> (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator>=(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator==(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator!=(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression is        (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression is_not    (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression in        (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression like      (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression glob      (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression regexp    (const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator&&(const expression& e1, const expression& e2);
MATCL_SQLITE_EXPORT_MACRO expression operator||(const expression& e1, const expression& e2);

// other

using expression_list = comma_list<expression>;

MATCL_SQLITE_EXPORT_MACRO 
expression call(const std::string& func, wildcard);

MATCL_SQLITE_EXPORT_MACRO 
expression call(const std::string& func);

MATCL_SQLITE_EXPORT_MACRO 
expression call(const std::string& func, const expression& e1); // 1 arg...

MATCL_SQLITE_EXPORT_MACRO 
expression call(const std::string& func, const expression& e1, const expression& e2); // 2 args...

MATCL_SQLITE_EXPORT_MACRO 
expression call(const std::string& func, const expression& e1, const expression& e2, const expression& e3); // 3 args...

MATCL_SQLITE_EXPORT_MACRO 
expression call(const std::string& func, const expression_list& e, bool distinct=false); // any number of args

MATCL_SQLITE_EXPORT_MACRO expression collate(const expression& e, collation c);

}}

