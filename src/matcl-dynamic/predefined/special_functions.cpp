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

#include "special_functions.h"

namespace matcl { namespace dynamic { namespace special_functions
{

function_name special_functions::convert_numeric_promotion_1()
{
    static function_name f("$convert_numeric_promotion_1");
    return f;
};
function_name special_functions::convert_numeric_promotion_2()
{
    static function_name f("$convert_numeric_promotion_2");
    return f;
};
function_name special_functions::convert_numeric_equivalent_1()
{
    static function_name f("$convert_numeric_equivalent_1");
    return f;
};
function_name special_functions::convert_numeric_equivalent_2()
{
    static function_name f("$convert_numeric_equivalent_2");
    return f;
};
function_name special_functions::convert_numeric_decay_1()
{
    static function_name f("$convert_numeric_decay_1");
    return f;
};
function_name special_functions::convert_numeric_decay_2()
{
    static function_name f("$convert_numeric_decay_2");
    return f;
};
function_name special_functions::convert_numeric_int_float()
{
    static function_name f("$convert_numeric_int_float");
    return f;
};

function_name special_functions::convert_numeric_explicit()
{
    static function_name f("$convert_numeric_explicit");
    return f;
};

function_name special_functions::convert_numeric_cast()
{
    static function_name f("$convert_numeric_cast");
    return f;
};

function_name special_functions::convert_any()
{
    static function_name f("$convert_any");
    return f;
};

function_name special_functions::convert_unit()
{
    static function_name f("$convert_unit");
    return f;
};

function_name special_functions::convert_id()
{
    static function_name f("$convert_id");
    return f;
};

function_name special_functions::convert_promotion()
{
    static function_name f("$convert_promotion");
    return f;
};
function_name special_functions::convert_equivalent()
{
    static function_name f("$convert_equivalent");
    return f;
};

function_name special_functions::convert_decay()
{
    static function_name f("$convert_decay");
    return f;
};

function_name special_functions::convert_explicit()
{
    static function_name f("$convert_explicit");
    return f;
};

function_name special_functions::convert_cast()
{
    static function_name f("$convert_cast");
    return f;
};

function_name special_functions::unifier()
{
    static function_name f("$unifier");
    return f;
};

function_name special_functions::assign()
{
    static function_name f("$assign");
    return f;
}

function_name special_functions::assign_id()
{
    static function_name f("$assign_id");
    return f;
}

bool special_functions::is_special(function_name func)
{
    return func.is_special();

}

};};};

