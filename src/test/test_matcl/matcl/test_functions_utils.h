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

#pragma once

#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

class utils_functions_list
{
    public:
        void	make();

        void operator()()
        {
            make();
        };

    private:
        void        test_thread();
        void        test_tuple();

    public:
        utils_functions_list();
};

class test_function_thread
{
    public:
        test_function_thread()	{};

        Real make();
};

class test_function_tuple
{
    public:
        test_function_tuple()	{};

        Real make();
};

};};
