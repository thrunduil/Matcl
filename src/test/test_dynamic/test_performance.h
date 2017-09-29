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

#include "matcl-dynamic/matcl_dynamic.h"

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

void test_performance();

class performance_tester
{
    public:
        void    make();

    private:
        template<class T>
        double  test_mult(T& res, Integer M, Integer K, Integer N, Integer prec);

        template<class T>
        double  test_mult_obj(T& res, Integer M, Integer K, Integer N, Integer prec);

        template<class T>
        void    init(T* data, Integer size, Integer prec);

        template<class T, class TR = T>
        void    mult(const T* A, const T* B, T* C, Integer M, Integer K, Integer N);

        void    disp_res(const std::string& type, double time, double otime, double otime2);
};

}};
