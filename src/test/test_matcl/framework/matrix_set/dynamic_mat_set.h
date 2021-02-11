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

#include "test/test_matcl/framework/test_functions/unary_function.h"
#include "test/test_matcl/framework/test_functions/bin_function.h"
#include "test/test_matcl/framework/data_struct/options.h"
#include "test/test_matcl/framework/matrix_set/matrix_rand.h"

#include <vector>
#include <map>

namespace matcl { namespace test
{

class dynamic_mat_set
{
    public:
        using container     = std::vector<Matrix>;
        using index         = std::pair<Integer,Integer>;
        using cont_mat      = std::map<index,container>;
        using scoped_lock   = std::unique_lock<std::mutex>;

    private:
        cont_mat		m_map;
        std::mutex	    m_mutex;
        rand_matrix_ptr m_rand;

        void            insert(Integer r, Integer c,container& mc);
        void            add_matrices_band(container& mat, Integer m, Integer n, Integer ld, Integer ud);
        void            add_matrices_dense(container& mat, Integer m, Integer n);
        void            add_matrices_sparse(container& mat, Integer m, Integer n, Real d);
        Matrix          rand_dense(Integer r, Integer c);
        Matrix          rand_sparse(Integer r, Integer c);
        Matrix          rand_band(Integer r, Integer c);

    public:
        dynamic_mat_set(rand_matrix_ptr rand)   : m_rand(rand) {};

        const container&	get(Integer s);
        const container&	get(Integer r,Integer c);
        Matrix              rand(Integer r, Integer c, Integer seed);
        Matrix              rand_dense(Integer r, Integer c, Integer seed);
        void				clear();
};

};};