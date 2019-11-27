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

#include <vector>
#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace test
{

class colon_set
{
    private:
        std::vector<colon>    m_colons;

        void    insert(Integer max);

    public:
        colon_set(Integer num, Integer max, Integer seed);

        int     size() const;
        colon   get(int pos) const;

    private:
        static Integer rand_int(Integer min, Integer max);
};

};};