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

#include "colon_set.h"
#include "test_options.h"

namespace matcl { namespace test
{

colon_set::colon_set(Integer num, Integer max, Integer seed)
{
    auto st = local_rand_state(seed + test_options::get_seed());

    for (Integer i = 0; i < num; ++i)
        insert(max);
};

int colon_set::size() const
{
    return (int)m_colons.size();
};

colon colon_set::get(int pos) const
{
    return m_colons[pos].copy();
};

void colon_set::insert(Integer max_val)
{
    int type    = abs(irand()) % 6;

    switch (type)
    {
        case 0:
        default:
        {
            m_colons.push_back(colon());
            return;
        }
        case 1:
        {
            int val = rand_int(1, max_val);
            m_colons.push_back(colon(val));
            return;
        }
        case 2:
        {
            int val1 = rand_int(1, max_val);
            int val2 = rand_int(1, max_val);
            m_colons.push_back(colon(val1,val2));
            return;
        }            
        case 3:
        {
            int val1 = rand_int(1, max_val);
            int val2 = rand_int(-max_val, max_val);
            int val3 = rand_int(1, max_val);
            m_colons.push_back(colon(val1,val2,val3));
            return;
        }
        case 4:
        {
            int len = rand_int(0,5);
            mat_row mr;

            for (int i = 0; i < len; ++i)
            {
                int val = rand_int(1, max_val);
                mr.add(val);
            };
            m_colons.push_back(colon(mr));
            return;
        }
        case 5:
        {
            int type2   = abs(irand()) % 7;

            switch(type2)
            {
                case 0:
                default:
                    m_colons.push_back(colon(end, -1, 1));
                    return;
                case 1:
                    m_colons.push_back(colon(end - 1, -1, 1));
                    return;
                case 2:
                    m_colons.push_back(colon(end, end));
                    return;
                case 3:
                    m_colons.push_back(colon(end, -1, end - 1));
                    return;
                case 4:
                    m_colons.push_back(colon(end - 2, end));
                    return;
                case 5:
                    m_colons.push_back(colon(end - 2, end - 1));
                    return;
                case 6:
                    m_colons.push_back(colon(end - 2, -1, 2));
                    return;
            }
            return;
        }
        case 6:
        {
            // empty colons - valid
            int type2   = abs(irand()) % 5;
            switch(type2)
            {
                case 0:
                default:
                    m_colons.push_back(colon(1, 0));
                    return;
                case 1:
                    m_colons.push_back(colon(1, 2, 0));
                    return;
                case 2:
                    m_colons.push_back(colon(0, -1, 1));
                    return;
                case 3:
                    m_colons.push_back(colon(end, end-1));
                    return;
                case 4:
                    m_colons.push_back(colon(end-1, -1, end));
                    return;
            }
            return;
        }
    };
};

Integer colon_set::rand_int(Integer min, Integer max)
{
    if (min > max)
        return min;

    Integer val = min + abs(irand()) % (max - min + 1);
    return val;    
};

};};