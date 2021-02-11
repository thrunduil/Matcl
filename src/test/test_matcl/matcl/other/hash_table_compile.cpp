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

#include "matcl-core/details/hash_table/hash_table.h"

namespace matcl { namespace test
{

template<class Type>
struct case_hasher
{
    size_t operator()(const Type* ptr) const
    {
        return ptr->hash_value();
    };
};
template<class Type>
struct case_equaler
{
    bool operator()(const Type* ptr1, size_t hash_value,const Type* ptr2) const
    {
        (void)hash_value;
        return ptr1->equal(*ptr2);
    };
};

template<class Type>
using hash_type = matcl::hash_table<Type, case_hasher<Type>, case_equaler<Type>>;

class Test
{
    private:
        int m_val;

    public:
        Test(int val)                           :m_val(val){};
        size_t  hash_value() const              { return m_val; };
        bool    equal(const Test& other) const  { return m_val == other.m_val; };
};

void hash_table_test()
{
    // only check compilation errors
    hash_type<Test> ht(10);

    ht.insert(new Test(1));
    ht.insert(new Test(2));
};

}};