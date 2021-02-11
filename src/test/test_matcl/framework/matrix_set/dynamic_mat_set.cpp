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

#include "dynamic_mat_set.h"
#include "matcl-core/IO/logger.h"

#pragma warning( push )
#pragma warning(disable:4702)	//  unreachable code
    #include <boost/lexical_cast.hpp>
#pragma warning( pop )


namespace matcl { namespace test
{

static const Integer max_int = 1000;

const dynamic_mat_set::container& dynamic_mat_set::get(Integer s)
{
    scoped_lock lock(m_mutex);

    using iterator      = cont_mat::iterator;
    using value_type    = cont_mat::value_type;

    iterator pos = m_map.find(index(s,-1));

    if (pos == m_map.end())
    {
        container mc;
        insert(s,-1,mc);
        m_map.insert(value_type(index(s,-1),mc));
        pos = m_map.find(index(s,-1));
        return pos->second;
    }
    else
    {
        return pos->second;
    };
};

const dynamic_mat_set::container& dynamic_mat_set::get(Integer r,Integer c)
{
    scoped_lock lock(m_mutex);

    using iterator      = cont_mat::iterator;
    using value_type    = cont_mat::value_type;

    iterator pos = m_map.find(index(r,c));
    if (pos == m_map.end())
    {
        container mc;
        insert(r,c,mc);
        m_map.insert(value_type(index(r,c),mc));
        pos = m_map.find(index(r,c));
        return pos->second;
    }
    else
    {
        return pos->second;
    };
};

Matrix dynamic_mat_set::rand_dense(Integer r, Integer c,Integer seed)
{
    rand_state rs = get_rand_state();
    init_genrand(seed);

    Matrix ret  = rand_dense(r, c);
    ret         = full(ret);

    set_rand_state(rs);
    return ret;
};

Matrix dynamic_mat_set::rand(Integer r, Integer c,Integer seed)
{
    rand_state rs = get_rand_state();
    init_genrand(seed);
    
    Integer t = abs(irand())%3;
    Matrix out;
    switch (t)
    {
        case 0:
            out = rand_dense(r,c);
            break;
        case 1:
            out = rand_sparse(r,c);
            break;
        case 2:
        default:
            out = rand_band(r,c);
            break;
    };

    set_rand_state(rs);
    return out;
};

Matrix dynamic_mat_set::rand_dense(Integer m, Integer n)
{
    Integer t = abs(irand()) % 5;

    switch (t)
    {
        case 0: //integer
        {
            Matrix tmp = m_rand->rand_dense_int(m,n);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_dense_int(m+5,n+5);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 1: //real
        {
            Matrix tmp = m_rand->rand_dense_real(m,n);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_dense_real(m+5,n+5);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 2: //complex
        {
            Matrix tmp = m_rand->rand_dense_compl(m,n);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_dense_compl(m+5,n+5);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 3: //float
        {
            Matrix tmp = m_rand->rand_dense_float(m,n);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_dense_float(m+5,n+5);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 4: //float complex
        {
            Matrix tmp = m_rand->rand_dense_fcompl(m,n);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_dense_fcompl(m+5,n+5);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
    };

    return m_rand->rand_dense_int(m,n);
};

Matrix dynamic_mat_set::rand_sparse(Integer m, Integer n)
{
    Real d = .2;

    Integer t = abs(irand() ) % 5;
    
    switch (t)
    {
        case 0: //integer
        {
            Matrix tmp = m_rand->rand_sparse_int(m,n,d);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_sparse_int(m,n+4,d);
            tmp = tmp(colon(),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand()) % 3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 1: //real
        {
            Matrix tmp = m_rand->rand_sparse_real(m,n,d);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_sparse_real(m,n+4,d);
            tmp = tmp(colon(),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 2: //complex
        {
            Matrix tmp = m_rand->rand_sparse_compl(m,n,d);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_sparse_compl(m,n+4,d);
            tmp = tmp(colon(),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 3: //float
        {
            Matrix tmp = m_rand->rand_sparse_float(m,n,d);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_sparse_float(m,n+4,d);
            tmp = tmp(colon(),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
        case 4: //float complex
        {
            Matrix tmp = m_rand->rand_sparse_fcompl(m,n,d);
            Matrix m1 = tmp;

            tmp.reserve(m+2,n+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_sparse_fcompl(m,n+4,d);
            tmp = tmp(colon(),colon(2,2+n-1));
            Matrix m3 = tmp;

            switch (abs(irand())%3)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
            };
        }
    };

    return m_rand->rand_sparse_compl(m,n,d);
};

Matrix dynamic_mat_set::rand_band(Integer m, Integer n)
{
    Integer t = abs(irand()) % 5;
    Integer ld = 1;
    Integer ud = 1;
    switch (t)
    {
        case 0: //integer
        {
            Matrix tmp = m_rand->rand_band_int(m,n,-ld,ud);
            Matrix m1 = tmp;

            tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_band_int(m+5,n+5,-ld,ud);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            tmp = m_rand->rand_band_int(m+5,n+5,-(ld+2),ud+2);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
            Matrix m4 = tmp;

            switch (abs(irand())%4)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
                case 3: return m4;
            };
        }
        case 1: //real
        {
            Matrix tmp = m_rand->rand_band_real(m,n,-ld,ud);
            Matrix m1 = tmp;

            tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_band_real(m+5,n+5,-ld,ud);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            tmp = m_rand->rand_band_real(m+5,n+5,-(ld+2),ud+2);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
            Matrix m4 = tmp;

            switch (abs(irand())%4)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
                case 3: return m4;
            };
        }
        case 2: //complex
        {
            Matrix tmp = m_rand->rand_band_compl(m,n,-ld,ud);
            Matrix m1 = tmp;

            tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_band_compl(m+5,n+5,-ld,ud);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            tmp = m_rand->rand_band_compl(m+5,n+5,-(ld+2),ud+2);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
            Matrix m4 = tmp;

            switch (abs(irand())%4)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
                case 3: return m4;
            };
        }
        case 3: //float
        {
            Matrix tmp = m_rand->rand_band_float(m,n,-ld,ud);
            Matrix m1 = tmp;

            tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_band_float(m+5,n+5,-ld,ud);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            tmp = m_rand->rand_band_float(m+5,n+5,-(ld+2),ud+2);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
            Matrix m4 = tmp;

            switch (abs(irand())%4)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
                case 3: return m4;
            };
        }
        case 4: //float complex
        {
            Matrix tmp = m_rand->rand_band_fcompl(m,n,-ld,ud);
            Matrix m1 = tmp;

            tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
            Matrix m2 = tmp;

            tmp = m_rand->rand_band_fcompl(m+5,n+5,-ld,ud);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            Matrix m3 = tmp;

            tmp = m_rand->rand_band_fcompl(m+5,n+5,-(ld+2),ud+2);
            tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
            tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
            Matrix m4 = tmp;

            switch (abs(irand())%4)
            {
                case 0: return m1;
                case 1: return m2;
                case 2: return m3;
                case 3: return m4;
            };
        }
    };

    return m_rand->rand_band_compl(m,n,-ld,ud);
};

void dynamic_mat_set::insert(Integer m, Integer n, container& mc)
{	
    if (n == -1)
    {
        insert(m,1,mc);
        insert(1,m,mc);

        if (m % 2 == 0 && m > 0)
        {
            insert(m/2,2,mc);
            insert(2,m/2,mc);
        };

        if (m % 3 == 0 && m > 0)
        {
            insert(m / 3, 3,mc);
            insert(3, m / 3,mc);
        };
        
        return;
    };

    if (m == 1 && n == 1)
    {
        mc.push_back(m_rand->rand_scalar_int());
        mc.push_back(m_rand->rand_scalar_real());
        mc.push_back(m_rand->rand_scalar_float());
        mc.push_back(m_rand->rand_scalar_compl());
        mc.push_back(m_rand->rand_scalar_fcompl());
    };

    add_matrices_dense(mc,m,n);
    add_matrices_sparse(mc,m,n,.2);
    add_matrices_band(mc,m,n,0,0);
    add_matrices_band(mc,m,n,1,1);
};

void dynamic_mat_set::add_matrices_dense(container& mc, Integer m, Integer n)
{
    //int
    {
        Matrix tmp = m_rand->rand_dense_int(m,n);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_int(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_dense_real(m,n);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_real(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };
    //complex
    {
        Matrix tmp = m_rand->rand_dense_compl(m,n);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_compl(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_dense_float(m,n);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_float(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };
    //float complex
    {
        Matrix tmp = m_rand->rand_dense_fcompl(m,n);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_dense_fcompl(m+5,n+5);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };
};

void dynamic_mat_set::add_matrices_sparse(container& mc, Integer m, Integer n, Real d)
{
    //int
    {
        Matrix tmp = m_rand->rand_sparse_int(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_int(m,n+4,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //real
    {
        Matrix tmp = m_rand->rand_sparse_real(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_real(m,n+4,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_sparse_compl(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_compl(m,n+4,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_sparse_float(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_float(m,n+4,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_sparse_fcompl(m,n,d);
        Matrix m1 = tmp;

        tmp.reserve(m+2,n+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_sparse_fcompl(m,n+4,d);
        tmp = tmp(colon(),colon(2,2+n-1));
        Matrix m3 = tmp;

        switch (abs(irand())%3)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
        };
    };
};

void dynamic_mat_set::add_matrices_band(container& mc, Integer m, Integer n, Integer ld, Integer ud)
{
    Integer diag_l = (m == 0)? 0 : 1;
    Integer diag_u = (n == 0)? 0 : 1;

    if (ld + diag_l > m)
        return;

    if (ud + diag_u > n)
        return;

    //real
    {
        Matrix tmp = m_rand->rand_band_real(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_real(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_real(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //complex
    {
        Matrix tmp = m_rand->rand_band_compl(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_compl(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_compl(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //float
    {
        Matrix tmp = m_rand->rand_band_float(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_float(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_float(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };

    //float complex
    {
        Matrix tmp = m_rand->rand_band_fcompl(m,n,-ld,ud);
        Matrix m1 = tmp;

        tmp.reserve_band(m+2,n+2,-(ld+2),ud+2);
        Matrix m2 = tmp;

        tmp = m_rand->rand_band_fcompl(m+5,n+5,-ld,ud);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        Matrix m3 = tmp;

        tmp = m_rand->rand_band_fcompl(m+5,n+5,-(ld+2),ud+2);
        tmp = tmp(colon(2,2+m-1),colon(2,2+n-1));
        tmp.resize_band(tmp.rows(), tmp.cols(),-ld,ud);
        Matrix m4 = tmp;

        switch (abs(irand())%4)
        {
            case 0: mc.push_back(m1); break;
            case 1: mc.push_back(m2); break;
            case 2: mc.push_back(m3); break;
            case 3: mc.push_back(m4); break;
        };
    };
};

void dynamic_mat_set::clear()
{
    scoped_lock lock(m_mutex);
    m_map.clear();
};

};};

