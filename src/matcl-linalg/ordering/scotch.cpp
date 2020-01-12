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

#include "matcl-linalg/graph/scotch.h"
#include "matcl-linalg/ordering/scotch_wrapper.h"

namespace matcl
{

//-----------------------------------------------------------------------
//                          scotch
//-----------------------------------------------------------------------
scotch::scotch()
    :m_impl(new details::scotch_graph_impl(1.0, false))
{}

scotch::scotch(const Matrix& A, bool use_weights)
    :m_impl(new details::scotch_graph_impl(A, use_weights))
{}

scotch::~scotch()
{};

void scotch::set_strategy(scotch_strategy strategy)
{
    m_impl->set_strategy(strategy);
};

void scotch::set_seed(Integer s)
{
    m_impl->set_seed(s);
};

Integer scotch::get_seed() const
{
    return m_impl->get_seed();
};

void scotch::set_node_weights(const Matrix& W)
{
    m_impl->set_node_weights(W);
};

Matrix scotch::assign(Integer n_part, Real load_balance)
{
    return m_impl->make_assign(n_part, load_balance);
};

Matrix scotch::assign(Integer n_part, const Matrix& work_share, Real load_balance)
{
    return m_impl->make_assign(n_part, work_share, load_balance);
};

Matrix scotch::coloring()
{
    return m_impl->make_coloring();
};

Matrix scotch::edge_separator(Integer n_part, Real lb)
{
    return m_impl->edge_separator(n_part, lb);
};

Matrix scotch::clustering(Integer n_part, Real lb)
{
    return m_impl->make_clustering(n_part, lb);
};

Matrix scotch::node_separator(Integer n_part, Real lb)
{
    return m_impl->node_separator(n_part, lb);
};

Matrix scotch::assign_fixed(Integer n_part, const Matrix& fixed, Real load_balance)
{
    return m_impl->make_assign_fixed(n_part, fixed, load_balance);
};

Matrix scotch::assign_fixed(Integer n_part, const Matrix& work_shares, 
                    const Matrix& fixed, Real load_balance)
{
    return m_impl->make_assign_fixed(n_part, work_shares, fixed, load_balance);
};

Matrix scotch::edge_separator_fixed(Integer n_part, const Matrix& fixed, Real lb)
{
    return m_impl->edge_separator_fixed(n_part, fixed, lb);
};

permvec scotch::ordering()
{
    permvec p;
    Integer nblocks;
    Matrix  blocks, tree;
    m_impl->make_ordering(false, p, nblocks, blocks, tree);
    return p;
};

scotch::order_ext_ret scotch::ordering_ext()
{
    permvec p;
    Integer nblocks;
    Matrix  blocks, tree;

    m_impl->make_ordering(true, p, nblocks, blocks, tree);
    return order_ext_ret(p, nblocks, blocks, tree);
};

//-----------------------------------------------------------------------
//                          scotch_mesh
//-----------------------------------------------------------------------
scotch_mesh::scotch_mesh()
    :m_impl(new details::scotch_mesh_impl(1.0))
{}

scotch_mesh::scotch_mesh(const Matrix& A)
    :m_impl(new details::scotch_mesh_impl(A))
{}

scotch_mesh::~scotch_mesh()
{};

void scotch_mesh::set_seed(Integer s)
{
    m_impl->set_seed(s);
};

Integer scotch_mesh::get_seed() const
{
    return m_impl->get_seed();
};

void scotch_mesh::ordering()
{
    m_impl->make_ordering();
    return;
};

};