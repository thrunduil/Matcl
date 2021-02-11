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

#include "matcl-linalg/linear_eq/linsolve_objects_ext.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/linalg_exception.h"

namespace matcl { namespace details
{

linsolve_extended_sol::linsolve_extended_sol(const linsolve_obj& lo, const options& opts)
    :linsolve_forwarding(lo), m_opts(opts)
{
    init_options(opts);
};

linsolve_extended_sol::~linsolve_extended_sol()
{};

void linsolve_extended_sol::init_options(const options& opts)
{
    m_do_refinement = false;
    m_test_sol      = 0;

    namespace ol = opt::linsolve;

    //iterative refinement
    if (base().is_direct() == true && base().is_modified() == false)
    {
        bool ir         = opts.get_option<bool>(ol::use_ir());
        m_do_refinement = ir;
    }

    //test sol
    m_test_sol  = opts.get_option<Integer>(ol::test_sol());
};

linsolve_extended_sol::data_ptr linsolve_extended_sol::convert(value_code nvc) const
{
    return data_ptr(new linsolve_extended_sol(base().convert(nvc), m_opts));
}

Matrix linsolve_extended_sol::solve(const Matrix& X, trans_type tA) const
{
    Matrix sol  = base().solve(X, tA);

    if (m_do_refinement)
        sol     = base().iterative_refinement(std::move(sol), X, tA, m_opts);

    if (m_test_sol != 0)
    {
        Matrix be, fe;
        tie(be, fe) = base().comp_error(sol, X, tA);
        check_solution(fe, m_test_sol);
    };

    return sol;
};
Matrix linsolve_extended_sol::solve(Matrix&& X, trans_type tA) const
{
    if (m_do_refinement == false && m_test_sol == 0)
        return base().solve(std::move(X), tA);

    return solve(X, tA);
};

Matrix linsolve_extended_sol::solve_rev(const Matrix& X, trans_type tA) const
{
    Matrix sol  = base().solve_rev(X, tA);

    if (m_do_refinement)
        sol     = base().iterative_refinement_rev(std::move(sol), X, tA, m_opts);

    if (m_test_sol != 0)
    {
        Matrix be, fe;
        tie(be, fe) = base().comp_error_rev(sol, X, tA);
        check_solution(fe, m_test_sol);
    };

    return sol;
};

Matrix linsolve_extended_sol::solve_rev(Matrix&& X, trans_type tA) const
{
    if (m_do_refinement == false && m_test_sol == 0)
        return base().solve_rev(std::move(X), tA);

    return solve_rev(X, tA);
};

void linsolve_extended_sol::check_solution(const Matrix& fe, Integer prec) const
{
    Matrix fe_log   = log10(fe);
    Integer prec2   = abs(prec);
    Matrix I        = find(fe_log > -prec2);

    if (I.length() == 0)
        return;

    Integer k       = I(1).get_scalar<Integer>();
    Real fe_k       = fe(k).get_scalar<Real>();

    if (prec < 0)
        throw error::inaccurate_solution(k, fe_k, prec2);
    else
        error::get_global_messanger_linalg()->warning_inaccurate_solution(k, fe_k, prec2);
};

};};