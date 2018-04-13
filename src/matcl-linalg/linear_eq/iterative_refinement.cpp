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

#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "iterative_refinement.h"
#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace details
{

void iterative_refinement::new_x_state(x_state& xs, Real norm_x, Real norm_dx_i, Real norm_dx_ip, 
                        Real rho_thresh, Real eps_w)
{
    Real prog   = div_1(norm_dx_ip,norm_dx_i);
    Real prog2  = div_0(norm_dx_ip,norm_x);

    if (norm_dx_i == constants::inf() && norm_dx_ip == constants::inf())
        prog    = 1.0;
    if (norm_dx_ip == constants::inf() && norm_x == constants::inf())
        prog2   = 1.0;

    if (xs == x_state::no_progress && prog <= rho_thresh)
	    xs  = x_state::working;

    if (xs == x_state::working)
    {
	    if (prog2 <= eps_w)
		    xs = x_state::converged;
        else if (prog > rho_thresh)
            xs = x_state::no_progress;
    };
}
    
void iterative_refinement::new_z_state(z_state& zs, Real dz_ip, Real dz_i, Real dz_thresh, 
                                       Real rho_thresh, Real eps_w)
{
    Real prog   = div_1(dz_ip,dz_i);

    if (dz_ip == constants::inf() && dz_i == constants::inf())
        prog    = 1.0;

    if (zs == z_state::unstable && dz_ip <= dz_thresh)
	    zs = z_state::working;
    if (zs == z_state::no_progress && prog <= rho_thresh)
	    zs = z_state::working;

    if (zs == z_state::working)
    {
	    if (dz_ip <= eps_w)
		    zs  = z_state::converged;
	    else if (dz_ip > dz_thresh)
		    zs = z_state::unstable;
	    else if (prog > rho_thresh )
		    zs = z_state::no_progress;
    };
};

iterative_refinement::iterative_refinement(const Matrix& X, trans_type tA, const linsolve_obj& lo,
                                           const options& opts)
    :m_X(X), m_tA(tA), m_lo(lo)
{
    n_iterations    = opts.get_option<Integer>(opt::linsolve::refinement_iter());
    m_rho_threshhold= opts.get_option<Real>(opt::linsolve::refinement_rho());
};

void iterative_refinement::eval(Matrix& Y_sol, bool rev)
{
    if (rev == false && m_X.rows() != m_lo.cols())
        throw error::invalid_size2(m_X.rows(), m_X.cols(), m_lo.rows(), m_X.cols());
    if (rev == true && m_X.cols() != m_lo.rows())
        throw error::invalid_size2(m_X.rows(), m_X.cols(), m_X.rows(), m_lo.rows());

    if (Y_sol.rows() != m_X.rows() || Y_sol.cols() != m_X.cols())
        throw error::invalid_size2(Y_sol.rows(), Y_sol.cols(), m_X.rows(), m_X.cols());

    if (n_iterations <= 0)
        return;

    if (rev == false && Y_sol.cols() == 1)
        return eval_1(Y_sol, rev);
    else if (rev == true && Y_sol.rows() == 1)
        return eval_1(Y_sol, rev);

    if (Y_sol.get_struct_code() == struct_code::struct_sparse)
        return eval_sparse(Y_sol, rev);
    else
        return eval_dense(Y_sol, rev);
};

Integer iterative_refinement::get_i_thresh()
{
    return n_iterations;
};
Real iterative_refinement::get_rho_thresh()
{
    //const Real rho_thresh   = 0.5;  //0.9
    return m_rho_threshhold;
};
Real iterative_refinement::get_dz_thresh()
{
    return 0.25;
};

void iterative_refinement::eval_1(Matrix& Y_sol, bool rev)
{
    Real eps            = constants::eps(Y_sol.get_value_code());

    Integer i_thresh    = get_i_thresh();
    Real rho_thresh     = get_rho_thresh();    
    Real dz_thresh      = get_dz_thresh();
    Real eps_w          = eps;

    x_state xs          = x_state::working;
    z_state zs          = z_state::unstable;        
    Real norm_dz_i      = constants::inf();
    Real norm_dx_i      = constants::inf();

    basic_vector_norm p = (rev == false) ? basic_vector_norm::norm_inf
                        : basic_vector_norm::norm_1;

    for (Integer i = 1; i <= i_thresh; ++i)
    {
        // Compute residual
        Matrix ri       = rev == false ? m_lo.resid(Y_sol, m_X, m_tA) 
                                       : m_lo.resid_rev(Y_sol, m_X, m_tA);

		// Compute correction to Y
        Matrix dY       = rev == false ? m_lo.solve(std::move(ri), m_tA)
                                       : m_lo.solve_rev(std::move(ri), m_tA);

		// Check error-related stopping criteria
        Real norm_x     = norm(Y_sol, p);
        Real norm_dx_ip = norm(dY, p);
        Real norm_dz_ip = norm(div_0(dY, Y_sol),p);

        new_x_state(xs, norm_x, norm_dx_i, norm_dx_ip, rho_thresh, eps_w);
        new_z_state(zs, norm_dz_ip, norm_dz_i, dz_thresh, rho_thresh, eps_w);

		// Update solution
		Y_sol           = std::move(Y_sol) - std::move(dY);

        if ((xs == x_state::working || zs == z_state::working) == false)
		    break;

        // update norms
        norm_dz_i       = norm_dz_ip;
        norm_dx_i       = norm_dx_ip;
    };

    return;
};

void iterative_refinement::initialize_state(Integer N)
{
    m_state.resize(N);
};

it_ref_state::it_ref_state()
{
    m_iter      = 0;
    m_finished  = false;
    m_norm_dz_i = constants::inf();
    m_norm_dx_i = constants::inf();
    m_xs        = x_state::working;
    m_zs        = z_state::unstable;
};

void iterative_refinement::get_working_columns(Integer done, Integer& f, Integer& l) const
{
    f           = done + 1;
    Integer N   = (int)m_state.size();

    for (Integer i = f; i <= N; ++i)
    {
        if (m_state[i-1].m_finished == true)
            ++f;
        else
            break;
    };

    l           = f;

    for (Integer i = f+1; i <= N; ++i)
    {
        if (m_state[i-1].m_finished == true)
            break;
        else
            ++l;
    };

    l           = std::min(l, N);
};

void iterative_refinement::new_x_state(state_vec& state, Integer f, Integer l, const Matrix& norm_x, 
                                       const Matrix& norm_dx_ip, Real rho_thresh, Real eps_w)
{
    const Real* ptr_norm_x      = norm_x.get_array<Real>();
    const Real* ptr_norm_dx_ip  = norm_dx_ip.get_array<Real>();

    for (Integer i = f-1, pos = 0; i < l; ++i, ++pos)
    {
        new_x_state(state[i].m_xs, ptr_norm_x[pos], state[i].m_norm_dx_i, ptr_norm_dx_ip[pos], rho_thresh, eps_w);
        state[i].m_norm_dx_i    = ptr_norm_dx_ip[pos];
    };
};

void iterative_refinement::new_z_state(state_vec& state, Integer f, Integer l, const Matrix& norm_dz_ip,
                                       Real dz_thresh, Real rho_thresh, Real eps_w)
{
    const Real* ptr_norm_dz_ip  = norm_dz_ip.get_array<Real>();

    for (Integer i = f-1, pos = 0; i < l; ++i, ++pos)
    {
        new_z_state(state[i].m_zs, ptr_norm_dz_ip[pos], state[i].m_norm_dz_i, dz_thresh, rho_thresh, eps_w);
        state[i].m_norm_dz_i = ptr_norm_dz_ip[pos];
    };
};

void iterative_refinement::update_state(state_vec& state, Integer max_it, Integer f, Integer l, Integer& done)
{
    for (Integer i = f-1, pos = 0; i < l; ++i, ++pos)
    {
        state[i].m_iter += 1;

        if ((state[i].m_xs == x_state::working || state[i].m_zs == z_state::working) == false)
		    state[i].m_finished = true;
        else if (state[i].m_iter >= max_it)
            state[i].m_finished = true;
    };

    done    = f - 1;

    for (Integer i = f-1, pos = 0; i < l; ++i, ++pos)
    {
        if (state[i].m_finished == true)
            ++done;
        else
            break;
    };
};

void iterative_refinement::eval_dense(Matrix& Y_sol, bool rev)
{
    basic_vector_norm p = basic_vector_norm::norm_inf;
    Integer d           = (rev == false) ? 1 : 2;
    Integer N           = (rev == false) ? Y_sol.cols() : Y_sol.rows();
    Real eps            = constants::eps(Y_sol.get_value_code());

    initialize_state(N);

    Real rho_thresh     = get_rho_thresh();    
    Real dz_thresh      = get_dz_thresh();
    Integer max_it      = get_i_thresh();
    Real eps_w          = eps;
    Integer done        = 0;

    for (;;)
    {
        //get working vectors
        Integer f, l;
        get_working_columns(done, f, l);

        if (f > l)
            break;

        Matrix X            = (rev == false) ? m_X(colon(), colon(f,l)) 
                                             : m_X(colon(f,l), colon());
        Matrix Y            = (rev == false) ? Y_sol(colon(), colon(f,l)) 
                                             : Y_sol(colon(f,l), colon());

        // Compute residual
        Matrix ri           = (rev == false) ? m_lo.resid(Y, X, m_tA) 
                                             : m_lo.resid_rev(Y, X, m_tA);

		// Compute correction to Y
        Matrix dY           = (rev == false) ? m_lo.solve(std::move(ri), m_tA) 
                                             : m_lo.solve_rev(std::move(ri), m_tA);

		// Check error-related stopping criteria
        Matrix norm_x       = norm_vec(Y, p, d);
        Matrix norm_dx_ip   = norm_vec(dY, p, d);
        Matrix norm_dz_ip   = norm_vec(div_0(dY, Y),p,d);

        new_x_state(m_state, f, l, norm_x, norm_dx_ip, rho_thresh, eps_w);
        new_z_state(m_state, f, l, norm_dz_ip, dz_thresh, rho_thresh, eps_w);
        update_state(m_state, max_it, f,l, done);

        // Update solution
        if (rev == false)
        {
            Y_sol(colon(), colon(f,l))
                            = std::move(Y) - std::move(dY);
        }
        else
        {
            Y_sol(colon(f,l), colon())
                            = std::move(Y) - std::move(dY);
        };
    };

    return;
};

final_sparse_sol::final_sparse_sol(bool rev)
    :m_rev(rev)
{};

void final_sparse_sol::add(const Matrix& Y)
{
    if (m_rev == false)
        m_mr.add(Y);
    else
        m_mc, Y;
};

Matrix final_sparse_sol::make() const
{
    if (m_rev == false)
        return m_mr.to_matrix();
    else
        return m_mc.to_matrix();
};


void iterative_refinement::eval_sparse(Matrix& Y_sol, bool rev)
{
    Integer N           = (rev == false) ? Y_sol.cols() : Y_sol.rows();

    initialize_state(N);

    final_sparse_sol fss(rev);
    eval_sparse(Y_sol, rev, 1, N, fss);

    Y_sol   = fss.make();
};

void iterative_refinement::eval_sparse(Matrix& Y, bool rev, Integer fv, Integer lv, 
                                       final_sparse_sol& fss)
{    
    basic_vector_norm p = basic_vector_norm::norm_inf;
    Integer d           = (rev == false) ? 1 : 2;
    Real eps            = constants::eps(Y.get_value_code());

    Real rho_thresh     = get_rho_thresh();    
    Real dz_thresh      = get_dz_thresh();
    Integer max_it      = get_i_thresh();
    Real eps_w          = eps;
    Integer done        = fv - 1;

    for (;;)
    {
        //get working vectors
        Integer f, l;
        get_working_columns(done, f, l);

        if (f > l || f > lv)
        {
            fss.add(Y);
            break;
        }

        if (f > fv)
        {
            Integer n_vec   = f - fv;
            Matrix Yf       = (rev == false) ? Y(colon(), colon(1,n_vec)) 
                                             : Y(colon(1,n_vec), colon());
            fss.add(std::move(Yf));

            Integer last    = l - fv + 1;
            Yf              = (rev == false) ? Y(colon(), colon(n_vec+1,last)) 
                                             : Y(colon(n_vec+1,last), colon());
            eval_sparse(Yf, rev, f, l, fss);
            
            Y               = (rev == false) ? Y(colon(), colon(last+1,end)) 
                                             : Y(colon(last+1,end), colon());
            fv              = l + 1;
            done            = l;
            continue;
        }
        else if (l < lv)
        {
            Integer last    = l - fv + 1;
            Matrix Yf       = (rev == false) ? Y(colon(), colon(1,last)) 
                                             : Y(colon(1,last), colon());
            eval_sparse(Yf, rev, f, l, fss);            
            Y               = (rev == false) ? Y(colon(), colon(last+1,end)) 
                                             : Y(colon(last+1,end), colon());
            fv              = l + 1;
            done            = l;
            continue;
        };

        Matrix X            = (rev == false) ? m_X(colon(), colon(f,l)) 
                                             : m_X(colon(f,l), colon());

        // Compute residual
        Matrix ri           = (rev == false) ? m_lo.resid(Y, X, m_tA) 
                                             : m_lo.resid_rev(Y, X, m_tA);

		// Compute correction to Y
        Matrix dY           = (rev == false) ? m_lo.solve(std::move(ri), m_tA) 
                                             : m_lo.solve_rev(std::move(ri), m_tA);

		// Check error-related stopping criteria
        Matrix norm_x       = norm_vec(Y, p, d);
        Matrix norm_dx_ip   = norm_vec(dY, p, d);
        Matrix norm_dz_ip   = norm_vec(div_0(dY, Y),p,d);

        new_x_state(m_state, f, l, norm_x, norm_dx_ip, rho_thresh, eps_w);
        new_z_state(m_state, f, l, norm_dz_ip, dz_thresh, rho_thresh, eps_w);
        update_state(m_state, max_it, f,l, done);

        // Update solution
        if (rev == false)
            Y               = std::move(Y) - std::move(dY);
        else
            Y               = std::move(Y) - std::move(dY);
    };

    return;
};

};};