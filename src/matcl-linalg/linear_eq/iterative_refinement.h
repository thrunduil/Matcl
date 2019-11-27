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
#include "matcl-matrep/matrix/matrix_concat.h"

namespace matcl { namespace details
{

enum class x_state    { working, no_progress, converged };
enum class z_state    { unstable, working, no_progress, converged };

struct it_ref_state
{
    Integer     m_iter;
    bool        m_finished;
    Real        m_norm_dz_i;
    Real        m_norm_dx_i;
    x_state     m_xs;
    z_state     m_zs;

    it_ref_state();
};

class final_sparse_sol
{
    private:
        bool    m_rev;
        mat_col m_mc;
        mat_row m_mr;

    public:
        final_sparse_sol(bool rev);

        void    add(const Matrix& Y);
        Matrix  make() const;
};

//algorith from "Error Bounds from Extra Precise Iterative Refinement",
//Demmel et.al., 2004.
class iterative_refinement
{
    private:
        using state_vec     = std::vector<it_ref_state>;

    private:
        const Matrix&       m_X;
        trans_type          m_tA;
        const linsolve_obj& m_lo;
        state_vec           m_state;
        Integer             n_iterations;
        Real                m_rho_threshhold;

    private:
        static void     new_x_state(x_state& xs, Real norm_x, Real norm_dx_i, Real norm_dx_ip, 
                                Real rho_thresh, Real eps_w);
    
        static void     new_z_state(z_state& zs, Real dz_ip, Real dz_i, Real dz_thresh, Real rho_thresh, Real eps_w);

        Integer         get_i_thresh();
        Real            get_rho_thresh();    
        static Real     get_dz_thresh();

        void            initialize_state(Integer N);
        void            get_working_columns(Integer done, Integer& f, Integer& l) const;
        static void     new_x_state(state_vec& state, Integer f, Integer l, const Matrix& norm_x, 
                                    const Matrix& norm_dx_ip, Real rho_thresh, Real eps_w);
        static void     new_z_state(state_vec& state, Integer f, Integer l, const Matrix& norm_dz_ip, 
                                    Real dz_thresh, Real rho_thresh, Real eps_w);
        static void     update_state(state_vec& state, Integer max_it, Integer f, Integer l, Integer& done);

        void            eval_1(Matrix& Y_sol, bool rev);
        void            eval_sparse(Matrix& Y_sol, bool rev);
        void            eval_dense(Matrix& Y_sol, bool rev);
        void            eval_sparse(Matrix& Y, bool rev, Integer fv, Integer lv, final_sparse_sol& fss);

    public:
        iterative_refinement(const Matrix& X, trans_type tA, const linsolve_obj& lo, const options& opts);

        void            eval(Matrix& Y_sol, bool rev);
};

};};