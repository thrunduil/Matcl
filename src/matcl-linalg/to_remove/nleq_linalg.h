#pragma once

#include "linalg/config_linalg.h"
#include "mmlib_basic/mmlib_header.h"

namespace forma_lib
{

/*
*  PURPOSE: find least square solution to linear equation AX = Y 
*           multiple times with different rhs, A must be square
*/
class MMLIB_LINALG_EXPORT linsolve_ls_obj
{
    public:
        typedef std::shared_ptr<linsolve_ls_obj> linsolve_ls_obj_ptr;
        typedef mmlib::tuple<mmlib::Matrix,bool>  matrix_bool;

    private:
        mmlib::Real         m_rank_atol;
        mmlib::Real         m_rank_rtol;
        mmlib::Real         m_tol;
        bool                m_do_scal_rows;
        bool                m_do_scal_cols;

        mmlib::Integer      m_rows;
        mmlib::Integer      m_cols;
        mmlib::Integer      m_rk;
        mmlib::Matrix       m_p;
        mmlib::Matrix       m_q;
        mmlib::Matrix       m_L;
        mmlib::Matrix       m_U;
        mmlib::Matrix       m_A;
        mmlib::Matrix       m_U_ind;
        mmlib::Matrix       m_U_val;
        mmlib::Matrix       m_scale_cols;
        mmlib::Matrix       m_scale_rows;
        mmlib::Real         m_scale_cols_norm;

        mmlib::Integer      m_r_rk;
        mmlib::Matrix       m_r_U_11t;
        mmlib::Matrix       m_r_U_12t;
        mmlib::Matrix       m_r_L_11t;
        mmlib::Matrix       m_r_null;
        mmlib::Matrix       m_r_null_l;
        mmlib::Matrix       m_r_pinv;
        mmlib::Matrix       m_r_q;
        mmlib::Real         m_r_norm_U;

    public:
        linsolve_ls_obj();

        matrix_bool          solve(const mmlib::Matrix& X) const;
        mmlib::Matrix        get_null_left() const;
        mmlib::Matrix        get_null_right() const;
        mmlib::Integer       get_rank() const; 
        mmlib::Matrix        project_on_left_null(const mmlib::Matrix& X) const;
        mmlib::Matrix        project_on_left_orth(const mmlib::Matrix& X) const;
        mmlib::Matrix        project_on_right_null(const mmlib::Matrix& X) const;
        mmlib::Matrix        project_on_right_orth(const mmlib::Matrix& X) const;
        void                 set_rank(mmlib::Integer rank);
        mmlib::Matrix        get_scal_rows() const;
        mmlib::Matrix        get_scal_cols() const;

        friend linsolve_ls_obj_ptr MMLIB_LINALG_EXPORT 
                             prepare_linsolve_ls(const mmlib::Matrix& A, const mmlib::Matrix& F);

    private:
        mmlib::Integer       get_rank(mmlib::Matrix& Ud, const mmlib::Matrix& Um);
        mmlib::Matrix        normalize_null(const mmlib::Matrix& N, int dim) const;
        mmlib::Matrix        scale_matrix(const mmlib::Matrix& A, const mmlib::Matrix& F);
        mmlib::Matrix        inv_scal(const mmlib::Matrix& scal) const;
};
typedef std::shared_ptr<linsolve_ls_obj> linsolve_ls_obj_ptr;

/*
*  PURPOSE: prepare linsolve_ls_obj allowing for finding least square solution
*           to linear equation AX = Y many times with different matrix Y
*
*  USAGE: out = prepare_linsolve_ls(A, Y)
*  INPUT: A                 matrix
*         Y                 rhs matrix, used only for internal scaling purpose
*                           use [] if not gived
*
*  OUTPUT: out              linolve_obj_ls object
*/
linsolve_ls_obj_ptr MMLIB_LINALG_EXPORT prepare_linsolve_ls(const mmlib::Matrix& A, const mmlib::Matrix& F = mmlib::zeros(0,0));

};