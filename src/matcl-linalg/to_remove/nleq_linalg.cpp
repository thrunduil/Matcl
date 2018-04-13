#include "mmlib_basic/general/exception.h"
#include "nleq_linalg.h"
#include "linalg/linalg.h"

#include <iostream> // DEBUG

#include <stdexcept>

//#include "lrem/lrem_linalg.h"
//#include "lrem/lrem.h"

using namespace mmlib;

namespace forma_lib
{

linsolve_ls_obj_ptr prepare_linsolve_ls(const mmlib::Matrix& A, const mmlib::Matrix& F)
{
    linsolve_ls_obj_ptr ret(new linsolve_ls_obj());

    ret->m_rows     = A.rows();
    ret->m_cols     = A.cols();
    ret->m_A        = A;

    Matrix A_sc     = ret->scale_matrix(A,F);

    Matrix L, U;
    permutation_vector p,q;
    lu_options opts;
    opts.pivot_type = lu_options::rook;

    std::tie(L,U,p,q)    = lu(trans(A_sc),opts);

    Matrix U_d      = full(get_diag(U));
    Matrix U_m      = full(max_abs(U));
    Integer rk      = ret->get_rank(U_d, U_m);

    {
        Matrix val, sel;
        tie(val,sel)    = sort2(abs(U_d),1,false);
        ret->m_U_ind    = sel(colon());
        ret->m_U_val    = val(colon());
    };

    ret->m_p        = p.to_matrix();    
    ret->m_q        = q.to_matrix();
    ret->m_L        = L;
    ret->m_U        = U;    
    ret->m_rk       = rk;
    ret->m_r_rk     = -1;    

    ret->set_rank(rk);

    return ret;  
};

void linsolve_ls_obj::set_rank(mmlib::Integer rank)
{
    if (std::min(m_rk, rank) == m_r_rk)
        return;

    bool rank_red   = m_r_rk > 0;
    this->m_r_rk    = std::min(m_rk, rank);

    //sorting makes U_11 matrix upper triangular
    Matrix I        = sort(m_U_ind(colon(1,m_r_rk)));
    Matrix J        = sort(m_U_ind(colon(m_r_rk+1,mmlib::end)));

    Matrix Jr       = (mat_col(), J, trans(irange(m_U_ind.length() + 1, m_cols)));
    Matrix Jc       = (mat_col(), J, trans(irange(m_U_ind.length() + 1, m_rows)));
    Matrix IJr      = (mat_col(), I, Jr);
    Matrix IJc      = (mat_col(), I, Jc);

    Matrix U_11     = m_U(colon(I),colon(I));
    Matrix U_12     = m_U(colon(I),colon(Jc));
    Matrix L_11     = m_L(colon(I),colon(I));

    this->m_r_U_11t = trans(U_11);
    this->m_r_U_12t = trans(U_12);
    this->m_r_L_11t = trans(L_11);

    Matrix IJc_i    = mmlib::invperm(IJc);
    Matrix IJr_i    = mmlib::invperm(IJr);

    this->m_r_q     = m_q(colon(IJc));
    Matrix p        = m_p(colon(IJr));

    this->m_r_pinv  = mmlib::invperm(p);

    if (m_r_rk == m_cols)
    {
        //full column rank, Nu is empty
        m_r_null    = spzeros(m_cols,0);        
    }
    else
    {
        Matrix L2   = m_L(colon(Jr),colon(I));
        m_r_null    = linsolve(m_r_L_11t, trans(L2));

        m_r_null    = (mat_col() , m_r_null , -speye(m_cols-m_r_rk));
        m_r_null    = m_r_null(colon(m_r_pinv),colon());

        m_r_null    = m_scale_cols * m_r_null;
        m_r_null    = normalize_null(m_r_null, 1);
    };

    if (m_r_rk == m_rows)
    {
        //full row rank, Nu_l is empty
        m_r_null_l  = spzeros(0,m_rows);  
    }
    else
    {
        m_r_null_l  = -mmlib::trans(linsolve(U_11, U_12));
        m_r_null_l  = (mat_row() , m_r_null_l , speye(m_rows - m_r_rk));
        m_r_null_l  = m_r_null_l(colon(), colon(mmlib::invperm(m_r_q)));

        m_r_null_l  = m_r_null_l * m_scale_rows;
        m_r_null_l  = normalize_null(m_r_null_l, 2);
    };

    m_r_norm_U      = mmlib::norm(m_r_U_12t,1);
};

linsolve_ls_obj::matrix_bool linsolve_ls_obj::solve(const mmlib::Matrix& X) const
{ 
    Integer K = X.cols();

    if (X.rows() != this->m_rows)
    {
        std::ostringstream os;
        os << "unable to solve linear equation due to nonconformant matrix sizes, lhs matrix size: "
            << this->m_rows << " x " << this->m_cols << ", rhs matrix size: "
            << X.rows() << " x " << X.cols() << ".";

        throw std::runtime_error(os.str());
    };

    if (this->m_rows == 0)
    {
        Matrix Y        = spzeros(this->m_cols,K);
        return matrix_bool(Y,true);
    };
    if (this->m_cols == 0)
    {
        Matrix Y        = spzeros(0,K);
        bool is_sol     = (mmlib::norm(X,1) == 0.);
        return matrix_bool(Y,is_sol);
    }

    Matrix Xs           = m_scale_rows * X;
    Matrix Y1           = Xs(colon(this->m_r_q),colon());
    Matrix Y2           = linsolve(this->m_r_U_11t, Y1(colon(1,this->m_r_rk),colon()));
    Y1                  = Y1(colon(this->m_r_rk+1, mmlib::end),colon());
    Matrix Y            = linsolve(this->m_r_L_11t,Y2);
    Y                   = (mat_col() , Y , spzeros(this->m_cols - this->m_r_rk,X.cols()));
    Y                   = Y(colon(this->m_r_pinv),colon());

    bool is_sol         = true;

    if (this->m_r_rk != this->m_rows)
    {
        //U12'*Y2       = Y1
        Real cond       = this->m_r_norm_U * norm(Y2,1) + norm(Y1,1);
        Real dif        = norm(this->m_r_U_12t*Y2-Y1,1);

        if (dif > cond * this->m_tol)
        {
            is_sol      = false;
        };
    };

    Y                   = m_scale_cols * Y;
    Matrix YP           = project_on_right_orth(Y);

    return matrix_bool(YP,is_sol);
};

mmlib::Matrix linsolve_ls_obj::project_on_left_orth(const mmlib::Matrix& X) const
{
    return X - project_on_left_null(X);
}
mmlib::Matrix linsolve_ls_obj::project_on_right_orth(const mmlib::Matrix& X) const
{
    return X - project_on_right_null(X);
}
mmlib::Matrix linsolve_ls_obj::project_on_right_null(const mmlib::Matrix& X) const
{
    if (X.rows() != this->m_cols)
    {
        std::ostringstream os;
        os << "invalid matrix size, expacting matrix with " << this->m_cols 
            << " columns, matrix size is: " << X.rows() << " x " << X.cols();

        throw std::runtime_error(os.str());
    };

    if (m_r_null.cols() == 0 || X.cols() == 0)
    {
        return spzeros(X.rows(), X.cols());
    };

    // m_r_null is orthogonal
    Matrix NT = trans(this->m_r_null);

    Matrix Y = this->m_r_null * (NT * X);
    return Y;
};
mmlib::Matrix linsolve_ls_obj::project_on_left_null(const mmlib::Matrix& X) const
{
    if (X.rows() != this->m_rows)
    {
        std::ostringstream os;
        os << "invalid matrix size, expacting matrix with " << this->m_rows 
            << " rows, matrix size is: " << X.rows() << " x " << X.cols();

        throw std::runtime_error(os.str());
    };

    if (m_r_null_l.rows() == 0 || X.cols() == 0)
    {
        return spzeros(X.rows(), X.cols());
    };

    // m_r_null is orthogonal
    Matrix NT = trans(this->m_r_null_l);

    Matrix Y = NT * (this->m_r_null_l * X);
    return Y;
};

Integer linsolve_ls_obj::get_rank(Matrix& Ud, const Matrix& Um)
{
    Integer N = Ud.length();

    if (N == 0)
        return 0;

    Real* ptr_d = Ud.get_array_unique<Real>();
    const Real* ptr_m = Um.get_array<Real>();

    Integer rk = 0;
    Real max_d = abs(ptr_d[0]);

    for (Integer i = 0; i < N; ++i)
    {
        max_d = std::max(max_d, (Real)abs(ptr_d[i]));
    };

    for (Integer i = 0; i < N; ++i)
    {
        if (abs(ptr_d[i]) > max_d * m_rank_atol && abs(ptr_d[i]) > ptr_m[i] * m_rank_rtol)
        {
            rk++;
        }
        else
        {
            ptr_d[i] = 0.;
        };
    };
    return rk;
}
linsolve_ls_obj::linsolve_ls_obj()
{
    m_rows  = 0;
    m_cols  = 0;
    m_rk    = 0;
    m_r_rk  = 0;

    m_r_norm_U  = 0;
    m_tol       = 0.;

    //TODO: below is the (commented out) setting of tolerance, based upon lrem.
    // For now we use random hardcoded values.

    m_tol       = 1.e-5;
    m_rank_atol = 1.e-5;
    m_rank_rtol = 1.e-5;

    //lrem::lrem_options lopts;
    //m_tol           = lopts.m_rtol_resid;
    //m_rank_atol     = 1e-1 * lopts.m_atol_singular;
    //m_rank_rtol     = 1e-1 * lopts.m_rtol_singular;

    m_do_scal_rows  = true;
    m_do_scal_cols  = true;
};

mmlib::Matrix linsolve_ls_obj::get_null_left() const
{
    return m_r_null_l;
};
mmlib::Matrix linsolve_ls_obj::get_null_right() const
{
    return m_r_null;
};
mmlib::Integer linsolve_ls_obj::get_rank() const
{
    return this->m_r_rk;
};

mmlib::Matrix linsolve_ls_obj::normalize_null(const mmlib::Matrix& N, int dim) const
{
    mmlib::Matrix Q, R;
    if (dim == 1)
    {
        tie(Q,R)    = mmlib::qr2(N,true);
    }
    else
    {
        tie(Q,R)    = mmlib::qr2(mmlib::trans(N),true);
        Q           = mmlib::trans(Q);
    };

    return Q;
}
mmlib::Matrix linsolve_ls_obj::scale_matrix(const mmlib::Matrix& A, const mmlib::Matrix& F)
{
    m_scale_cols        = mmlib::ones(A.cols(), 1);
    m_scale_rows        = mmlib::ones(A.rows(), 1);

    if (m_do_scal_cols == false && m_do_scal_rows == false)
    {
        return A;
    };

    mmlib::Matrix sF    = 0.;

    mmlib::Matrix As    = mmlib::abs(A);
    mmlib::Matrix Fs    = mmlib::abs(F);

    if (m_do_scal_rows )
    {
        if (F.is_empty() == false)
        {
            sF          = Fs * mmlib::ones(Fs.cols(),1);
            std::cerr << sF << std::endl; // DEBUG
        };
        m_scale_rows    = As * m_scale_cols + sF;
        m_scale_rows    = inv_scal(m_scale_rows);
    }
    if (m_do_scal_cols)
    {
        m_scale_cols    = mmlib::trans(mmlib::trans(m_scale_rows) * As);
        m_scale_cols    = inv_scal(m_scale_cols);
    };

    m_scale_cols_norm   = mmlib::norm(m_scale_cols,2) / mmlib::sqrt_nc(m_scale_cols.length());
    m_scale_rows        = mmlib::spdiag(m_scale_rows);
    m_scale_cols        = mmlib::spdiag(m_scale_cols);    

    mmlib::Matrix As2   = m_scale_rows * A * m_scale_cols;   
    return As2;
};

struct test_zero : mmlib::test_function_templ<test_zero>
{
    template<class T>
    bool eval_templ(const T& val) const
    {
        return val == 0.;
    };
};

mmlib::Matrix linsolve_ls_obj::inv_scal(const mmlib::Matrix& scal) const
{
    mmlib::Matrix sc    = scal;
    mmlib::Matrix I     = mmlib::find(scal,test_zero());
    sc(I)               = 1.;

    return mmlib::div(1., sc);
};
mmlib::Matrix linsolve_ls_obj::get_scal_rows() const
{
    return m_scale_rows * m_scale_cols_norm;
};
mmlib::Matrix linsolve_ls_obj::get_scal_cols() const
{
    return div(m_scale_cols, m_scale_cols_norm);
};

};
