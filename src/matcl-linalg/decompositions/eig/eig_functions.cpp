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
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/gschur.h"
#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/decompositions/eig/schur_range.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/utils/linalg_utils.h"

namespace matcl
{

mat_tup_2 matcl::schur(const Matrix& A, schur_sym_alg alg)
{
    schur_decomposition obj(A, alg);
    return  mat_tup_2(obj.U(), obj.TA());
}
mat_tup_2 matcl::schur(Matrix&& A, schur_sym_alg alg)
{
    schur_decomposition obj(std::move(A), alg);
    return  mat_tup_2(obj.U(), obj.TA());
}

mat_tup_2 matcl::schur_compl(const Matrix& A, schur_sym_alg alg)
{
    value_code vc   = A.get_value_code();
    value_code vcc  = matrix_traits::complex_value_type(vc);

    if (vc == vcc)
    {
        schur_decomposition obj(A, alg);
        return  mat_tup_2(obj.U(), obj.TA());
    }
    else
    {
        mat_code mc = matrix_traits::get_matrix_type(vcc, A.get_struct_code());
        Matrix Ac = convert(A, mc);
        schur_decomposition obj(std::move(Ac), alg);
        return  mat_tup_2(obj.U(), obj.TA());
    };
}
mat_tup_2 matcl::schur_compl(Matrix&& A, schur_sym_alg alg)
{
    value_code vc   = A.get_value_code();
    value_code vcc  = matrix_traits::complex_value_type(vc);

    if (vc == vcc)
    {
        schur_decomposition obj(std::move(A), alg);
        return  mat_tup_2(obj.U(), obj.TA());
    }
    else
    {
        mat_code mc = matrix_traits::get_matrix_type(vcc, A.get_struct_code());
        Matrix Ac = convert(A, mc);
        schur_decomposition obj(std::move(Ac), alg);
        return  mat_tup_2(obj.U(), obj.TA());
    };
}

mat_tup_4 matcl::gschur(const Matrix& A, const Matrix& B)
{
    gschur_decomposition obj(A, B);
    return  mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
}
mat_tup_4 matcl::gschur(Matrix&& A, const Matrix& B)
{
    gschur_decomposition obj(std::move(A), B);
    return  mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
}
mat_tup_4 matcl::gschur(const Matrix& A, Matrix&& B)
{
    gschur_decomposition obj(A, std::move(B));
    return  mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
}
mat_tup_4 matcl::gschur(Matrix&& A, Matrix&& B)
{
    gschur_decomposition obj(std::move(A), std::move(B));
    return  mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
}

template<class Mat1, class Mat2>
static mat_tup_4 gschur_compl_impl(Mat1&& A, Mat2&& B)
{
    value_code vc_A = A.get_value_code();
    value_code vc_B = B.get_value_code();

    bool compl  = matrix_traits::is_float_complex(vc_A) || matrix_traits::is_float_complex(vc_B);

    if (compl == true)
    {
        gschur_decomposition obj(std::forward<Mat1>(A), std::forward<Mat2>(B));
        return  mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
    }
    else
    {
        value_code vcc_A    = matrix_traits::complex_value_type(vc_A);
        value_code vcc_B    = matrix_traits::complex_value_type(vc_B);

        mat_code mc_A       = matrix_traits::get_matrix_type(vcc_A, A.get_struct_code());
        mat_code mc_B       = matrix_traits::get_matrix_type(vcc_B, B.get_struct_code());

        Matrix Ac = convert(A,mc_A);
        Matrix Bc = convert(B,mc_B);

        gschur_decomposition obj(std::move(Ac), std::move(Bc));
        return  mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
    };
}

mat_tup_4 matcl::gschur_compl(const Matrix& A, const Matrix& B)
{
    return gschur_compl_impl(A,B);
}
mat_tup_4 matcl::gschur_compl(Matrix&& A, const Matrix& B)
{
    return gschur_compl_impl(std::move(A),B);
}
mat_tup_4 matcl::gschur_compl(const Matrix& A, Matrix&& B)
{
    return gschur_compl_impl(A,std::move(B));
}
mat_tup_4 matcl::gschur_compl(Matrix&& A, Matrix&& B)
{
    return gschur_compl_impl(std::move(A),std::move(B));
}

mat_tup_2 matcl::ordschur(const Matrix& U, const Matrix& T, const Matrix& select)
{
    schur_decomposition obj(U, T);
    obj.select(select);
    return  mat_tup_2(obj.U(), obj.TA());
}

mat_tup_4 matcl::ordgschur(const Matrix& Q, const Matrix& Z,
                         const Matrix& S, const Matrix& T, const Matrix& select)
{
    gschur_decomposition obj(Q,Z,S,T);
    obj.select(select);
    return mat_tup_4(obj.Q(), obj.Z(), obj.TA(), obj.TB() );
}

Matrix matcl::eig(const Matrix& A, schur_sym_alg alg)
{
    schur_decomposition obj(A, alg, false);
    return obj.eig();
};
Matrix matcl::eig(Matrix&& A, schur_sym_alg alg)
{
    schur_decomposition obj(std::move(A), alg, false);
    return obj.eig();
};

Matrix matcl::eig(const Matrix& A, const Matrix& B)
{
    gschur_decomposition obj(A, B, false);
    return obj.eig();
};
Matrix matcl::eig(Matrix&& A, const Matrix& B)
{
    gschur_decomposition obj(std::move(A), B, false);
    return obj.eig();
};
Matrix matcl::eig(const Matrix& A, Matrix&& B)
{
    gschur_decomposition obj(A, std::move(B), false);
    return obj.eig();

};
Matrix matcl::eig(Matrix&& A, Matrix&& B)
{
    gschur_decomposition obj(std::move(A), std::move(B), false);
    return obj.eig();
};

Matrix matcl::eig_tridiag(const Matrix &A, const Matrix &B)
{
    check_tridiag(A, B);

    Integer N           = A.length();
    value_code vc_A     = A.get_value_code();
    value_code vc_B     = B.get_value_code();
    value_code vc       = matrix_traits::unify_value_types(vc_A, vc_B);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (N == 0)
    {
        vc              = matrix_traits::real_value_type(vc);

        Matrix S        = zeros(0,0,vc);
        return S;
    }

    if (N == 1)
        return A;

    Matrix AB           = make_band_noinit(N,N,-1,1,vc);
    AB.diag(0)          = A;
    AB.diag(-1)         = B(colon(1,N-1));
    AB.diag(1)          = conj(B(colon(1,N-1)));    

    AB.add_struct(predefined_struct_type::her);

    return eig(std::move(AB));
};

mat_tup_2 matcl::eig_tridiag2(const Matrix &A, const Matrix &B)
{
    check_tridiag(A, B);

    Integer N           = A.length();
    value_code vc_A     = A.get_value_code();
    value_code vc_B     = B.get_value_code();
    value_code vc       = matrix_traits::unify_value_types(vc_A, vc_B);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (N == 0)
    {
        vc              = matrix_traits::real_value_type(vc);

        Matrix S        = zeros(0,0,vc);
        return mat_tup_2(S,S);
    }

    if (N == 1)
    {
        Matrix I        = eye(1,1,vc);
        return mat_tup_2(A,I);
    }

    Matrix AB           = make_band_noinit(N,N,-1,1,vc);
    AB.diag(0)          = A;
    AB.diag(-1)         = B(colon(1,N-1));
    AB.diag(1)          = conj(B(colon(1,N-1)));    

    AB.add_struct(predefined_struct_type::her);

    schur_decomposition sd(AB);
    Matrix E            = sd.eig();
    Matrix V            = sd.U();
    return mat_tup_2(E,V);
};

Matrix matcl::eigsel_tridiag_range(const Matrix& diag, const Matrix& subdiag, Real VL, Real VU)
{
    Matrix E;
    schur_range::eval_tridiag_range(E, diag, subdiag, VL, VU);
    return E;
};
Matrix matcl::eigsel_tridiag_index(const Matrix& diag, const Matrix& subdiag, Integer IF, Integer IL)
{
    Matrix E;
    schur_range::eval_tridiag_index(E, diag, subdiag, IF, IL);
    return E;
};

mat_tup_2 matcl::eigsel_tridiag_range2(const Matrix& diag, const Matrix& subdiag, Real VL, Real VU)
{
    mat_tup_2 ret;
    schur_range::eval_tridiag2_range(ret, diag, subdiag, VL, VU);
    return ret;
}
mat_tup_2 matcl::eigsel_tridiag_index2(const Matrix& diag, const Matrix& subdiag, Integer IF, Integer IL)
{
    mat_tup_2 ret;
    schur_range::eval_tridiag2_index(ret, diag, subdiag, IF, IL);
    return ret;
}

Matrix matcl::eigsel_range(const Matrix& A, Real VL, Real VU)
{
    Matrix E;
    schur_range::eval_range(E, A, VL, VU);
    return E;
}
Matrix matcl::eigsel_index(const Matrix& A, Integer IF, Integer IL)
{
    Matrix E;
    schur_range::eval_index(E, A, IF, IL);
    return E;
}
mat_tup_2 matcl::eigsel_range2(const Matrix& A, Real VL, Real VU)
{
    mat_tup_2 ret;
    schur_range::eval2_range(ret, A, VL, VU);
    return ret;
}
mat_tup_2 matcl::eigsel_index2(const Matrix& A, Integer IF, Integer IL)
{
    mat_tup_2 ret;
    schur_range::eval2_index(ret, A, IF, IL);
    return ret;
}

mat_tup_2 matcl::gen_sym_eigen(const Matrix& A, const Matrix& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(A,B,type,true);
    return mat_tup_2(gs.V(), gs.D());
};
mat_tup_2 matcl::gen_sym_eigen(Matrix&& A, const Matrix& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(std::move(A),B,type,true);
    return mat_tup_2(gs.V(), gs.D());
};
mat_tup_2 matcl::gen_sym_eigen(const Matrix& A, Matrix&& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(A,std::move(B),type,true);
    return mat_tup_2(gs.V(), gs.D());
};
mat_tup_2 matcl::gen_sym_eigen(Matrix&& A, Matrix&& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(std::move(A),std::move(B),type,true);
    return mat_tup_2(gs.V(), gs.D());
};

Matrix matcl::eig_sym(const Matrix& A, const Matrix& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(A,B,type,true);
    return gs.eig();
};
Matrix matcl::eig_sym(Matrix&& A, const Matrix& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(std::move(A),B,type,true);
    return gs.eig();
};
Matrix matcl::eig_sym(const Matrix& A, Matrix&& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(A,std::move(B),type,true);
    return gs.eig();
};

Matrix matcl::eig_sym(Matrix&& A, Matrix&& B, gschur_sym_type type)
{
    gschur_sym_decomposition gs(std::move(A),std::move(B),type,true);
    return gs.eig();
};

pschur_return matcl::pschur(const linear_operator& A, const Integer k,  cluster_type ec, const options& opts)
{
    pschur_decomposition sd(A, k, ec, opts);
    Matrix U    = sd.U();
    Matrix T    = sd.TA();
    bool conv   = sd.converged();

    return pschur_return(U,T,conv);
};
pschur_return matcl::pbschur(const linear_operator& A, const linear_operator& B, bool hermitian, const Integer k,  
                             cluster_type ec, const options& opts)
{
    pbschur_decomposition sd(A, B, hermitian, k, ec, opts);
    Matrix U    = sd.U();
    Matrix T    = sd.TA();
    bool conv   = sd.converged();

    return pschur_return(U,T,conv);
};

eigs_return matcl::eigs(const linear_operator& A, const Integer k,  cluster_type ec, const options& opts)
{
    pschur_decomposition sd(A, k, ec, opts);
    Matrix E    = sd.eig();
    bool conv   = sd.converged();

    return eigs_return(E,conv);
};

eigs_return matcl::eigs(const linear_operator& A, const linear_operator& B, bool herm, const Integer k, cluster_type ec, 
                        const options& opts)
{
    pbschur_decomposition sd(A, B, herm, k, ec, opts);
    Matrix E    = sd.eig();
    bool conv   = sd.converged();

    return eigs_return(E,conv);
};

eigs_return matcl::svds1(const linear_operator& A, const Integer k, cluster_type ec, const options& opts)
{
    Integer M       = A.rows(); 
    Integer N       = A.cols();
    bool trans      = (M < N)? false : true;
    bool use_her    = A.is_hermitian();

    switch(ec)
    {
        case cluster_type::LM:
        case cluster_type::SM:
            break;
        default:
            use_her = false;
    };

    if (use_her)
    {
        pschur_decomposition sd(A, k, ec, opts);

        Matrix E    = sd.eig();
        bool conv   = sd.converged();
        E           = abs(std::move(E));
        E           = sort(std::move(E), 1, false);
        return eigs_return(E,conv);
    }
    else
    {
        linear_operator AAt = linop_symprod(A, trans);
        pschur_decomposition sd(AAt, k, ec, opts);

        Matrix E    = sd.eig();
        bool conv   = sd.converged();
        E           = sqrt(abs(std::move(E)));
        E           = sort(std::move(E), 1, false);
        return eigs_return(E,conv);
    };
};

svds_return matcl::svds(const linear_operator& A, const Integer k, cluster_type ec, const options& opts)
{
    Integer M       = A.rows(); 
    Integer N       = A.cols();
    bool trans      = (M < N)? false : true;
    bool use_her    = A.is_hermitian();

    switch(ec)
    {
        case cluster_type::LM:
        case cluster_type::SM:
            break;
        default:
            use_her = false;
    };

    if (use_her)
    {
        pschur_decomposition sd(A, k, ec, opts);

        Matrix E    = sd.eig();
        Matrix U    = sd.U();
        Matrix V    = U;
        bool conv   = sd.converged();

        Matrix I    = find(E < 0.0f);
        Matrix C    = ones(U.cols(), 1, U.get_value_code());
        C(I)        = -1.0f;
        V           = scale_cols(std::move(V), C);

        E           = abs(std::move(E));
        tie(E,I)    = sort2(std::move(E), 1, false);
        U           = U(colon(), I);        
        V           = V(colon(), I);  

        E           = bdiag(E);
        return svds_return(U,E,V,conv);
    }
    else
    {
        linear_operator AAt = linop_symprod(A, trans);
        pschur_decomposition sd(AAt, k, ec, opts);

        Matrix E    = sd.eig();
        Matrix U2   = sd.U();
        bool conv   = sd.converged();
        Matrix U, V;
        
        if (trans == false)
        {
            //A = U * S * V'
            //A * A' = U * S2 * U', 
            //V*S = A'*U

            U           = U2;
            E           = sqrt(abs(std::move(E)));
            Matrix I;
            tie(E,I)    = sort2(std::move(E), 1, false);
            U           = U(colon(), I);
            V           = A.mmul_right(U, trans_type::conj_trans);

            unitary_matrix Q;
            Matrix R;
            tie(Q,R)    = qr2(std::move(V), true);
            V           = Q.to_matrix();
            E           = bdiag(E);

            // correct signs
            R           = R.diag(0);
            I           = find(R < 0.0f);
            Matrix C    = ones(U.cols(), 1, U.get_value_code());
            C(I)        = -1.0f;
            V           = scale_cols(std::move(V), C);
        }
        else
        {
            //A = U * S * V'
            //A' * A = V * S2 * V'
            //U = A * V * S^-1

            V           = U2;
            E           = sqrt(abs(std::move(E)));
            Matrix I;
            tie(E,I)    = sort2(std::move(E), 1, false);
            V           = V(colon(), I);
            U           = A.mmul_right(V, trans_type::no_trans);
            //U         = scale_cols(std::move(U), div(1.0f, E));

            unitary_matrix Q;
            Matrix R;
            tie(Q,R)    = qr2(std::move(U), true);
            U           = Q.to_matrix();
            E           = bdiag(E);

            // correct signs
            R           = R.diag(0);
            I           = find(R < 0.0f);
            Matrix C    = ones(U.cols(), 1, U.get_value_code());
            C(I)        = -1.0f;
            U           = scale_cols(std::move(U), C);
        };
        
        return svds_return(U,E,V,conv);
    };
};

};