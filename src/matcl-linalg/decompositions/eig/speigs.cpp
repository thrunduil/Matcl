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

#include "matcl-linalg/decompositions/speigs.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "arpack_wrapper.h"
#include "matcl-linalg/options/options_speigs.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/func/converter.h"

namespace matcl { namespace details
{

namespace md = matcl :: details;

class pschur_impl
{
    public:
        using impl_type = std::shared_ptr<details::pschur_impl>;

    public:
        static void set_impl(pschur_decomposition& ret, const impl_type& impl)
        {
            ret.m_impl = impl;
        };
        static void set_impl(pbschur_decomposition& ret, const impl_type& impl)
        {
            ret.m_impl = impl;
        };

        virtual void    set_x0(const Matrix& x0) = 0;
        virtual void    set_cluster(cluster_type ec) = 0;
        virtual void    compute(Integer k, bool return_nonconverged) = 0;
        virtual Matrix  get_U() const = 0;
        virtual Matrix  get_TA() const = 0;
        virtual Matrix  get_eig() const = 0;
        virtual bool    get_converged() const = 0;
        virtual Integer size() const = 0;
        virtual Integer number_converged_eigenvalues() const = 0;
        virtual Integer number_Arnoldi_iterations() const = 0;
        virtual Integer number_operations_Ax() const = 0;
        virtual Integer number_operations_Bx() const = 0;
        virtual Integer number_reorthogonalizations() const = 0;
};

template<class Val>
class pschur_impl_val : public pschur_impl
{
    private:
        using raw_impl_type = md::arpack_wrapper<Val>;

    private:
        raw_impl_type   m_impl;

    public:
        pschur_impl_val(const linear_operator& A, const options& opts)
            :m_impl(A, arpack_output::schur, get_max_iter(opts), get_num_arnoldi(opts), get_tol(opts))
        {};
        pschur_impl_val(const linear_operator& A, const linear_operator& B, bool hermitian, const options& opts)
            :m_impl(A, B, hermitian, arpack_output::schur, get_max_iter(opts), get_num_arnoldi(opts), get_tol(opts))
        {};

        static void construct(const linear_operator& A, const options& opts, pschur_decomposition& sd)
        {            
            impl_type impl(new pschur_impl_val(A, opts));
            set_impl(sd, impl);
        };
        static void construct_g(const linear_operator& A, const linear_operator& B, bool hermitian,
                                const options& opts, pbschur_decomposition& sd)
        {            
            impl_type impl(new pschur_impl_val(A, B, hermitian, opts));
            set_impl(sd, impl);
        };

        virtual void set_x0(const Matrix& x0) override
        {
            Integer M       = x0.rows();
            Integer N       = x0.cols();

            if (m_impl.size() != M || N != 1)
                throw error::invalid_size2(M, N, m_impl.size(), std::max(N,1));

            value_code vc   = m_impl.get_value_code();

            if (x0.get_value_code() != vc)
            {
                mat_code mc = matrix_traits::get_matrix_type(vc, struct_code::struct_dense);
                Matrix x    = convert(x0, mc);
                m_impl.set_x0(x);
            }
            else
            {
                m_impl.set_x0(x0);
            };
        };

        virtual void set_cluster(cluster_type ec) override
        {
            m_impl.set_cluster(ec);
        };

        virtual void compute(Integer k, bool return_nonconverged) override
        {
            m_impl.calculate(k, return_nonconverged);
        };

        virtual Matrix get_U() const override
        {
            return m_impl.get_U();
        };
        virtual Matrix get_TA() const override
        {
            return m_impl.get_T();
        };
        virtual Matrix get_eig() const override
        {
            return m_impl.get_D();
        };
        virtual bool get_converged() const override
        {
            return m_impl.get_converged();
        };
        virtual Integer size() const override
        {
            return m_impl.size();
        };
        virtual Integer number_converged_eigenvalues() const override
        {
            return m_impl.number_converged_eigenvalues();
        };
        virtual Integer number_Arnoldi_iterations() const override
        {
            return m_impl.number_Arnoldi_iterations();
        };
        virtual Integer number_operations_Ax() const override
        {
            return m_impl.number_operations_Ax();
        };
        virtual Integer number_operations_Bx() const override
        {
            return m_impl.number_operations_Bx();
        };
        virtual Integer number_reorthogonalizations() const override
        {
            return m_impl.number_reorthogonalizations();
        };
    private:
        Integer get_max_iter(const options& opts) const
        {
            Integer max_it  = opts.get_option<Integer>(opt::speigs::maxit());
            return std::max(max_it, 1);
        };

        Real get_tol(const options& opts) const
        {
            Real tol = opts.get_option<Real>(opt::speigs::tol());
            return tol;
        };

        Integer get_num_arnoldi(const options& opts) const
        {
            Integer num_arnoldi  = opts.get_option<Integer>(opt::speigs::n_arnoldi());
            return num_arnoldi;
        };
};

template<>
struct pschur_impl_val<Integer>
{
    static void construct(const linear_operator& A, const options& opts, pschur_decomposition& sd)
    {
        linear_operator Ac  = A.convert(value_code::v_real);
        return pschur_impl_val<Real>::construct(Ac,opts,sd);
    };
    static void construct_g(const linear_operator& A, const linear_operator& B, bool hermitian, 
                            const options& opts, pbschur_decomposition& sd)
    {
        linear_operator Ac  = A.convert(value_code::v_real);
        linear_operator Bc  = B.convert(value_code::v_real);
        return pschur_impl_val<Real>::construct_g(Ac,Bc,hermitian, opts,sd);
    };
};

template<>
struct pschur_impl_val<Object>
{
    static void construct(const linear_operator&, const options&, pschur_decomposition&)
    {
        throw error::object_value_type_not_allowed("speigs");
    };
    static void construct_g(const linear_operator&, const linear_operator&, bool, const options&, 
                            pbschur_decomposition&)
    {
        throw error::object_value_type_not_allowed("speigs");
    };
};

struct vis_speigs : public extract_type_switch<void,vis_speigs,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, const linear_operator& A, const options& opts, 
                     pschur_decomposition& sd)
    {
        using V = typename T::value_type;
        return pschur_impl_val<V>::construct(A, opts, sd);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const linear_operator& A, const options& opts, 
                     pschur_decomposition& sd)
    {
        return pschur_impl_val<T>::construct(A, opts, sd);
    };
};

struct vis_speigs_g : public extract_type_switch<void,vis_speigs_g,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, const linear_operator& A, const linear_operator& B, 
                     bool hermitian, const options& opts, pbschur_decomposition& sd)
    {
        using V = typename T::value_type;
        return pschur_impl_val<V>::construct_g(A, B, hermitian, opts, sd);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const linear_operator& A, const linear_operator& B, 
                            bool hermitian, const options& opts, pbschur_decomposition& sd)
    {
        return pschur_impl_val<T>::construct_g(A, B, hermitian, opts, sd);
    };
};

}};

namespace matcl
{

pschur_decomposition::pschur_decomposition()
    :m_initialized(false)
{};

pschur_decomposition::pschur_decomposition(const linear_operator& A, Integer k, cluster_type ec,
                                           const options& opts)
{
    initialize(A,opts);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);
};

pschur_decomposition::pschur_decomposition(const linear_operator& A, Integer k, const Matrix& x0, 
                                           cluster_type ec, const options& opts)
{
    initialize(A,opts);
    set_x0(x0);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);
};

pschur_decomposition::pschur_decomposition(const Matrix& A, Real sig, bool sing,
        Integer k, cluster_type ec, const options& opts)
    :pschur_decomposition(linsolve_shift_invert(A, sig, sing, opts), k, ec, opts)
{
    transform(A);
};

pschur_decomposition::pschur_decomposition(const Matrix& A, Real sig, bool sing, 
                Integer k, const Matrix& x0, cluster_type ec, const options& opts)
    :pschur_decomposition(linsolve_shift_invert(A, sig, sing, opts), k, x0, ec, opts)
{
    transform(A);
};

pschur_decomposition::~pschur_decomposition()
{};

pschur_decomposition& pschur_decomposition::operator()(const linear_operator& A, Integer k, 
                                                       cluster_type ec, const options& opts)
{
    initialize(A,opts);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);

    return *this;
};
pschur_decomposition& pschur_decomposition::operator()(const linear_operator& A, Integer k, const Matrix& x0, 
                                    cluster_type ec, const options& opts)
{
    initialize(A,opts);
    set_x0(x0);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);

    return *this;
};

pschur_decomposition& pschur_decomposition::operator()(const Matrix& A, Real sig, bool sing, 
                            Integer k, cluster_type ec, const options& opts)
{
    this->operator()(linsolve_shift_invert(A, sig, sing, opts), k, ec, opts);
    this->transform(A);
    return *this;
};
pschur_decomposition& pschur_decomposition::operator()(const Matrix& A, Real sig, bool sing, 
                            Integer k, const Matrix& x0, cluster_type ec, const options& opts)
{
    this->operator()(linsolve_shift_invert(A, sig, sing, opts), k, x0, ec, opts);
    this->transform(A);
    return *this;
};

void pschur_decomposition::initialize(const linear_operator& A, const options& opts)
{
    if (A.rows() != A.cols())
        throw error::error_size_eig(A.rows(), A.cols());

    value_code vc_A = A.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_A, value_code::v_float);

    linear_operator Ac  = (vc == vc_A)? A : A.convert(vc);

    m_initialized   = true;
    Matrix dum      = zeros(0,0,vc);
    return details::vis_speigs::make<const Matrix&>(dum, Ac, opts, *this);
};

void pschur_decomposition::compute(Integer k, bool return_nonconverged)
{
    if (k < 1 || k > m_impl->size() - 3)
        throw error::invalid_speigs_k(k, m_impl->size() - 3);
        
    m_impl->compute(k,return_nonconverged);

    Matrix U    = m_impl->get_U();
    Matrix TA   = m_impl->get_TA();

    base_type::set_factors(U, TA);
};
void pschur_decomposition::transform(const Matrix& A)
{
    Matrix U    = this->U();
    Matrix UAU  = ctrans(U) * A * U;

    bool is_re  = matrix_traits::is_float_real(A.get_value_code());
    bool herm   = A.get_struct().is_hermitian(A.is_square(), is_re);

    if (herm == true)
        A.add_struct(predefined_struct_type::her);

    Matrix TA;
    tie(U,TA)   = schur(UAU);
    base_type::set_factors(U, TA);
};

void pschur_decomposition::set_x0(const Matrix& x0)
{
    m_impl->set_x0(x0);
};
void pschur_decomposition::set_cluster(cluster_type ec)
{
    m_impl->set_cluster(ec);
};

bool pschur_decomposition::converged() const
{
    check();
    return m_impl->get_converged();
};

Integer pschur_decomposition::number_converged_eigenvalues() const
{
    check();
    return m_impl->number_converged_eigenvalues();
};

Integer pschur_decomposition::number_Arnoldi_iterations() const
{
    check();
    return m_impl->number_Arnoldi_iterations();
};

Integer pschur_decomposition::number_operations_Ax() const
{
    check();
    return m_impl->number_operations_Ax();
};
Integer pschur_decomposition::number_reorthogonalizations() const
{
    check();
    return m_impl->number_reorthogonalizations();
};

void pschur_decomposition::check() const
{
    if (!m_impl)
        throw error::uninitialized_object_used("speigs");
};

pbschur_decomposition::pbschur_decomposition()
{};

pbschur_decomposition::~pbschur_decomposition()
{};

pbschur_decomposition::pbschur_decomposition(const linear_operator& A, const linear_operator& B, 
                        bool hermitian, Integer k, cluster_type ec, const options& opts)
{
    initialize_G(A,B,hermitian,opts);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);
};

void pbschur_decomposition::initialize_G(const linear_operator& A, const linear_operator& B, 
                                         bool hermitian, const options& opts)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    if (B.rows() != B.cols())
        throw error::square_matrix_required(B.rows(), B.cols());

    if (A.rows() != B.rows())
        throw error::invalid_size2(B.rows(), B.cols(), A.rows(), A.rows());

    value_code vc_A = A.get_value_code();
    value_code vc_B = B.get_value_code();
    value_code vc0  = matrix_traits::unify_value_types(vc_A, vc_B);
    value_code vc   = matrix_traits::unify_value_types(vc0, value_code::v_float);

    linear_operator Ac  = (vc == vc_A)? A : A.convert(vc);
    linear_operator Bc  = (vc == vc_B)? B : B.convert(vc);

    m_initialized   = true;
    Matrix dum      = zeros(0,0,vc);

    return details::vis_speigs_g::make<const Matrix&>(dum, Ac, Bc, hermitian, opts, *this);
};

pbschur_decomposition::pbschur_decomposition(const linear_operator& A, const linear_operator& B, 
                        bool hermitian, Integer k, const Matrix& x0, cluster_type ec, const options& opts)
{
    initialize_G(A,B,hermitian,opts);
    set_x0(x0);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);
};
pbschur_decomposition& pbschur_decomposition::operator()(const linear_operator& A, const linear_operator& B, 
                            bool hermitian, Integer k, cluster_type ec, const options& opts)
{
    initialize_G(A,B,hermitian, opts);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);

    return *this;
};
pbschur_decomposition& pbschur_decomposition::operator()(const linear_operator& A, const linear_operator& B, 
                            bool hermitian, Integer k, const Matrix& x0, cluster_type ec, const options& opts)
{
    initialize_G(A,B,hermitian,opts);
    set_x0(x0);
    set_cluster(ec);

    bool return_nonconvergent   = opts.get_option<bool>(opt::speigs::return_nonconvergent());
    compute(k,return_nonconvergent);

    return *this;
};

Integer pbschur_decomposition::number_operations_Bx() const
{
    check();
    return m_impl->number_operations_Bx();
};

Matrix pbschur_decomposition::U() const
{
    return base_type::U();
};

linsolve_obj matcl::linsolve_shift_invert(const Matrix& A, Real sig, bool sing, const options& opts)
{
    if (A.is_square() == false)
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix B    = A;

    if (sig != 0.0)
    {
        Matrix sig_m    = sig;
        sig_m           = details::convert_value(sig_m, A.get_value_code());
        B               = A - sig_m * speye(A.rows(), A.cols(), A.get_value_code());
    };

    namespace ol    = opt::linsolve;
    options opts2   = opts;

    if (sing == true)
    {
        //rank revealing factorization must be used
        opts2.set(ol::use_rr(true));

        //iterative refinement can only produce worse solution
        opts2.set(ol::use_ir(false));

        //solution cannot be valid
        opts2.set(ol::test_sol(false));

        //perform balancing
        if (opts2.has_option(ol::do_balancing_rr()) == false)
            opts2.set(ol::do_balancing_rr(true));
        
        //small pivots must be perturbed
        if (opts2.has_option(ol::tol_sing()) == false)
        {
            opts2.set(ol::tol_sing(0.7));
        }
        else
        {
            Real tol    = opts2.get_option<Real>(ol::tol_sing());
            if (tol == 0.0)
                opts2.set(ol::tol_sing(0.7));
        };
    }
    else
    {
        opts2.set(ol::tol_sing(0.0));

        //other options are OK
    };

    return make_linsolve_obj(B, opts2);
};

};