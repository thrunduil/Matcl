/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017-2018
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

#include "matcl-scalar/lib_functions/sequence/limest.h"
#include "matcl-scalar/lib_functions/sequence/sequence.h"
#include "matcl-scalar/lib_functions/func_binary.h"

#pragma warning(push)
#pragma warning(disable: 4127) //conditional expression is constant

namespace matcl { namespace details
{

template<class T>
struct has_huge_exponent
{
    static const bool value     = false;
};

template<>
struct has_huge_exponent<mp_float>
{
    static const bool value     = true;
};

template<class Float>
class scaling
{
    private:
        Float           m_scaling;
        Float           m_scaling2;
        Float           m_small_x;

        // average convergence rate
        Float           m_mean_conv_rate;

        // number of iterations with slow (linear) convergence
        int             m_num_slow_iters;

        // number of iteration without improvement of the limit estimation
        int             m_num_failed_iters;

    public:
        void            set_epsilon(const Float& epsilon);

        void            initialize(const Float& initial_scaling);
        void            update(bool better, const Float& x);

        Float           apply(const Float& x) const;

        void            report_progress(const Float& loc_err, const Float& best_err,
                            int iter);
        void            report_no_progress();

    private:
        void            report_fast_progress();
        void            report_slow_progress();
};

template<class Float>
class limit_estimator_impl
{
    public:
        using function_type     = std::function<Float (const Float&)>;
        using extrapolator_ptr  = matcl::extrapolator_ptr<Float>;
        using scaling_type      = scaling<Float>;

    private:
        Float               m_epsilon;
        Float               m_huge_value;
        Float               m_tiny_value;

        int                 m_precision;
        extrapolator_ptr    m_extrapolator;
        int                 m_max_iter;
        Float               m_initial_scaling;
        Float               m_point_scaling;
        
        int                 m_iterations;
        scaling_type        m_scaling;

    public:
        limit_estimator_impl(const extrapolator_ptr& extr, int prec);

        Float           eval(const function_type& f, const Float& x, 
                            limit_type lim_type, Float& abs_err);

        void            set_extrapolator(const extrapolator_ptr& extr);
        extrapolator_ptr get_extrapolator() const;

        void            set_max_iterations(int max_iter);
        int             get_max_iterations() const;

        void            set_initial_scaling(const Float& scaling);
        Float           get_initial_scaling() const;

        void            set_initial_point_scaling(const Float& scaling);
        Float           get_initial_point_scaling() const;

        int             number_iterations() const;
        void            clear(int precision);

    private:
        int             make(const function_type& f0, Float x0, 
                            limit_type lim_type, Float& lim, Float& abs_err);

        bool            get_step_params(const function_type& f, const Float& x0, 
                            limit_type lim_type, Float& dx, Float& scale);

        function_type   get_inverse_function(const function_type& f) const;
        void            eval_function(const function_type& f, const Float& x, Float& val);

        int             default_max_iterations() const;
        Float           default_scaling() const;
        Float           default_point_scaling() const;

        void            error_invalid_starting_point(const Float& x) const;
        void            error_initial_point_not_found(const Float& x) const;
};

template<class Float>
limit_estimator_impl<Float>::limit_estimator_impl(const extrapolator_ptr& extr, int prec)
    :m_extrapolator(extr)
{
    clear(prec);

    m_max_iter          = default_max_iterations();
    m_initial_scaling   = default_scaling();
    m_point_scaling     = default_point_scaling();
}

template<class Float>
void limit_estimator_impl<Float>::clear(int prec)
{
    m_precision     = prec;
    m_huge_value    = seq_helpers::max_value<Float>::value(prec) * Float(0.1);
    m_tiny_value    = seq_helpers::min_value<Float>::value(prec) * Float(10.0);
    m_epsilon       = seq_helpers::epsilon<Float>::value(prec);
    m_iterations    = 0;

    m_scaling.set_epsilon(m_epsilon);
}

template<class Float>
int limit_estimator_impl<Float>::default_max_iterations() const
{
    return 100;
}

template<class Float>
Float limit_estimator_impl<Float>::default_scaling() const
{
    return Float(0.5);
}

template<class Float>
Float limit_estimator_impl<Float>::default_point_scaling() const
{
    return Float(1.0);
}

template<class Float>
void limit_estimator_impl<Float>::set_extrapolator(const extrapolator_ptr& extr)
{
    m_extrapolator = extr;
}

template<class Float>
void limit_estimator_impl<Float>::error_invalid_starting_point(const Float& x) const
{
    (void)x;
    throw error::limest_error_nonfinite_starting_point();
}

template<class Float>
void limit_estimator_impl<Float>::error_initial_point_not_found(const Float& x) const
{
    throw error::limest_error_valid_initial_point_not_found
                    (seq_helpers::cast_double<Float>::eval(x));
}

void details::seq_error_omega_is_zero()
{
    throw error::seq_error_omega_is_zero();
}

void details::seq_error_two_equal_points(double x)
{
    throw error::seq_error_two_equal_points(x);
}

template<class Float>
typename limit_estimator_impl<Float>::extrapolator_ptr
limit_estimator_impl<Float>::get_extrapolator() const
{
    return m_extrapolator;
}

template<class Float>
void limit_estimator_impl<Float>::set_max_iterations(int max_iter)
{
    if (max_iter <= 0)
        m_max_iter  = default_max_iterations();
    else
        m_max_iter  = max_iter;
}

template<class Float>
int limit_estimator_impl<Float>::get_max_iterations() const
{
    return m_max_iter;
}

template<class Float>
void limit_estimator_impl<Float>::set_initial_scaling(const Float& scaling)
{
    if (scaling <= Float(0.0))
        m_initial_scaling   = default_scaling();
    else
        m_initial_scaling   = scaling;
}

template<class Float>
Float limit_estimator_impl<Float>::get_initial_scaling() const
{
    return m_initial_scaling;
}

template<class Float>
void limit_estimator_impl<Float>::set_initial_point_scaling(const Float& scaling)
{
    if (scaling <= Float(0.0))
        m_point_scaling = default_point_scaling();
    else
        m_point_scaling = scaling;
}

template<class Float>
Float limit_estimator_impl<Float>::get_initial_point_scaling() const
{
    return m_point_scaling;
}

template<class Float>
int limit_estimator_impl<Float>::number_iterations() const
{
    return m_iterations;
}

template<class Float>
bool limit_estimator_impl<Float>::get_step_params(const function_type& f, 
            const Float& x0, limit_type lim_type, Float& dx, Float& scale)
{
    Float abs_x         = matcl::abs(x0);
    Float scale_x       = Float(0.5);

    Float initial_step  = m_point_scaling;

    scale               = (lim_type == limit_type::both) ? -m_initial_scaling 
                                : m_initial_scaling;
    int max_trials      = 50;

    dx                  = (abs_x == Float(0.0)) ? Float(1.0) : abs_x;
    dx                  = dx * initial_step;

    if (lim_type == limit_type::left)
        dx              = -dx;

    for (int i = 0; i < max_trials; ++i)
    {
        Float val       = f(x0 + dx);

        if (matcl::is_finite(val) == false)
        {
            dx          = dx * scale_x;
            scale_x     = scale_x * scale_x;
            continue;
        }

        return true;
    };

    return false;
}

template<class Float>
void limit_estimator_impl<Float>::eval_function(const function_type& f, 
            const Float& x, Float& val)
{
    val         = f(x);
}

template<class Float>
typename limit_estimator_impl<Float>::function_type 
limit_estimator_impl<Float>::get_inverse_function(const function_type& f) const
{
    return [f](const Float& x) -> Float { return f(Float(1.0) / x); };
}

template<class Float>
Float limit_estimator_impl<Float>::eval(const function_type& f, const Float& x0, 
                limit_type lim_type, Float& abs_err)
{
    Float x     = seq_helpers::prepare_value<Float>::eval(x0, m_precision);

    m_extrapolator->clear(m_precision);

    Float lim;
    m_iterations = make(f, x, lim_type, lim, abs_err);

    return lim;
}

template<class Float>
int limit_estimator_impl<Float>::make(const function_type& f0, Float x0, 
                limit_type lim_type, Float& lim, Float& abs_err)
{       
    //-----------------------------------------------------------------------
    //                   INITIALIZE PARAMETERS
    //-----------------------------------------------------------------------
    
    // maximum allowed number of iterations with the same estimation of
    // the limit; if this number is reached, then convergence is assumed
    const int max_no_progress_it    = 3;

    // relative tolerance used to detect convergence; if relative error of
    // the estimated limit is less than this number, then convergence is 
    // assumed
    const Float tol_conv            = m_epsilon * Float(100.0);

    // maximum number of iterations allowed
    const int max_it                = m_max_iter;

    // it is assumed, that estimated limit is close to zero if 
    // |lim| <= |abs_err| * tol_zero; tol_zero should be > 1 since error
    // estimator is not precise
    const Float tol_zero            = 3.0;

    // minimum number of iterations required; if this number is not reached
    // then convergence tests are not performed
    const int min_iter              = 8;

    // minimum required decrease of the distance between value of function
    // and estimated limit per iteration
    const Float scale_dif_decrease  = Float(0.7);

    // tolerance used to test divergence
    const Float tol_divergence      = Float(100.0);

    const Float mean_req_accuracy   = matcl::sqrt(m_epsilon);

    //-----------------------------------------------------------------------
    //                   TEST ARGUMENTS
    //-----------------------------------------------------------------------

    function_type f         = f0;

    if (matcl::is_finite(x0) == false)
    {
        if (matcl::is_inf(x0) == true)
        {
            // instead of lim_{x -> Inf} f(x)
            // compute lim_{x -> 0} f(1/x)

            if (x0 > 0)
                lim_type    = limit_type::right;
            else
                lim_type    = limit_type::left;

            f               = get_inverse_function(f);
            x0              = 0.0;
        }
        else
        {
            error_invalid_starting_point(x0);
            return 0;
        }
    };

    //-----------------------------------------------------------------------
    //                   INITIALIZE LOCAL VARIABLES
    //-----------------------------------------------------------------------

    const Float zero        = Float(0.0);
    const Float one         = Float(1.0);

    // get initial point and scaling
    Float x_init;
    Float scaling_init;

    bool ok                 = get_step_params(f, x0, lim_type, x_init, scaling_init);

    if (ok == false)
    {
        // initial point could not be found; exit
        error_initial_point_not_found(x0);
        lim                 = seq_helpers::nan_value<Float>::value(m_precision);
        abs_err             = seq_helpers::inf_value<Float>::value(m_precision);
        return 0;
    }

    m_scaling.initialize(scaling_init);

    Float x                 = x_init;

    // initialize output arguments
    lim                     = m_huge_value;
    abs_err                 = m_huge_value;

    // best relative error of the limit estimator so far
    Float best_err_rel      = m_huge_value;

    // value of function from previous iteration
    Float val_prev          = m_huge_value;

    // value of dx from previous iteration
    Float x_prev            = m_huge_value;

    // |val - val_prev| from previous iteration
    Float best_val_lim_dif  = m_huge_value;

    // estimation of the limit from previous iteration
    Float lim_est_prev      = m_huge_value;

    // true if the limit estimated in previous iteration is close
    // to zero
    bool is_zero_prev       = false;

    // number of iterations with the same estimation of the limit
    int no_progress_it      = 0;

    // interation number with lowest error
    int best_iter           = 0;    

    // minimum required decrease of the distance between value of function
    // and estimated limit
    Float req_dif_decrease  = Float(1.0);

    // value of the function at current point
    Float val               = m_huge_value;

    Float small_x           = matcl::pow(m_epsilon, Float(0.33));

    //-----------------------------------------------------------------------
    //                   MAIN LOOP
    //-----------------------------------------------------------------------
    int iter;

    for (iter = 1; iter <= max_it; ++iter)
    {
        // evaluate the function
        Float x_loc         = x0 + x;

        if (x_loc - x0 == 0.0)
        {
            // x is too small; just quit
            break;
        }
        
        eval_function(f, x_loc, val);

        // function value is invalid; just return
        if (matcl::is_finite(val) == false)
            break;
                
        // update limit estimation

        Float loc_lim, loc_err;
        try
        {
            m_extrapolator->eval(x, val, loc_lim, loc_err);
        }
        catch(std::exception&)
        {
            break;
        }

        Float abs_lim       = matcl::abs(loc_lim);
        bool is_zero        = abs_lim <= loc_err * tol_zero;
        Float rel_err       = (is_zero == false) ? (loc_err / abs_lim) : m_huge_value;        

        // select best estimation
        bool better         = false;
        
        // prefer relative error test; however if limit is close to zero
        // we must use absolute error estimator
        if (is_zero == false && is_zero_prev == false && rel_err < best_err_rel)
            better          = true;
        else if ((is_zero == true || is_zero_prev == true) && loc_err < abs_err)
            better          = true;

        Float val_lim_dif   = matcl::abs(val - loc_lim) + loc_err;
        bool func_conv      = val_lim_dif < req_dif_decrease * best_val_lim_dif;

        (void)small_x;

        // limit estimation error must be lower and function value must be
        // closer to current limit estimation
        if (better == true && func_conv == true)
        {   
            m_scaling.report_progress(loc_err, abs_err, iter);

            // update estimation
            no_progress_it      = 0;
            lim                 = loc_lim;
            abs_err             = std::max(loc_err, m_tiny_value);
            best_err_rel        = rel_err;
            best_iter           = iter;

            best_val_lim_dif    = val_lim_dif + loc_err;
            req_dif_decrease    = one;
        }
        else
        {
            m_scaling.report_no_progress();
        }

        //-----------------------------------------------------------------------
        //                   CONVERGENCE TESTS
        //-----------------------------------------------------------------------

        // we need a few iterations to get meaningful estimations from 
        // the extrapolator
        if (iter >= min_iter)
        {
            // error is very small; assume convergence
            if (is_zero == false && loc_err <= abs_lim * tol_conv)
            {
                break;
            }
            else if (abs_err <= m_tiny_value)
            {
                break;
            }
        
            // value of function is the same in two consecutive iterations;
            // assume convergence
            if (val == val_prev)
                break;

            // limit estimation is the same in two consecutive iterations;
            // we need few additional iterations to detect spurious convergence
            // or decrease estimation error
            if (loc_lim == lim_est_prev)
            {
                no_progress_it  += 1;

                if (no_progress_it >= max_no_progress_it)
                    break;
            }
            else
            {
                no_progress_it  = 0;
            }
        };

        // update scaling
        m_scaling.update(better, x);

        // prepare next step
        is_zero_prev        = is_zero;
        val_prev            = val;
        lim_est_prev        = loc_lim;
        req_dif_decrease    = req_dif_decrease * scale_dif_decrease;
        x_prev              = x;

        x                   = m_scaling.apply(x);
    };

    return iter;
}

//-----------------------------------------------------------------------
//                   SCALING
//-----------------------------------------------------------------------
template<class Float>
void scaling<Float>::set_epsilon(const Float& epsilon)
{
    using seq_helpers::sqr;

    m_small_x   = sqr(sqr(epsilon));
}

template<class Float>
void scaling<Float>::initialize(const Float& initial_scaling)
{
    m_scaling           = initial_scaling;
    m_scaling2          = Float(0.5);
    m_mean_conv_rate    = Float(1.0);
    m_num_slow_iters    = 0;
    m_num_failed_iters  = 0;
}

template<class Float>
Float scaling<Float>::apply(const Float& x) const
{
    return x * m_scaling;
}

template<class Float>
void scaling<Float>::report_progress(const Float& loc_err, const Float& best_err,
                                    int iter)
{
    // convergence rate must be at least mean_conv_rate * min_conv_decrease
    // to classify convergence in given iteration as fast
    const Float min_conv_decrease   = Float(0.3);

    // minimum number of iterations required; if this number is not reached
    // then convergence rate is not updated
    const int min_iter              = 6;

    Float conv_rate         = loc_err / best_err;

    // if estimation error decrease is not large enough, then
    // increase slow convergence iteration counter
    if (conv_rate < m_mean_conv_rate * min_conv_decrease)
        this->report_fast_progress();
    else
        this->report_slow_progress();

    // update average convergence rate
    if (iter > min_iter)
        m_mean_conv_rate    = matcl::sqrt(m_mean_conv_rate * conv_rate);
}

template<class Float>
void scaling<Float>::report_fast_progress()
{
    m_num_slow_iters    = 0;
    m_num_failed_iters  = 0;
}

template<class Float>
void scaling<Float>::report_slow_progress()
{
    m_num_slow_iters    += 1;
    m_num_failed_iters  = 0;
}

template<class Float>
void scaling<Float>::report_no_progress()
{
    m_num_slow_iters    = 0;
    m_num_failed_iters  += 1;
}

template<class Float>
void scaling<Float>::update(bool better, const Float& x)
{
    // maximum allowed number of iterations with slow convergence; 
    // if this number is reached, then scaling factor is increased
    static const 
    int max_slow_iters      = 4;

    // maximum number of iterations without improvement of the
    // estimated limit; if this number is reached, then scaling factor
    // is increased
    static const 
    int max_failed_iters    = 3;

    // test if value can have very large exponent
    static const
    bool huge_exponent      = has_huge_exponent<Float>::value;

    // update scaling
    if (better == true && m_num_slow_iters >= max_slow_iters)
    {
        // we have convergence, but convergence is slow; 
        // faster convergence x -> 0 may speed up limit convergence
        m_scaling       = m_scaling2 * m_scaling;

        if (huge_exponent == true)
        {
            m_scaling2  = m_scaling2 * m_scaling2;

            if (matcl::abs(x) < m_small_x)
                m_scaling2  = m_scaling2 * m_scaling2;
        }
        else
        {
            m_scaling2  = m_scaling2 * matcl::sqrt(m_scaling2);
        };
    }
    else if (m_num_failed_iters >= 2 * max_failed_iters)
    {
        // no progress for a long time; decrease x in order to safe
        // iterations
        m_scaling       = m_scaling2 * m_scaling;

        if (huge_exponent == true)
        {
            m_scaling2  = m_scaling2 * m_scaling2;

            if (matcl::abs(x) < m_small_x)
                m_scaling2  = m_scaling2 * m_scaling2;
        }
        else
        {
            m_scaling2  = m_scaling2 * matcl::sqrt(m_scaling2);
        };
    }
    else if (m_num_failed_iters >= max_failed_iters)
    {
        // no progress and function is well defined for very small
        // aguments; speed up convergence of x to zero in order to
        // safe iterations and possibly improve converge of the limit
        m_scaling       = m_scaling2 * m_scaling;

        if (huge_exponent == true && matcl::abs(x) < m_small_x)
            m_scaling2  = m_scaling2 * m_scaling2;
    }
}

}}

namespace matcl
{

template<class Float>
limit_estimator<Float>::limit_estimator(int precision)
    :m_impl(new details::limit_estimator_impl<Float>
            (matcl::make_extrapolator_wynn_epsilon<Float>(precision), precision))
{}

template<class Float>
limit_estimator<Float>::~limit_estimator()
{};

template<class Float>
Float limit_estimator<Float>::eval(const function_type& f, const Float& x0, limit_type lim_type,
            Float& abs_err)
{
    return m_impl->eval(f, x0, lim_type, abs_err);
}

template<class Float>
void limit_estimator<Float>::set_extrapolator(const extrapolator_type& extr)
{
    m_impl->set_extrapolator(extr);
};

template<class Float>
typename limit_estimator<Float>::extrapolator_type
limit_estimator<Float>::get_extrapolator() const
{
    return m_impl->get_extrapolator();
};

template<class Float>
void limit_estimator<Float>::set_max_iterations(int max_iter)
{
    return m_impl->set_max_iterations(max_iter);
}

template<class Float>
int limit_estimator<Float>::get_max_iterations() const
{
    return m_impl->get_max_iterations();
}

template<class Float>
void limit_estimator<Float>::set_initial_scaling(const Float& scaling)
{
    return m_impl->set_initial_scaling(scaling);
}

template<class Float>
Float limit_estimator<Float>::get_initial_scaling() const
{
    return m_impl->get_initial_scaling();
}

template<class Float>
void limit_estimator<Float>::set_initial_point_scaling(const Float& scaling)
{
    return m_impl->set_initial_point_scaling(scaling);
}

template<class Float>
Float limit_estimator<Float>::get_initial_point_scaling() const
{
    return m_impl->get_initial_point_scaling();
}

template<class Float>
int limit_estimator<Float>::number_iterations() const
{
    return m_impl->number_iterations();
};

template<class Float>
void limit_estimator<Float>::clear(int precision)
{
    return m_impl->clear(precision);
}

template limit_estimator<mp_float>;
template limit_estimator<double>;
template limit_estimator<float>;

}

#pragma warning(pop)