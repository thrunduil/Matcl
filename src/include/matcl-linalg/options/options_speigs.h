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

#include "matcl-core/options/matcl_options.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include <sstream>

namespace matcl { namespace opt { namespace speigs
{

// options for speigs_decomposition

// return nonconvergent eigenvalues?
class return_nonconvergent : public option_base<bool, return_nonconvergent>
{
    private:
        using base_type             = option_base<bool, return_nonconvergent>;
        using opt_type              = optional<bool>;

    public:
        return_nonconvergent()              : base_type() {};
        return_nonconvergent(opt_type x)    : base_type(x) {};

        static void config()
        {
            m_description       = "if true, return also eigenvalues which have not converged";
            m_default_value     = false;
        };
};

// tolerance used to test convergence of eigenvalues;
// the relative accuracy of the i-th Ritz value l is considered acceptable if
// bounds <= tol * abs( l ), where bounds is the estimated precision of l
// default value: tol = eps; set negative tol to get default tolerance
class tol : public option_base<Real, tol>
{
    private:
        using base_type         = option_base<Real, tol>;
        using opt_type          = optional<Real>;

    public:
        tol()               : base_type() {};
        tol(opt_type x)     : base_type(x) {};

        static void config()
        {
            m_description       = "tolerance used to test convergence of eigenvalues";
            m_default_value     = Real(-1.0);
            m_validator         = validator_finite<Real>();
        };
};

// maximum number of iterations
class maxit : public option_base<Integer, maxit>
{
    private:
        using base_type         = option_base<Integer, maxit>;
        using opt_type          = optional<Integer>;

    public:
        maxit()                 : base_type() {};
        maxit(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description       = "maximum number of iterations";
            m_default_value     = 1000;
            m_validator         = validator_positive<Integer>();
        };
};

// how many Arnoldi vectors are generated at each iteration. There are always
// generated at least k + 3 Arnoldi vectors, where k is the number of eigenvalues
// to find. In order to obtain faster convergence, higher number is recommended,
// for example 2 * k + 1, however optimal number of Arnoldi vectors can only be
// determined empirically; if this value is negative than default value is used
class n_arnoldi : public option_base<Integer, n_arnoldi>
{
    private:
        using base_type         = option_base<Integer, n_arnoldi>;
        using opt_type          = optional<Integer>;

    public:
        n_arnoldi()             : base_type() {};
        n_arnoldi(opt_type x)   : base_type(x) {};

        static void config()
        {
            m_description       = "number of Arnoldi vectors generated in each iteration";
            m_default_value     = -1;
        };
};

};};};
