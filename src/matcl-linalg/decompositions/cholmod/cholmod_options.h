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

#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace details
{

/// represent options used by the cholmod function
class cholmod_options
{
    public:
        /// correction algorithm used
        enum correction_alg
        {
            E_NONE,     /// no correction, factorization stops if the
                        /// matrix appeared to be indefinite. 
            E_GMW,      /// Gill, Murray, and Wright (GMW) algorithm with
                        /// 2 phase strategy as in Schnabel and Eskow, proposed by Feng.
            E_SE        /// Schnabel and Eskow (SE99) algorithm, default
        };

        /// type of correction
        enum correction_type
        {
            TYPE_I,     /// elements d_k on the diagonal are increased to max (|d_k|, delta)
            TYPE_II     /// elements d_k on the diagonal are increased to max (d_k, delta); 
                        /// default
        };

        /// pivoting strategy in the correction phase (second phase)
        enum pivot_type
        {
            PIV_NONE,           /// no pivoting
            PIV_DIAG,           /// diagonal element with maximal value is selected
            PIV_ABS_DIAG,       /// diagonal element with maximal absolute value is selected
            PIV_GERSHGORIN      /// diagonal element with maximal lower Gerschgorin bound;
                                /// default, available only if E_SE correction is used
        };

    private:
        using self_type = cholmod_options;

    private:
        Real                    m_tol;
        Real                    m_mu;
        correction_alg          m_corr_alg;
        correction_type         m_corr_type;        
        bool                    m_nondecreasing_strategy;
        pivot_type              m_piv_type;

    public:
        /// create default options
        cholmod_options();

        /// get tolerance
        Real                   tol() const                     { return m_tol; };

        /// get mu
        Real                    mu() const                      { return m_mu; };

        /// get correction type
        correction_type         corr_type() const               { return m_corr_type; };

        /// get correction algorithm
        correction_alg          corr_alg() const                { return m_corr_alg; };        

        /// true if nondecreasing strategy is used
        bool                    nondecreasing_strategy() const  { return m_nondecreasing_strategy; };

        /// get pivoting strategy
        pivot_type              piv_type() const                { return m_piv_type; };

        /// set tolerance; correction algorithm is performed if maximal
        /// element on diagonal of the submatrix not factorized yet is lower
        /// than val; if val <= 0, then default tolerance is used
        void                    tol(Real val)                   { m_tol = val; };

        /// set mu (???)
        void                    mu(Real val);

        /// set correction type
        void                    corr_type(correction_type val)  { m_corr_type = val; };

        /// set correction algorithm
        void                    corr_alg(correction_alg val)    { m_corr_alg = val; };

        /// if v is true then nondecreasing strategy is used;
        /// in this strategy diagonal corrections never decrease; default
        void                    nondecreasing_strategy(bool v)  { m_nondecreasing_strategy = v; };

        /// set pivoting strategy
        void                    piv_type(pivot_type piv)        { m_piv_type = piv; };

        /// check consistency of options; throw exception if options
        /// are not consistent
        void                    check() const;
};


};};