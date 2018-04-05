/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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
#include <sstream>

namespace matcl { namespace opt { namespace disp
{

// terminal width used in printing
class terminal_width : public option_base<Integer, terminal_width>
{
    private:
        using base_type             = option_base<Integer, terminal_width>;
        using opt_type              = optional<Integer>;

    public:
        terminal_width()            : base_type() {};
        terminal_width(opt_type x)  : base_type(x) {};

        static void config()
        {
            m_description       = "terminal width used in printing; set -1 "
                                  "to take width from output stream"; 
            m_default_value     = -1;
            m_validator         = [](opt_type x)->opt_type
                                {
                                    if (x)
                                    {
                                        if (x.value() <= -1)
                                            return opt_type(-1);
                                        else if (x.value() <= 20)
                                            return opt_type(20);
                                    }
                                    return x;
                                };
        };
};

// maximum number of columns printed
class max_cols : public option_base<Integer, max_cols>
{
    private:
        using base_type             = option_base<Integer, max_cols>;
        using opt_type              = optional<Integer>;

    public:
        max_cols()                  : base_type() {};
        max_cols(opt_type x)        : base_type(x) {};

        static void config()
        {
            m_description           = "maximum number of columns printed";
            m_default_value         = 20;
            m_validator             = validator_positive<Integer>();
        };
};

// maximum number of rows printed
class max_rows : public option_base<Integer, max_rows>
{
    private:
        using base_type             = option_base<Integer, max_rows>;
        using opt_type              = optional<Integer>;

    public:
        max_rows()                  : base_type() {};
        max_rows(opt_type x)        : base_type(x) {};

        static void config()
        {
            m_description           = "maximum number of rows printed";
            m_default_value         = 100;
            m_validator             = validator_positive<Integer>();
        };
};

// maximum number of nonzero elements printed
class max_nnz : public option_base<Integer, max_nnz>
{
    private:
        using base_type             = option_base<Integer, max_nnz>;
        using opt_type              = optional<Integer>;

    public:
        max_nnz()                   : base_type() {};
        max_nnz(opt_type x)         : base_type(x) {};

        static void config()
        {
            m_description           = "maximum number of nonzero elements printed";
            m_default_value         = 200;
            m_validator             = validator_positive<Integer>();
        };
};

// number of significant digits printed for floating point values
class precision : public option_base<Integer, precision>
{
    private:
        using base_type             = option_base<Integer, precision>;
        using opt_type              = optional<Integer>;

    public:
        precision()                 : base_type() {};
        precision(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description           = "number of significant digits printed for floating point values";
            m_default_value         = 5;
            m_validator             = validator_nonnegative<Integer>();
        };
};

// do add limit on number of rows and columns printed for sparse matrices?
class restrict_sparse_matrix_size : public option_base<bool, restrict_sparse_matrix_size>
{
    private:
        using base_type         = option_base<bool, restrict_sparse_matrix_size>;

    public:
        restrict_sparse_matrix_size()                   : base_type() {};
        restrict_sparse_matrix_size(optional<bool> x)   : base_type(x) {};

        static void config()
        {
            m_description       = "do add limit on number of rows and columns printed "
                                  "for sparse matrices?";
            m_default_value     = false;
        };
};

// do display zero elements?
class display_zero : public option_base<bool, display_zero>
{
    private:
        using base_type         = option_base<bool, display_zero>;

    public:
        display_zero()                   : base_type() {};
        display_zero(optional<bool> x)   : base_type(x) {};

        static void config()
        {
            m_description       = "do display zero elements?";
            m_default_value     = false;
        };
};

// display lower traingle for symmetric and hermitian matrices?
class ignore_lower_triangle : public option_base<bool, ignore_lower_triangle>
{
    private:
        using base_type         = option_base<bool, ignore_lower_triangle>;

    public:
        ignore_lower_triangle()                   : base_type() {};
        ignore_lower_triangle(optional<bool> x)   : base_type(x) {};

        static void config()
        {
            m_description       = "display lower traingle for symmetric and hermitian matrices?";
            m_default_value     = true;
        };
};

// how to display scalars and sparse matrices?; see disp_mode
class disp_mode : public option_base<Integer, disp_mode>
{
    private:
        using base_type             = option_base<Integer, disp_mode>;
        using opt_type              = optional<Integer>;

    public:
        disp_mode()                 : base_type() {};
        disp_mode(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description           = "how to display scalars and sparse matrices?";
            m_default_value         = (Integer)matcl::disp_mode::standard;
            m_validator             = validator_enum((Integer)matcl::disp_mode::size);
        };
};

// do display only matrix header?
class header_only : public option_base<bool, header_only>
{
    private:
        using base_type         = option_base<bool, header_only>;

    public:
        header_only()                   : base_type() {};
        header_only(optional<bool> x)   : base_type(x) {};

        static void config()
        {
            m_description       = "display only matrix header?";
            m_default_value     = false;
        };
};

// distance between row separators for dense and band matrices;
// negative value or zero implies, that row separators are not
// displayed
class row_block : public option_base<Integer, row_block>
{
    private:
        using base_type             = option_base<Integer, row_block>;
        using opt_type              = optional<Integer>;

    public:
        row_block()                 : base_type() {};
        row_block(opt_type x)       : base_type(x) {};

        static void config()
        {
            m_description           = "distance between row separators";
            m_default_value         = 5;
        };
};

};};};
