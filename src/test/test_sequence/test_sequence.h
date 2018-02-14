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

#include "matcl-core/config.h"
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-scalar/lib_functions/sequence/limest.h"

#include <string>

namespace matcl { namespace test
{

class test_function
{
    public:
        virtual double      eval(double x) const = 0;
        virtual double      eval_omega(double x) const { return x; };
        virtual double      get_limit() const = 0;
        virtual double      get_offset() const = 0;
        virtual std::string get_info() const = 0;
        virtual bool        is_rational() const = 0;
};

class seq_generator
{
    private:
        double  m_x0;
        double  m_scal;
        double  m_dx;
        bool    m_altern;

    public:
        seq_generator(double x0, double dx, double scal, bool altern);

        double      get(int iter) const;
        std::string get_info() const;
};

class test_sequence
{
    private:
        template<class T>
        using limit_function    = std::function<T (const T&)>;

    public:
        void        test_lim();
        void        test_extrapolation();
        void        test_extrapolation_inf();
        void        test_extrapolation_gen();
        void        test_compile();

    private:
        void        test_extrapol_function(const test_function& f, 
                        const seq_generator& seq);
        void        test_extrapol_function(const test_function& f);

        void        test_extrapol_function_inf(const test_function& f);
        void        test_extrapol_function_inf(const test_function& f, 
                        const seq_generator& seq);

        void        test_extrapol_function_gen(const test_function& f);
        void        test_extrapol_function_gen(const test_function& f, 
                        const seq_generator& seq);

        template<class T>
        void        test_lim_type(int prec);

        template<class T>
        void        test_limit_function(formatted_disp& fd, const limit_function<T>& f, 
                        const T& x, limit_type lim_type, const T& lim, const std::string& func,
                        bool has_limit, int prec);

        template<class Float>
        std::string get_limest_status(bool has_limit, const Float& diff_err, double tol_err);

        std::string get_label(int it);
        int         get_max_iter();
};

class Func_expm1 : public test_function
{
    public:
        double      off;
        precision   prec;

        Func_expm1(precision p, double off_)
            :prec(p), off(off_)
        {};

        double eval(double x) const override
        {
            mp_float xp     = mp_float(x, prec);
            mp_float res    = expm1(xp) / xp - off;
            return res.cast_float();
        }
        
        double get_limit() const override
        {
            return 1.0;
        }
        
        double get_offset() const override
        {
            return off;
        }

        bool is_rational() const override
        {
            return true;
        };

        std::string get_info() const override
        {
            std::ostringstream os;
            os << "expm1(x)/x - " << off;

            return os.str();
        }
};

class Func_rat : public test_function
{
    public:
        double      off;
        precision   prec;

        Func_rat(precision p, double off_)
            :prec(p), off(off_)
        {};

        double eval(double x) const override
        {
            mp_float xp     = mp_float(x, prec);
            mp_float res    = (xp + 2.0) / (xp + 1.0) - off;
            return res.cast_float();
        }
        
        double get_limit() const override
        {
            return 2.0;
        }
        
        double get_offset() const override
        {
            return off;
        }

        bool is_rational() const override
        {
            return true;
        };

        std::string get_info() const override
        {
            std::ostringstream os;
            os << "(x+2)/(x+1) - " << off;

            return os.str();
        }
};

class Func_sqrt : public test_function
{
    public:
        double      off;
        precision   prec;

        Func_sqrt(precision p, double off_)
            :prec(p), off(off_)
        {};

        double eval(double x) const override
        {
            mp_float xp     = mp_float(x, prec);
            mp_float res    = sqrt(abs(xp)) - off;
            return res.cast_float();
        }
        
        double get_limit() const override
        {
            return 0.0;
        }
        
        double get_offset() const override
        {
            return off;
        }

        std::string get_info() const override
        {
            std::ostringstream os;
            os << "sqrt(x) - " << off;

            return os.str();
        }

        bool is_rational() const override
        {
            return false;
        };
};

class Func_rat_inf : public test_function
{
    public:
        double      off;
        precision   prec;
        double      mult_sqrt;

        Func_rat_inf(precision p, double off_, double mult_sqrt_)
            :prec(p), off(off_), mult_sqrt(mult_sqrt_)
        {};

        double eval(double x) const override
        {
            mp_float xp     = mp_float(x, prec);
            mp_float xp2    = pow(xp, 2.0);
            mp_float xsr    = sqrt(xp);
            mp_float res    = (xp2 + 1.0 + mult_sqrt * xsr) / xp2 - off;

            return res.cast_float();
        }
        
        double get_limit() const override
        {
            return 1.0;
        }
        
        double get_offset() const override
        {
            return off;
        }

        bool is_rational() const override
        {
            return true;
        };

        std::string get_info() const override
        {
            std::ostringstream os;
            os << "(x^2 + 1 + " << mult_sqrt << " * sqrt(x))/x^2 - " << off;
            return os.str();
        }
};

class Func_rat_sqrt : public test_function
{
    public:
        double      off;
        precision   prec;

        Func_rat_sqrt(precision p, double off_)
            :prec(p), off(off_)
        {};

        double eval(double x) const override
        {
            mp_float xp     = mp_float(x, prec);
            mp_float xsq    = sqrt(abs(xp));
            mp_float res    = xsq * (xp + 2.0) / (xp + 1.0) - off;

            return res.cast_float();
        }
        
        double eval_omega(double x) const 
        { 
            mp_float xp     = mp_float(x, prec);
            mp_float xsq    = sqrt(abs(xp));
            mp_float res    = xsq;

            return res.cast_float();
        };

        double get_limit() const override
        {
            return 0.0;
        }
        
        double get_offset() const override
        {
            return off;
        }

        bool is_rational() const override
        {
            return false;
        };

        std::string get_info() const override
        {
            std::ostringstream os;
            os << "sqrt(|x|) * (x+2)/(x+1) - " << off;
            return os.str();
        }
};

}}
