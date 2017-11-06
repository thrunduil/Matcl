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

#include "test_simd_config.h"
#include "matcl-core/config.h"
#include "matcl-scalar/IO/formatted_disp.h"

#include <string>

namespace matcl { namespace test
{

void test_performance_real();
void test_values_real();

class test_simd
{
    private:
        std::string m_instr_tag;
        bool        m_test_values;

    public:
        test_simd(bool test_values);

        void    make_ternary();
        void    make_binary();
        void    make_unary();

    private:
        int     get_N() const;
        int     get_M() const;

        template<class T>
        bool    test_equal(int size, const T* res, const T* res_gen, double max_dist, double& dist);

        template<class T>
        bool    test_equal(const T& res, const T& res_gen, double max_dist, double& dist);

        template<class T>
        void    test_functions();

        template<class T>
        void    test_functions_bin();

        template<class T>
        void    test_functions_3();

        template<class T, class Func>
        void    test_function(formatted_disp& fd, int size, int n_rep, const T* in, 
                    T* out, T* out_gen);

        template<class T, class Func>
        void    test_function_bin(formatted_disp& fd, int size, int n_rep, const T* in_1, 
                    const T* in_2, T* out, T* out_gen);

        template<class T, class Func>
        void    test_function_3(formatted_disp& fd, int size, int n_rep, const T* in_1, 
                    const T* in_2, const T* in_3, T* out, T* out_gen);

        template<class T, class Func>
        void    test_function_block(formatted_disp& fd, int size, int n_rep, const T* in, 
                    T* out, T* out_gen);

        template<class T, class Simd_type, class Func>
        double  test_function_simd(int size, int n_rep, const T* in, T* out);

        template<class T, class Simd_type, class Func>
        double  test_function_bin_simd(int size, int n_rep, const T* in1, const T* in2, T* out);

        template<class T, class Simd_type, class Func>
        double  test_function_3_simd(int size, int n_rep, const T* in1, const T* in2, const T* in3, T* out);

        template<class T, class Func>
        double  test_function_generic(int size, int n_rep, const T* in, T* out);
};

// missing scalar functions
template<class T>
inline T reverse(const T& x)    { return x; };

struct Func_reverse
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return reverse(x); 
    }

    static std::string name()
    { 
        return "reverse"; 
    };
};

struct Func_uminus
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return -x; 
    }

    static std::string name()
    { 
        return "uminus"; 
    };
};

struct Func_any_inf
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(any_inf(x)); 
    }

    static std::string name()
    { 
        return "any_inf"; 
    };
};

struct Func_any_nan
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(any_nan(x)); 
    }

    static std::string name()
    { 
        return "any_nan"; 
    };
};

struct Func_conj
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return conj(x); 
    }

    static std::string name()
    { 
        return "conj"; 
    };
};

struct Func_abs
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return abs(x); 
    }

    static std::string name()
    { 
        return "abs"; 
    };
};

struct Func_sum_all
{
    template<class T>    
    force_inline static T eval(const T& x) 
    { 
        return T(sum_all(x)); 
    }

    static std::string name()
    { 
        return "sum_all"; 
    };
};

struct Func_sqrt
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return sqrt(x); 
    }

    static std::string name()
    { 
        return "sqrt"; 
    };
};

struct Func_round
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return round(x); 
    }

    static std::string name()
    { 
        return "round"; 
    };
};

struct Func_floor
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return floor(x); 
    }
    
    static std::string name()
    { 
        return "floor"; 
    };
};

struct Func_ceil
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return ceil(x); 
    }
    
    static std::string name() 
    { 
        return "ceil"; 
    };
};

struct Func_trunc
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return trunc(x); 
    }
    
    static std::string name()
    { 
        return "trunc"; 
    };
};

struct Func_mult
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 * x2; 
    }

    static std::string name()
    { 
        return "mult"; 
    };
};

struct Func_mult_RC
{
    template<class T1, class T2>    
    force_inline static T2 eval(const T1& x1, const T2& x2)
    { 
        return x1 * x2; 
    }

    static std::string name()
    { 
        return "mult rc"; 
    };
};

struct Func_mult_CR
{
    template<class T1, class T2>    
    force_inline static T1 eval(const T1& x1, const T2& x2)
    { 
        return x1 * x2; 
    }

    static std::string name()
    { 
        return "mult cr"; 
    };
};

struct Func_div
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 / x2; 
    }

    static std::string name()
    { 
        return "div"; 
    };
};

struct Func_div_RC
{
    template<class T1, class T2>    
    force_inline static T2 eval(const T1& x1, const T2& x2)
    { 
        return x1 / x2; 
    }

    static std::string name()
    { 
        return "div rc"; 
    };
};

struct Func_div_CR
{
    template<class T1, class T2>    
    force_inline static T1 eval(const T1& x1, const T2& x2)
    { 
        return x1 / x2; 
    }

    static std::string name()
    { 
        return "div cr"; 
    };
};

struct Func_plus
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 + x2; 
    }
    
    static std::string name()
    { 
        return "plus"; 
    };
};

struct Func_minus
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 - x2; 
    }

    static std::string name()
    { 
        return "minus"; 
    };
};

struct Func_sub_add
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return sub_add(x1, x2); 
    }
    
    static std::string name()
    { 
        return "sub_add"; 
    };
};

struct Func_max
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return max(x1, x2); 
    }
    
    static std::string name()
    { 
        return "max"; 
    };
};

struct Func_min
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return min(x1, x2); 
    }

    static std::string name()
    { 
        return "min"; 
    };
};

struct Func_eeq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return eeq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "eeq"; 
    };
};

struct Func_neq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return neq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "neq"; 
    };
};

struct Func_leq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return leq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "leq"; 
    };
};

struct Func_geq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return geq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "geq"; 
    };
};

struct Func_lt
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return lt(x1, x2); 
    }

    static std::string name()
    { 
        return "lt"; 
    };
};

struct Func_gt
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return gt(x1, x2); 
    }

    static std::string name()
    { 
        return "gt"; 
    };
};

struct Func_fma
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fma(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fma"; 
    };
};

struct Func_fms
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fms(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fms"; 
    };
};

}}
