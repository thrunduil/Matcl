/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas-lapack/level1/level1_axpby.h"

namespace matcl { namespace level1
{

template<class TY, class TX, class TA, class TB, Integer Rows>
struct axpby_test<TY, TX, TA, TB, Rows, details::true_t>
{
    static void eval(TY* Y1, const TX* X1, Integer rows, const TA& a, const TB& b)
    {        
        if (b == TB(0.0))
        {
            if (a == TA(0.0))
            {
                //Y = 0
                return set_val<TY,Rows>::eval(Y1, rows, TY(0.0));
            }
            else if ( a == TA(1.0) )
            {
                //Y = X
                return copy<TY, TX, Rows>::eval(Y1, X1, rows);
            }
            else if ( a == TA(-1.0) )
            {
                //Y = -X
                return mx<TY, TX, Rows>::eval(Y1, X1, rows);
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X
                return ax<TY, TX, TAR, Rows>::eval(Y1, X1, rows, details::eval_real<TA>::eval(a));
            }
            else
            {
                //Y = a*X
                return ax<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
            };
        }
        else if ( b == TB(1.0) )
        {            
            if (a == TA(0.0))
            {
                //Y = Y
                return;
            }
            else if ( a == TA(1.0) )
            {
                //Y = X + Y
                return ypx<TY, TX, Rows>::eval(Y1, X1, rows);
            }
            else if (a == TA(-1.0))
            {
                //Y = Y-X
                return ymx<TY, TX, Rows>::eval(Y1, X1, rows);
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X + Y
                return axpy<TY, TX, TAR, Rows>::eval(Y1, X1, rows, details::eval_real<TA>::eval(a));
            }
            else
            {
                //Y = a*X + Y
                return axpy<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
            };
        }
        else if (b == TB(-1.0))
        {
            if (a == TA(0.0))
            {
                //Y = -Y
                return my<TY, Rows>::eval(Y1, rows);
            }
            else if ( a == TA(1.0) )
            {
                //Y = X-Y
                return xmy<TY, TX, Rows>::eval(Y1, X1, rows);
            }
            else if (a == TA(-1.0))
            {
                //Y = -X-Y
                return ypxm<TY, TX, Rows>::eval(Y1, X1, rows);
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X - Y
                return axmy<TY, TX, TAR, Rows>::eval(Y1, X1, rows, details::eval_real<TA>::eval(a));
            }
            else
            {
                //Y = a*X - Y
                return axmy<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
            };
        }
        else if (details::is_real<TB>::eval(b))
        {
            using TBR   = typename details::real_type<TB>::type;
            TBR br      = details::eval_real<TB>::eval(b);

            if (a == TA(0.0))
            {
                //Y = bY
                return ay<TY, TBR, Rows>::eval(Y1, rows, br);
            }
            else if ( a == TA(1.0) )
            {
                //Y = X+bY
                return xpby<TY, TX, TBR, Rows>::eval(Y1, X1, rows, br);
            }
            else if (a == TA(-1.0))
            {
                //Y = -X+bY
                return mxpby<TY, TX, TBR, Rows>::eval(Y1, X1, rows, br);
            }
            else
            {
                if (a == br)
                {
                    //Y = a*(X + Y)
                    return xpya<TY, TX, TBR, Rows>::eval(Y1, X1, rows, br);
                }
                else if (a == -br)
                {
                    //Y = a*(X - Y)
                    return xmya<TY, TX, TBR, Rows>::eval(Y1, X1, rows, br);
                }
                else if (details::is_real<TA>::eval(a))
                {
                    using TAR   = typename details::real_type<TA>::type;

                    //Y = a*X + bY
                    return axpby<TY, TX, TAR, TBR, Rows>::eval(Y1, X1, rows, 
                                                    details::eval_real<TA>::eval(a), br);
                }
                else
                {
                    //Y = a*X + bY
                    return axpby<TY, TX, TA, TBR, Rows>::eval(Y1, X1, rows, a, br);
                };
            };
        }
        else
        {
            if (a == TA(0.0))
            {
                //Y = bY
                return ay<TY, TB, Rows>::eval(Y1, rows, b);
            }
            else if ( a == TA(1.0) )
            {
                //Y = X+bY
                return xpby<TY, TX, TB, Rows>::eval(Y1, X1, rows, b);
            }
            else if (a == TA(-1.0))
            {
                //Y = -X+bY
                return mxpby<TY, TX, TB, Rows>::eval(Y1, X1, rows, b);
            }
            else
            {               
                if (a == b)
                {
                    //Y = a*(X + Y)
                    return xpya<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
                }
                else if (a == -b)
                {
                    //Y = a*(X - Y)
                    return xmya<TY, TX, TA, Rows>::eval(Y1, X1, rows, a);
                }
                else if (details::is_real<TA>::eval(a))
                {
                    using TAR   = typename details::real_type<TA>::type;

                    //Y = a*X + bY
                    return axpby<TY, TX, TAR, TB, Rows>::eval(Y1, X1, rows, 
                                                    details::eval_real<TA>::eval(a), b);
                }
                else
                {
                    //Y = a*X + bY
                    return axpby<TY, TX, TA, TB, Rows>::eval(Y1, X1, rows, a, b);
                };
            };
        };
    };

    static void eval(TY* Y1, Integer Y_step, const TX* X1, Integer X_step, Integer rows, 
                     const TA& a, const TB& b)
    {        
        if (b == TB(0.0))
        {
            if (a == TA(0.0))
            {
                //Y = 0
                return set_val<TY,Rows>::eval(Y1, Y_step, rows, TY(0.0));
            }
            else if ( a == TA(1.0) )
            {
                //Y = X
                return copy<TY, TX, Rows>::eval(Y1, Y_step, X1, X_step, rows);
            }
            else if ( a == TA(-1.0) )
            {
                //Y = -X
                return mx<TY, TX, Rows>::eval(Y1, Y_step, X1, X_step, rows);
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X
                return ax<TY, TX, TAR, Rows>::eval(Y1, Y_step, X1, X_step, rows, 
                                                   details::eval_real<TA>::eval(a));
            }
            else
            {
                //Y = a*X
                return ax<TY, TX, TA, Rows>::eval(Y1, Y_step, X1, X_step, rows, a);
            };
        }
        else if ( b == TB(1.0) )
        {            
            if (a == TA(0.0))
            {
                //Y = Y
                return;
            }
            else if ( a == TA(1.0) )
            {
                //Y = X + Y
                return ypx<TY, TX, Rows>::eval(Y1, Y_step, X1, X_step, rows);
            }
            else if (a == TA(-1.0))
            {
                //Y = Y-X
                return ymx<TY, TX, Rows>::eval(Y1, Y_step, X1, X_step, rows);
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X + Y
                return axpy<TY, TX, TAR, Rows>::eval(Y1, Y_step, X1, X_step, rows, 
                                                     details::eval_real<TA>::eval(a));
            }
            else
            {
                //Y = a*X + Y
                return axpy<TY, TX, TA, Rows>::eval(Y1, Y_step, X1, X_step, rows, a);
            };
        }
        else if (b == TB(-1.0))
        {
            if (a == TA(0.0))
            {
                //Y = -Y
                return my<TY, Rows>::eval(Y1, Y_step, rows);
            }
            else if ( a == TA(1.0) )
            {
                //Y = X-Y
                return xmy<TY, TX, Rows>::eval(Y1, Y_step, X1, X_step, rows);
            }
            else if (a == TA(-1.0))
            {
                //Y = -X-Y
                return ypxm<TY, TX, Rows>::eval(Y1, Y_step, X1, X_step, rows);
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X - Y
                return axmy<TY, TX, TAR, Rows>::eval(Y1, Y_step, X1, X_step, rows, 
                                                     details::eval_real<TA>::eval(a));
            }
            else
            {
                //Y = a*X - Y
                return axmy<TY, TX, TA, Rows>::eval(Y1, Y_step, X1, X_step, rows, a);
            };
        }
        else if (details::is_real<TB>::eval(b))
        {
            using TBR   = typename details::real_type<TB>::type;
            TBR br      = details::eval_real<TB>::eval(b);

            if (a == TA(0.0))
            {
                //Y = bY
                return ay<TY, TBR, Rows>::eval(Y1, Y_step, rows, br);
            }
            else if ( a == TA(1.0) )
            {
                //Y = X+bY
                return xpby<TY, TX, TBR, Rows>::eval(Y1, Y_step, X1, X_step, rows, br);
            }
            else if (a == TA(-1.0))
            {
                //Y = -X+bY
                return mxpby<TY, TX, TBR, Rows>::eval(Y1, Y_step, X1, X_step, rows, br);
            }
            else
            {
                if (a == br)
                {
                    //Y = a*(X + Y)
                    return xpya<TY, TX, TBR, Rows>::eval(Y1, Y_step, X1, X_step, rows, br);
                }
                else if (a == -br)
                {
                    //Y = a*(X - Y)
                    return xmya<TY, TX, TBR, Rows>::eval(Y1, Y_step, X1, X_step, rows, br);
                }
                else if (details::is_real<TA>::eval(a))
                {
                    using TAR   = typename details::real_type<TA>::type;

                    //Y = a*X + bY
                    return axpby<TY, TX, TAR, TBR, Rows>::eval(Y1, Y_step, X1, X_step, rows, 
                                                    details::eval_real<TA>::eval(a), br);
                }
                else
                {
                    //Y = a*X + bY
                    return axpby<TY, TX, TA, TBR, Rows>::eval(Y1, Y_step, X1, X_step, rows, a, br);
                };
            };
        }
        else
        {
            if (a == TA(0.0))
            {
                //Y = bY
                return ay<TY, TB, Rows>::eval(Y1, Y_step, rows, b);
            }
            else if ( a == TA(1.0) )
            {
                //Y = X+bY
                return xpby<TY, TX, TB, Rows>::eval(Y1, Y_step, X1, X_step, rows, b);
            }
            else if (a == TA(-1.0))
            {
                //Y = -X+bY
                return mxpby<TY, TX, TB, Rows>::eval(Y1, Y_step, X1, X_step, rows, b);
            }
            else
            {               
                if (a == b)
                {
                    //Y = a*(X + Y)
                    return xpya<TY, TX, TA, Rows>::eval(Y1, Y_step, X1, X_step, rows, a);
                }
                else if (a == -b)
                {
                    //Y = a*(X - Y)
                    return xmya<TY, TX, TA, Rows>::eval(Y1, Y_step, X1, X_step, rows, a);
                }
                else if (details::is_real<TA>::eval(a))
                {
                    using TAR   = typename details::real_type<TA>::type;

                    //Y = a*X + bY
                    return axpby<TY, TX, TAR, TB, Rows>::eval(Y1, Y_step, X1, X_step, rows, 
                                                    details::eval_real<TA>::eval(a), b);
                }
                else
                {
                    //Y = a*X + bY
                    return axpby<TY, TX, TA, TB, Rows>::eval(Y1, Y_step, X1, X_step, rows, a, b);
                };
            };
        };
    };
};

//-----------------------------------------------------------------------
//                  				axpby_test_mat
//-----------------------------------------------------------------------
template<bool Select, class TY, class TX, class TA, class TB, Integer Rows, 
        Integer Cols, Integer Continuous>
struct axpby_test_mat<Select, TY, TX, TA, TB, Rows, Cols, Continuous, details::true_t>
{
    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols, 
                     const TA& a, const TB& b)
    {
        if (b == TB(0.0))
        {
            if (a == TA(0.0))
            {
                return eval_mat_func_Y<Select, TY, Rows, Cols, Continuous, details::func_set_val_zero<TY>>
                            ::eval(Y, Y_ld, rows, cols, details::func_set_val_zero<TY>());
            }
            else if ( a == TA(1.0) )
            {
                //Y = X
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_copy<TY,TX>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_copy<TY,TX>());
            }
            else if ( a == TA(-1.0) )
            {
                //Y = -X
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_mx<TY,TX>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_mx<TY,TX>());
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ax<TY,TX,TAR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                                   details::func_ax<TY,TX,TAR>(details::eval_real<TA>::eval(a)));
            }
            else
            {
                //Y = a*X
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ax<TY,TX,TA>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ax<TY,TX,TA>(a));
            };
        }        
        else if ( b == TB(1.0) )
        {
            if (a == TA(0.0))
            {
                //Y = Y
                return;
            }
            else if ( a == TA(1.0) )
            {
                //Y = X + Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ypx<TY,TX>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ypx<TY,TX>());
            }
            else if (a == TA(-1.0))
            {
                //Y = Y-X
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ymx<TY,TX>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ymx<TY,TX>());                
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X + Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpy<TY,TX,TAR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                                   details::func_axpy<TY,TX,TAR>(details::eval_real<TA>::eval(a)));
            }
            else
            {
                //Y = a*X + Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpy<TY,TX,TA>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_axpy<TY,TX,TA>(a));
            };
        }
        else if (b == TB(-1.0))
        {
            if (a == TA(0.0))
            {
                //Y = -Y
                return eval_mat_func_Y<Select, TY, Rows, Cols, Continuous, details::func_my<TY>>
                            ::eval(Y, Y_ld, rows, cols, details::func_my<TY>());
            }
            else if ( a == TA(1.0) )
            {
                //Y = X-Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xmy<TY,TX>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xmy<TY,TX>());
            }
            else if (a == TA(-1.0))
            {
                //Y = -X-Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_ypxm<TY,TX>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_ypxm<TY,TX>());
            }
            else if (details::is_real<TA>::eval(a))
            {
                using TAR   = typename details::real_type<TA>::type;
                //Y = real(a)*X - Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axmy<TY,TX,TAR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                                   details::func_axmy<TY,TX,TAR>(details::eval_real<TA>::eval(a)));
            }
            else
            {
                //Y = a*X - Y
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axmy<TY,TX,TA>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_axmy<TY,TX,TA>(a));
            };
        }
        else if (details::is_real<TB>::eval(b))
        {
            using TBR   = typename details::real_type<TB>::type;
            TBR br      = details::eval_real<TB>::eval(b);

            if (a == TA(0.0))
            {
                //Y = bY
                return eval_mat_func_Y<Select, TY, Rows, Cols, Continuous, details::func_ay<TY,TBR>>
                            ::eval(Y, Y_ld, rows, cols, details::func_ay<TY,TBR>(br));
            }
            else if ( a == TA(1.0) )
            {
                //Y = X+bY
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xpby<TY,TX,TBR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xpby<TY,TX,TBR>(br));
            }
            else if (a == TA(-1.0))
            {
                //Y = -X+bY
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_mxpby<TY,TX,TBR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_mxpby<TY,TX,TBR>(br));
            }
            else
            {
                if (a == br)
                {
                    //Y = a*(X + Y)
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xpya<TY,TX, TBR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xpya<TY,TX,TBR>(br));
                }
                else if (a == -br)
                {
                    //Y = a*(X - Y)
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xmya<TY,TX,TBR>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xmya<TY,TX,TBR>(br));                    
                }
                else if (details::is_real<TA>::eval(a))
                {
                    using TAR   = typename details::real_type<TA>::type;

                    //Y = a*X + bY
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpby<TY,TX,TAR,TBR>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                               details::func_axpby<TY,TX,TAR,TBR>(details::eval_real<TA>::eval(a),br));
                }
                else
                {
                    //Y = a*X + bY
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpby<TY,TX,TA,TBR>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_axpby<TY,TX,TA,TBR>(a,br));
                };
            };
        }
        else
        {
            if (a == TA(0.0))
            {
                //Y = bY
                return eval_mat_func_Y<Select, TY, Rows, Cols, Continuous, details::func_ay<TY, TB>>
                            ::eval(Y, Y_ld, rows, cols, details::func_ay<TY,TB>(b));
            }
            else if ( a == TA(1.0) )
            {
                //Y = X+bY
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xpby<TY,TX,TB>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xpby<TY,TX,TB>(b));
            }
            else if (a == TA(-1.0))
            {
                //Y = -X+bY
                return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_mxpby<TY,TX,TB>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_mxpby<TY,TX,TB>(b));
            }
            else
            {
                if (a == b)
                {
                    //Y = a*(X + Y)
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xpya<TY,TX,TA>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xpya<TY,TX,TA>(a));
                }
                else if (a == -b)
                {
                    //Y = a*(X - Y)
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_xmya<TY,TX,TA>>
                            ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_xmya<TY,TX,TA>(a));
                }
                else if (details::is_real<TA>::eval(a))
                {
                    using TAR   = typename details::real_type<TA>::type;

                    //Y = a*X + bY
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpby<TY,TX,TAR,TB>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, 
                               details::func_axpby<TY,TX,TAR,TB>(details::eval_real<TA>::eval(a),b));
                }
                else
                {
                    //Y = a*X + bY
                    return eval_mat_func_XY<Select, TY, TX, Rows, Cols, Continuous, details::func_axpby<TY,TX,TA,TB>>
                        ::eval(Y, Y_ld, X, X_ld, rows, cols, details::func_axpby<TY,TX,TA,TB>(a,b));
                };
            };
        };
    };
};

}}