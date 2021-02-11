/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-simd/details/math/impl/simd_sincos_double.h"
#include "matcl-simd/details/math/impl/simd_sincos_float.h"

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                              DOUBLE
//-----------------------------------------------------------------------

const double simd_sincos_table_double_data::poly_sin[] = 
{
    -1.666666666666665916382629555125560e-01,
    8.333333333328349062140114618895757e-03,
    -1.984126983386195863616397631859518e-04,
    2.755731488022942163032964623184407e-06,
    -2.505091491904481811209312506444584e-08,
    1.590445674389578325979670578185607e-10,
};


const double simd_sincos_table_double_data::poly_cos[] = 
{                                    
    4.166666666666659896134445436786686e-02,
    -1.388888888887487817166513568166090e-03,
    2.480158728990971958442638966674033e-05,
    -2.755731443473966320824186766381266e-07,
    2.087572839454679266204771377713960e-09,
    -1.135961908867298976458649289359009e-11,
};

const double simd_sincos_table_double_data::poly_sincos[] =
{
    -1.666666666666665916382629555125560e-01,   4.166666666666659896134445436786686e-02,
    8.333333333328349062140114618895757e-03,   -1.388888888887487817166513568166090e-03,
    -1.984126983386195863616397631859518e-04,   2.480158728990971958442638966674033e-05,
    2.755731488022942163032964623184407e-06,    -2.755731443473966320824186766381266e-07,
    -2.505091491904481811209312506444584e-08,   2.087572839454679266204771377713960e-09,
    1.590445674389578325979670578185607e-10,    -1.135961908867298976458649289359009e-11,
};

//-----------------------------------------------------------------------
//                              FLOAT
//-----------------------------------------------------------------------
const float simd_sincos_table_float_data::poly_sin[] = 
{
    -1.666665997462497498532501453535245e-01f,
    8.332358159410846135491186125998975e-03f,
    -1.953172998286545971259259264194477e-04f,
};

const float simd_sincos_table_float_data::poly_cos[] = 
{                                    
    4.166665997182391716567229450482557e-02f,
    -1.388791273368356551816308637378690e-03f,
    2.449173792393513916022597242670008e-05f,
};

const float simd_sincos_table_float_data::poly_sincos[] =
{
    -1.666665997462497498532501453535245e-01f,  -1.666665997462497498532501453535245e-01f,
    4.166665997182391716567229450482557e-02f,   4.166665997182391716567229450482557e-02f,
    8.332358159410846135491186125998975e-03f,   8.332358159410846135491186125998975e-03f,
    -1.388791273368356551816308637378690e-03f,  -1.388791273368356551816308637378690e-03f,
    -1.953172998286545971259259264194477e-04f,  -1.953172998286545971259259264194477e-04f,
    2.449173792393513916022597242670008e-05f,   2.449173792393513916022597242670008e-05f,
};

const double simd_sincos_table_float_data::poly_sin_double[] = 
{
    -1.666665997462497498532501453535245e-01,
    8.332358159410846135491186125998975e-03,
    -1.953172998286545971259259264194477e-04,
};

const double simd_sincos_table_float_data::poly_cos_double[] = 
{                                    
    4.166665997182391716567229450482557e-02,
    -1.388791273368356551816308637378690e-03,
    2.449173792393513916022597242670008e-05,
};

const double simd_sincos_table_float_data::poly_sincos_double[] =
{
    -1.666665997462497498532501453535245e-01,  4.166665997182391716567229450482557e-02,
    8.332358159410846135491186125998975e-03,   -1.388791273368356551816308637378690e-03,
    -1.953172998286545971259259264194477e-04,  2.449173792393513916022597242670008e-05,
};

}}}
