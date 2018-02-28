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

#include "matcl-matrep/details/fwd_decls.h"

namespace matcl { namespace raw
{
    using matcl::Integer;
    using matcl::Real;
    using matcl::Float;
    using matcl::Complex;
    using matcl::Float_complex;
    using matcl::Object;

};};

#define MACRO_EVAL_SCALAR_TYPES(macro,arg,arg2)     \
macro(Integer,arg,arg2)                             \
macro(Real,arg,arg2)                                \
macro(Float,arg,arg2)                               \
macro(Complex,arg,arg2)                             \
macro(Float_complex,arg,arg2)                       \
macro(Object,arg,arg2)

#define MACRO_EVAL_SCALAR_TYPES_2(macro,arg,arg2)   \
macro(Integer,arg,arg2)                             \
macro(Real,arg,arg2)                                \
macro(Float,arg,arg2)                               \
macro(Complex,arg,arg2)                             \
macro(Float_complex,arg,arg2)                       \
macro(Object,arg,arg2)

#define MACRO_EVAL_GEN_TYPES(macro,arg,arg2)        \
macro(integer_dense,arg,arg2)                       \
macro(float_dense,arg,arg2)                         \
macro(real_dense,arg,arg2)                          \
macro(float_complex_dense,arg,arg2)                 \
macro(complex_dense,arg,arg2)                       \
macro(object_dense,arg,arg2)                        \
macro(integer_sparse,arg,arg2)                      \
macro(float_sparse,arg,arg2)                        \
macro(real_sparse,arg,arg2)                         \
macro(float_complex_sparse,arg,arg2)                \
macro(complex_sparse,arg,arg2)                      \
macro(object_sparse,arg,arg2)                    

#define MACRO_EVAL_DENSE_TYPES(macro,arg,arg2)      \
macro(integer_dense,arg,arg2)                       \
macro(float_dense,arg,arg2)                         \
macro(real_dense,arg,arg2)                          \
macro(float_complex_dense,arg,arg2)                 \
macro(complex_dense,arg,arg2)                       \
macro(object_dense,arg,arg2)

#define MACRO_EVAL_SPARSE_TYPES(macro,arg,arg2)     \
macro(integer_sparse,arg,arg2)                      \
macro(float_sparse,arg,arg2)                        \
macro(real_sparse,arg,arg2)                         \
macro(float_complex_sparse,arg,arg2)                \
macro(complex_sparse,arg,arg2)                      \
macro(object_sparse,arg,arg2)                    

#define MACRO_EVAL_BAND_TYPES(macro,arg,arg2)       \
macro(integer_band,arg,arg2)                        \
macro(float_band,arg,arg2)                          \
macro(real_band,arg,arg2)                           \
macro(float_complex_band,arg,arg2)                  \
macro(complex_band,arg,arg2)                        \
macro(object_band,arg,arg2)

#define MACRO_EVAL_GEN_TYPES_2(macro,arg,arg2)      \
macro(integer_dense,arg,arg2)                       \
macro(float_dense,arg,arg2)                         \
macro(real_dense,arg,arg2)                          \
macro(float_complex_dense,arg,arg2)                 \
macro(complex_dense,arg,arg2)                       \
macro(object_dense,arg,arg2)                        \
macro(integer_sparse,arg,arg2)                      \
macro(float_sparse,arg,arg2)                        \
macro(real_sparse,arg,arg2)                         \
macro(float_complex_sparse,arg,arg2)                \
macro(complex_sparse,arg,arg2)                      \
macro(object_sparse,arg,arg2)                    

#define MACRO_EVAL_VEC_TYPES_2(macro,arg,arg2)      \
macro(integer_dense,arg,arg2)                       \
macro(real_dense,arg,arg2)                          \
macro(float_dense,arg,arg2)                         \
macro(complex_dense,arg,arg2)                       \
macro(float_complex_dense,arg,arg2)                 \
macro(object_dense,arg,arg2)                        \
macro(integer_sparse,arg,arg2)                      \
macro(float_sparse,arg,arg2)                        \
macro(real_sparse,arg,arg2)                         \
macro(complex_sparse,arg,arg2)                      \
macro(float_complex_sparse,arg,arg2)                \
macro(object_sparse,arg,arg2)                    

#define MACRO_EVAL_STR_TYPES(macro,arg,arg2)        \
macro(integer_band,arg,arg2)                        \
macro(float_band,arg,arg2)                          \
macro(real_band,arg,arg2)                           \
macro(float_complex_band,arg,arg2)                  \
macro(complex_band,arg,arg2)                        \
macro(object_band,arg,arg2)


#define MACRO_EVAL_STR_TYPES_2(macro,arg,arg2)      \
macro(integer_band,arg,arg2)                        \
macro(float_band,arg,arg2)                          \
macro(real_band,arg,arg2)                           \
macro(float_complex_band,arg,arg2)                  \
macro(complex_band,arg,arg2)                        \
macro(object_band,arg,arg2)                   

#define MACRO_EVAL_STRUCT_TYPES(macro,arg,arg2)     \
macro(struct_banded,arg,arg2)


#define MACRO_INSTANTIATE_21_IMPL(type_name2,type_name1,cl_name)    \
    template struct cl_name<matcl::raw::type_name1,matcl::raw::type_name2>;

#define MACRO_INSTANTIATE_1_IMPL(type_name,cl_name,arg2)            \
    template struct cl_name<matcl::raw::type_name>;


#define MACRO_INSTANTIATE_2_G_IMPL(type_name,cl_name,arg2)          \
    MACRO_EVAL_GEN_TYPES_2(MACRO_INSTANTIATE_21_IMPL,type_name,cl_name)

#define MACRO_INSTANTIATE_2_S_IMPL(type_name,cl_name,arg2)          \
    MACRO_EVAL_SCALAR_TYPES_2(MACRO_INSTANTIATE_21_IMPL,type_name,cl_name)

#define MACRO_INSTANTIATE_2_ST_IMPL(type_name,cl_name,arg2)         \
    MACRO_EVAL_STR_TYPES_2(MACRO_INSTANTIATE_21_IMPL,type_name,cl_name)

#define MACRO_INSTANTIATE_SCAL_1(cl_name)                           \
MACRO_EVAL_SCALAR_TYPES(MACRO_INSTANTIATE_1_IMPL,cl_name,)                

#define MACRO_INSTANTIATE_G_1(cl_name)                              \
MACRO_EVAL_GEN_TYPES(MACRO_INSTANTIATE_1_IMPL,cl_name,)                

#define MACRO_INSTANTIATE_DENSE_1(cl_name)                          \
MACRO_EVAL_DENSE_TYPES(MACRO_INSTANTIATE_1_IMPL,cl_name,)

#define MACRO_INSTANTIATE_BAND_1(cl_name)                           \
MACRO_EVAL_BAND_TYPES(MACRO_INSTANTIATE_1_IMPL,cl_name,)

#define MACRO_INSTANTIATE_SPARSE_1(cl_name)                         \
MACRO_EVAL_SPARSE_TYPES(MACRO_INSTANTIATE_1_IMPL,cl_name,)

#define MACRO_INSTANTIATE_S_1(cl_name)                              \
MACRO_EVAL_STR_TYPES(MACRO_INSTANTIATE_1_IMPL,cl_name,)

#define MACRO_INSTANTIATE_SS_2_F(cl_name)                           \
MACRO_EVAL_SCALAR_TYPES(MACRO_INSTANTIATE_2_S_IMPL,cl_name,)

#define MACRO_INSTANTIATE_SG_2_F(cl_name)                           \
MACRO_EVAL_SCALAR_TYPES(MACRO_INSTANTIATE_2_G_IMPL,cl_name,)

#define MACRO_INSTANTIATE_SST_2_F(cl_name)                          \
MACRO_EVAL_SCALAR_TYPES(MACRO_INSTANTIATE_2_ST_IMPL,cl_name,)

#define MACRO_INSTANTIATE_GS_2_F(cl_name)                           \
MACRO_EVAL_GEN_TYPES(MACRO_INSTANTIATE_2_S_IMPL,cl_name,)

#define MACRO_INSTANTIATE_GG_2_F(cl_name)                           \
MACRO_EVAL_GEN_TYPES(MACRO_INSTANTIATE_2_G_IMPL,cl_name,)

#define MACRO_INSTANTIATE_GG_2_ALL_F(cl_name)                       \
MACRO_EVAL_GEN_TYPES_ALL(MACRO_INSTANTIATE_2_G_ALL_IMPL,cl_name,)

#define MACRO_INSTANTIATE_GST_2_F(cl_name)                          \
MACRO_EVAL_GEN_TYPES(MACRO_INSTANTIATE_2_ST_IMPL,cl_name,)

#define MACRO_INSTANTIATE_GST_2_ALL_F(cl_name)                      \
MACRO_EVAL_GEN_TYPES_ALL(MACRO_INSTANTIATE_2_ST_ALL_IMPL,cl_name,)

#define MACRO_INSTANTIATE_STG_2_F(cl_name)                          \
MACRO_EVAL_STR_TYPES(MACRO_INSTANTIATE_2_G_IMPL,cl_name,)

#define MACRO_INSTANTIATE_STG_2_ALL_F(cl_name)                      \
MACRO_EVAL_STR_TYPES_ALL(MACRO_INSTANTIATE_2_G_ALL_IMPL,cl_name,)

#define MACRO_INSTANTIATE_STS_2_F(cl_name)                          \
MACRO_EVAL_STR_TYPES(MACRO_INSTANTIATE_2_S_IMPL,cl_name,)

#define MACRO_INSTANTIATE_STST_2_F(cl_name)                         \
MACRO_EVAL_STR_TYPES(MACRO_INSTANTIATE_2_ST_IMPL,cl_name,)

#define MACRO_INSTANTIATE_SS_2(cl_name)                             \
MACRO_INSTANTIATE_21_IMPL(Integer,Integer,cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(Real,Real,cl_name)                        \
MACRO_INSTANTIATE_21_IMPL(Float,Float,cl_name)                      \
MACRO_INSTANTIATE_21_IMPL(Complex,Real,cl_name)                     \
MACRO_INSTANTIATE_21_IMPL(Real,Complex,cl_name)                     \
MACRO_INSTANTIATE_21_IMPL(Float_complex,Float,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(Float,Float_complex,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(Float_complex,Float_complex,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(Complex,Complex,cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(Object,Object,cl_name)                    \

#define MACRO_INSTANTIATE_GS_2(cl_name)                             \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_dense,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_sparse,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(Real,real_dense,cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(Real,real_sparse,cl_name)                 \
MACRO_INSTANTIATE_21_IMPL(Real,complex_dense,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(Real,complex_sparse,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(Float,float_dense,cl_name)                \
MACRO_INSTANTIATE_21_IMPL(Float,float_sparse,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(Float,float_complex_dense,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(Float,float_complex_sparse,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(Complex,real_dense,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(Complex,real_sparse,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_dense,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_sparse,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_dense,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_sparse,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_dense,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Object,object_dense,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(Object,object_sparse,cl_name)

#define MACRO_INSTANTIATE_GS_DIAG_2(cl_name)                        \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_dense,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_sparse,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(Real,real_dense,cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(Real,real_sparse,cl_name)                 \
MACRO_INSTANTIATE_21_IMPL(Float,float_dense,cl_name)                \
MACRO_INSTANTIATE_21_IMPL(Float,float_sparse,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_dense,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_sparse,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_dense,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Object,object_dense,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(Object,object_sparse,cl_name)

#define MACRO_INSTANTIATE_GS_DIAG_DENSE_2(cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_dense,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Real,real_dense,cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(Float,float_dense,cl_name)                \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_dense,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_dense,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Object,object_dense,cl_name)              \

#define MACRO_INSTANTIATE_STS_2(cl_name)                        \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_band,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(Real,real_band,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(Real,complex_band,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Float,float_band,cl_name)             \
MACRO_INSTANTIATE_21_IMPL(Float,float_complex_band,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(Complex,real_band,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_band,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_band,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_band,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Object,object_band,cl_name)           \

#define MACRO_INSTANTIATE_STS_DIAG_2(cl_name)                   \
MACRO_INSTANTIATE_21_IMPL(Integer,integer_band,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(Real,real_band,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(Float,float_band,cl_name)             \
MACRO_INSTANTIATE_21_IMPL(Complex,complex_band,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(Float_complex,float_complex_band,cl_name)\
MACRO_INSTANTIATE_21_IMPL(Object,object_band,cl_name)           \

#define MACRO_INSTANTIATE_SG_2(cl_name)                         \
MACRO_INSTANTIATE_21_IMPL(integer_dense,Integer,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,Integer,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(real_dense,Real,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(real_sparse,Real,cl_name)             \
MACRO_INSTANTIATE_21_IMPL(float_dense,Float,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(float_sparse,Float,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(complex_dense,Real,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,Real,cl_name)          \
MACRO_INSTANTIATE_21_IMPL(real_dense,Complex,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(real_sparse,Complex,cl_name)          \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,Float,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,Float,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_dense,Float_complex,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(float_sparse,Float_complex,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(complex_dense,Complex,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,Complex,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_dense,Object,cl_name)          \
MACRO_INSTANTIATE_21_IMPL(object_sparse,Object,cl_name)         \

#define MACRO_INSTANTIATE_SG_DIAG_2(cl_name)                    \
MACRO_INSTANTIATE_21_IMPL(integer_dense,Integer,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,Integer,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(real_dense,Real,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(real_sparse,Real,cl_name)             \
MACRO_INSTANTIATE_21_IMPL(float_dense,Float,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(float_sparse,Float,cl_name)           \
MACRO_INSTANTIATE_21_IMPL(complex_dense,Complex,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,Complex,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_dense,Object,cl_name)          \
MACRO_INSTANTIATE_21_IMPL(object_sparse,Object,cl_name)         \

#define MACRO_INSTANTIATE_SG_DIAG_DENSE_2(cl_name)              \
MACRO_INSTANTIATE_21_IMPL(integer_dense,Integer,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(real_dense,Real,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(float_dense,Float,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(complex_dense,Complex,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_dense,Object,cl_name)          \

#define MACRO_INSTANTIATE_GG_2(cl_name)                         \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,integer_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,integer_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(real_dense,real_dense,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(real_dense,real_sparse,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(real_sparse,real_dense,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(real_sparse,real_sparse,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_sparse,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_dense,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_sparse,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(complex_dense,real_dense,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(complex_dense,real_sparse,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,real_dense,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,real_sparse,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(real_dense,complex_dense,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(real_dense,complex_sparse,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(real_sparse,complex_dense,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(real_sparse,complex_sparse,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_complex_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_complex_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,complex_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_complex_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_dense,object_dense,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(object_dense,object_sparse,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(object_sparse,object_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(object_sparse,object_sparse,cl_name)  \

#define MACRO_INSTANTIATE_GG_DIAG_2(cl_name)                    \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,integer_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,integer_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(real_dense,real_dense,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(real_dense,real_sparse,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(real_sparse,real_dense,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(real_sparse,real_sparse,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_sparse,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_dense,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_sparse,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,complex_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_sparse,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_complex_dense,cl_name) \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_complex_sparse,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_dense,object_dense,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(object_dense,object_sparse,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(object_sparse,object_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(object_sparse,object_sparse,cl_name)  \

#define MACRO_INSTANTIATE_STG_2(cl_name)                        \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_band,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,integer_band,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(real_dense,real_band,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(real_sparse,real_band,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_band,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_band,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(complex_dense,real_band,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,real_band,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(real_dense,complex_band,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(real_sparse,complex_band,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_band,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_band,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_complex_band,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_complex_band,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_band,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,complex_band,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_band,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_complex_band,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(object_dense,object_band,cl_name)                 \
MACRO_INSTANTIATE_21_IMPL(object_sparse,object_band,cl_name)

#define MACRO_INSTANTIATE_STG_DIAG_2(cl_name)                   \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_band,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(integer_sparse,integer_band,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(real_dense,real_band,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(real_sparse,real_band,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_band,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_sparse,float_band,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_band,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(complex_sparse,complex_band,cl_name)              \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_band,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_complex_sparse,float_complex_band,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(object_dense,object_band,cl_name)                 \
MACRO_INSTANTIATE_21_IMPL(object_sparse,object_band,cl_name)

#define MACRO_INSTANTIATE_SST_2(cl_name)                        \
MACRO_INSTANTIATE_21_IMPL(integer_band,Integer,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(real_band,Real,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(float_band,Float,cl_name)             \
MACRO_INSTANTIATE_21_IMPL(complex_band,Real,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(real_band,Complex,cl_name)            \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,Float,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_band,Float_complex,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(complex_band,Complex,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_band,Object,cl_name)           \

#define MACRO_INSTANTIATE_SST_DIAG_2(cl_name)                   \
MACRO_INSTANTIATE_21_IMPL(integer_band,Integer,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(real_band,Real,cl_name)               \
MACRO_INSTANTIATE_21_IMPL(float_band,Float,cl_name)             \
MACRO_INSTANTIATE_21_IMPL(complex_band,Complex,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,Float_complex,cl_name)\
MACRO_INSTANTIATE_21_IMPL(object_band,Object,cl_name)           \

#define MACRO_INSTANTIATE_GST_2(cl_name)                        \
MACRO_INSTANTIATE_21_IMPL(integer_band,integer_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(integer_band,integer_sparse,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(real_band,real_dense,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(real_band,real_sparse,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(float_band,float_dense,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_band,float_sparse,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(real_band,complex_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(real_band,complex_sparse,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_band,float_complex_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_band,float_complex_sparse,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(complex_band,real_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(complex_band,real_sparse,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_sparse,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(complex_band,complex_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(complex_band,complex_sparse,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_complex_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_complex_sparse,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(object_band,object_dense,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(object_band,object_sparse,cl_name)    \

#define MACRO_INSTANTIATE_GST_DIAG_2(cl_name)                   \
MACRO_INSTANTIATE_21_IMPL(integer_band,integer_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(integer_band,integer_sparse,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(real_band,real_dense,cl_name)         \
MACRO_INSTANTIATE_21_IMPL(real_band,real_sparse,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(float_band,float_dense,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_band,float_sparse,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(complex_band,complex_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(complex_band,complex_sparse,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_complex_dense,cl_name)   \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_complex_sparse,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(object_band,object_dense,cl_name)     \
MACRO_INSTANTIATE_21_IMPL(object_band,object_sparse,cl_name)    \

#define MACRO_INSTANTIATE_STST_2(cl_name)                       \
MACRO_INSTANTIATE_21_IMPL(integer_band,integer_band,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(real_band,real_band,cl_name)          \
MACRO_INSTANTIATE_21_IMPL(float_band,float_band,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(real_band,complex_band,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(complex_band,real_band,cl_name)       \
MACRO_INSTANTIATE_21_IMPL(float_band,float_complex_band,cl_name)\
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_band,cl_name)\
MACRO_INSTANTIATE_21_IMPL(complex_band,complex_band,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_complex_band,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(object_band,object_band,cl_name)      \

#define MACRO_INSTANTIATE_STST_DIAG_2(cl_name)                  \
MACRO_INSTANTIATE_21_IMPL(integer_band,integer_band,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(real_band,real_band,cl_name)          \
MACRO_INSTANTIATE_21_IMPL(float_band,float_band,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(complex_band,complex_band,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(float_complex_band,float_complex_band,cl_name)    \
MACRO_INSTANTIATE_21_IMPL(object_band,object_band,cl_name)      \

#define MACRO_INSTANTIATE_BIN_ALL(cl_name)  \
MACRO_INSTANTIATE_GG_2(cl_name)             \
MACRO_INSTANTIATE_GST_2(cl_name)            \
MACRO_INSTANTIATE_STG_2(cl_name)            \
MACRO_INSTANTIATE_STST_2(cl_name)

#define MACRO_INSTANTIATE_BIN_DIAG(cl_name) \
MACRO_INSTANTIATE_GG_DIAG_2(cl_name)        \
MACRO_INSTANTIATE_GST_DIAG_2(cl_name)       \
MACRO_INSTANTIATE_STG_DIAG_2(cl_name)       \
MACRO_INSTANTIATE_STST_DIAG_2(cl_name)

#define MACRO_INSTANTIATE_BIN_DIAG_DENSE(cl_name) \
MACRO_INSTANTIATE_21_IMPL(integer_dense,integer_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(real_dense,real_dense,cl_name)        \
MACRO_INSTANTIATE_21_IMPL(float_dense,float_dense,cl_name)      \
MACRO_INSTANTIATE_21_IMPL(complex_dense,complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(float_complex_dense,float_complex_dense,cl_name)  \
MACRO_INSTANTIATE_21_IMPL(object_dense,object_dense,cl_name)    \
