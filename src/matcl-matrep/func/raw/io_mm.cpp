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

#include "matcl-matrep/func/raw/io.h"
#include <vector>
#include "matcl-core/error/exception_classes.h"
#include "matcl-matrep/func/raw/mvgen.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/integer.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/lib_functions/func_unary.h"

#include <iomanip>
#include "boost/io/ios_state.hpp"

namespace matcl { namespace raw 
{
 
namespace mrd = matcl::raw::details;

//--------------------------------------------------------------------
//                      mm_helper
//--------------------------------------------------------------------
static inline 
char tolower_impl(char c)
{
    return (char)::tolower(c);
};

struct mm_helper
{
    enum class struct_type
    {
        general, sym, her, skew_sym
    };

    static void save_preambule(std::ostream& os)
    {
        os << "%%MatrixMarket" << " ";
    };

    static void load_preambule(std::istream& is)
    {
        if (!is.good())
            throw error::unable_to_read_mm_matrix("invalid stream");

        std::string header;
        
        if (stream_helpers::read(is, header) == false)
            throw error::unable_to_read_mm_matrix("unable to read header");

        std::transform(header.begin(), header.end(), header.begin(), &tolower_impl);

        if (header != "%%matrixmarket")
            throw error::unable_to_read_mm_matrix("not a valid MatrixMarket header");
    };

    static void load_matrix_type(std::istream& is, std::string& object, std::string& storage, 
                                 std::string& value, std::string& struct_type)
    {
        if (stream_helpers::read(is, object) == false)
            throw error::unable_to_read_mm_matrix("unrecognized header format");
        if (stream_helpers::read(is, storage) == false)
            throw error::unable_to_read_mm_matrix("unrecognized header format");
        if (stream_helpers::read(is, value) == false)
            throw error::unable_to_read_mm_matrix("unrecognized header format");
        if (stream_helpers::read(is, struct_type) == false)
            throw error::unable_to_read_mm_matrix("unrecognized header format");
    };

    static void save_matrix_type(std::ostream& os, const std::string& object, const std::string& storage, 
                                 const std::string& value, const std::string& st)
    {
        os << object << " " << storage << " " << value << " " << st << "\n";
    };

    static void load_header_finalize(std::istream& is)
    {
        if (stream_helpers::check_nl(is) == false)
            throw error::unable_to_read_mm_matrix("unrecognized header format");
    };

    static void load_comments(std::istream& is, std::string& comments)
    {
        stream_helpers::load_comments(is, comments);
    };

    static void save_comments(std::ostream& os, const std::string& comments)
    {
        stream_helpers::save_comments(os, comments);
    };

    static void load_dense_matrix_dims(std::istream& is, Integer& M, Integer& N)
    {
        is >> M;
        is >> N;

        stream_helpers::skip_white(is);

        if (stream_helpers::check_nl(is) == false || !is.good())
            throw error::unable_to_read_mm_matrix("unrecognized format of dense matrix size");
    };

    static void load_sparse_matrix_dims(std::istream& is, Integer& M, Integer& N, Integer& NZ)
    {
        is >> M;
        is >> N;
        is >> NZ;

        stream_helpers::skip_white(is);

        if (stream_helpers::check_nl(is) == false || !is.good())
            throw error::unable_to_read_mm_matrix("unrecognized format of sparse matrix size");
    };

    template<class V>
    static void make_sym(Integer& NZ, Integer* ptr_c, Integer* ptr_r, V* ptr_x, struct_type str)
    {
        if (str == struct_type::general)
            return;

        Integer NZ0 = NZ;

        if (str == struct_type::sym)
        {
            for (Integer i = 0; i < NZ0; ++i)
            {
                if (ptr_r[i] != ptr_c[i])
                {
                    ptr_r[NZ]   = ptr_c[i];
                    ptr_c[NZ]   = ptr_r[i];
                    ptr_x[NZ]   = ptr_x[i];
                    ++NZ;
                };
            };
        }
        else if (str == struct_type::her)
        {
            for (Integer i = 0; i < NZ0; ++i)
            {
                if (ptr_r[i] != ptr_c[i])
                {
                    ptr_r[NZ]   = ptr_c[i];
                    ptr_c[NZ]   = ptr_r[i];
                    ptr_x[NZ]   = conj(ptr_x[i]);
                    ++NZ;
                };
            };
        }
        else
        {       
            for (Integer i = 0; i < NZ0; ++i)
            {
                if (ptr_r[i] != ptr_c[i])
                {
                    ptr_r[NZ]   = ptr_c[i];
                    ptr_c[NZ]   = ptr_r[i];
                    ptr_x[NZ]   = -ptr_x[i];
                    ++NZ;
                };
            };
        }
    };

    static void make_sym_pattern(Integer& NZ, Integer* ptr_c, Integer* ptr_r, struct_type str)
    {
        if (str == struct_type::general)
            return;

        Integer NZ0 = NZ;

        for (Integer i = 0; i < NZ0; ++i)
        {
            if (ptr_r[i] != ptr_c[i])
            {
                ptr_r[NZ]   = ptr_c[i];
                ptr_c[NZ]   = ptr_r[i];
                ++NZ;
            };
        };
    };

    template<class V>
    static void load_sparse_data(std::istream& is, Integer M, Integer N, Integer NZ, Integer* ptr_c, 
                                 Integer* ptr_r, V* ptr_x)
    {
        Integer i;
        for(i = 0; i < NZ; ++i)
        {
            Integer I, J;
            V val;

	        if (stream_helpers::read(is, I) == false)
                goto lab_error;
	        if (stream_helpers::read(is, J) == false)
                goto lab_error;
	        if (stream_helpers::read(is, val) == false)
                goto lab_error;

            ptr_r[i]    = I;
            ptr_c[i]    = J;
            ptr_x[i]    = val;

            if (I < 0 || I > M)
                throw error::unable_to_read_mm_matrix("row index out of range");
            if (J < 0 || J > N)
                throw error::unable_to_read_mm_matrix("column index out of range");

            stream_helpers::skip_white(is);

            if (stream_helpers::check_nl(is) == false || !is.good())
                goto lab_error;
        };

        return;

      lab_error:
        std::ostringstream msg;
        msg << "unable to read data for element " << i + 1;
        throw error::unable_to_read_mm_matrix(msg.str());
    };

    static void load_sparse_pattern(std::istream& is, Integer M, Integer N, Integer NZ, Integer* ptr_c, Integer* ptr_r)
    {
        Integer i;
        for(i = 0; i < NZ; ++i)
        {
            Integer I, J;

	        if (stream_helpers::read(is, I) == false)
                goto lab_error;
	        if (stream_helpers::read(is, J) == false)
                goto lab_error;

            ptr_r[i]    = I;
            ptr_c[i]    = J;

            if (I < 0 || I > M)
                throw error::unable_to_read_mm_matrix("row index out of range");
            if (J < 0 || J > N)
                throw error::unable_to_read_mm_matrix("column index out of range");

            stream_helpers::skip_white(is);

            if (stream_helpers::check_nl(is) == false || !is.good())
                goto lab_error;
        };

        return;

      lab_error:
        std::ostringstream msg;
        msg << "unable to read data for element " << i + 1;
        throw error::unable_to_read_mm_matrix(msg.str());
    };

    template<class V>
    static void load_dense_data(std::istream& is, Integer M, Integer N, V* ptr, struct_type str)
    {
        Integer LD  = M;

        Integer i   = 0;
        Integer NZ  = M * N;

        if (str == struct_type::general)
        {
            for(; i < NZ; ++i)
            {
                V val;

	            if (stream_helpers::read(is, val) == false)
                    goto lab_error;

                ptr[i]      = val;

                stream_helpers::skip_white(is);

                if (stream_helpers::check_nl(is) == false || !is.good())
                    goto lab_error;
            };
        }
        else
        {
            Integer off = (str == struct_type::skew_sym) ? 1 : 0;
            V* ptr_0    = ptr;

            //only lower triangle is stored
            for(Integer c = 0; c < N; ++c)
            {
                if (str == struct_type::sym)
                {
                    for (Integer r = 0; r < c; ++r)
                        ptr[r]  = ptr_0[c + r * LD];
                }
                else if (str == struct_type::her)
                {
                    for (Integer r = 0; r < c; ++r)
                        ptr[r]  = matcl::conj(ptr_0[c + r * LD]);
                }
                else
                {
                    //skew-sym
                    for (Integer r = 0; r < c; ++r)
                        ptr[r]  = -ptr_0[c + r * LD];

                    ptr[c]      = V(0);
                }

                for (Integer r = c + off; r < M; ++r, ++i)
                {
                    V val;

	                if (stream_helpers::read(is, val) == false)
                        goto lab_error;

                    ptr[r]      = val;

                    stream_helpers::skip_white(is);

                    if (stream_helpers::check_nl(is) == false || !is.good())
                        goto lab_error;
                }

                ptr += LD;
            };
        }

        return;

      lab_error:
        std::ostringstream msg;
        msg << "unable to read data for element " << i + 1;
        throw error::unable_to_read_mm_matrix(msg.str());
    };    
};

template<class V, class S>
struct save_matrix_data
{};

template<class V>
struct save_matrix_data<V,struct_dense>
{
    using Mat   = raw::Matrix<V,struct_dense>;
    static void eval(std::ostream& os, const Mat& A, bool is_symher)
    {
        boost::io::ios_flags_saver old_flags(os);
        boost::io::ios_precision_saver old_prec(os);

        int precision = mr::get_stream_precision<V>::eval();

        if (precision > 0)
            os << std::setprecision(precision) << std::scientific;

        Integer M       = A.rows();
        Integer N       = A.cols();

        const V* ptr    = A.ptr();
        Integer A_ld    = A.ld();

        os << M << " " << N << "\n";

        if (is_symher == true)
        {
            for (Integer c = 0; c < N; ++c)
            {
                for (Integer r = c; r < M; ++r)
                {
                    matcl::raw::stream_helpers::write(os, ptr[r]);
                    os << "\n";
                }

                ptr         += A_ld;
            };
        }
        else
        {
            for (Integer c = 0; c < N; ++c)
            {
                for (Integer r = 0; r < M; ++r)
                {
                    matcl::raw::stream_helpers::write(os, ptr[r]);
                    os << "\n";
                }

                ptr         += A_ld;
            };
        };
    };
};

template<class V>
struct save_matrix_data<V,struct_sparse>
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    static void eval(std::ostream& os, const Mat& A, bool is_symher)
    {
        boost::io::ios_flags_saver old_flags(os);
        boost::io::ios_precision_saver old_prec(os);

        int precision = mr::get_stream_precision<V>::eval();

        if (precision > 0)
            os << std::setprecision(precision) << std::scientific;

        Integer M       = A.rows();
        Integer N       = A.cols();

        const mrd::sparse_ccs<V>& rep   = A.rep();

        const V* ptr_x          = rep.ptr_x();
        const Integer* ptr_r    = rep.ptr_r();
        const Integer* ptr_c    = rep.ptr_c();        

        if (is_symher == true)
        {
            //calculate NZ
            Integer NZ          = 0;

            if (ptr_c != nullptr)
            {
                for (Integer c = 0; c < N; ++c)
                {
                    Integer J           = c + 1;

                    for (Integer k = ptr_c[c]; k < ptr_c[c+1]; ++k)
                    {
                        Integer I       = ptr_r[k] + 1;                        

                        if (I >= J)
                            ++NZ;
                    }
                };
            };

            os << M << " " << N << " " << NZ << "\n";

            if (ptr_c != nullptr)
            {
                for (Integer c = 0; c < N; ++c)
                {
                    Integer J           = c + 1;

                    for (Integer k = ptr_c[c]; k < ptr_c[c+1]; ++k)
                    {
                        Integer I       = ptr_r[k] + 1;                        
                        const V& val    = ptr_x[k];

                        if (I >= J)
                        {
                            os << I << " " << J << " ";
                            matcl::raw::stream_helpers::write(os, val);
                            os << "\n";
                        };
                    }
                };
            };
        }
        else
        {
            Integer NZ          = rep.nnz();

            os << M << " " << N << " " << NZ << "\n";

            if (ptr_c != nullptr)
            {
                for (Integer c = 0; c < N; ++c)
                {
                    Integer J           = c + 1;

                    for (Integer k = ptr_c[c]; k < ptr_c[c+1]; ++k)
                    {
                        Integer I       = ptr_r[k] + 1;                    
                        const V& val    = ptr_x[k];

                        os << I << " " << J << " ";
                        matcl::raw::stream_helpers::write(os, val);
                        os << "\n";
                    }
                };
            };
        };
    };
};

template<class V>
struct save_matrix_data<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;
    static void eval(std::ostream& os, const Mat& A, bool is_symher)
    {
        boost::io::ios_flags_saver old_flags(os);
        boost::io::ios_precision_saver old_prec(os);

        int precision = mr::get_stream_precision<V>::eval();

        if (precision > 0)
            os << std::setprecision(precision) << std::scientific;

        Integer M       = A.rows();
        Integer N       = A.cols();        
        Integer A_ld    = A.ld();

        const V* ptr    = A.rep_ptr();

        if (is_symher == true)
        {
            Integer NZ  = 0;

            for (Integer c = 0; c < N; ++c)
            {
                Integer fr  = std::max(c,A.first_row(c));
                Integer lr  = A.last_row(c);
                NZ          += std::max(lr - fr + 1, 0);
            };

            os << M << " " << N << " " << NZ << "\n";

            for (Integer c = 0; c < N; ++c)
            {
                Integer fr  = A.first_row(c);
                Integer lr  = A.last_row(c);
                Integer pos = A.first_elem_pos(c);

                for (Integer k = fr; k <= lr; ++k, ++pos)
                {
                    Integer I       = k + 1;
                    Integer J       = c + 1;
                    const V& val    = ptr[pos];

                    if (I >= J)
                    {
                        os << I << " " << J << " ";
                        matcl::raw::stream_helpers::write(os, val);
                        os << "\n";
                    };
                }

                ptr         += A_ld;
            };
        }
        else
        {
            Integer NZ  = A.nnz();
            os << M << " " << N << " " << NZ << "\n";

            for (Integer c = 0; c < N; ++c)
            {
                Integer fr  = A.first_row(c);
                Integer lr  = A.last_row(c);
                Integer pos = A.first_elem_pos(c);

                for (Integer k = fr; k <= lr; ++k, ++pos)
                {
                    Integer I       = k + 1;
                    Integer J       = c + 1;
                    const V& val    = ptr[pos];

                    os << I << " " << J << " ";
                    matcl::raw::stream_helpers::write(os, val);
                    os << "\n";
                }

                ptr         += A_ld;
            };
        };
    };
};

template<class Mat>
void mm_helper_save<Mat>::eval(std::ostream& os, const Mat& A, const std::string& comments)
{
    using S         = typename Mat::struct_type;
    using V         = typename Mat::value_type;

    mm_helper::save_preambule(os);
    
    std::string object, storage, value, struct_type;

    bool is_dense   = std::is_same<S, struct_dense>::value;
    bool is_int     = std::is_same<V, Integer>::value;
    bool is_compl   = md::is_complex<V>::value;
    bool is_real    = (is_compl == false);
    bool is_A_sq    = A.rows() == A.cols();

    bool is_her0    = A.get_struct().is_hermitian(is_A_sq, is_real);
    bool is_sym0    = A.get_struct().is_symmetric(is_A_sq, is_real);
    bool is_gen     = (is_her0 == false) && (is_sym0 == false);
    bool is_sym     = (is_her0 == true && is_compl == false) || is_sym0;
    bool is_symher  = is_gen == false;

    object          = "matrix";
    storage         = is_dense ? "array" : "coordinate";
    value           = is_int? "integer" : (is_compl ? "complex" : "real");
    struct_type     = is_gen? "general" : (is_sym ? "symmetric" : "hermitian");

    mm_helper::save_matrix_type(os, object, storage, value, struct_type);
    mm_helper::save_comments(os, comments);

    save_matrix_data<V,S>::eval(os, A, is_symher);
};

void mm_helper_load::eval(std::istream& is, matcl::Matrix& ret, std::string& comments)
{
    mm_helper::load_preambule(is);

    std::string object, storage, value, struct_type;

    mm_helper::load_matrix_type(is, object, storage, value, struct_type);

    std::transform(object.begin(), object.end(), object.begin(), &tolower_impl);
    std::transform(storage.begin(), storage.end(), storage.begin(), &tolower_impl);
    std::transform(value.begin(), value.end(), value.begin(), &tolower_impl);
    std::transform(struct_type.begin(), struct_type.end(), struct_type.begin(), &tolower_impl);

    // decode header
    if (object != "matrix")
        throw error::unable_to_read_mm_matrix("unrecognized object type: " + object);

    bool is_coordinate  = false;
    bool is_dense       = false;

    if (storage == "coordinate")
        is_coordinate   = true;
    else if (storage == "array")
        is_dense        = true;

    if (is_coordinate == false && is_dense == false)
        throw error::unable_to_read_mm_matrix("unrecognized storage type: " + storage);

    bool is_real        = false;
    bool is_complex     = false;
    bool is_int         = false;
    bool is_pattern     = false;

    if (value == "real")
        is_real         = true;
    else if (value == "complex")
        is_complex      = true;
    else if (value == "integer")
        is_int          = true;
    else if (value == "pattern")
        is_pattern      = true;

    if (is_real == false && is_complex == false && is_int == false && is_pattern == false)
        throw error::unable_to_read_mm_matrix("unrecognized value type: " + value);

    if (is_pattern == true && is_coordinate == false)
        throw error::unable_to_read_mm_matrix("invalid value type: pattern is allowed only for"
                                              " coordinate storage");

    bool is_general     = false;
    bool is_sym         = false;
    bool is_her         = false;
    bool is_skew_sym    = false;

    mm_helper::struct_type str = mm_helper::struct_type::general;

    if (struct_type == "general")
    {
        is_general      = true;
        str             = mm_helper::struct_type::general;
    }
    else if (struct_type == "symmetric")
    {
        is_sym          = true;
        str             = mm_helper::struct_type::sym;
    }
    else if (struct_type == "hermitian")
    {
        is_her          = true;
        str             = mm_helper::struct_type::her;
    }
    else if (struct_type == "skew-symmetric")
    {
        is_skew_sym     = true;
        str             = mm_helper::struct_type::skew_sym;
    }

    if (is_general == false && is_sym == false && is_her == false && is_skew_sym == false)
        throw error::unable_to_read_mm_matrix("unrecognized structure: " + struct_type);

    if (is_her == true && is_complex == false)
        throw error::unable_to_read_mm_matrix("invalid structure: hermitian is allowed only for"
                                              " complex value");
    if (is_skew_sym == true && is_pattern == true)
        throw error::unable_to_read_mm_matrix("invalid structure: skew-symmetric is not allowed"
                                              " for pattern storage");

    mm_helper::load_header_finalize(is);
    mm_helper::load_comments(is, comments);

    Integer M, N;
    Integer NZ = 0;

    if (is_dense == true)
        mm_helper::load_dense_matrix_dims(is, M, N);
    else
        mm_helper::load_sparse_matrix_dims(is, M, N, NZ);    

    bool is_gen = (str == mm_helper::struct_type::general);

    if (is_gen == false && M != N)
        throw error::unable_to_read_mm_matrix("nonsymmetric matrix with symmetric/hermitian/skew-symmetric structure");

    struct_flag sf;
    if (str == mm_helper::struct_type::sym)
        sf  = predefined_struct_type::sym;
    else if (str == mm_helper::struct_type::her)
        sf  = predefined_struct_type::her;

    //no skew-symmetric struct flag in matcl

    if (is_dense)
    {
        if (is_real == true)
        {
            using Mat   = raw::Matrix<Real,struct_dense>;
            Mat mat(ti::ti_empty(), M, N);

            Real* ptr   = mat.ptr();
            mm_helper::load_dense_data<Real>(is, M, N, ptr, str);

            mat.add_struct(sf);
            ret = matcl::Matrix(mat,true);
            return;
        }
        else if (is_complex == true)
        {
            using Mat   = raw::Matrix<Complex,struct_dense>;
            Mat mat(ti::ti_empty(), M, N);

            Complex* ptr   = mat.ptr();
            mm_helper::load_dense_data<Complex>(is, M, N, ptr, str);

            mat.add_struct(sf);
            ret = matcl::Matrix(mat,true);
            return;
        }
        else
        {
            //int
            using Mat   = raw::Matrix<Integer,struct_dense>;
            Mat mat(ti::ti_empty(), M, N);

            Integer* ptr   = mat.ptr();
            mm_helper::load_dense_data<Integer>(is, M, N, ptr, str);

            mat.add_struct(sf);
            ret = matcl::Matrix(mat,true);
            return;
        };
    }
    else if (is_pattern == false)
    {
        using Mat_I     = raw::Matrix<Integer,struct_dense>;

        Integer mult    = is_gen ? 1 : 2;

        Mat_I irows(ti::ti_empty(), NZ*mult, 1);
        Mat_I icols(ti::ti_empty(), NZ*mult, 1);

        Integer* ptr_r  = irows.ptr();
        Integer* ptr_c  = icols.ptr();

        if (is_real == true)
        {        
            using V     = Real;
            using Mat_D = raw::Matrix<V,struct_dense>;
            using Mat_S = raw::Matrix<V,struct_sparse>;

            Mat_D x(ti::ti_empty(), NZ*mult, 1);

            V* ptr_x    = x.ptr();

            mm_helper::load_sparse_data<V>(is, M, N, NZ, ptr_c, ptr_r, ptr_x);
            mm_helper::make_sym(NZ, ptr_c, ptr_r, ptr_x, str); 

            Mat_S mat(ti::ti_empty(), ptr_r, ptr_c, ptr_x, M, N, NZ, NZ);

            mat.add_struct(sf);
            ret = matcl::Matrix(mat,true);
            return;
        }
        else if (is_complex == true)
        {
            using V     = Complex;
            using Mat_D = raw::Matrix<V,struct_dense>;
            using Mat_S = raw::Matrix<V,struct_sparse>;

            Mat_D x(ti::ti_empty(), NZ*mult, 1);

            V* ptr_x    = x.ptr();

            mm_helper::load_sparse_data<V>(is, M, N, NZ, ptr_c, ptr_r, ptr_x);
            mm_helper::make_sym(NZ, ptr_c, ptr_r, ptr_x, str);

            Mat_S mat(ti::ti_empty(), ptr_r, ptr_c, ptr_x, M, N, NZ, NZ);

            mat.add_struct(sf);
            ret = matcl::Matrix(mat,true);
            return;
        }
        else
        {
            //int
            using V     = Integer;
            using Mat_D = raw::Matrix<V,struct_dense>;
            using Mat_S = raw::Matrix<V,struct_sparse>;

            Mat_D x(ti::ti_empty(), NZ*mult, 1);

            V* ptr_x    = x.ptr();

            mm_helper::load_sparse_data<V>(is, M, N, NZ, ptr_c, ptr_r, ptr_x);
            mm_helper::make_sym(NZ, ptr_c, ptr_r, ptr_x, str);

            Mat_S mat(ti::ti_empty(), ptr_r, ptr_c, ptr_x, M, N, NZ, NZ);

            mat.add_struct(sf);
            ret = matcl::Matrix(mat,true);
            return;
        };
    }
    else
    {
        Integer mult    = is_gen ? 1 : 2;

        using Mat_I     = raw::Matrix<Integer,struct_dense>;

        Mat_I irows(ti::ti_empty(), NZ*mult, 1);
        Mat_I icols(ti::ti_empty(), NZ*mult, 1);

        Integer* ptr_r  = irows.ptr();
        Integer* ptr_c  = icols.ptr();

        using V     = Real;
        using Mat_D = raw::Matrix<V,struct_dense>;
        using Mat_S = raw::Matrix<V,struct_sparse>;

        Mat_D x(ti::ti_empty(), NZ*mult, 1);

        V* ptr_x    = x.ptr();

        mm_helper::load_sparse_pattern(is, M, N, NZ, ptr_c, ptr_r);
        mm_helper::make_sym_pattern(NZ, ptr_c, ptr_r, str);

        for (Integer i = 0; i < NZ; ++i)
            ptr_x[i]    = 1.0;

        Mat_S mat(ti::ti_empty(), ptr_r, ptr_c, ptr_x, M, N, NZ, NZ);

        mat.add_struct(sf);
        ret = matcl::Matrix(mat,true);
        return;
    };
};

//--------------------------------------------------------------------
//                      INSTANTIATIONS
//--------------------------------------------------------------------
template struct mm_helper_save<Matrix<Integer,struct_dense>>;
template struct mm_helper_save<Matrix<Integer,struct_sparse>>;
template struct mm_helper_save<Matrix<Integer,struct_banded>>;
template struct mm_helper_save<Matrix<Float,struct_dense>>;
template struct mm_helper_save<Matrix<Float,struct_sparse>>;
template struct mm_helper_save<Matrix<Float,struct_banded>>;
template struct mm_helper_save<Matrix<Real,struct_dense>>;
template struct mm_helper_save<Matrix<Real,struct_sparse>>;
template struct mm_helper_save<Matrix<Real,struct_banded>>;
template struct mm_helper_save<Matrix<Complex,struct_dense>>;
template struct mm_helper_save<Matrix<Complex,struct_sparse>>;
template struct mm_helper_save<Matrix<Complex,struct_banded>>;
template struct mm_helper_save<Matrix<Float_complex,struct_dense>>;
template struct mm_helper_save<Matrix<Float_complex,struct_sparse>>;
template struct mm_helper_save<Matrix<Float_complex,struct_banded>>;

};};

