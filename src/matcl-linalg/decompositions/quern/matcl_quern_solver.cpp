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

#include "matcl-linalg/decompositions/quern/matcl_quern_solver.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/linalg_exception.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-linalg/decompositions/quern/include/quern.h"
#include "matcl-linalg/decompositions/quern/quern_unitary_mat_mult.h"
#include "matcl-internals/error/error_check_basic.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/utils/linalg_utils.h"

#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/norms_error/norm.h"

namespace matcl { namespace details
{

//------------------------------------------------------------------------
//                          quern_matrix
//------------------------------------------------------------------------
template<class Val>
class quern_matrix
{
    public:
        int*    row_start;
        int*    column_index;
        Val*    value;
        Integer cols;

    public:
        quern_matrix()
            : row_start(nullptr), column_index(nullptr), value(nullptr), cols(0)
        {}

        quern_matrix(int * row_start, int * column_index, Val * value, Integer c)
            : row_start(row_start), column_index(column_index), value(value), cols(c)
        {}

        quern_matrix(quern_matrix&& other)
            :row_start(other.row_start), column_index(other.column_index), value(other.value)
            ,cols(other.cols)
        {
            other.row_start     = nullptr;
            other.column_index  = nullptr;
            other.value         = nullptr;
        };

        ~quern_matrix()
        {
            std::free(row_start);
            std::free(column_index);
            std::free(value);
        }

        void    save_matrix(std::ostream& os) const;
        void    save_matrix(oarchive& os) const;
        void    load_matrix(iarchive& is);
        void    load_matrix(std::istream& is);

        bool    all_finite() const;

        template<class T>
        quern_matrix<T> convert_impl() const;

    public:
        void    allocate_rows(Integer m);
        void    allocate_cols_values(Integer nz, Integer nz_mult);

        quern_matrix(const quern_matrix&) = delete;
        quern_matrix& operator=(const quern_matrix&) = delete;
};

template<class Val>
void quern_matrix<Val>::allocate_rows(Integer m)
{
    int* Q_row_start = (int*)std::malloc((m+1)*sizeof(int));

    if (!Q_row_start)
        throw error::alloc((m+1)*sizeof(int));

    row_start = Q_row_start;
};

template<class Val>
void quern_matrix<Val>::allocate_cols_values(Integer nz, Integer nz_mult)
{
    int* Q_column_index = (int*)std::malloc(nz*sizeof(int));

    if (!Q_column_index)
        throw error::alloc(nz*sizeof(int));

    column_index    = Q_column_index;

    Val* Q_value    = (Val*)std::malloc(nz*sizeof(Val)*nz_mult);

    if (!Q_value)
        throw error::alloc(nz*sizeof(Val)*nz_mult);

    value           = Q_value;
};

template<class Val>
void quern_matrix<Val>::save_matrix(oarchive& os) const
{
    os << cols;

    //row indices
    for (Integer i = 0; i <= cols; ++i)
        os << row_start[i];

    Integer nz = row_start[cols];

    //column indices
    for (Integer i = 0; i < nz; ++i)
        os << column_index[i];

    Integer nz2 = nz * 2;

    //values
    for (Integer i = 0; i < nz2; ++i)
        os << value[i];
};

template<class Val>
void quern_matrix<Val>::load_matrix(iarchive& is)
{
    is >> cols;

    allocate_rows(cols);

    //row indices
    for (Integer i = 0; i <= cols; ++i)
        is >> row_start[i];

    Integer nz = row_start[cols];

    allocate_cols_values(nz, 2);

    //column indices
    for (Integer i = 0; i < nz; ++i)
        is >> column_index[i];

    Integer nz2 = nz * 2;

    //values
    for (Integer i = 0; i < nz2; ++i)
        is >> value[i];
};

template<class Val>
template<class T>
quern_matrix<T> quern_matrix<Val>::convert_impl() const
{
    quern_matrix<T> ret;

    Integer nz = row_start[cols];

    ret.cols = cols;
    ret.allocate_rows(cols);
    ret.allocate_cols_values(nz, 2);

    for (Integer i = 0; i <= cols; ++i)
        ret.row_start[i]    = this->row_start[i];

    for (Integer i = 0; i < nz; ++i)
        ret.column_index[i] = this->column_index[i];

    Integer nz2 = nz * 2;

    for (Integer i = 0; i < nz2; ++i)
        ret.value[i]    = raw::converter<T,Val>::eval(this->value[i]);

    return ret;
};

template<class Val>
bool quern_matrix<Val>::all_finite() const
{
    Integer nz  = row_start[cols];
    Integer nz2 = nz * 2;

    for (Integer i = 0; i < nz2; ++i)
        if (mrd::isfinite_helper<Val>::eval(this->value[i]) == false)
            return false;

    return true;
};

template<class Val>
void quern_matrix<Val>::save_matrix(std::ostream& os) const
{
    os << " ";

    os << cols << " ";

    //row indices
    os << "rows" << " ";
    for (Integer i = 0; i <= cols; ++i)
        os << row_start[i] << " ";

    Integer nz = row_start[cols];

    //column indices
    os << "columns" << " ";
    for (Integer i = 0; i < nz; ++i)
        os << column_index[i] << " ";

    Integer nz2 = nz * 2;

    //values
    os << "values" << " ";
    for (Integer i = 0; i < nz2; ++i)
        os << value[i] << " ";
};

template<class Val>
void quern_matrix<Val>::load_matrix(std::istream& is)
{
    is >> cols;

    allocate_rows(cols);

    std::string sep;

    //row indices
    is >> sep;
    for (Integer i = 0; i <= cols; ++i)
        is >> row_start[i];

    Integer nz = row_start[cols];

    allocate_cols_values(nz, 2);

    //column indices
    is >> sep;
    for (Integer i = 0; i < nz; ++i)
        is >> column_index[i];

    Integer nz2 = nz * 2;

    //values
    is >> sep;
    for (Integer i = 0; i < nz2; ++i)
        is >> value[i];
};

//------------------------------------------------------------------------
//                          quern_Q_matrix
//------------------------------------------------------------------------

template<class Val>
class quern_Q_matrix : public quern_matrix<Val>, public unitary_matrix_data
{
    private:
        using base_type = quern_matrix<Val>;

    private:
        Integer     m_rows;
        Integer     m_cols;

    public:
        quern_Q_matrix(base_type&& base, Integer rows, Integer cols)
            :base_type(std::move(base)), m_rows(rows), m_cols(cols)
        {};

    public:
        quern_Q_matrix()
            : quern_matrix<Val>(nullptr, nullptr, nullptr, 0), m_rows(0), m_cols(0)
        {}

        quern_Q_matrix(int * row_start, int * column_index, Val * value, Integer rows, Integer cols)
            : quern_matrix<Val>(row_start, column_index, value, rows), m_rows(rows), m_cols(cols)
        {}

        virtual Integer         rows() const override           { return m_rows; }
        virtual Integer         cols() const override           { return m_cols; }
        virtual value_code      get_value_code() const override { return matrix_traits::value_code<Val>::value; }
        virtual ti::ti_object   get_type() const override         { return ti::ti_object_type<Val>(); };
        virtual bool            all_finite() const override     { return true; };

        virtual void            to_matrix(matcl::Matrix& ret) const override;
        virtual data_ptr        convert(value_code new_val_code) const override;
        virtual void            mult_right(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                        trans_type t_unitary) const override;
        virtual void            mult_left(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                        trans_type t_unitary) const override;

        virtual serialization_helper<unitary_matrix_data>*
                                get_serialization_helper() const override;
        virtual void            save(oarchive& os) const override;
        virtual void            save(std::ostream& os) const override;

        static unitary_matrix_data* load(std::istream& is);
        static unitary_matrix_data* load(iarchive& ar);

        template<class T>
        data_ptr                convert_impl() const;
};

template<class Val>
serialization_helper<unitary_matrix_data>*
quern_Q_matrix<Val>::get_serialization_helper() const
{
    std::ostringstream msg;
    msg << "quern_Q_matrix" << "<" << get_type_name<Val>::eval() << ">";

    return serialization_helper<unitary_matrix_data>
            ::get<quern_Q_matrix<Val>>(msg.str());
};

template<class Val>
void quern_Q_matrix<Val>::mult_right(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                     trans_type t_unitary) const
{
    if (mat.is_scalar())
    {
        this->to_matrix(ret);
        ret = ret * mat;
        return;
    };

    error::check_mul(m_rows, m_cols, mat.rows(), mat.cols(), t_unitary, trans_type::no_trans);

    if (mat.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(mat.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        Integer M;

        if (t_unitary == trans_type::no_trans)
            M   = this->m_rows;
        else
            M   = this->m_cols;

        ret = spzeros(M, mat.cols(), 0, vt);
        return;
    };

    quern_unitary_mat_mult<Val> qmult(this->m_rows, this->m_cols, this->row_start, this->column_index,
                                      this->value);

    return qmult.mult_right(ret, mat, t_unitary);
};
template<class Val>
void quern_Q_matrix<Val>::mult_left(matcl::Matrix& ret, const matcl::Matrix& mat, 
                                     trans_type t_unitary) const
{
    if (mat.is_scalar())
    {
        this->to_matrix(ret);
        ret = ret * mat;
        return;
    };

    error::check_mul(mat.rows(), mat.cols(), m_rows, m_cols, trans_type::no_trans, t_unitary);

    if (mat.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(mat.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        Integer N;

        if (t_unitary == trans_type::no_trans)
            N   = this->m_cols;
        else
            N   = this->m_rows;

        ret = spzeros(mat.rows(), N, 0, vt);
        return;
    };

    quern_unitary_mat_mult<Val> qmult(this->m_rows, this->m_cols, this->row_start, this->column_index,
                                      this->value);

    return qmult.mult_left(ret, mat, t_unitary);
};

template<class Val>
void quern_Q_matrix<Val>::to_matrix(matcl::Matrix& ret) const
{
    Matrix I = speye(m_cols,m_cols,this->get_value_code());
    return mult_right(ret, I, trans_type::no_trans);
};

template<class Val>
struct quern_Q_convert_vis : public extract_type_switch<void, quern_Q_convert_vis<Val>,true>
{
    using umatrix   = quern_Q_matrix<Val>;
    using data_ptr  = typename unitary_matrix_data::data_ptr;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        using VM    = typename T::value_type;

        ret = data.convert_impl<VM>();
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        ret = data.convert_impl<T>();
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
};

template<class Val>
typename quern_Q_matrix<Val>::data_ptr 
quern_Q_matrix<Val>::convert(value_code new_val_code) const
{
    value_code vc   = matrix_traits::unify_value_types(new_val_code, value_code::v_float);
    Matrix v        = zeros(0,0,vc);

    data_ptr ret;
    quern_Q_convert_vis<Val>::make<const Matrix&>(v,*this, ret);
    return ret;
};

template<class Val>
template<class T>
typename quern_Q_matrix<Val>::data_ptr
quern_Q_matrix<Val>::convert_impl() const
{
    using TR        = typename details::unify_types<T, Float>::type;

    return data_ptr(new quern_Q_matrix<TR>(base_type::convert_impl<TR>(), m_rows, m_cols));
};

template<class Val>
void quern_Q_matrix<Val>::save(oarchive& os) const
{
    //general parameters
    os << this->m_rows;
    os << this->m_cols;

    quern_matrix<Val>::save_matrix(os);
};

template<class Val>
unitary_matrix_data* quern_Q_matrix<Val>::load(std::istream& is)
{
    using ptr_type = std::unique_ptr<quern_Q_matrix<Val>>;
    ptr_type ret = ptr_type(new quern_Q_matrix<Val>());

    //general parameters
    is >> ret->m_rows;
    is >> ret->m_cols;

    ret->quern_matrix<Val>::load_matrix(is);

    return ret.release();
};

template<class Val>
void quern_Q_matrix<Val>::save(std::ostream& os) const
{
    os << " ";

    //general parameters
    os << this->m_rows << " ";
    os << this->m_cols << " ";
    
    quern_matrix<Val>::save_matrix(os);
};

template<class Val>
unitary_matrix_data* quern_Q_matrix<Val>::load(iarchive& ar)
{
    using ptr_type  = std::unique_ptr<quern_Q_matrix<Val>>;
    ptr_type ret    = ptr_type(new quern_Q_matrix<Val>());

    //general parameters
    ar >> ret->m_rows;
    ar >> ret->m_cols;

    ret->quern_matrix<Val>::load_matrix(ar);

    return ret.release();
};

//------------------------------------------------------------------------
//                          quern_solver_impl
//------------------------------------------------------------------------
class quern_solver_impl
{
    public:
        using umatrix = unitary_matrix;

    public:
        Integer     nnz_A;
        Integer     m_rows_At;
        Integer     m_cols_At;
        bool        m_all_results_nan;
        bool        from_empty_matrix;

    public:
        quern_solver_impl()
            :nnz_A(0), m_all_results_nan(false), from_empty_matrix(false), 
            m_rows_At(0), m_cols_At(0)
        {};

        virtual ~quern_solver_impl(){};

        virtual Matrix  get_r_trans() const = 0;
        virtual umatrix get_unitary_matrix() const = 0;
};

template<class Val>
class quern_solver_impl_val : public quern_solver_impl
{
    private:
        using Mat           = raw::Matrix<Val,struct_sparse>;
        using umatrix       = unitary_matrix;
        using quern_impl    = std::shared_ptr<quern_Q_matrix<Val>>;

    private:
        umatrix         m_Q;
        matcl::Matrix   m_RT;
        bool            m_economy;

    public:
        quern_solver_impl_val(const Matrix& A, bool with_q, bool economy);

        ~quern_solver_impl_val()
        {}

        virtual Matrix  get_r_trans() const override;
        virtual umatrix get_unitary_matrix() const override;

    private:
        void            make_r_trans(const quern_matrix<Val>& r, matcl::Matrix& ret) const;
        void            make_q(const quern_impl& Q_impl, unitary_matrix& Q) const;
};

template<class Val>
quern_solver_impl_val<Val>::quern_solver_impl_val(const Matrix& As, bool with_q, bool economy)
    :m_economy(economy)
{
    // NOTE the cols/rows swap
    m_cols_At   = As.cols();
    m_rows_At   = As.rows();
    nnz_A       = As.structural_nnz();
    Integer K   = (economy == false)? m_cols_At : std::min(m_rows_At,m_cols_At);

    if (m_cols_At == 0 || m_rows_At == 0) 
        from_empty_matrix = true;    
            
    bool isv    = As.all_finite();

    if (isv == false) 
        m_all_results_nan = true;

    if(nnz_A == 0 || m_all_results_nan || from_empty_matrix) 
    {
        details::quern_matrix<Val> R_impl;
        make_r_trans(R_impl, m_RT);

        quern_impl Q_impl;
        make_q(Q_impl, m_Q);

        return;
    }

    // NOTE we are effectively transposing A by treating CSC as CSR
    Mat As_impl = As.get_impl<Mat>();    

    // will be allocated by quern
    int* ptr_Q_row_start    = nullptr;
    int* ptr_Q_column_index = nullptr;
    int* ptr_R_row_start    = nullptr;
    int* ptr_R_column_index = nullptr;
    Val* ptr_Q_value        = nullptr;
    Val* ptr_R_value        = nullptr;

    // this routine prevents leaks in case of trouble
    int quern_return_code;

    if (with_q == true)
    {
        quern_return_code 
            = quern::QUERN_compute_qr(m_cols_At ,m_rows_At, As_impl.rep().ptr_c(), // column start as row start
                                 As_impl.rep().ptr_r(), // row index as colummn index
                                 As_impl.rep().ptr_x(),
                                 nullptr, // row order
                                 &ptr_Q_row_start, &ptr_Q_column_index, &ptr_Q_value,
                                 &ptr_R_row_start, &ptr_R_column_index, &ptr_R_value);
    }
    else
    {
        quern_return_code 
            = quern::QUERN_compute_qr_without_q(m_cols_At ,m_rows_At, As_impl.rep().ptr_c(), // column start as row start
                                 As_impl.rep().ptr_r(), // row index as colummn index
                                 As_impl.rep().ptr_x(),
                                 nullptr, // row order
                                 &ptr_R_row_start, &ptr_R_column_index, &ptr_R_value);
    };

    if (quern_return_code != QUERN_OK) 
    {
        if (quern_return_code == QUERN_OUT_OF_MEMORY)
            throw error::alloc();
        else
            throw error::error_qrs(quern_return_code);
    };

    if (with_q == true)
    {
        quern_impl Q_impl(new details::quern_Q_matrix<Val>(ptr_Q_row_start, ptr_Q_column_index, 
                        ptr_Q_value, m_cols_At, K));
        make_q(Q_impl, m_Q);
    };

    details::quern_matrix<Val> R_impl(ptr_R_row_start, ptr_R_column_index, ptr_R_value, m_rows_At);
    make_r_trans(R_impl, m_RT);

    return;
};

template<class Val>
void quern_solver_impl_val<Val>::make_r_trans(const details::quern_matrix<Val>& m_R, matcl::Matrix& ret) const
{
    Integer K0          = std::min(m_rows_At, m_cols_At);
    Integer K           = (m_economy == true)? K0 : m_cols_At;

    if (from_empty_matrix || nnz_A == 0)
    {
        using VR        = typename md::real_type<Val>::type;
        value_code v    = matrix_traits::value_code<VR>::value;
       
        ret             = spzeros(m_rows_At,K,0,v);
        return;
    };

    // NOTE we are effectively transposing R by treating CSC as CSR
    if (m_all_results_nan) 
    {
        ret = details::make_nan_matrix<Val>(m_rows_At, K);
        return;
    };

    Integer nnz = m_R.row_start[K0];

    // NOTE the cols/rows swap
    Mat ret_impl(ti::ti_type<Val>(),m_rows_At,K,nnz);

    Integer * ptr_r = ret_impl.rep().ptr_r();
    Integer * ptr_c = ret_impl.rep().ptr_c();
    Val * ptr_x     = ret_impl.rep().ptr_x();

    for (int i = 0; i < K0; i++)
        ptr_c[i] = m_R.row_start[i];

    // complete the part of column-starts' which quern does not keep
    for (int i = K0; i < K; i++)
        ptr_c[i] = nnz;

    for (int i = 0; i < nnz; i++)
    {
        ptr_r[i] = m_R.column_index[i];
        ptr_x[i] = m_R.value[i];
    }

    ptr_c[K] = nnz;
    ret_impl.rep().sort();
    ret_impl.add_struct(predefined_struct_type::tril);
    
    ret = Matrix(ret_impl,true);
    return;
}

template<class Val>
Matrix quern_solver_impl_val<Val>::get_r_trans() const
{
    return m_RT;
}

template<class Val>
void quern_solver_impl_val<Val>::make_q(const quern_impl& Q_impl, unitary_matrix& Q) const
{
    Integer K0          = std::min(m_rows_At, m_cols_At);
    Integer K           = (m_economy == true)? K0 : m_cols_At;

    if (from_empty_matrix || nnz_A == 0)
    {
        using VR        = typename md::real_type<Val>::type;
        value_code v    = matrix_traits::value_code<VR>::value;        
        Matrix I        = speye(m_cols_At,K,v);
        Q               = unitary_matrix(I,false);
        return;
    };

    if (m_all_results_nan) 
    {
        Q               = unitary_matrix::from_nan(m_cols_At, K, matrix_traits::value_code<Val>::value);
        return;
    };

    Q = unitary_matrix(Q_impl);
    return;
}

template<class Val>
unitary_matrix quern_solver_impl_val<Val>::get_unitary_matrix() const
{
    return m_Q;
};

}};

//------------------------------------------------------------------------
//                          quern_solver
//------------------------------------------------------------------------
namespace matcl
{

quern_solver::quern_solver(const Matrix& A0, bool with_q, bool economy, bool trans_in)
    :m_economy(economy)
{
    matcl::value_code v  = A0.get_value_code();

    switch (v)
    {
        case value_code::v_integer:
        {
            Matrix A = convert(A0, mat_code::real_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Real>(A, with_q, economy));
            return;
        }
        case value_code::v_float:
        {
            Matrix A(A0);
            A = convert(A, mat_code::float_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Float>(A, with_q, economy));
            return;
        }
        case value_code::v_real:
        {
            Matrix A(A0);
            A = convert(A, mat_code::real_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Real>(A, with_q, economy));
            return;
        }
        case value_code::v_float_complex:
        {
            Matrix A(A0);
            A = convert(A, mat_code::float_complex_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Float_complex>(A, with_q, economy));
            return;
        }
        case value_code::v_complex:
        {
            Matrix A(A0);
            A = convert(A, mat_code::complex_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Complex>(A, with_q, economy));
            return;
        }
        case value_code::v_object:
            throw error::object_value_type_not_allowed("qr");
        default:
            throw error::error_general("impossible type case in qrs");
    }
}
quern_solver::quern_solver(Matrix&& A0, bool with_q, bool economy, bool trans_in)
    :m_economy(economy)
{
    matcl::value_code v  = A0.get_value_code();

    switch (v)
    {
        case value_code::v_integer:
        {
            Matrix As = convert(std::move(A0), mat_code::real_sparse);
            As = trans_in ? As : trans(As);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Real>(As,with_q,economy));
            return;
            return;
        }
        case value_code::v_float:
        {
            Matrix A(std::move(A0));
            A = convert(A, mat_code::float_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Float>(A,with_q,economy));
            return;
        }
        case value_code::v_real:
        {
            Matrix A(std::move(A0));
            A = convert(A, mat_code::real_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Real>(A,with_q,economy));
            return;
        }
        case value_code::v_float_complex:
        {
            Matrix A(std::move(A0));
            A = convert(A, mat_code::float_complex_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Float_complex>(A,with_q,economy));
            return;
        }
        case value_code::v_complex:
        {
            Matrix A(std::move(A0));
            A = convert(A, mat_code::complex_sparse);
            A = trans_in ? A : trans(A);
            m_impl = impl_ptr(new details::quern_solver_impl_val<Complex>(A,with_q,economy));
            return;
        }
        case value_code::v_object:
            throw error::object_value_type_not_allowed("qr");
        default:
            throw error::error_general("impossible type case in qrs");
    }
}

Matrix quern_solver::get_r() const
{
    return trans(get_r_trans());
}

Matrix quern_solver::get_r_trans() const
{
    return m_impl->get_r_trans();
};

unitary_matrix quern_solver::get_unitary_matrix() const
{
    return m_impl->get_unitary_matrix();
};

/*
MATCL_LINALG_EXPORT void matcl::test_qr_sparse()
{
    Integer M   = 10;
    Integer N   = 10;
    Integer K   = std::min(M,N);
    Matrix A    = csprandn(M,N,1.5/M);

    quern_solver qs(A,true,true,false);
    Matrix R    = qs.get_r();

    Matrix Q    = mmul(qs.get_unitary_matrix(),eye(K));
    disp(Q);

    disp(full(R));

    Matrix QR   = mmul(qs.get_unitary_matrix(),full(R));
    disp(full(QR));
    disp(full(A-QR));

    {
        Matrix X    = crandn(N,2);
        Matrix QX   = mmul(qs.get_unitary_matrix(), X, trans_type::no_trans);
        disp(full(QX));

        Matrix QX2  = Q*X;
        disp(full(QX2-QX));
    };
    {
        Matrix X    = crandn(M,2);
        Matrix QtX  = mmul(qs.get_unitary_matrix(), X, trans_type::trans);
        disp(full(QtX));

        Matrix QtX2 = trans(Q)*X;
        disp(full(QtX2-QtX));
    };
    {
        Matrix X    = crandn(M,2);
        Matrix QtX  = mmul(qs.get_unitary_matrix(), X, trans_type::conj_trans);
        disp(full(QtX));

        Matrix QtX2 = ctrans(Q)*X;
        disp(full(QtX2-QtX));
    };

    Matrix Q  = qs.get_unitary_matrix().to_matrix();

    Matrix QR = Q*R;
    
    disp(full(Q));    
    disp(A);
    disp(QR);

    Matrix dif = A - QR;
    Real n = norm(dif);

    disp(dif);
    disp(n);

    Matrix X    = randn(M,2);
    Matrix QX   = mmul(qs.get_unitary_matrix(),X);
    Matrix QTX  = mmul(ctrans(qs.get_unitary_matrix()),X);

    Matrix QX2      = Q*X;
    Matrix QTX2     = ctrans(Q)*X;

    Real n1 = norm(QX - QX2);
    Real n2 = norm(QTX - QTX2);

    disp(n1);
    disp(n2);
};
*/
}