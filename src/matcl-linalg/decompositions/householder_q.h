/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl { namespace details
{

//------------------------------------------------------------------------
//                          householder_q
//------------------------------------------------------------------------

template<class Val>
class householder_q : public unitary_matrix_data
{
    private:
        using Mat = raw::Matrix<Val,struct_dense>;

    public:
        Integer                 m_mat_cols;
        Integer                 m_ldiags;
        Integer                 m_offset;
        Mat                     m_reflectors;
        Mat                     m_tau_vec;

    private:
        householder_q();

    public:
        householder_q(Integer N, const Mat& Qc, const Mat& tau, Integer ld, Integer offset);

        virtual ~householder_q();

        Integer                 rep_rows() const                { return m_reflectors.rows(); };
        Integer                 number_reflectors() const       { return std::min(m_mat_cols, m_tau_vec.length()); };
        Integer                 reflector_length() const        { return rep_rows(); };
        virtual Integer         rows() const override           { return rep_rows(); }
        virtual Integer         cols() const override           { return m_mat_cols; }
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

//------------------------------------------------------------------------
//                          householder_band_q
//------------------------------------------------------------------------

template<class Val>
class householder_band_q : public unitary_matrix_data
{
    private:
        using Mat_D = raw::Matrix<Val,struct_dense>;
        using Mat_B = raw::Matrix<Val,struct_banded>;

    public:
        Integer                 m_mat_cols;
        Mat_B                   m_reflectors;
        Mat_D                   m_tau_vec;

    private:
        householder_band_q();

    public:
        householder_band_q(Integer N, const Mat_B& Qc, const Mat_D& tau);

        virtual ~householder_band_q();

        Integer                 rep_rows() const                { return m_reflectors.rows(); };
        Integer                 number_reflectors() const       { return std::min(m_mat_cols, m_tau_vec.length()); };
        Integer                 reflector_length() const        { return rep_rows(); };

        virtual Integer         rows() const override           { return rep_rows(); }
        virtual Integer         cols() const override           { return m_mat_cols; }
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

// if X is unique, than inplace version is used
void rand_unitary_impl(unitary_matrix& ret, const matcl::Matrix& X);

}};