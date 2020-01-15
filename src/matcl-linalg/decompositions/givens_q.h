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

#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/special_matrices/unitary_matrix.h"

namespace matcl { namespace details
{

//------------------------------------------------------------------------
//                          givens_q
//------------------------------------------------------------------------

template<class Val>
class givens_q : public unitary_matrix_data
{
    private:
        using VR    = typename md::real_type<Val>::type;

        using Mat   = raw::Matrix<Val,struct_dense>;
        using Mat_R = raw::Matrix<VR,struct_dense>;
        using Mat_I = raw::Matrix<Integer,struct_dense>;

    public:
        Integer                 m_mat_size;
        mr::const_matrix<Mat_R> m_C;
        mr::const_matrix<Mat>   m_S;
        mr::const_matrix<Mat_I> m_I;
        bool                    m_from_left;

    private:
        givens_q();

    public:
        givens_q(Integer N, const Mat_R& C, const Mat& S, const Mat_I& I, bool from_left);

        virtual ~givens_q();

        Integer                 seq_size() const                { return m_I.get().rows(); };

        virtual Integer         rows() const override           { return m_mat_size; }
        virtual Integer         cols() const override           { return m_mat_size; }
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

}};