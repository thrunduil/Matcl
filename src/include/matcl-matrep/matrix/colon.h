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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/matrix/permvec.h"

namespace matcl
{

namespace details
{
    struct make_index_1_helper;
    struct make_index_2_helper;
};

// class representing end subscript
class colonend 
{ 
    private:
        // offset from last element
        Integer     offset;        

        friend class colon;

    public:
        // construct end with zero offset
        colonend()  : offset(0) {};

        // get position in array of length s
        Integer     apply(Integer s) const { return s + offset; };

        // operator + and - defined for end representing subscript
        // end + k, end - k
        friend colonend operator+(colonend old, Integer val);
        friend colonend operator-(colonend old, Integer val);
        friend colonend operator+(Integer val, colonend old);
        friend colonend operator-(Integer val, colonend old);
};

namespace 
{
    // value representing colon end
    colonend end;
};

// class representing matrix subscript
class MATCL_MATREP_EXPORT colon
{
    //internal use
    public:    
        using cend              = colonend;

        struct matrix_size_type
        {
            Integer rows, cols;
        };

        enum colon_type { t_all, t_one, t_range_simple, t_range_mat, t_end, t_rev_end, 
                          t_end_end, t_last, t_matrix1, t_matrix2 };
        
        Integer                 m_s;
        Integer                 m_i;
        Integer                 m_e;
        colon_type              m_flag;

        union mat_or_size
        {
            // integer dense matrix
            Matrix*             m_imat_12;
            matrix_size_type*   m_mat_size;
        };

        mat_or_size             m_mat_size;

        friend details::make_index_1_helper;
        friend details::make_index_2_helper;

    public:
        // represent A(:), A(:,i)
        colon()                                 :m_mat_size{nullptr},m_s(0),m_i(0),m_e(0),m_flag(t_all){};

        // represent A(s), A(s,i)
        colon(Integer s)                        :m_mat_size{nullptr},m_s(s),m_i(1),m_e(s),m_flag(t_one){};

        // represent A(end), A(end, i)
        colon(colonend e)                       :m_mat_size{nullptr},m_s(0),m_i(0),m_e(e.offset),m_flag(t_last){};

        // represent A(s:e), A(s:e, i)
        colon(Integer s, Integer e)             :m_mat_size{nullptr},m_s(s),m_i(1),m_e(e),m_flag(t_range_simple){};

        // represent A(s:end), A(s:end,i)
        colon(Integer s, colonend e)            :m_mat_size{nullptr},m_s(s),m_i(1),m_e(e.offset),m_flag(t_end){};

        // represent A(s:i:e), A(s:i:e,j)
        colon(Integer s, Integer i, Integer e)  :m_mat_size{nullptr},m_s(s),m_i(i),m_e(e),m_flag(t_range_simple){};

        // represent A(s:i:end), A(s:i:end,j)
        colon(Integer s, Integer i, cend e)     :m_mat_size{nullptr},m_s(s),m_i(i),m_e(e.offset),m_flag(t_end){};

        // represent A(end:i:e), A(end:i:e,j)
        colon(cend s, Integer i, Integer e)     :m_mat_size{nullptr},m_s(s.offset),m_i(i),m_e(e),m_flag(t_rev_end){};

        // represent A(end:e), A(end:e,j)
        colon(cend s, Integer e)                :m_mat_size{nullptr},m_s(s.offset),m_i(1),m_e(e),m_flag(t_rev_end){};

        // represent A(end:i:end), A(end:i:end,j), for example  A(end-2:1:end)
        colon(cend s, Integer i,cend e)         :m_mat_size{nullptr},m_s(s.offset),m_i(i),m_e(e.offset)
                                                ,m_flag(t_end_end){};

        // represent A(end:end), A(end:end,j), for example  A(end-2:end)
        colon(cend s, cend e)                   :m_mat_size{nullptr},m_s(s.offset),m_i(1),m_e(e.offset)
                                                ,m_flag(t_end_end){};

        // represent A(mat), A(mat,i), where mat is any integer matrix 
        // containing element positions
        colon(const Matrix& mat);
        colon(Matrix&& mat);

        // matrix substript from permutation vector
        colon(const permvec& mat);
        colon(permvec&& mat);

        // create colon from any type convertible to Matrix (perfect forwarding)
        template<class T, class Enable = typename details::enable_convertible_to_matrix<T,void>::type>
        inline colon(T&& mat);

        // standard destructor
        ~colon();

        // move constructor
        colon(colon&&);

        // move assignment
        colon& operator=(colon&&);

        // copyt costructor and assignment are deleted
        colon(const colon&)             = delete;
        colon& operator=(const colon&)  = delete;

        // create indenendent copy
        colon   copy() const;

    private:
        // TODO: move to public when finished

        // represent A(mat_1 + (mat_2 - 1) * rows) without forming single index;
        // mat_1, mat_2 are integer matrices of equal size
        colon(const Matrix& mat_1, const Matrix& mat_2);

    private:
        void    create(const Matrix&);
        void    create(Matrix&&);
        void    create(const Matrix& mat1, const Matrix& mat2);
        Integer get_matrix_rows() const;
        Integer get_matrix_cols() const;
};

// --------------------------------------------------------------------
//                     inline functions
// --------------------------------------------------------------------
inline colonend operator+(colonend old, Integer val)
{
    colonend out;
    out.offset = old.offset + val;
    return out;
};

inline colonend operator+(Integer val, colonend old)
{
    colonend out;
    out.offset = old.offset + val;
    return out;
};

inline colonend operator-(colonend old, Integer val)
{
    colonend out;
    out.offset = old.offset - val;
    return out;
};

inline colonend operator-(Integer val, colonend old)
{
    colonend out;
    out.offset = val - old.offset;
    return out;
};

template<class T, class Enable>
inline colon::colon(T&& mat)
    :colon(Matrix(std::forward<T>(mat)))
{};

};
