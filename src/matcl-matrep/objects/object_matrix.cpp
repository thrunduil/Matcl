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

#include "matcl-matrep/objects/object_matrix.h"
#include "matcl-core/IO/archive.h"

//#include "matcl-matrep/matcl_header.h"

namespace matcl
{

template<class T>
T operator/(const T&, Integer);

template<class T>
void test_compile(const object_matrix<T>& v)
{
    disp(v);

    {
        std::ostringstream os;
        std::istringstream is;

        os << v;

        object_matrix<T> tmp;
        is >> tmp;

        oarchive oa(os);
        save(oa, v);

        iarchive ia(is);
        load(ia, tmp);
    };

    using vect  = std::vector<object_matrix<T>>;

    delrows(v, colon());
    delcols(v,colon());
    delrowscols(v,colon(), colon());
    horzcat(v,v);
    vertcat(v,v);
    horzcat({v,v});
    vertcat({v,v});
    blkdiag({v,v});

    horzcat(vect{v,v});
    vertcat(vect{v,v});
    blkdiag(vect{v,v});

    repmat(v,1,1);

    full(v);
    sparse(v);
    band(v);

    clone(v);
    trans(v);
    ctrans(v);
    trans(v,trans_type::trans);
    trans(v,trans_type_ext::trans);
    vec(v);
    tril(v);
    triu(v);
    flipud(v);
    fliplr(v);
    reshape(v,1,1);
    get_diag(v);
    rot90(v);
    find(v);
    find(v, *(test_function*)nullptr);
    find(v, *(test_object_function<T>*)nullptr);
    find2(v);
    find2(v,*(test_function*)nullptr);
    find2(v, *(test_object_function<T>*)nullptr);
    find3(v);
    find3(v,*(test_function*)nullptr);
    find3(v, *(test_object_function<T>*)nullptr);
    sort(v);
    sort2(v);
    sortrows(v);
    sortrows2(v);
    sortrows(v,Matrix());
    sortrows2(v,Matrix());
    sortcols(v);
    sortcols2(v);
    sortcols(v,Matrix());
    sortcols2(v,Matrix());
    issorted(v);
    drop_sparse(v,0);
    convert_object<Float>(v);

    nnz(v, 1);
    all(v);
    all(v,*(test_function*)nullptr);
    all(v, *(test_object_function<T>*)nullptr);
    any(v);
    any(v,*(test_function*)nullptr);
    any(v, *(test_object_function<T>*)nullptr);
    sum(v);
    prod(v);
    cumsum(v);
    cumprod(v);
    min_d(v);
    max_d(v);
    min2(v);
    max2(v);
    min(v);
    max(v);

    using ret = object_matrix<T>;
    ret tmp;

    using T2 = decltype(decltype(abs(T()))()/Integer());
    using T3 = decltype(abs(T()));

    tmp = mean(v);
    
    object_matrix<T2> tmp2 = matcl::std<T>(v);
    object_matrix<T3> tmp3 = min_abs_d(v);
    object_matrix<T3> tmp4 = max_abs_d(v);
    object_matrix<T3> tmp5 = min_abs2(v).get<1>();
    object_matrix<T3> tmp6 = max_abs2(v).get<1>();

    object_matrix<T3> tmp7 = min_abs(v);
    object_matrix<T3> tmp8 = max_abs(v);

    nnz_vec(v);
    all_vec(v);
    all_vec(v,*(test_function*)nullptr);
    all_vec(v,*(test_object_function<T>*)nullptr);
    any_vec(v);
    any_vec(v,*(test_function*)nullptr);
    any_vec(v, *(test_object_function<T>*)nullptr);
    sum_vec(v);
    prod_vec(v);
    min_vec(v);
    max_vec(v);
    min2_vec(v);
    max2_vec(v);
    tmp = mean_vec(v);
    object_matrix<T2> tmp9 = std_vec(v);
    object_matrix<T3> tmp10 = min_abs_vec(v);
    object_matrix<T3> tmp11 = max_abs_vec(v);
    object_matrix<T3> tmp12 = min_abs2_vec(v).get<1>();
    object_matrix<T3> tmp13 = max_abs2_vec(v).get<1>();
};

template class matcl::object_matrix<Integer>;
template class matcl::object_matrix<Real>;
template class matcl::object_matrix<Float>;
template class matcl::object_matrix<Complex>;
template class matcl::object_matrix<Float_complex>;
//template class MATCL_MATREP_EXPORT matcl::object_matrix<Object>;

template class matcl::sub_object_matrix<Integer>;
template class matcl::sub_object_matrix<Real>;
template class matcl::sub_object_matrix<Float>;
template class matcl::sub_object_matrix<Complex>;
template class matcl::sub_object_matrix<Float_complex>;
//template class MATCL_MATREP_EXPORT matcl::sub_object_matrix<Object>;

template class matcl::sub_object_matrix_1<Integer>;
template class matcl::sub_object_matrix_1<Real>;
template class matcl::sub_object_matrix_1<Float>;
template class matcl::sub_object_matrix_1<Complex>;
template class matcl::sub_object_matrix_1<Float_complex>;
//template class MATCL_MATREP_EXPORT matcl::sub_object_matrix_1<Object>;

template class matcl::sub_object_matrix_2<Integer>;
template class matcl::sub_object_matrix_2<Real>;
template class matcl::sub_object_matrix_2<Float>;
template class matcl::sub_object_matrix_2<Complex>;
template class matcl::sub_object_matrix_2<Float_complex>;
//template class MATCL_MATREP_EXPORT matcl::sub_object_matrix_2<Object>;

template class matcl::object_row<Integer>;
template class matcl::object_row<Real>;
template class matcl::object_row<Float>;
template class matcl::object_row<Complex>;
template class matcl::object_row<Float_complex>;

template class matcl::object_col<Integer>;
template class matcl::object_col<Real>;
template class matcl::object_col<Float>;
template class matcl::object_col<Complex>;
template class matcl::object_col<Float_complex>;

template void test_compile(const object_matrix<Integer>& v);
template void test_compile(const object_matrix<Real>& v);
template void test_compile(const object_matrix<Float>& v);
template void test_compile(const object_matrix<Complex>& v);
template void test_compile(const object_matrix<Float_complex>& v);

};
