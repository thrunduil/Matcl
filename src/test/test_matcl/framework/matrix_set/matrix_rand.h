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

#pragma once

#include <boost/shared_ptr.hpp>
#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace test
{

class rand_matrix
{
    public:
        virtual ~rand_matrix() {};

        virtual Integer         rand_scalar_int() = 0;
        virtual Real            rand_scalar_real() = 0;
        virtual Float           rand_scalar_float() = 0;
        virtual Complex         rand_scalar_compl() = 0;
        virtual Float_complex   rand_scalar_fcompl() = 0;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) = 0;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) = 0;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) = 0;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) = 0;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) = 0;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) = 0;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) = 0;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) = 0;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) = 0;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) = 0;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) = 0;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) = 0;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) = 0;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) = 0;        
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) = 0;

        virtual bool            is_nan_allowed() const = 0;

        virtual matcl::Matrix   rand_any_matrix()                   { return 0.; };
        virtual bool            generate_random_matrices() const    { return false; };
};

using rand_matrix_ptr = std::shared_ptr<rand_matrix>;


class rand_matrix_1 : public rand_matrix
{
    public:
        virtual Integer         rand_scalar_int() override;
        virtual Real            rand_scalar_real() override;
        virtual Float           rand_scalar_float() override;
        virtual Complex         rand_scalar_compl() override;
        virtual Float_complex   rand_scalar_fcompl() override;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) override;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) override;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) override;

        virtual bool is_nan_allowed() const override
        {
            return false;
        };
};

//rand only sparse matrix, struct and size is ignored
class rand_matrix_sparse : public rand_matrix
{
    private:
        rand_matrix_ptr         m_base_generator;

    public:
        rand_matrix_sparse(const rand_matrix_ptr& base_gen);

        virtual Integer         rand_scalar_int() override;
        virtual Real            rand_scalar_real() override;
        virtual Float           rand_scalar_float() override;
        virtual Complex         rand_scalar_compl() override;
        virtual Float_complex   rand_scalar_fcompl() override;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) override;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) override;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) override;

        virtual bool            is_nan_allowed() const override           { return false; };
        virtual bool            generate_random_matrices() const override { return true; };
        virtual matcl::Matrix   rand_any_matrix() override;
};

//rand only dense matrices, struct and size is ignored
class rand_matrix_dense : public rand_matrix
{
    private:
        rand_matrix_ptr         m_base_generator;

    public:
        rand_matrix_dense(const rand_matrix_ptr& base_gen);

        virtual Integer         rand_scalar_int() override;
        virtual Real            rand_scalar_real() override;
        virtual Float           rand_scalar_float() override;
        virtual Complex         rand_scalar_compl() override;
        virtual Float_complex   rand_scalar_fcompl() override;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) override;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) override;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) override;

        virtual bool            is_nan_allowed() const override           { return false; };
        virtual bool            generate_random_matrices() const override { return true; };
        virtual matcl::Matrix   rand_any_matrix();
};

class rand_matrix_obj : public rand_matrix
{
    public:
        virtual Integer         rand_scalar_int() override;
        virtual Real            rand_scalar_real() override;
        virtual Float           rand_scalar_float() override;
        virtual Complex         rand_scalar_compl() override;
        virtual Float_complex   rand_scalar_fcompl() override;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) override;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) override;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) override;

        virtual bool is_nan_allowed() const override
        {
            return false;
        };
};

class rand_matrix_1_sp : public rand_matrix
{
    private:
        bool with_nan;

    public:
        rand_matrix_1_sp(bool with_nan);

        virtual Integer         rand_scalar_int() override;
        virtual Real            rand_scalar_real() override;
        virtual Float           rand_scalar_float() override;
        virtual Complex         rand_scalar_compl() override;
        virtual Float_complex   rand_scalar_fcompl() override;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) override;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) override;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) override;

        virtual bool is_nan_allowed() const override
        {
            return with_nan;
        };
};

class rand_matrix_1_str : public rand_matrix
{
    public:
        rand_matrix_1_str();

        virtual Integer         rand_scalar_int() override;
        virtual Real            rand_scalar_real() override;
        virtual Float           rand_scalar_float() override;
        virtual Complex         rand_scalar_compl() override;
        virtual Float_complex   rand_scalar_fcompl() override;

        virtual matcl::Matrix   rand_dense_int(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_real(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_float(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_compl(Integer m,Integer n) override;
        virtual matcl::Matrix   rand_dense_fcompl(Integer m,Integer n) override;

        virtual matcl::Matrix   rand_sparse_int(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_real(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_float(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_compl(Integer m,Integer n, Real d) override;
        virtual matcl::Matrix   rand_sparse_fcompl(Integer m,Integer n, Real d) override;

        virtual matcl::Matrix   rand_band_int(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_real(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_float(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_compl(Integer m,Integer n, Integer fd, Integer ld) override;
        virtual matcl::Matrix   rand_band_fcompl(Integer m,Integer n, Integer fd, Integer ld) override;

        virtual bool is_nan_allowed() const override
        {
            return false;
        };
};

};};

