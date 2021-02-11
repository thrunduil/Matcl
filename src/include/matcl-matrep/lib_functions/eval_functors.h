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
#include "matcl-scalar/object.h"

namespace matcl
{

// test functions used in all, any, find functions
class MATCL_MATREP_EXPORT test_function
{
    public:
        virtual ~test_function() {};

        // test integer scalar
        virtual bool eval(Integer) const = 0;

        // test double scalar
        virtual bool eval(Real) const = 0;

        // test float scalar
        virtual bool eval(Float) const = 0;

        // test double complex scalar
        virtual bool eval(const Complex&) const = 0;

        // test float complex scalar
        virtual bool eval(const Float_complex&) const = 0;

        // test object scalar
        virtual bool eval(const Object&) const = 0;
};

// template version of test_function
// derived class must implement function
//     template<class Type> bool eval_templ(const Type& t) const
template<class derived>
class test_function_templ : public test_function
{
    private:
        virtual bool eval(Integer val) const override;
        virtual bool eval(Real val) const override;
        virtual bool eval(Float val) const override;
        virtual bool eval(const Complex& val) const override;
        virtual bool eval(const Float_complex& val) const override;
        virtual bool eval(const Object&) const override;

    private:
        template<class Type>
        bool eval_impl(const Type& t) const
        {
            return static_cast<const derived*>(this)->eval_templ(t);
        };
};

// test_function used in functions defined for typed matrices (i.e. dense_matrix,
// sparse_matrix, band_matrix, object_matrix); for these matrices types of stored
// elements are known and only test for this type must be implemented.
template<class T>
class test_type_function : private test_function_templ<test_type_function<T>>
{
    public:
        using base_type = test_function_templ<test_type_function<T>>;

    public:
        // check if given object passes the test;
        // this function must be implemented in derived class
        virtual bool        eval_type(const T& val) const = 0;

        // conversion to base class; this conversion is explicit, this
        // test function cannot be used in functions that expect 
        // test_function class unless it is guaranteed, that stored elements
        // have type T
        const base_type&    to_test_function() const        { return *this; };

    private:
        bool                eval_templ(const T& val) const  { return this->eval_type(val); }
        template<class S>
        bool                eval_templ(const S&) const      { return false; };

        friend base_type;
};

// test_function used in functions defined for object_matrix; for these matrices
// types of stored elements are known and only test for this type must be implemented.
template<class T>
class test_object_function : private test_function_templ<test_object_function<T>>
{
    public:
        using base_type = test_function_templ<test_object_function<T>>;

    public:
        // check if given object passes the test;
        // this function must be implemented in derived class
        virtual bool        eval_type(const dynamic::object_type<T>& val) const = 0;

        // conversion to base class; this conversion is explicit, this
        // test function cannot be used in functions that expect 
        // test_function class unless it is guaranteed, that stored elements
        // have type T
        const base_type&    to_test_function() const        { return *this; };

    private:
        bool eval_templ(const Object& val) const
        { 
            const matcl::dynamic::object& tmp = val;
            return this->eval_type(reinterpret_cast<const object_type<T>&>(tmp)); 
        }

        template<class S>
        bool eval_templ(const S&) const
        { 
            return false; 
        };

        friend base_type;
};

};

#include "matcl-matrep/details/eval_functors.inl"
