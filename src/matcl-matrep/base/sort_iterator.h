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

#include <boost/iterator/iterator_facade.hpp>
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-scalar/details/matfunc_helpers.h"

namespace matcl { namespace details
{

namespace mrd = matcl::raw::details;

template<class T>
struct less_nan
{
    static bool eval(const T& A, const T& B)
    {
        if (mrd::isnan_helper<T>::eval(A))
            return false;

        if (mrd::isnan_helper<T>::eval(B))
            return true;

        return (bool)mrd::lt_helper<T,T>::eval(A,B);
    };
};

template<class T, bool asceding>
struct compare_nan
{};

template<class T>
struct compare_nan<T,true>
{
    static bool eval(const T& A, const T& B)
    {
        return less_nan<T>::eval(A,B);
    };
};

template<class T>
struct compare_nan<T,false>
{
    static bool eval(const T& A, const T& B)
    {
        return less_nan<T>::eval(B,A);
    };
};

template<class T>
struct leq_nan
{
    static bool eval(const T& A, const T& B)
    {
        if (mrd::isnan_helper<T>::eval(A))
        {
            if (mrd::isnan_helper<T>::eval(B))
                return true;
            else
                return false;
        };

        if (mrd::isnan_helper<T>::eval(B))
            return true;

        return mrd::leq_helper<T,T>::eval(A,B);
    };
};

template<class T,class S> class value_ind_ref;
template<class T,class S>
class value_ind
{
    public:
        T			x;
        S			ind;		

    public:
        value_ind(const T& ptr, const S& index) : x(ptr),ind(index) {}

        value_ind(const value_ind_ref<T,S>& other) : x(other.get_x()),ind(other.get_ind()) {}

        bool operator<(const value_ind& other)
        {
            return less_nan<T>::eval(x,other.x);
        }

        value_ind& operator=(const value_ind& other)
        {
            mrd::reset_helper(x,other.x);
            mrd::reset_helper(ind,other.ind);
            return *this;
        }

        value_ind& operator=(const value_ind_ref<T,S>& other)
        {
            x	= other.get_x();
            ind = other.get_ind();
            return *this;
        }

};
template<class T,class S>
class value_ind_ref
{
    private:
        T*			x;
        S*			ind;		

    public:
        value_ind_ref(T& ptr, S& index) : x(&ptr),ind(&index) {};
        value_ind_ref(value_ind<T,S>& other) : x(&other.x),ind(&other.ind) {};

        bool operator<(const value_ind_ref& other) const
        {
            return less_nan<T>::eval(get_x(),other.get_x());
        }

        const value_ind_ref& operator=(const value_ind_ref& other) const
        {
            mrd::reset_helper(*x,*other.x);
            mrd::reset_helper(*ind,*other.ind);
            return *this;
        }

        const value_ind_ref& operator=(const value_ind<T,S>& other) const
        {
            mrd::reset_helper(*x,other.x);
            mrd::reset_helper(*ind,other.ind);
            return *this;
        }

        T& get_x() const { return *x; }
        S& get_ind() const { return *ind; }
};

template <class T,class S, bool asceding> 
struct value_ind_compare 
{ 
    bool operator()(const  value_ind_ref<T,S>& t1, const value_ind_ref<T,S>& t2) 
    { 
        return compare_nan<T,asceding>::eval(t1.get_x(),t2.get_x());
    } 
};

template<class T,class S>
class iterator_helper_2 : public boost::iterator_facade
                                <	iterator_helper_2<T,S>,
                                    value_ind<T,S>, 
                                    std::random_access_iterator_tag,
                                    value_ind_ref<T,S> const
                                >
{
    private:
        T*				x_ptr;
        S*				ind_ptr;
        Integer			ld;

        using base_type = boost::iterator_facade
                        <	iterator_helper_2<T,S>,
                            value_ind<T,S>,
                            std::random_access_iterator_tag,
                            value_ind_ref<T,S> const
                        >;

    public:
        using value_type        = typename base_type::value_type;
        using reference         = typename base_type::reference;
        using difference_type   = typename base_type::difference_type;
        using pointer           = typename base_type::pointer;
        using iterator_category = typename base_type::iterator_category;

    public:
        iterator_helper_2() : x_ptr(0),ind_ptr(0),ld(0) {};
        iterator_helper_2(T* x_ptr, S* ind_ptr, Integer ld) 
            :x_ptr(x_ptr),ind_ptr(ind_ptr), ld(ld)
        {};

    private:
        friend class boost::iterator_core_access;

        void increment() 
        { 
            x_ptr	+= ld; 
            ind_ptr	+= ld; 
        };

        void decrement() 
        { 
            x_ptr	-= ld; 
            ind_ptr	-= ld; 
        };

        bool equal(const iterator_helper_2 & other) const
        {
            return this->x_ptr == other.x_ptr;
        }

        reference dereference() const
        { 
            return value_ind_ref<T,S>(*x_ptr,*ind_ptr); 
        }

        difference_type distance_to(const iterator_helper_2& j) const
        {
            return (j.x_ptr - this->x_ptr)/ld;
        }

        iterator_helper_2 & advance(difference_type n)
        {
            this->x_ptr	    += imult(cast_int64(static_cast<Integer_64>(n)),ld);
            this->ind_ptr   += imult(cast_int64(static_cast<Integer_64>(n)),ld);

            return *this;
        };
};

};};