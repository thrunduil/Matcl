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

//TODO: is it needed?
namespace data_struct
{

template<class T, class initializer>
class vector
{
    private:
        using block = char[sizeof(T)];

    private:
        T*          m_ptr;
        initializer m_init;
        size_t      m_size;

        inline void delete_table()
        {
            for(size_t i = 0; i < m_size; ++i)
                m_ptr[i].~T();

            delete[] reinterpret_cast<block*>(m_ptr);
            m_ptr = nullptr;
        }

    public:
        vector(const initializer& init)
            :m_ptr(nullptr), m_init(init), m_size(0)
        {};

        vector(const initializer& init,size_t size)
            :m_init(init),m_size(size)
        {
            m_ptr = reinterpret_cast<T*>(new block[size]);
            for (size_t i = 0; i < size; ++i)
            {
                m_init.initialize(&m_ptr[i]);
            };
        };

        ~vector()
        {
            delete_table();
        };

        vector(const vector& other)
            :m_init(other.m_init),m_size(other.m_size)
        {
            m_ptr = reinterpret_cast<T*>(new block[m_size]);
            for (size_t i = 0; i < m_size; ++i)
            {
                new(&m_ptr[i]) T(other.m_ptr[i]);
            };
        };

        vector& operator=(const vector& other)
        {
            if (m_ptr == other.m_ptr)
                return *this;

            delete_table();

            m_size = other.m_size;
            m_init = other.m_init;

            m_ptr = reinterpret_cast<T*>(new block[m_size]);
            for (size_t i = 0; i < m_size; ++i)
                new(&m_ptr[i]) T(other.m_ptr[i]);

            return *this;
        };

        T& operator[](size_t pos)
        {
            return m_ptr[pos];
        };

        const T& operator[](size_t pos) const
        {
            return m_ptr[pos];
        };

        void clear()
        {
            delete_table();
        };

        void resize(size_t new_size)
        {
            if (new_size == m_size)
                return;

            T* tmp_ptr = reinterpret_cast<T*>(new block[new_size]);
            for (size_t i = 0; i < m_size; ++i)
                new(&tmp_ptr[i]) T(m_ptr[i]);
            for (size_t i = m_size; i < new_size; ++i)
                m_init.initialize(&tmp_ptr[i]);

            delete_table();
            m_size = new_size;
            m_ptr= tmp_ptr;
        };

        void resize(size_t new_size, const T& val)
        {
            if (new_size == m_size)
                return;

            T* tmp_ptr = reinterpret_cast<T*>(new block[new_size]);
            for (size_t i = 0; i < m_size; ++i)
                new(&tmp_ptr[i]) T(m_ptr[i]);

            for (size_t i = m_size; i < new_size; ++i)
                new(&tmp_ptr[i]) T(val);

            delete_table();
            m_size = new_size;
            m_ptr= tmp_ptr;
        };
};

};