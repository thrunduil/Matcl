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

#include <list>
#include <string>
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-file/config_matcl_file.h"

#pragma warning( push )
#pragma warning(disable:4251)    // warning C4251: class 'X' needs to have dll-interface to be used by clients of class Y

namespace matcl
{

class MATCL_FILE_EXPORT blob_data
{
    private:
        void*           m_ptr;
        size_t          m_bytes;

    public:
        blob_data(void* ptr, size_t bytes);
        explicit blob_data(size_t bytes);

        ~blob_data()                            { delete[] m_ptr; };

        void*           pointer()               { return m_ptr; };
        const void*     pointer() const         { return m_ptr; };
        size_t          bytes() const           { return m_bytes; };

    private:
        blob_data(const blob_data& bl) = delete;
        blob_data& operator=(const blob_data& bl) = delete;
};

using blob_ptr = std::shared_ptr<blob_data>;

namespace details
{
    class matcl_file_data;
    using matcl_file_data_ptr = std::shared_ptr<matcl_file_data>;
};

class MATCL_FILE_EXPORT matrix_info
{
    private:
        Integer                m_flags;
        Integer                m_rows;
        Integer                m_cols;
        Integer                m_ldiags;
        Integer                m_udiags;
        Integer                m_nnz;

        matcl::value_code      m_value_type;
        matcl::struct_code     m_struct_type;
        matcl::mat_code        m_matrix_type;

        std::string            m_name;

    public:
        Integer                flags() const            { return m_flags; };
        Integer                rows() const             { return m_rows; };
        Integer                cols() const             { return m_cols; };
        Integer                structural_ldiags() const{ return m_ldiags; };
        Integer                structural_udiags() const{ return m_udiags; };
        Integer                structural_nnz() const   { return m_nnz; };

        matcl::value_code      get_value_code() const   { return m_value_type; };
        matcl::struct_code     get_struct_code() const  { return m_struct_type; };
        matcl::mat_code        get_matrix_code() const  { return m_matrix_type; };

        const std::string      name() const             { return m_name; };

        friend class matcl_file;
};

class MATCL_FILE_EXPORT matcl_file
{
    public:
        using matrix_list       = std::list<Matrix>;
        using data_list         = std::list<blob_ptr>;
        using matrix_info_list  = std::list<matrix_info>;
        using string_list       = std::list<std::string>;

    private:
        using matcl_file_data_ptr   = details::matcl_file_data_ptr;
    
    public:
        matcl_file(const std::string& file_name, open_mode om, thread_mode = thread_mode::multi_thread);
        ~matcl_file();

        void                timeout_limit(int msec);

        bool                exist(const std::string& mat_name);
        bool                exist_data(const std::string& data_name);

        Matrix              load(const std::string& mat_name);
        Matrix              load(const std::string& mat_name,std::string& mat_string);
        blob_ptr            load_data(const std::string& data_name);
        blob_ptr            load_data(const std::string& mat_name,std::string& data_string);

        std::string         load_mat_string(const std::string& mat_name);
        std::string         load_data_string(const std::string& data_name);
        matrix_info         load_info(const std::string& mat_name);

        matrix_list         load_all();
        data_list           load_data_all();
        void                load_all_info(matrix_info_list& ml);
        void                load_all_mat_string(string_list& sl);
        void                load_data_all_mat_string(string_list& sl);

        void                save(const Matrix& mat, const std::string& mat_name, 
                                const std::string& mat_string = "", bool allow_replace = true);
        void                save_data(const void* data, size_t bytes, const std::string& data_name, 
                                const std::string& data_string = "", bool allow_replace = true);

        void                remove(const std::string& mat_name);
        void                remove_data(const std::string& data_name);

        void                add_to_save_list(const Matrix& mat, const std::string& mat_name, 
                                const std::string& mat_string = "", bool allow_replace = true);
                            
                            //pointers must be valid, when save_list is called
        void                add_data_to_save_list(const void* data, size_t bytes, 
                                const std::string& data_name, const std::string& data_string = "", 
                                bool allow_replace = true);
        void                save_list();

        void                modify_mat_string(const std::string& mat_name, 
                                const std::string& mat_string);
        void                modify_data_mat_string(const std::string& data_name, 
                                const std::string& mat_string);

    private:
        matcl_file_data_ptr    m_data;

        void                insert_main_table_if_not_exist();
        void                insert_data_table_if_not_exist();
        void                insert_main_table();
        void                insert_data_table();
        bool                exist_data_table() const;
        matrix_info         get_mat_info(void* ptr);
};

};

#pragma warning( pop )
