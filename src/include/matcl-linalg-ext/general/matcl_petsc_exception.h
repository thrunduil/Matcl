/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
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

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/general/linalg_exception.h"

#include <stdexcept>
#include <memory>

#pragma warning( push )
#pragma warning( disable: 4251 ) // class 'shared_ptr<>' needs to have dll-interface to be used by clients
#pragma warning( disable: 4275 ) // 

//TODO:

namespace matcl { namespace error
{

class exception_message_petsc;
using exception_message_petsc_ptr                  = std::shared_ptr<exception_message_petsc>;

MATCL_LINALG_EXPORT void                            set_global_messanger_petsc(exception_message_petsc_ptr msg);
MATCL_LINALG_EXPORT exception_message_petsc_ptr     get_global_messanger_petsc();

class MATCL_LINALG_EXPORT exception_message_petsc
{
    public:
        virtual ~exception_message_petsc(){};        

        virtual const char* multithreading_not_allowed() = 0;
        virtual const char* petsc_size_error(Integer expected, Integer got) = 0;
        virtual const char* internal_petsc_error(Integer err_code) = 0;
        virtual const char* internal_petsc_error(const std::string& err_msg) = 0;
        virtual const char* unrecognized_petsc_solver(const std::string& solver) = 0;
        virtual const char* unrecognized_petsc_preconditioner(const std::string& precond) = 0;
        virtual const char* invalid_argument() = 0;
        virtual const char* petsc_solution_std_not_available() = 0;
        virtual const char* unable_partition_external_preconditioner() = 0;
        virtual const char* preconditioner_is_not_partitioned() = 0;
        virtual const char* invalid_block_index(Integer block, Integer total_blocks) = 0;

        virtual const char* null_petsc_matrix() = 0;
        virtual const char* unable_convert_to_mmlib_shell() = 0;
        virtual const char* unable_convert_to_mmlib_unsupported_format(const std::string& format_name) = 0;
        virtual const char* unable_to_create_petsc_object(const std::string& obj_name) = 0;
        virtual const char* invalid_number_of_blocks(Integer nb, Integer nb_exp) = 0;
        virtual const char* invalid_mg_set_smoother(Integer level) = 0;
        virtual const char* invalid_multigrid_levels(Integer nlev) = 0;
        virtual const char* invalid_mg_set_solver(Integer current_level) = 0;
        virtual const char* invalid_mg_level(Integer lev) = 0;
        virtual const char* multigrid_solver_not_set() = 0;
        virtual const char* multigrid_smoother_not_set() = 0;
};

class MATCL_LINALG_EXPORT default_exception_message_petsc : public exception_message_petsc
{
    private:
        safe_string_message m_message;

    public:
        virtual ~default_exception_message_petsc(){};        

        virtual const char* multithreading_not_allowed() override;
        virtual const char* petsc_size_error(Integer expected, Integer got) override;
        virtual const char* internal_petsc_error(Integer err_code) override;
        virtual const char* internal_petsc_error(const std::string& err_msg) override;
        virtual const char* unrecognized_petsc_solver(const std::string& solver) override;
        virtual const char* unrecognized_petsc_preconditioner(const std::string& precond) override;
        virtual const char* invalid_argument() override;
        virtual const char* petsc_solution_std_not_available() override;
        virtual const char* unable_partition_external_preconditioner() override;
        virtual const char* preconditioner_is_not_partitioned() override;
        virtual const char* invalid_block_index(Integer block, Integer total_blocks) override;
        virtual const char* null_petsc_matrix() override;
        virtual const char* unable_convert_to_mmlib_shell() override;
        virtual const char* unable_convert_to_mmlib_unsupported_format(const std::string& format_name) override;
        virtual const char* unable_to_create_petsc_object(const std::string& obj_name) override;
        virtual const char* invalid_number_of_blocks(Integer nb, Integer nb_exp) override;
        virtual const char* invalid_mg_set_smoother(Integer level) override;
        virtual const char* invalid_multigrid_levels(Integer nlev) override;
        virtual const char* invalid_mg_set_solver(Integer current_level) override;
        virtual const char* invalid_mg_level(Integer lev) override;
        virtual const char* multigrid_solver_not_set() override;
        virtual const char* multigrid_smoother_not_set() override;
};

class MATCL_LINALG_EXPORT mmlib_petsc_error : public matcl_exception
{
    public:
        virtual const char* what() const throw() override;

    private:
        virtual const char* what(exception_message& em) const override;
        virtual const char* what(exception_message_petsc& em) const = 0;
};

class MATCL_LINALG_EXPORT multithreading_not_allowed : public mmlib_petsc_error
{
    public:
        multithreading_not_allowed(){};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT petsc_size_error : public mmlib_petsc_error
{
    private:
        Integer m_expected;
        Integer m_got;

    public:
        petsc_size_error(int expected, const int got)
            : m_expected(expected), m_got(got)
        {}

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_argument : public mmlib_petsc_error
{
    public:
        invalid_argument(){};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT internal_petsc_error : public mmlib_petsc_error
{
    private:
        Integer         m_err_code;
        std::string     m_err_msg;
        bool            m_from_code;

    public:
        internal_petsc_error(const int ierr)
            : m_err_code(ierr), m_from_code(true)
        {}
        internal_petsc_error(const std::string& mess)
            : m_err_msg(mess), m_from_code(false)
        {}

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT unrecognized_petsc_solver : public mmlib_petsc_error
{
    private:
        std::string m_solver;

    public:
        unrecognized_petsc_solver(const std::string& solver)
            : m_solver(solver)
        {}

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT unrecognized_petsc_preconditioner : public mmlib_petsc_error
{
    private:
        std::string     m_precond;

    public:
        unrecognized_petsc_preconditioner(const std::string& precond)
            : m_precond(precond)
        {}

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT petsc_solution_std_not_available : public mmlib_petsc_error
{
    public:
        petsc_solution_std_not_available(){};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT unable_partition_external_preconditioner : public mmlib_petsc_error
{
    public:
        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT preconditioner_is_not_partitioned : public mmlib_petsc_error
{
    public:
        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_block_index : public mmlib_petsc_error
{
    private:
        Integer     m_block;
        Integer     m_total_blocks;

    public:
        invalid_block_index(Integer block, Integer n_blocks)
            :m_block(block), m_total_blocks(n_blocks)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT null_petsc_matrix : public mmlib_petsc_error
{
    public:
        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT unable_convert_to_mmlib_shell : public mmlib_petsc_error
{
    public:
        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT unable_convert_to_mmlib_unsupported_format : public mmlib_petsc_error
{
    private:
        std::string m_format_name;

    public:
        unable_convert_to_mmlib_unsupported_format(const std::string& format)
            :m_format_name(format)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT unable_to_create_petsc_object : public mmlib_petsc_error
{
    private:
        std::string m_object_name;

    public:
        unable_to_create_petsc_object(const std::string& object_name)
            :m_object_name(object_name)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_number_of_blocks : public mmlib_petsc_error
{
    private:
        Integer n_blocks;
        Integer n_blocks_exp;

    public:
        invalid_number_of_blocks(Integer nb, Integer nb_exp)
            :n_blocks(nb), n_blocks_exp(nb_exp)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_mg_set_smoother : public mmlib_petsc_error
{
    private:
        Integer m_level;

    public:
        invalid_mg_set_smoother(Integer level)
            :m_level(level)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_multigrid_levels : public mmlib_petsc_error
{
    private:
        Integer m_levels;

    public:
        invalid_multigrid_levels(Integer nlev)
            :m_levels(nlev)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_mg_set_solver : public mmlib_petsc_error
{
    private:
        Integer m_current_level;

    public:
        invalid_mg_set_solver(Integer current_level)
            :m_current_level(current_level)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT invalid_mg_level : public mmlib_petsc_error
{
    private:
        Integer m_level;

    public:
        invalid_mg_level(Integer lev)
            :m_level(lev)
        {};

        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT multigrid_solver_not_set : public mmlib_petsc_error
{
    public:
        virtual const char* what(exception_message_petsc& em) const;
};

class MATCL_LINALG_EXPORT multigrid_smoother_not_set : public mmlib_petsc_error
{
    public:
        virtual const char* what(exception_message_petsc& em) const;
};

}}

#pragma warning( pop )

