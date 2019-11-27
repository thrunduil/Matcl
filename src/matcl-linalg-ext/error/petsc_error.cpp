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

#include "matcl-linalg/general/mmlib_petsc_exception.h"
#include <sstream>

namespace matcl { namespace error
{

static exception_message_petsc_ptr messanger(new default_exception_message_petsc());

void set_global_messanger_petsc(exception_message_petsc_ptr msg)
{
    if (msg)
    {
        messanger = msg;
    };
};
exception_message_petsc_ptr get_global_messanger_petsc()
{
    return messanger;
};

const char* mmlib_petsc_error::what() const throw()
{
    m_message = what(*get_global_messanger_petsc());
    return m_message.c_str();
};

const char* mmlib_petsc_error::what(exception_message&) const
{
    m_message = "invalid exception_message";
    return m_message.c_str();
};

const char* default_exception_message_petsc::multithreading_not_allowed()
{
    m_message = "multithreaded use of petsc tools is forbidden";
    return m_message.c_str();
};

const char* default_exception_message_petsc::petsc_size_error(Integer expected, Integer got)
{
    std::ostringstream os;
    os << "wrong size of petsc object; expected: " << expected << " got: " << got;

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::internal_petsc_error(Integer err_code)
{
    std::ostringstream os;
    os << "petsc returned code " << err_code;

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::internal_petsc_error(const std::string& err_msg)
{
    m_message = err_msg;
    return m_message.c_str();
};

const char* default_exception_message_petsc::unrecognized_petsc_solver(const std::string& solver)
{
    std::ostringstream os;
    os << "unrecognized petsc solver: " << solver;

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::unrecognized_petsc_preconditioner(const std::string& precond)
{
    std::ostringstream os;
    os << "unrecognized petsc preconditioner: " << precond;

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::invalid_argument()
{
    m_message = "invalid argument for petsc operation";
    return m_message.c_str();
};

const char* default_exception_message_petsc::petsc_solution_std_not_available()
{
    m_message = "estimation of standard error of solution to a least squares problem is not available; "
            "this estimation was not requested; set option opts::petsc::lsqr_std_est to calculate "
            "standard error";
    return m_message.c_str();
};

const char* default_exception_message_petsc::unable_partition_external_preconditioner()
{
    m_message = "preconditioner given by external operator cannot be partitioned";
    return m_message.c_str();
};

const char* default_exception_message_petsc::preconditioner_is_not_partitioned()
{
    m_message = "preconditioner is not divided on blocks";
    return m_message.c_str();
};
const char* default_exception_message_petsc::invalid_block_index(Integer block, Integer total_blocks)
{
    std::ostringstream os;
    os << "invalid block index: " << block << ", total number of blocks is: " << total_blocks;

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::null_petsc_matrix()
{
    m_message = "unable to convert Petsc matrix to mmlib: null matrix supplied";
    return m_message.c_str();
};

const char* default_exception_message_petsc::unable_convert_to_mmlib_shell()
{
    m_message = "unable to convert Petsc matrix to mmlib: shell matrix supplied";
    return m_message.c_str();
};

const char* default_exception_message_petsc
            ::unable_convert_to_mmlib_unsupported_format(const std::string& format_name)
{
    m_message = "unable to convert Petsc matrix to mmlib: unsupported format " + format_name;
    return m_message.c_str();
};

const char* default_exception_message_petsc::unable_to_create_petsc_object(const std::string& obj_name)
{
    m_message = "unable to create Petsc object: " + obj_name;
    return m_message.c_str();
};

const char* default_exception_message_petsc::invalid_number_of_blocks(Integer nb, Integer nb_exp)
{
    std::ostringstream os;
    os  << "invalid number of block: matrix is divided on " << nb << " blocks, expecting " << nb_exp
        << " blocks";

    m_message = os.str();
    return m_message.c_str();
};
const char* default_exception_message_petsc::invalid_mg_set_smoother(Integer level)
{
    std::ostringstream os;
    os  << "unable to set smoother on level: " << level << "; level must be greater than zero";

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::invalid_multigrid_levels(Integer nlev)
{
    std::ostringstream os;
    os  << "invalid number of multigrid levels: " << nlev << "; expecting lev >= 2";

    m_message = os.str();
    return m_message.c_str();
};

const char* default_exception_message_petsc::invalid_mg_set_solver(Integer current_level)
{
    std::ostringstream os;
    os  << "unable to set solver on level " << current_level << "; level must be zero";

    m_message = os.str();
    return m_message.c_str();
};
const char* default_exception_message_petsc::invalid_mg_level(Integer lev)
{
    std::ostringstream os;
    os  << "invalid current multigrid level: " << lev;

    m_message = os.str();
    return m_message.c_str();
};
const char* default_exception_message_petsc::multigrid_solver_not_set()
{
    std::ostringstream os;
    os  << "multigrid building not finished: coarsest solver is not set";

    m_message = os.str();
    return m_message.c_str();
};
const char* default_exception_message_petsc::multigrid_smoother_not_set()
{
    std::ostringstream os;
    os  << "multigrid building not finished: smoother is not set";

    m_message = os.str();
    return m_message.c_str();
};


const char* multithreading_not_allowed::what(exception_message_petsc& em) const
{
    return em.multithreading_not_allowed();
};

const char* petsc_size_error::what(exception_message_petsc& em) const
{
    return em.petsc_size_error(m_expected, m_got);
}

const char* invalid_argument::what(exception_message_petsc& em) const
{
    return em.invalid_argument();
};

const char* internal_petsc_error::what(exception_message_petsc& em) const
{
    if (m_from_code)
        return em.internal_petsc_error(m_err_code);
    else
        return em.internal_petsc_error(m_err_msg);
}

const char* unrecognized_petsc_solver::what(exception_message_petsc& em) const
{
    return em.unrecognized_petsc_solver(m_solver);
}

const char* unrecognized_petsc_preconditioner::what(exception_message_petsc& em) const
{
    return em.unrecognized_petsc_preconditioner(m_precond);
}

const char* petsc_solution_std_not_available::what(exception_message_petsc& em) const
{
    return em.petsc_solution_std_not_available();
}

const char* unable_partition_external_preconditioner::what(exception_message_petsc& em) const
{
    return em.unable_partition_external_preconditioner();
}

const char* preconditioner_is_not_partitioned::what(exception_message_petsc& em) const
{
    return em.preconditioner_is_not_partitioned();
}

const char* invalid_block_index::what(exception_message_petsc& em) const
{
    return em.invalid_block_index(m_block, m_total_blocks);
}

const char* null_petsc_matrix::what(exception_message_petsc& em) const
{
    return em.null_petsc_matrix();
}

const char* unable_convert_to_mmlib_shell::what(exception_message_petsc& em) const
{
    return em.unable_convert_to_mmlib_shell();
}

const char* unable_convert_to_mmlib_unsupported_format::what(exception_message_petsc& em) const
{
    return em.unable_convert_to_mmlib_unsupported_format(m_format_name);
}
const char* unable_to_create_petsc_object::what(exception_message_petsc& em) const
{
    return em.unable_to_create_petsc_object(m_object_name);
}
const char* invalid_number_of_blocks::what(exception_message_petsc& em) const
{
    return em.invalid_number_of_blocks(n_blocks, n_blocks_exp);
}

const char* invalid_mg_set_smoother::what(exception_message_petsc& em) const
{
    return em.invalid_mg_set_smoother(m_level);
}

const char* invalid_multigrid_levels::what(exception_message_petsc& em) const
{
    return em.invalid_multigrid_levels(m_levels);
}

const char* invalid_mg_set_solver::what(exception_message_petsc& em) const
{
    return em.invalid_mg_set_solver(m_current_level);
}

const char* invalid_mg_level::what(exception_message_petsc& em) const
{
    return em.invalid_mg_level(m_level);
}

const char* multigrid_solver_not_set::what(exception_message_petsc& em) const
{
    return em.multigrid_solver_not_set();
}

const char* multigrid_smoother_not_set::what(exception_message_petsc& em) const
{
    return em.multigrid_smoother_not_set();
}

}}
