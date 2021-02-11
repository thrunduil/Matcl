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

#include "matcl-blas-lapack/blas/blas.h"

namespace matcl { namespace lapack
{

class optim_params
{
    private:
        static const i_type     m_gghrd_block               = 32;
        static const i_type     m_gghrd_mult_block          = 64;

        //static const i_type     m_gghrd_block               = 2;
        //static const i_type     m_gghrd_mult_block          = 2;

    public:
        static const i_type     dtgsen3_max_eig_in_window   = 15;
        static const i_type     dtgsen3_max_window_size     = 2*dtgsen3_max_eig_in_window;
        static const i_type     dtrsen3_max_eig_in_window   = 30;
        static const i_type     dtrsen3_max_window_size     = 2*dtrsen3_max_eig_in_window;
        static const i_type     getrfr_block                = 64;
        static const i_type     potfp3_block                = 64;
        
        static const i_type     rotseq_min_length_to_acc    = 8;

        static const i_type     rotseq_panel_col_size_max   = 256;
        static const i_type     rotseq_panel_col_size_min   = 128;
        //static const i_type     rotseq_panel_col_size_max   = 32;
        //static const i_type     rotseq_panel_col_size_min   = 16;

        static const i_type     gghrd_block()               { return m_gghrd_block; };
        static const i_type     gghrd_mult_block()          { return m_gghrd_mult_block; };
};

};}
