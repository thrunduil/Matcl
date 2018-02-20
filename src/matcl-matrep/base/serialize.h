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

#include "matcl-core/IO/serialize.h"

namespace matcl { namespace details
{

template<class Archive>
void serialize_save(Archive&, ti::ti_empty, const unsigned int)
{};

template<class Archive>
void serialize_load(Archive&, ti::ti_empty&, const unsigned int)
{};

template<class Archive>
void serialize_save(Archive& ar, ti::ti_object ti, const unsigned int)
{
   ar << ti;
};

template<class Archive>
void serialize_load(Archive& ar, ti::ti_object& ret_ti, const unsigned int)
{
    ar >> ret_ti;
};

};};

