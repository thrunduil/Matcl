#pragma once

#include <iosfwd>

#include "mkgen/mkgen_fwd.h"
#include "mkgen/TODO/utils/utils.h"
#include "mkgen/TODO/expression/colon.h"

namespace matcl { namespace mkgen
{

}};

namespace matcl { namespace mkgen { namespace details
{

//--------------------------------------------------------
//              print priorities
//--------------------------------------------------------
static const int prior_start    = 0;
static const int prior_assign   = 10;
static const int prior_plus     = 20;
static const int prior_mult     = 30;


}}}
