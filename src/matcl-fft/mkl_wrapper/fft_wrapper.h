#pragma once

#include <complex>

#include "matcl-fft/fft_context.h"

namespace matcl { namespace fft { namespace details { namespace mkl_fft
{

using fft_int       = int;
using fft_cplx      = std::complex<double>;
using fft_real      = double;
using fft_bool      = bool;
using fft_status    = long;

class fft_wrapper
{
    public:
        void        error_string(fft_status info, std::ostream& os);
        void*       make_mkl_context(Integer length, fft_context::fft_type type);
        void        destroy_mkl_context(void* context);

        fft_status  eval_inv_conj_even (const fft_cplx* in, fft_int in_ld, fft_real* out, fft_int out_ld, 
                             fft_int rows, fft_int cols, fft_int dim, fft_real scale, fft_context& cont);

        fft_status  eval_comp(const fft_cplx* in, fft_int in_ld, fft_cplx* out, fft_int out_ld, 
                             fft_int rows, fft_int cols, fft_int dim, fft_real scale, bool inverse,
                             bool inplace, fft_context& cont);

        fft_status  eval_real_pack(const fft_real* in, fft_int in_ld, fft_real* out, fft_int out_ld, 
                             fft_int rows, fft_int cols, fft_int dim, fft_real scale, bool inverse, 
                             bool inplace, fft_context& cont);

        fft_status  eval_real(const fft_real* in, fft_int in_ld, fft_cplx* out, fft_int out_ld, 
                             fft_int rows, fft_int cols, fft_int dim, fft_real scale, fft_context& cont);

        fft_status  eval_real_inv(fft_cplx* out, fft_int out_ld, fft_int rows, fft_int cols, fft_int dim,
                             fft_real scale, fft_context& cont);

        //inplace transform, out as input in packed format
        fft_status  eval_real_inv_conj_even(fft_real* out, fft_int out_ld, fft_int rows, fft_int cols,
                            fft_int dim, fft_real scale, fft_context& cont); 
};

}}}}
