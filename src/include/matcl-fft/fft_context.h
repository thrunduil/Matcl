#pragma once

#include "matcl-fft/matcl_fft_config.h"
#include "matcl-fft/matcl_fft_exception.h"
#include "matcl-core/matrix/scalar_types.h"

#pragma warning(push)
#pragma warning(disable:4251)   //needs to have dll-interface to be used by clients of class 

namespace matcl { namespace fft
{

namespace details
{
    class fft_context_impl;
};

// Stores internal data used by fft and dct routines. The same context
// can be used for computation of different type. Then data required by
// given computation are appended to this context. Memory is released, when
// given context is no longer in use. One should not use the same context
// for many computations of different types, since this results in large
// memory consumption. Problem type is determined by its length and type.
// Objects of this type cannot be shared between threads.
class MATCL_FFT_EXPORT fft_context
{
    private:
        std::shared_ptr<details::fft_context_impl>  m_impl;

    public:
        // problem type
        enum fft_type { dct_1, dct_23_1, dct_23_2, fft_complex, fft_real};

    public:
        // create empty context
        fft_context();

        // destroy all data stored
        ~fft_context();

        // append this context with data for a problem of given length and type
        void        initialize(Integer length, fft_type type);

    // internal use only
    public:
        const Real* get_table(Integer length, fft_type type) const;
        void*       initialize_fft(Integer length, fft_type type);
};

}}

#pragma warning(pop)