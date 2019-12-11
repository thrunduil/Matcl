#pragma once

#include "matcl-fft/matcl_fft_config.h"
#include "matcl-core/error/exception.h"

#pragma warning(push)
#pragma warning(disable: 4251)  //needs to have dll-interface to be used by clients

namespace matcl { namespace error
{

// base class for exceptions thrown in matcl-fft library
class MATCL_FFT_EXPORT matcl_fft_exception : public matcl_exception
{
    private:
        std::string     m_msg;

    public:
        matcl_fft_exception(const std::string& msg)  
            : m_msg(msg)
        {};

        // get exception message
        virtual const char* what() const throw() override;
        virtual const char* what(exception_message& em) const override;
};

}}

#pragma warning(pop)