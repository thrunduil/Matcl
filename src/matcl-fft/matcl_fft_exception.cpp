#include "matcl-fft/matcl_fft_exception.h"
#include "matcl-core\error\exception_message.h"

namespace matcl { namespace error
{

const char* matcl_fft_exception::what() const throw()
{
    return m_msg.c_str();   
};

const char* matcl_fft_exception::what(exception_message& em) const
{
    return em.error_general(m_msg);
};

}}
