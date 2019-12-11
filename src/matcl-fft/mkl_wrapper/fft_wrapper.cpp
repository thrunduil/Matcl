#include "fft_wrapper.h"
#include <mkl_dfti.h>

namespace matcl { namespace fft { namespace details { namespace mkl_fft
{

struct DFTI_desc_mmlib;

inline DFTI_Descriptor_struct* _rc(DFTI_desc_mmlib* in)
{
    return reinterpret_cast<DFTI_Descriptor_struct*>(in);
}

inline DFTI_Descriptor_struct** _rc(DFTI_desc_mmlib** in)
{
    return reinterpret_cast<DFTI_Descriptor_struct**>(in);
}

static fft_status common(DFTI_desc_mmlib* desc, bool inplace, 
                         fft_int ld_in, fft_int stride_in, 
                         fft_int ld_out, fft_int stride_out, 
                         fft_int number_of_transforms, 
                         fft_real forward_scale, fft_real backward_scale)
{
    fft_status status = 0;

    if(status == 0 ) 
    {
        if (inplace)
            status = DftiSetValue( _rc(desc), DFTI_PLACEMENT, DFTI_INPLACE);
        else
            status = DftiSetValue( _rc(desc), DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    }

    fft_int strides_in[2]   = {0, stride_in};
    fft_int strides_out[2]  = {0, stride_out};

    if(status == 0 ) 
        status = DftiSetValue( _rc(desc), DFTI_INPUT_DISTANCE, ld_in);
    
    if(status == 0 ) 
        status = DftiSetValue( _rc(desc), DFTI_OUTPUT_DISTANCE, ld_out);

    if(status == 0 ) 
        status = DftiSetValue( _rc(desc), DFTI_INPUT_STRIDES, strides_in);

    if(status == 0 ) 
        status = DftiSetValue( _rc(desc), DFTI_OUTPUT_STRIDES, strides_out);

    if(status == 0 )
        status = DftiSetValue( _rc(desc), DFTI_NUMBER_OF_TRANSFORMS, number_of_transforms);

    if(status == 0 )
        status = DftiSetValue( _rc(desc), DFTI_FORWARD_SCALE, forward_scale);

    if(status == 0 )
        status = DftiSetValue( _rc(desc), DFTI_BACKWARD_SCALE, backward_scale);

    if(status == 0 )
        status = DftiCommitDescriptor(_rc(desc));
    
    return status;
}
static fft_status common2(DFTI_desc_mmlib* desc, bool inplace, fft_int ld_in0, fft_int ld_out0, 
                          fft_int rows, fft_int cols, fft_int dim, fft_real scale, bool inverse)
{
    fft_int ld_in;
    fft_int stride_in;
    fft_int ld_out;
    fft_int stride_out;
    fft_int length;
    fft_int number_of_transforms;

    fft_real scale_forward;
    fft_real scale_backward;

    if (inverse == true)
    {
        scale_forward  = 1.0;
        scale_backward = scale;
    }
    else
    {
        scale_forward  = scale;
        scale_backward = 1.0;
    };

    if (dim == 1)
    {
        ld_in               = ld_in0;
        stride_in           = 1;
        ld_out              = ld_out0;    
        stride_out          = 1;
        length              = rows;
        number_of_transforms= cols;
    }
    else
    {
        ld_in               = 1;
        stride_in           = ld_in0;
        ld_out              = 1;    
        stride_out          = ld_out0;
        length              = cols;
        number_of_transforms= rows;
    };
 
    return common(desc, inplace, ld_in, stride_in, ld_out, stride_out, number_of_transforms, 
                    scale_forward, scale_backward);
};

struct mkl_context
{
    DFTI_desc_mmlib*    m_descriptor;
    bool                m_inplace;
    fft_int             m_ld_in0;
    fft_int             m_ld_out0;
    fft_int             m_rows;
    fft_int             m_cols;
    fft_int             m_dim;
    fft_real            m_scale;
    bool                m_inverse;
    int                 m_conj_storage;
    int                 m_packed_format;

    mkl_context();

    bool    need_reset(int conj_storage, bool inplace, fft_int ld_in0, fft_int ld_out0, 
                       fft_int rows, fft_int cols, fft_int dim, fft_real scale, bool inverse,
                       int packed_format = DFTI_CCS_FORMAT);
};
mkl_context::mkl_context()
{
    memset(this, 0, sizeof(mkl_context));
};

bool mkl_context::need_reset(int conj_storage, bool inplace, fft_int ld_in0, fft_int ld_out0, fft_int rows, fft_int cols, 
                    fft_int dim, fft_real scale, bool inverse, int packed_format)
{
    bool need = 
        false
        || m_inplace != inplace
        || m_ld_in0 != ld_in0
        || m_ld_out0 != ld_out0
        || m_rows != rows
        || m_cols != cols
        || m_dim != dim
        || m_scale != scale
        || m_inverse != inverse
        || m_conj_storage != conj_storage
        || m_packed_format != packed_format;
        ;

    if (need == true)
    {
        m_inplace           = inplace;
        m_ld_in0            = ld_in0;
        m_ld_out0           = ld_out0;
        m_rows              = rows;
        m_cols              = cols;
        m_dim               = dim;
        m_scale             = scale;
        m_inverse           = inverse;
        m_conj_storage      = conj_storage;
        m_packed_format     = packed_format;
    };

    return need;
};

void* fft_wrapper::make_mkl_context(Integer length, fft_context::fft_type type)
{    
    DFTI_desc_mmlib* desc   = nullptr;
    fft_status status       = 0;
 
    switch(type)
    {
        case fft_context::fft_real:
            status = DftiCreateDescriptor(_rc(&desc), DFTI_DOUBLE, DFTI_REAL, 1, length);
            break;
        case fft_context::fft_complex:
            status = DftiCreateDescriptor(_rc(&desc), DFTI_DOUBLE, DFTI_COMPLEX, 1, length);
            break;
        default:
            throw std::runtime_error("invalid fft context type");
    };

    if (status != 0)
    {
        std::ostringstream msg;
        msg << "unable to initialize fft library, reason: ";
        error_string(status, msg);

        throw error::matcl_fft_exception(msg.str());
    };

    mkl_context* cont       = new mkl_context();
    cont->m_descriptor      = desc;

    return cont;
};

void  fft_wrapper::destroy_mkl_context(void* context)
{
    mkl_context* cont       = reinterpret_cast<mkl_context*>(context);
    DFTI_desc_mmlib* desc   = (DFTI_desc_mmlib*) cont->m_descriptor;
    DftiFreeDescriptor(_rc(&desc));

    delete cont;
};

fft_status fft_wrapper::eval_comp(const fft_cplx* in, fft_int ld_in0, fft_cplx* out, fft_int ld_out0,
                                  fft_int rows, fft_int cols, fft_int dim, fft_real scale, bool inverse, 
                                  bool inplace, fft_context& cont)
{   
    fft_status status = 0;
    fft_int length;

    if (dim == 1)
        length              = rows;
    else
        length              = cols;

    mkl_context* ctx        = (mkl_context*)cont.initialize_fft(length, fft_context::fft_complex);
    DFTI_desc_mmlib* desc   = ctx->m_descriptor;

    bool need_reset         = ctx->need_reset(DFTI_COMPLEX_REAL, 
                                 inplace, ld_in0, ld_out0, rows, cols, dim, scale, inverse);

    if (need_reset == true)
    {
        if(status == 0) 
            status = common2(desc, inplace, ld_in0, ld_out0, rows, cols, dim, scale, inverse);
    };

    if(status == 0 ) 
    {
        if (inplace == false)
            status = inverse ? DftiComputeBackward( _rc(desc), const_cast<fft_cplx*>(in), out )
                             : DftiComputeForward(  _rc(desc), const_cast<fft_cplx*>(in), out );
        else
            status = inverse ? DftiComputeBackward( _rc(desc), out)
                             : DftiComputeForward(  _rc(desc), out);
    }

    return status;
}
fft_status fft_wrapper::eval_real(const fft_real* in, fft_int ld_in0, fft_cplx* out, fft_int ld_out0, 
                                  fft_int rows, fft_int cols, fft_int dim, fft_real scale, fft_context& cont)
{    
    bool inplace            = false;

    fft_int length;

    if (dim == 1)
        length              = rows;
    else
        length              = cols;

    //DFTI_CONJUGATE_EVEN_STORAGE
    mkl_context* ctx        = (mkl_context*)cont.initialize_fft(length, fft_context::fft_real);
    DFTI_desc_mmlib* desc   = ctx->m_descriptor;

    bool need_reset         = ctx->need_reset(DFTI_COMPLEX_COMPLEX, inplace, ld_in0, ld_out0, 
                                                 rows, cols, dim, scale, false);
    fft_status status = 0;       

    if (need_reset == true)
    {
        if(status == 0 )
            status = DftiSetValue( _rc(desc), DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);

        if(status == 0 )
            status = common2(desc, inplace, ld_in0, ld_out0, rows, cols, dim, scale, false);
    };

    if(status == 0 )
        status = DftiComputeForward(_rc(desc), const_cast<fft_real*>(in), out );

    return status;
}
fft_status fft_wrapper::eval_real_inv(fft_cplx* out, fft_int ld_out0, fft_int rows, fft_int cols,
                                      fft_int dim, fft_real scale, fft_context& cont)
{
    bool inplace            = true;
    fft_int length;

    if (dim == 1)
        length              = rows;
    else
        length              = cols;

    mkl_context* ctx        = (mkl_context*)cont.initialize_fft(length, fft_context::fft_complex);
    DFTI_desc_mmlib* desc   = ctx->m_descriptor;

    bool need_reset         = ctx->need_reset(DFTI_COMPLEX_REAL, 
                                inplace, ld_out0, ld_out0, rows, cols, dim, scale, true);

    fft_status status = 0;
    
    if (need_reset == true)
    {
        if(status == 0 )
            status = common2(desc, inplace, ld_out0, ld_out0, rows, cols, dim, scale, true);
    };

    if(status == 0 )
        status = DftiComputeBackward(_rc(desc), out );

    return status;
}
fft_status fft_wrapper::eval_real_inv_conj_even(fft_real* out, fft_int ld_out0, fft_int rows, fft_int cols,
                                      fft_int dim, fft_real scale, fft_context& cont)
{
    bool inplace            = true;
    fft_int length;

    if (dim == 1)
        length              = rows;
    else
        length              = cols;

    mkl_context* ctx        = (mkl_context*)cont.initialize_fft(length, fft_context::fft_real);
    DFTI_desc_mmlib* desc   = ctx->m_descriptor;

    bool need_reset         = ctx->need_reset(DFTI_COMPLEX_REAL, inplace, ld_out0, ld_out0, rows, cols,
                                              dim, scale, true, DFTI_PACK_FORMAT);

    fft_status status = 0;
    
    if (need_reset == true)
    {
        if(status == 0 )
            status = DftiSetValue( _rc(desc), DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL);

        if(status == 0 )
            status = DftiSetValue( _rc(desc), DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT);

        if(status == 0 )
            status = common2(desc, inplace, ld_out0, ld_out0, rows, cols, dim, scale, true);
    };

    if(status == 0 )
        status = DftiComputeBackward(_rc(desc), out );

    return status;
}

fft_status fft_wrapper::eval_real_pack(const fft_real* in, fft_int ld_in0, fft_real* out, fft_int ld_out0, 
                        fft_int rows, fft_int cols, fft_int dim, fft_real scale, bool inverse, 
                        bool inplace0, fft_context& cont)
{
    bool inplace            = inplace0;
    fft_int length;

    if (dim == 1)
        length              = rows;
    else
        length              = cols;

    mkl_context* ctx        = (mkl_context*)cont.initialize_fft(length, fft_context::fft_real);
    DFTI_desc_mmlib* desc   = ctx->m_descriptor;

    bool need_reset         = ctx->need_reset(DFTI_COMPLEX_REAL, inplace, ld_in0, ld_out0, rows, 
                                              cols, dim, scale, inverse, DFTI_PACK_FORMAT);

    fft_status status = 0;
    
    if (need_reset == true)
    {
        if(status == 0 )
            status = DftiSetValue( _rc(desc), DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL);

        if(status == 0 )
            status = DftiSetValue( _rc(desc), DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT);

        if(status == 0 )
            status = common2(desc, inplace, ld_in0, ld_out0, rows, cols, dim, scale, inverse);
    };

    if(status == 0 )
    {
        if (inplace == false)
        {
            if (inverse == false)
                status = DftiComputeForward(_rc(desc), const_cast<fft_real*>(in), out );
            else
                status = DftiComputeBackward(_rc(desc), const_cast<fft_real*>(in), out );
        }
        else
        {
            if (inverse == false)
                status = DftiComputeForward(_rc(desc), out );
            else
                status = DftiComputeBackward(_rc(desc), out );
        };
    }

    return status;
};
fft_status fft_wrapper::eval_inv_conj_even(const fft_cplx* in, fft_int ld_in0, fft_real* out, fft_int ld_out0,
                                           fft_int rows, fft_int cols, fft_int dim, fft_real scale, 
                                           fft_context& cont)
{
    fft_status status       = 0;

    bool inplace            = false;
    fft_int length;
    if (dim == 1)
        length              = rows;
    else
        length              = cols;

    mkl_context* ctx        = (mkl_context*)cont.initialize_fft(length, fft_context::fft_real);
    DFTI_desc_mmlib* desc   = ctx->m_descriptor;

    bool need_reset         = ctx->need_reset(DFTI_COMPLEX_COMPLEX,
                                    inplace, ld_in0, ld_out0, rows, cols, dim, scale, true);

    if (need_reset == true)
    {
        if(status == 0 )
            status = DftiSetValue( _rc(desc), DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);

        if(status == 0 )
            status = common2(desc, inplace, ld_in0, ld_out0, rows, cols, dim, scale, true);
    };

    if(status == 0 )
        status = DftiComputeBackward( _rc(desc), const_cast<fft_cplx*>(in), out );

    return status;
}

void fft_wrapper::error_string(fft_status info, std::ostream& os)
{
    switch(info)
    {
        case DFTI_MEMORY_ERROR:
            os << "memory error";
            break;
        case DFTI_INVALID_CONFIGURATION:
            os << "invalid configuration";
            break;
        case DFTI_INCONSISTENT_CONFIGURATION:
            os << "inconsistent configuration";
            break;
        case DFTI_NUMBER_OF_THREADS_ERROR:
            os << "number of threads error";
            break;
        case DFTI_MULTITHREADED_ERROR:
            os << "thread error";
            break;
        case DFTI_BAD_DESCRIPTOR:
            os << "bad destriptor";
            break;
        case DFTI_UNIMPLEMENTED:
            os << "not implemented";
            break;
        case DFTI_MKL_INTERNAL_ERROR:
            os << "internal error";
            break;
        default:
            os << "unknown error";
            break;
    };
}

}}}}
