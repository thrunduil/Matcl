#include "matcl-fft/matcl_dct.h"
#include "matcl-fft/matcl_fft.h"

//#include "mmlib/details/raw_fwd.h"
#include "matcl-internals/container/sparse_ccs.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/func/converter.h"

#include "matcl-blas-lapack/level1/level1.h"

namespace matcl { namespace fft
{

using matcl::level1::get_int;

struct prepare_dct1_dim1
{
    static void eval(Real* Y, Integer Y_ld, const Real* ptr_x,
                     Integer x_ld, Integer rows, Integer cols, const Real* cs);  
};

template<Integer Rows, Integer Y_LD, Integer X_LD>
struct prepare_dct1_dim2
{
    static void eval(Real* Y, Integer Y_ld, const Real* ptr_x,
                     Integer x_ld, Integer rows, Integer cols, const Real* cs);  
};

struct prepare_dct1_inpl_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);  
};

template<Integer Rows, Integer Y_LD>
struct prepare_dct1_inpl_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);  
};

struct prepare_dct2_dim1
{
    static void eval(Real* Y, Integer Y_ld, const Real* ptr_x,
                     Integer x_ld, Integer rows, Integer cols, const Real* cs);  
};

struct prepare_dct2_inpl_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);  
};

template<Integer Rows, Integer Y_LD, Integer X_LD>
struct prepare_dct2_dim2
{
    static void eval(Real* Y, Integer Y_ld, const Real* ptr_x,
                     Integer x_ld, Integer rows, Integer cols, const Real* cs);  
};

template<Integer Rows, Integer Y_LD>
struct prepare_dct2_inpl_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);  
};

struct prepare_dct3_dim1
{
    static void eval(Real* Y, Integer Y_ld, const Real* ptr_x,
                     Integer x_ld, Integer rows, Integer cols, const Real* cs);  
};
struct prepare_dct3_inpl_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);  
};

template<Integer Rows, Integer Y_LD, Integer X_LD>
struct prepare_dct3_dim2
{
    static void eval(Real* Y, Integer Y_ld, const Real* ptr_x,
                     Integer x_ld, Integer rows, Integer cols, const Real* cs);  
};
template<Integer Rows, Integer Y_LD>
struct prepare_dct3_inpl_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);  
};

struct recurrence_dct1_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, Real scale);
};
template<Integer Rows, Integer Y_LD>
struct recurrence_dct1_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, Real scale);
};

struct recurrence_dct2_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);
};
template<Integer Rows, Integer Y_LD>
struct recurrence_dct2_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);
};

struct recurrence_dct3_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);
};
template<Integer Rows, Integer Y_LD>
struct recurrence_dct3_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);
};


struct make_permute_dct2_inpl_dim1
{
    //work length at least rows0
    static void eval(Real* Y, Integer Y_ld, Integer rows0, Integer cols, Real* work);
};
template<Integer Rows, Integer LD_Y>
struct make_permute_dct2_inpl_dim2
{
    //work length at least cols
    static void eval(Real* Y, Integer Y_ld, Integer rows0, Integer cols, Real* work);
};

struct extract_all_from_pack_with_mult_dim1
{
    //work must be at least of size rows
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, 
                    const Real* cs_table, Real* work);
};

template<Integer Rows, Integer Y_LD>
struct extract_all_from_pack_with_mult_dim2
{
    //work must be at least of size cols
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, 
                    const Real* cs_table, Real* work);
};

struct make_permute_dct2_dim1
{
    static void eval(Real* Y, Integer Y_ld, const Real* X, Integer X_ld, Integer rows0, Integer cols);
};
template<Integer Rows, Integer LD_Y, Integer LD_X>
struct make_permute_dct2_dim2
{
    static void eval(Real* Y, Integer Y_ld, const Real* X, Integer X_ld, Integer rows0, Integer cols);
};

struct make_pack_dct3_inpl_dim1
{
    //work size at least rows
    static void eval(Real* ptr_Y, Integer Y_ld, Integer rows, Integer cols, 
                     const Real* cs, Real* work);
};

template<Integer Rows, Integer LD_Y>
struct make_pack_dct3_inpl_dim2
{
    //work size at least cols
    static void eval(Real* ptr_Y, Integer Y_ld, Integer rows0, Integer cols, const Real* cs,
                     Real* work);
};

struct make_permute_dct3_dim1
{
    //work must be at least of size rows
    static void eval(Real* ptr_Y, Integer Y_ld, Integer rows, Integer cols, Real* work);
};
template<Integer Rows, Integer LD_Y>
struct make_permute_dct3_dim2
{
    //work must be at least of size cols
    static void eval(Real* ptr_Y, Integer Y_ld, Integer rows, Integer cols, Real* work);
};

struct make_pack_dct3_dim1
{
    static void eval(Real* ptr_Y, Integer Y_ld, const Real* ptr_x, Integer x_ld, 
                      Integer rows, Integer cols, const Real* cs);
};

template<Integer Rows, Integer LD_Y, Integer LD_X>
struct make_pack_dct3_dim2
{
    static void eval(Real* ptr_Y, Integer Y_ld, const Real* ptr_x, Integer x_ld, 
                      Integer rows0, Integer cols, const Real* cs);
};

struct mult_by_cos_inpl_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);
};

template<Integer Rows, Integer Y_LD>
struct mult_by_cos_inpl_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs);
};

struct recursion_dct4_dim1
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols);
};

template<Integer Rows, Integer Y_LD>
struct recursion_dct4_dim2
{
    static void eval(Real* Y, Integer Y_ld, Integer rows, Integer cols);
};

struct mult_by_cos_dim1
{
    static void eval(Real* Y, Integer Y_ld, const Real* X, Integer X_ld,
                     Integer rows, Integer cols, const Real* cs);
};

template<Integer Rows, Integer Y_LD, Integer LD_X>
struct mult_by_cos_dim2
{
    static void eval(Real* Y, Integer Y_ld, const Real* X, Integer X_ld,
                     Integer rows, Integer cols, const Real* cs);
};

Real* create_dct_table(Integer length, fft_context::fft_type type);
void  destroy_dct_table(Real* table);

void check_work(Integer supplied, Integer required);

}}
