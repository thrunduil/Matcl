#include "dct_helpers.h"
#include "matcl-blas-lapack/level1/level1.h"

namespace matcl { namespace fft
{

void fft::check_work(Integer supplied, Integer required)
{
    if (supplied < required)
    {
        std::ostringstream msg;
        msg << "insufficient workspace, " << "required size: " << required
            << ", supplied: " << supplied;

        throw error::matcl_fft_exception(msg.str());
    };
};

static Real* make_table_dct23_1(Integer length)
{
    Integer M           = length;
    Real omega          = - matcl::constants::pi() / (2.0 * M);
    Integer N1          = M /2;
    Integer pos         = 0;
    Integer i           = 1;

    Real* table         = new Real[2*N1];

    for (i = 1; i <= N1; ++i, pos += 2)
    {
        table[pos]      = cos(omega * i);
        table[pos+1]    = -sin(omega * i);
    };

    return table;
};
static Real* make_table_dct23_2(Integer length)
{
    Integer M           = length;
    Real omega          = matcl::constants::pi() / (M);
    Integer N1          = M /2;
    Integer pos         = 0;
    Integer i           = 1;

    Real* table         = new Real[2*N1 + 2*N1];

    for (i = 1; i <= N1; ++i, pos += 2)
    {
        table[pos]      = cos(omega * i);
        table[pos+1]    = sin(omega * i);
    };

    Real* table2        = table + 2*N1;
    for (i = 0; i < N1; ++i)
    {        
        table2[i]       = 2.0 * sin(omega * (i + 0.5));
    };
    for (Integer k = 0; i < 2*N1; ++i, ++k)
    {        
        table2[i]       = 1.0 / table2[k];
    };
    return table;
};

static Real* make_table_dct1(Integer length)
{
    Integer M           = length - 1;
    Real omega          = matcl::constants::pi() / M;
    Integer N           = M/2;
    Integer pos         = 0;
    Integer i           = 1;

    Real* table         = new Real[2*N];

    for (i = 0; i < N; ++i, ++pos)
    {
        table[i]        = 2.0 * sin(omega * (i+1));
        table[N+i]      = 2.0 * cos(omega * (i+1));
    };

    return table;
};
Real* create_dct_table(Integer length, fft_context::fft_type type)
{
    switch(type)
    {
        case fft_context::dct_1:
            return make_table_dct1(length);
        case fft_context::dct_23_1:
            return make_table_dct23_1(length);
        case fft_context::dct_23_2:
            return make_table_dct23_2(length);
        default:
            throw std::runtime_error("invalid dct problem type");
    };
};
void fft::destroy_dct_table(Real* table)
{
    delete[] table;
}

void fft::extract_all_from_pack_with_mult_dim1
            ::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs, Real* work)
{
    --cs;

    bool even       = (rows % 2 == 0);
    Integer last    = even  ? rows - 1 : rows;

    for (Integer j = 0; j < cols; ++j)
    {
        Integer i, k;

        for (i = 1, k = 1; i < last; i += 2, ++k)
        {
            Real re = Y[i] * cs[i]   + Y[i+1] * cs[i+1];
            Real im = Y[i] * cs[i+1] - Y[i+1] * cs[i];

            Y[k]    = re;
            work[k] = im;
        }

        Integer pos = k - 1;

        if (even == true)
        {
            Real re = Y[i] * cs[i];
            Y[k]    = re;
            ++k;
        };            

        for (; k < rows; ++k, --pos)
        {
            Y[k]    = work[pos];
        };

        Y           += Y_ld;
    };

    return;
}

template<Integer Rows, Integer Y_LD>
void fft::extract_all_from_pack_with_mult_dim2<Rows, Y_LD>
            ::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols, const Real* cs, Real* work)
{
    //row by row version. Column by column version is possible in two stages:
    //first create matrix with required values but permuted, then permute columns.
    //This may be faster in case of matrices with many rows.

    Integer Y_ld    = get_int<Y_LD>::eval(Y_ld0);

    --cs;

    bool even       = (cols % 2 == 0);
    Integer last    = even  ? cols - 1 : cols;

    for (Integer j = 0; j < get_int<Rows>::eval(rows0); ++j)
    {
        Integer i, k;

        Real* Y_1   = Y + 1 * Y_ld;
        Real* Y_2   = Y + 1 * Y_ld;
        Real* Y_3   = Y + 2 * Y_ld;

        for (i = 1, k = 1; i < last; i += 2, ++k)
        {
            Real re = Y_2[j] * cs[i]   + Y_3[j] * cs[i+1];
            Real im = Y_2[j] * cs[i+1] - Y_3[j] * cs[i];

            Y_1[j]  = re;
            work[k] = im;

            Y_1     += 1 * Y_ld;
            Y_2     += 2 * Y_ld;
            Y_3     += 2 * Y_ld;
        }

        Integer pos = k - 1;

        if (even == true)
        {
            Real re = Y_2[j] * cs[i];
            Y_1[j]  = re;
            Y_1     += Y_ld;
            ++k;
        };            

        for (; k < cols; ++k, --pos)
        {
            Y_1[j]  = work[pos];
            Y_1     += 1 * Y_ld;
        };
    };

    return;
}

void fft::make_permute_dct2_dim1::eval(Real* Y, Integer Y_ld, const Real* X, Integer X_ld, 
                                        Integer rows, Integer cols)
{
    bool even       = (rows % 2 == 0);
    Integer last    = even  ? rows - 1 : rows;

    for (Integer i = 0; i < cols; ++i)
    {
        Y[0]        = X[0];

        Integer ps  = 1;
        Integer pe  = rows - 1;

        for (Integer j = 2; j <= last; j += 2, ++ps)
        {
            Y[ps]   = X[j];
        };
        for (Integer j = 1; j < last; j += 2, --pe)
        {
            Y[pe]   = X[j];
        };

        if (even == true)
            Y[pe]   = X[rows - 1];

        Y           += Y_ld;
        X           += X_ld;
    };
};

template<Integer Rows, Integer LD_Y, Integer LD_X>
void fft::make_permute_dct2_dim2<Rows, LD_Y, LD_X>::eval(Real* Y, Integer Y_ld0, const Real* X, Integer X_ld0, 
                                        Integer rows0, Integer cols)
{
    Integer Y_ld    = get_int<LD_Y>::eval(Y_ld0);
    Integer X_ld    = get_int<LD_X>::eval(X_ld0);

    bool even       = (cols % 2 == 0);
    Integer last    = even  ? cols - 1 : cols;

    Real* Y2        = Y + (cols - 1) * Y_ld;

    matcl::level1::copy_mat<false,Real,Real,Rows,1,1>::eval(Y, Y_ld, X, X_ld, rows0, cols);
    Y               += Y_ld;
    X               += X_ld;        

    for (Integer i = 1; i < last; i += 2)
    {
        matcl::level1::copy_mat<false,Real,Real,Rows,1,1>::eval(Y2, Y_ld, X, X_ld, rows0, cols);
        Y2          -= Y_ld;
        X           += X_ld;

        matcl::level1::copy_mat<false,Real,Real,Rows,1,1>::eval(Y, Y_ld, X, X_ld, rows0, cols);
        Y           += Y_ld;
        X           += X_ld;
    };

    if (even == true)
    {
        matcl::level1::copy_mat<false,Real,Real,Rows,1,1>::eval(Y, Y_ld, X, X_ld, rows0, cols);
    };
};

void fft::make_permute_dct2_inpl_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, Real* work)
{
    if (rows <= 2)
        return;

    bool even           = (rows % 2 == 0);
    Integer ncopy       = (rows - 2)/3 + 1 + (rows % 3 == 0? 1 : 0);

    for (Integer i = 0; i < cols; ++i)
    {
        for (Integer k = rows - ncopy; k < rows; ++k)
        {
            work[k]     = Y[k];
        };

        Integer j; 
        Integer ps      = 1;
        Integer pe      = rows - 1;        

        for (j = 1; j + 1 < pe; j += 2, pe -= 1, ps += 1)
        {
            Y[pe]       = Y[j];
            Y[ps]       = Y[j+1];
        };
        
        for (; j < rows - 1; j += 2, --pe, ++ps)
        {
            Y[pe]       = work[j];
            Y[ps]       = work[j+1];
        };

        if (even == true)
            Y[pe]       = work[j];

        Y               += Y_ld;
    };
};
template<Integer Rows, Integer LD_Y>
void fft::make_permute_dct2_inpl_dim2<Rows, LD_Y>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols,
                                                        Real* work)
{
    if (cols <= 2)
        return;

    Integer rows        = get_int<Rows>::eval(rows0);
    Integer Y_ld        = get_int<LD_Y>::eval(Y_ld0);

    bool even           = (cols % 2 == 0);
    Integer ncopy       = (cols - 2)/3 + 1 + (cols % 3 == 0? 1 : 0);

    for (Integer i = 0; i < rows; ++i)
    {
        for (Integer k = cols - ncopy; k < cols; ++k)
        {
            work[k]     = Y[i+k*Y_ld];
        };

        Integer j; 
        Integer ps      = 1;
        Integer pe      = cols - 1;        
     
        for (j = 1; j + 1 < pe; j += 2, pe -= 1, ps += 1)
        {
            Y[i+pe*Y_ld]    = Y[i+j*Y_ld];
            Y[i+ps*Y_ld]    = Y[i+(j+1)*Y_ld];
        };

        for (; j < cols - 1; j += 2, --pe, ++ps)
        {
            Y[i+pe*Y_ld]    = work[j];
            Y[i+ps*Y_ld]    = work[j+1];            
        };

        if (even == true)
            Y[i+pe*Y_ld]    = work[j];
    };
};

void make_pack_dct3_dim1::eval(Real* ptr_Y, Integer Y_ld, const Real* ptr_x, Integer x_ld, 
                      Integer rows, Integer cols, const Real* cs)
{
    --cs;

    bool even           = (rows % 2 == 0);
    Integer last        = even  ? rows - 1 : rows;

    Integer k1          = rows - 1;
    for (Integer i = 0; i < cols; ++i)
    {
        ptr_Y[0]        = ptr_x[0];

        Integer k       = k1;
        Integer m       = 1;
        Integer j;

        for (j = 1; m < last; ++j, --k, m += 2)
        {
            //Real val  = (ptr_x[j] - I * ptr_x[k]) * (cs[m] + I * cs[m+1]);
            Real re     = ptr_x[j] * cs[m]   + ptr_x[k] * cs[m+1];
            Real im     = ptr_x[j] * cs[m+1] - ptr_x[k] * cs[m];

            ptr_Y[m]    = re;
            ptr_Y[m+1]  = im;
        };

        if (even == true)
        {
            Real re     = 2.0 * ptr_x[j] * cs[m];
            ptr_Y[m]    = re;
        }

        ptr_Y           += Y_ld;
        ptr_x           += x_ld;
    };
};

template<Integer Rows,Integer LD_Y, Integer LD_X>
void make_pack_dct3_dim2<Rows,LD_Y,LD_X>::eval(Real* Y, Integer Y_ld0, const Real* X, Integer X_ld0, 
                      Integer rows0, Integer cols, const Real* cs)
{
    --cs;

    Integer Y_ld        = get_int<LD_Y>::eval(Y_ld0);
    Integer X_ld        = get_int<LD_X>::eval(X_ld0);

    bool even           = (cols % 2 == 0);
    Integer last        = cols  ? cols - 1 : cols;

    matcl::level1::copy_mat<false, Real, Real, Rows, 1, 1>
            ::eval(Y, Y_ld, X, X_ld, rows0, 1);    

    Real* restricted  Y_1       = Y + Y_ld;
    Real* restricted  Y_2       = Y + 2*Y_ld;
    const Real* restricted  X_1 = X + X_ld;
    const Real* restricted  X_2 = X + (cols - 1) * X_ld;

    Integer m           = 1;

    for (; m < last; m += 2)
    {
        Real cs1        = cs[m];
        Real cs2        = cs[m+1];

        //This should be vectorized by compiler (checked)
        for (Integer j = 0; j < get_int<Rows>::eval(rows0); ++j)
        {            
            Real re     = X_1[j] * cs1 + X_2[j] * cs2;
            Real im     = X_1[j] * cs2 - X_2[j] * cs1;

            Y_1[j]      = re;
            Y_2[j]      = im;
        };

        X_1             += X_ld;
        X_2             -= X_ld;
        Y_1             += 2*Y_ld;
        Y_2             += 2*Y_ld;
    };

    if (even == true)
    {
        Real scal       = 2.0 * cs[m];
        matcl::level1::ax<Real, Real, Real, Rows>::eval(Y_1, X_1, rows0, scal);
    };
};

void make_pack_dct3_inpl_dim1::eval(Real* ptr_Y, Integer Y_ld, Integer rows, Integer cols, 
                                       const Real* cs, Real* work)
{
    --cs;

    if (rows == 2)
    {
        for (Integer i = 0; i < cols; ++i)
        {
            Real re     = 2.0 * ptr_Y[1] * cs[1];
            ptr_Y[1]    = re;

            ptr_Y       += Y_ld;
        };
        return;
    };

    bool even           = (rows % 2 == 0);
    Integer last        = even  ? rows - 1 : rows;

    Integer k1          = rows - 1;
    for (Integer i = 0; i < cols; ++i)
    {
        Integer k       = k1;
        Integer m       = 1;
        Integer j;

        work[1]         = ptr_Y[1];

        for (j = 1; m <= k; ++j, --k, m += 2)
        {
            Real re     = work[j] * cs[m]   + ptr_Y[k] * cs[m+1];
            Real im     = work[j] * cs[m+1] - ptr_Y[k] * cs[m];

            work[m]     = ptr_Y[m];
            work[m+1]   = ptr_Y[m+1];

            ptr_Y[m]    = re;
            ptr_Y[m+1]  = im;
        };
        for (; m < last; ++j, --k, m += 2)
        {
            Real re     = work[j] * cs[m]   + work[k] * cs[m+1];
            Real im     = work[j] * cs[m+1] - work[k] * cs[m];

            ptr_Y[m]    = re;
            ptr_Y[m+1]  = im;
        };
      
        if (even == true)
        {
            Real re     = 2.0 * work[j] * cs[m];
            ptr_Y[m]    = re;
        }

        ptr_Y           += Y_ld;
    };
};

template<Integer Rows,Integer LD_Y>
void make_pack_dct3_inpl_dim2<Rows,LD_Y>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols, 
                                                  const Real* cs, Real* work)
{
    --cs;

    Integer rows        = get_int<Rows>::eval(rows0);
    Integer Y_ld        = get_int<LD_Y>::eval(Y_ld0);

    if (cols == 2)
    {
        Y               += Y_ld;
        Real cs1        = 2.0 * cs[1];

        matcl::level1::ay<Real, Real, Rows>::eval(Y, rows0, cs1);

        return;
    };

    bool even           = (cols % 2 == 0);
    Integer last        = even  ? cols - 1 : cols;

    Integer k1          = cols - 1;
    for (Integer i = 0; i < rows; ++i)
    {
        Integer k       = k1;
        Integer m       = 1;
        Integer j;

        work[1]         = Y[i+1*Y_ld];

        for (j = 1; m <= k; ++j, --k, m += 2)
        {
            Real re     = work[j] * cs[m]   + Y[i+k*Y_ld] * cs[m+1];
            Real im     = work[j] * cs[m+1] - Y[i+k*Y_ld] * cs[m];

            work[m]     = Y[i+m*Y_ld];
            work[m+1]   = Y[i+(m+1)*Y_ld];

            Y[i+m*Y_ld]     = re;
            Y[i+(m+1)*Y_ld] = im;
        };

        for (; m < last; ++j, --k, m += 2)
        {
            //Real val  = (ptr_x[j] - I * ptr_x[k]) * (cs[m] + I * cs[m+1]);
            Real re     = work[j] * cs[m]   + work[k] * cs[m+1];
            Real im     = work[j] * cs[m+1] - work[k] * cs[m];

            Y[i+m*Y_ld]     = re;
            Y[i+(m+1)*Y_ld] = im;
        };
      
        if (even == true)
        {
            Real re     = 2.0 * work[j] * cs[m];
            Y[i+m*Y_ld] = re;
        }
    };
};

void make_permute_dct3_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, Real* work)
{
    if (rows <= 2)
        return;

    bool even       = (rows % 2 == 0);
    Integer last    = even  ? rows - 1 : rows;

    for (Integer i = 0; i < cols; ++i)
    {
        Integer j; 

        Integer pe      = rows - 1;
        Integer ps      = 1;

        for (j = 1; j <= pe; j += 2, pe -= 1, ps += 1)
        {
            work[j]     = Y[j];
            work[j+1]   = Y[j+1];

            Y[j]        = Y[pe];
            Y[j+1]      = work[ps];            
        };
        for (; j < last; j += 2, --pe, ++ps)
        {
            Y[j]        = work[pe];
            Y[j+1]      = work[ps];            
        };

        if (even == true)
            Y[j]        = work[pe];

        Y               += Y_ld;
    };
};

template<Integer Rows, Integer Y_LD>
void make_permute_dct3_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols,
                                             Real* work)
{
    if (cols <= 2)
        return;

    //row by row version
    bool even       = (cols % 2 == 0);
    Integer last    = even  ? cols - 1 : cols;

    Integer rows    = get_int<Rows>::eval(rows0);
    Integer Y_ld    = get_int<Y_LD>::eval(Y_ld0);

    for (Integer i = 0; i < rows; ++i)
    {
        Integer j; 

        Integer pe      = cols - 1;
        Integer ps      = 1;

        for (j = 1; j <= pe; j += 2, pe -= 1, ps += 1)
        {
            work[j]     = Y[i+(j) * Y_ld];
            work[j+1]   = Y[i+(j+1) * Y_ld];

            Y[i+(j)*Y_ld]   = Y[i + pe * Y_ld];
            Y[i+(j+1)*Y_ld] = work[ps];            
        };

        for (; j < last; j += 2, --pe, ++ps)
        {
            Y[i+j*Y_ld]     = work[pe];
            Y[i+(j+1)*Y_ld] = work[ps];            
        };

        if (even == true)
            Y[i+j*Y_ld]     = work[pe];
    };
};

void mult_by_cos_dim1::eval(Real* Y, Integer Y_ld, const Real* X, Integer X_ld,
                    Integer rows, Integer cols, const Real* cs)
{
    for (Integer i = 0; i < cols; ++i)
    {
        for (Integer j = 0; j < rows; ++j)
        {
            Y[j]            = X[j] * cs[j];
        };

        Y                   += Y_ld;
        X                   += X_ld;
    };
};

template<Integer Rows, Integer Y_LD, Integer X_LD>
void mult_by_cos_dim2<Rows,Y_LD,X_LD>::eval(Real* Y, Integer Y_ld0, const Real* X, Integer X_ld0,
                     Integer rows0, Integer cols, const Real* cs)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);
    Integer X_ld            = get_int<X_LD>::eval(X_ld0);

    (void)rows;

    for (Integer i = 0; i < cols; ++i)
    {
        matcl::level1::ax<Real, Real, Real, Rows>::eval(Y, X, rows0, cs[i]);

        Y                   += Y_ld;
        X                   += X_ld;
    };
};

void mult_by_cos_inpl_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs)
{
    for (Integer i = 0; i < cols; ++i)
    {
        for (Integer j = 0; j < rows; ++j)
        {
            Y[j]            = Y[j] * cs[j];
        };

        Y                   += Y_ld;
    };
};

template<Integer Rows, Integer Y_LD>
void mult_by_cos_inpl_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols, 
                                       const Real* cs)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    (void)rows;

    for (Integer i = 0; i < cols; ++i)
    {
        matcl::level1::ay<Real, Real, Rows>::eval(Y, rows0, cs[i]);

        Y                   += Y_ld;
    };
};


void recursion_dct4_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols)
{
    for (Integer i = 0; i < cols; ++i)
    {
        Y[0]                = Y[0] * 0.5;
        for (Integer j = 1; j < rows; ++j)
        {
            Y[j]            = Y[j] - Y[j-1];
        };

        Y                   += Y_ld;
    }
};

template<Integer Rows, Integer Y_LD>
void recursion_dct4_dim2<Rows, Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    (void)rows;

    matcl::level1::ay<Real, Real, Rows>::eval(Y, rows0, 0.5);

    Real* Y_prev            = Y;
    Y                       += Y_ld;

    for (Integer i = 1; i < cols; ++i)
    {
        matcl::level1::ymx<Real, Real,Rows>::eval(Y, Y_prev, rows0);

        Y                   += Y_ld;
        Y_prev              += Y_ld;
    }
};

void prepare_dct1_dim1::eval(Real* Y, Integer Y_ld, const Real* X, 
                             Integer x_ld, Integer rows, Integer cols, const Real* cs)
{
    Integer N               = (rows-1)/2;
    const Real* st          = cs;
    const Real* ct          = cs + N;

    --st;
    --ct;

    bool even               = (rows % 2 == 0);
    Integer last_row        = rows/2;

    for (Integer i = 0; i < cols; ++i)
    {
        Real sum            = X[0] - X[rows - 1];
        Y[0]                = X[0] + X[rows - 1];

        Integer k           = rows - 2;
        for (Integer j = 1; j < last_row; ++j, --k)
        {
            Real tmp1       = X[j] + X[k];
            Real tmp2       = X[j] - X[k];
            Real tmp3       = tmp2 * st[j];
            Real tmp4       = tmp2 * ct[j];

            Y[j]            = tmp1 - tmp3;
            Y[k]            = tmp1 + tmp3;

            sum             += tmp4;
        };

        if (even == false)
        {
            Y[k]            = 2.0 * X[k];
        };

        Y[rows-1]           = sum;

        Y                   += Y_ld;
        X                   += x_ld;
    };    
};
template<Integer Rows, Integer Y_LD, Integer X_LD>
void prepare_dct1_dim2<Rows,Y_LD,X_LD>::eval(Real* Y, Integer Y_ld0, const Real* X, 
                             Integer X_ld0, Integer rows0, Integer cols, const Real* cs)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);
    Integer X_ld            = get_int<X_LD>::eval(X_ld0);

    Integer N               = (cols-1)/2;
    const Real* st          = cs;
    const Real* ct          = cs + N;

    --st;
    --ct;

    bool even               = (cols % 2 == 0);
    Integer last_col        = cols/2;

    const Real* restricted  X_last  = X + (cols -  1) * X_ld;
    Real* restricted  Y_last        = Y + (cols -  1) * Y_ld;
    Real* Y_sum                     = Y_last;

    //This should be vectorized by compiler (checked)
    for (Integer i = 0; i < rows; ++i)
    {
        Real sum            = X[i] - X_last[i];
        Y[i]                = X[i] + X_last[i];
        Y_sum[i]            = sum;
    };

    X                       += X_ld;
    Y                       += Y_ld;
    X_last                  -= X_ld;
    Y_last                  -= Y_ld;

    Integer k               = cols - 2;

    for (Integer j = 1; j < last_col; ++j, --k)
    {
        Real stj            = st[j];
        Real ctj            = ct[j];

        //This should be vectorized by compiler (checked)
        for (Integer i = 0; i < rows; ++i)
        {
            Real tmp1       = X[i] + X_last[i];
            Real tmp2       = X[i] - X_last[i];
            Real tmp3       = tmp2 * stj;

            Y_sum[i]        += tmp2 * ctj;

            Y[i]            = tmp1 - tmp3;
            Y_last[i]       = tmp1 + tmp3;
        };

        X                   += X_ld;
        Y                   += Y_ld;
        X_last              -= X_ld;
        Y_last              -= Y_ld;
    };

    if (even == false)
    {
        matcl::level1::ax<Real,Real, Real, Rows>::eval(Y, X, rows, 2.0);
    };
};
template<Integer Rows, Integer Y_LD>
void prepare_dct1_inpl_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, 
                                              Integer cols, const Real* cs)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    Integer N               = (cols-1)/2;
    bool even               = (cols % 2 == 0);
    const Real* st          = cs;
    const Real* ct          = cs + N;

    --st;
    --ct;

    Integer last_col        = cols/2;

    Real* restricted Y_last = Y + (cols -  1) * Y_ld;
    Real* restricted Y_sum  = Y_last;

    //This should be vectorized by compiler (checked)
    for (Integer i = 0; i < rows; ++i)
    {
        Real sum            = Y[i] - Y_last[i];
        Y[i]                = Y[i] + Y_last[i];
        Y_sum[i]            = sum;
    };

    Y                       += Y_ld;
    Y_last                  -= Y_ld;

    Integer k               = cols - 2;

    for (Integer j = 1; j < last_col; ++j, --k)
    {
        Real stj            = st[j];
        Real ctj            = ct[j];

        //This should be vectorized by compiler (checked)
        for (Integer i = 0; i < rows; ++i)
        {
            Real tmp1       = Y[i] + Y_last[i];
            Real tmp2       = Y[i] - Y_last[i];
            Real tmp3       = tmp2 * stj;

            Y_sum[i]        += tmp2 * ctj;

            Y[i]            = tmp1 - tmp3;
            Y_last[i]       = tmp1 + tmp3;
        };

        Y                   += Y_ld;
        Y_last              -= Y_ld;
    };

    if (even == false)
    {
        matcl::level1::ay<Real,Real, Rows>::eval(Y, rows, 2.0);
    };
};

void prepare_dct1_inpl_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs)
{
    Integer N               = (rows-1)/2;
    bool even               = (rows % 2 == 0);
    const Real* st          = cs;
    const Real* ct          = cs + N;

    --st;
    --ct;

    Integer last_row        = rows/2;

    for (Integer i = 0; i < cols; ++i)
    {
        Real sum            = Y[0] - Y[rows - 1];
        Y[0]                = Y[0] + Y[rows - 1];

        Integer k           = rows - 2;
        for (Integer j = 1; j < last_row; ++j, --k)
        {
            Real tmp1       = Y[j] + Y[k];
            Real tmp2       = Y[j] - Y[k];
            Real tmp3       = tmp2 * st[j];

            sum             += tmp2 * ct[j];

            Y[j]            = tmp1 - tmp3;
            Y[k]            = tmp1 + tmp3;
        };

        if (even == false)
        {
            Y[k]            = 2.0 * Y[k];
        };

        Y[rows-1]           = sum;

        Y                   += Y_ld;
    };    
};
void prepare_dct2_dim1::eval(Real* Y, Integer Y_ld, const Real* X, 
                             Integer x_ld, Integer rows, Integer cols, const Real* cs)
{
    bool even               = (rows % 2 == 0);
    Integer last_row        = rows/2;

    cs                      += (rows/2)*2;

    for (Integer i = 0; i < cols; ++i)
    {
        Integer k           = rows - 1;
        for (Integer j = 0; j < last_row; ++j, --k)
        {
            Real sj         = cs[j];
            Real tmp1       = X[j] + X[k];
            Real tmp2       = X[j] - X[k];
            Real tmp3       = tmp2 * sj;

            Y[j]            = tmp1 + tmp3;
            Y[k]            = tmp1 - tmp3;
        };

        if (even == false)
        {
            Y[k]            = 2.0 * X[k];
        };

        Y                   += Y_ld;
        X                   += x_ld;
    };    
};
void prepare_dct2_inpl_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs)
{
    cs                      += (rows/2)*2;
    Integer last_row        = rows/2;
    bool even               = (rows % 2 == 0);

    for (Integer i = 0; i < cols; ++i)
    {
        Integer k           = rows - 1;
        for (Integer j = 0; j < last_row; ++j, --k)
        {
            Real tmp1       = Y[j] + Y[k];
            Real tmp2       = Y[j] - Y[k];
            Real tmp3       = tmp2 * cs[j];

            Y[j]            = tmp1 + tmp3;
            Y[k]            = tmp1 - tmp3;
        };

        if (even == false)
            Y[k]            = 2.0 * Y[k];

        Y                   += Y_ld;
    };    
};

template<Integer Rows, Integer Y_LD, Integer X_LD>
void prepare_dct2_dim2<Rows,Y_LD,X_LD>::eval(Real* Y, Integer Y_ld0, const Real* X, 
                             Integer x_ld0, Integer rows0, Integer cols, const Real* cs)
{
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);
    Integer X_ld            = get_int<X_LD>::eval(x_ld0);
    Integer rows            = get_int<Rows>::eval(rows0);

    cs                      += (cols/2)*2;

    bool even               = (cols % 2 == 0);
    Integer last_col        = cols/2;

    const Real* restricted X_j  = X;
    const Real* restricted X_k  = X + (cols - 1) * X_ld;
    Real* restricted Y_j        = Y;
    Real* restricted Y_k        = Y + (cols - 1) * Y_ld;

    for (Integer j = 0; j < last_col; ++j)
    {
        //This should be vectorized by compiler (checked)
        for (Integer i = 0; i < rows; ++i)
        {
            Real tmp1       = X_j[i] + X_k[i];
            Real tmp2       = X_j[i] - X_k[i];
            Real tmp3       = tmp2 * cs[j];

            Y_j[i]          = tmp1 + tmp3;
            Y_k[i]          = tmp1 - tmp3;
        };

        X_j                 += X_ld;
        Y_j                 += Y_ld;
        X_k                 -= X_ld;
        Y_k                 -= Y_ld;
    };

    if (even == false)
    {
        matcl::level1::ax<Real,Real, Real, Rows>::eval(Y_k, X_k, rows0, 2.0);
    };
};
template<Integer Rows, Integer Y_LD>
void prepare_dct2_inpl_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols,
                                        const Real* cs)
{
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);
    Integer rows            = get_int<Rows>::eval(rows0);

    cs                      += (cols/2)*2;

    bool even               = (cols % 2 == 0);
    Integer last_col        = cols/2;

    Real* restricted Y_j    = Y;
    Real* restricted Y_k    = Y + (cols - 1) * Y_ld;

    for (Integer j = 0; j < last_col; ++j)
    {
        //This should be vectorized by compiler (checked)
        for (Integer i = 0; i < rows; ++i)
        {
            Real tmp1       = Y_j[i] + Y_k[i];
            Real tmp2       = Y_j[i] - Y_k[i];
            Real tmp3       = tmp2 * cs[j];

            Y_j[i]          = tmp1 + tmp3;
            Y_k[i]          = tmp1 - tmp3;
        };

        Y_j                 += Y_ld;
        Y_k                 -= Y_ld;
    };

    if (even == false)
    {
        matcl::level1::ay<Real, Real, Rows>::eval(Y_k, rows0, 2.0);
    };

};

void prepare_dct3_dim1::eval(Real* Y, Integer Y_ld, const Real* X,
                     Integer X_ld, Integer rows, Integer cols, const Real* cs)
{
    Integer N               = rows;

    bool even               = (rows % 2 == 0);
    Integer last            = even ? N / 2 : N / 2 + 1;

    --cs;

    for (Integer j = 0; j < cols; ++j)
    {
        Y[0]                = X[0];
        if (even == true)
        {            
            for (Integer k = 1; k <= N-3; k += 2)
            {
                Y[k]        = X[k] - X[k+2];
                Y[k+1]      = X[k+1];
            };
            Y[rows-1]       = X[rows-1] * 2.0;
        }
        else
        {
            Integer k       = 1;
            for (; k <= N-4; k += 2)
            {
                Y[k]        = X[k] - X[k+2];
                Y[k+1]      = X[k+1];
            };

            for(; k < rows; ++k)
                Y[k]        = X[k];
        };

        Integer pos         = 1;

        for (Integer i = 1; i < last; ++i, pos += 2)
        {
            Real re         = Y[pos];
            Real im         = Y[pos + 1];

            Real ci         = cs[pos];
            Real si         = cs[pos+1];

            Real re1        = ci* im;
            Real re2        = si* im;
            Real im1        = si * re;
            Real im2        = ci* re;

            Y[pos]          = re1 + im1;
            Y[pos+1]        = re2 - im2;
        };

        Y                   += Y_ld;
        X                   += X_ld;
    };
};
void prepare_dct3_inpl_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs)
{
    Integer N               = rows;

    bool even               = (rows % 2 == 0);
    Integer last            = even ? N / 2 : N / 2 + 1;

    --cs;

    for (Integer j = 0; j < cols; ++j)
    {
        if (even == true)
        {            
            for (Integer k = 1; k <= N-3; k += 2)
            {
                Y[k]        = Y[k] - Y[k+2];
            };
            Y[rows-1]       = Y[rows-1] * 2.0;
        }
        else
        {
            Integer k       = 1;
            for (; k <= N-4; k += 2)
            {
                Y[k]        = Y[k] - Y[k+2];
            };
        };

        Integer pos         = 1;

        for (Integer i = 1; i < last; ++i, pos += 2)
        {
            Real re         = Y[pos];
            Real im         = Y[pos + 1];

            Real ci         = cs[pos];
            Real si         = cs[pos+1];

            Real re1        = ci* im;
            Real re2        = si* im;
            Real im1        = si * re;
            Real im2        = ci* re;

            Y[pos]          = re1 + im1;
            Y[pos+1]        = re2 - im2;
        };

        Y                   += Y_ld;
    };

};

template<Integer Rows, Integer Y_LD, Integer X_LD>
void prepare_dct3_dim2<Rows,Y_LD,X_LD>::eval(Real* Y, Integer Y_ld0, const Real* X,
                    Integer X_ld0, Integer rows0, Integer cols, const Real* cs)
{
    Integer N               = cols;
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);
    Integer X_ld            = get_int<X_LD>::eval(X_ld0);

    bool even               = (cols % 2 == 0);
    Integer last            = even ? N / 2 : N / 2 + 1;

    --cs;

    matcl::level1::copy_mat<false,Real,Real,Rows,1,1>::eval(Y, Y_ld, X, X_ld, rows0, cols);

    X                       += X_ld;
    Y                       += Y_ld;

    if (even == true)
    {   
        Real* restricted Y_k            = Y;
        const Real* restricted X_k      = X;
        const Real* restricted X_kp2    = X + 2 * X_ld;

        for (Integer k = 1; k <= N-3; k += 2)
        {
            //This should be vectorized by compiler (checked)
            for (Integer i = 0; i < rows; ++i)
            {
                Y_k[i]      = X_k[i] - X_kp2[i];
            };

            X_k             += X_ld;
            Y_k             += Y_ld;

            matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
                ::eval(Y_k, Y_ld, X_k, X_ld, rows0, cols);
            X_k             += X_ld;
            Y_k             += Y_ld;

            X_kp2           += 2*X_ld;
        };

        matcl::level1::ax<Real,Real,Real, Rows>::eval(Y_k, X_k, rows0, 2.0);
    }
    else
    {
        Real* restricted Y_k            = Y;
        const Real* restricted X_k      = X;
        const Real* restricted X_kp2    = X + 2 * X_ld;

        Integer k       = 1;
        for (; k <= N-4; k += 2)
        {
            //This should be vectorized by compiler (checked)
            for (Integer i = 0; i < rows; ++i)
            {
                Y_k[i]      = X_k[i] - X_kp2[i];
            };

            X_k             += X_ld;
            Y_k             += Y_ld;

            matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
                ::eval(Y_k, Y_ld, X_k, X_ld, rows0, cols);
            X_k             += X_ld;
            Y_k             += Y_ld;

            X_kp2           += 2*X_ld;
        };

        matcl::level1::copy_mat<false,Real,Real,Rows,2,0>
            ::eval(Y_k, Y_ld, X_k, X_ld, rows0, cols);
    };

    Real* restricted Y_k    = Y;
    Real* restricted Y_kp1  = Y + Y_ld;

    for (Integer i = 1, k = 1; i < last; ++i, k += 2)
    {
        Real ci             = cs[k];
        Real si             = cs[k+1];

        //This should be vectorized by compiler (checked)
        for (Integer j = 0; j < rows; ++j)
        {
            Real re         = Y_k[j];
            Real im         = Y_kp1[j];

            Real re1        = ci* im;
            Real re2        = si* im;
            Real im1        = si * re;
            Real im2        = ci* re;

            Y_k[j]          = re1 + im1;
            Y_kp1[j]        = re2 - im2;
        };

        Y_k                 += 2*Y_ld;
        Y_kp1               += 2*Y_ld;
    };
};
template<Integer Rows, Integer Y_LD>
void prepare_dct3_inpl_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols,
                                             const Real* cs)
{
    Integer N               = cols;
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    bool even               = (cols % 2 == 0);
    Integer last            = even ? N / 2 : N / 2 + 1;

    --cs;

    Y                       += Y_ld;

    if (even == true)
    {   
        Real* restricted Y_k            = Y;
        const Real* restricted Y_kp2    = Y + 2 * Y_ld;

        for (Integer k = 1; k <= N-3; k += 2)
        {
            //This should be vectorized by compiler (checked)
            for (Integer i = 0; i < rows; ++i)
            {
                Y_k[i]      = Y_k[i] - Y_kp2[i];
            };

            Y_k             += 2 * Y_ld;
            Y_kp2           += 2 * Y_ld;
        };

        matcl::level1::ay<Real, Real, Rows>::eval(Y_k, rows0, 2.0);
    }
    else
    {
        Real* restricted Y_k            = Y;
        const Real* restricted Y_kp2    = Y + 2 * Y_ld;

        Integer k       = 1;
        for (; k <= N-4; k += 2)
        {
            //This should be vectorized by compiler (checked)
            for (Integer i = 0; i < rows; ++i)
            {
                Y_k[i]      = Y_k[i] - Y_kp2[i];
            };

            Y_k             += 2*Y_ld;
            Y_kp2           += 2*Y_ld;
        };
    };

    Real* restricted Y_k    = Y;
    Real* restricted Y_kp1  = Y + Y_ld;

    for (Integer i = 1, k = 1; i < last; ++i, k += 2)
    {
        Real ci             = cs[k];
        Real si             = cs[k+1];

        //This should be vectorized by compiler (checked)
        for (Integer j = 0; j < rows; ++j)
        {
            Real re         = Y_k[j];
            Real im         = Y_kp1[j];

            Real re1        = ci* im;
            Real re2        = si* im;
            Real im1        = si * re;
            Real im2        = ci* re;

            Y_k[j]          = re1 + im1;
            Y_kp1[j]        = re2 - im2;
        };

        Y_k                 += 2*Y_ld;
        Y_kp1               += 2*Y_ld;
    };
};

void recurrence_dct1_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, Real mult)
{
    for (Integer i = 0; i < cols; ++i)
    {
        Real sum            = Y[rows - 1] * mult;

        for (Integer k = rows-1; k >= 2; --k)
        {
            Y[k]            = Y[k-1];
        };

        Y[1]                = sum;

        for (Integer k = 3; k < rows; k += 2)
        {
            Y[k]            = Y[k-2] - Y[k];
        };

        Y                   += Y_ld;
    };
};
template<Integer Y_LD>
struct recurrence_dct1_dim2<0,Y_LD>
{
    static const Integer Rows   = 0;
    static const Integer nb     = 32;

    static void eval(Real* Y, Integer Y_ld0, Integer rows, Integer cols, Real scal)
    {
        Integer n_iter          = rows/nb;

        for (Integer i = 0; i < n_iter; ++i)
        {
            recurrence_dct1_dim2<nb,Y_LD>::eval(Y, Y_ld0, nb, cols, scal);
            Y                   += nb;
        };

        Integer rem             = rows - n_iter * nb;

        if (rem > 0)
            eval_tail(Y,Y_ld0, rem, cols, scal);
    };
    static void eval_tail(Real* Y, Integer Y_ld0, Integer rows0, Integer cols, Real scal)
    {       
        Integer rows            = get_int<Rows>::eval(rows0);
        Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

        Real* Y_last            = Y + (cols - 1) * Y_ld;    

        Real buf[nb];
        Real* Y_sum             = buf;

        if (scal == 1.0)
        {
            matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
                ::eval(Y_sum, nb, Y_last, Y_ld, rows, cols);
        }
        else
        {
            matcl::level1::ax<Real,Real,Real, Rows>::eval(Y_sum, Y_last, rows, scal);
        };

        for (Integer k = cols - 1; k >= 2; --k)
        {
            matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
                    ::eval(Y + k * Y_ld, 1, Y + (k-1)*Y_ld, Y_ld, rows, cols);
        };

        matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
                ::eval(Y + 1 * Y_ld, 1, Y_sum, nb, rows, cols);

        for (Integer k = 3; k < cols; k += 2)
        {
            matcl::level1::xmy<Real,Real,Rows>
                    ::eval(Y + k * Y_ld, Y + (k-2)*Y_ld, rows);
        };
    };
};
template<Integer Rows, Integer Y_LD>
void recurrence_dct1_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols, Real scal)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    Real* Y_last            = Y + (cols - 1) * Y_ld;    

    Real buf[Rows];

    Real* Y_sum             = buf;

    if (scal == 1.0)
    {
        matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
            ::eval(Y_sum, Rows, Y_last, Y_ld, rows, cols);
    }
    else
    {
        matcl::level1::ax<Real,Real,Real, Rows>
            ::eval(Y_sum, Y_last, rows, scal);
    };

    for (Integer k = cols - 1; k >= 2; --k)
    {
        matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
                ::eval(Y + k * Y_ld, 1, Y + (k-1)*Y_ld, Y_ld, rows, cols);
    };

    matcl::level1::copy_mat<false,Real,Real,Rows,1,1>
            ::eval(Y + 1 * Y_ld, 1, Y_sum, Rows, rows, cols);

    for (Integer k = 3; k < cols; k += 2)
    {
        matcl::level1::xmy<Real,Real,Rows>
                ::eval(Y + k * Y_ld, Y + (k-2)*Y_ld, rows);
    };
};

void recurrence_dct2_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs)
{
    Integer N               = rows;

    bool even               = (rows % 2 == 0);
    Integer last            = even ? N / 2 : N / 2 + 1;

    --cs;

    for (Integer j = 0; j < cols; ++j)
    {
        Integer pos         = 1;

        for (Integer i = 1; i < last; ++i, pos += 2)
        {
            Real re         = Y[pos];
            Real im         = Y[pos + 1];

            Real ci         = cs[pos];
            Real si         = cs[pos+1];

            Real re1        = re * ci;
            Real re2        = re * si;
            Real im1        = im * ci;
            Real im2        = im * si;

            Y[pos]          = re2 - im1;
            Y[pos+1]        = re1 + im2;
        };

        if (even == true)
        {
            Y[rows-1]       = Y[rows-1] * 0.5;
            for (Integer k = N-3; k >= 1; k -= 2)
            {
                Y[k]        = Y[k] + Y[k+2];
            };
        }
        else
        {
            for (Integer k = N-4; k >= 1; k -= 2)
            {
                Y[k]        = Y[k] + Y[k+2];
            };
        };

        Y                   += Y_ld;
    };
};

template<Integer Rows, Integer Y_LD>
void recurrence_dct2_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols,
                                           const Real* cs)
{
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    Integer N               = cols;

    bool even               = (cols % 2 == 0);
    Integer last            = even ? N / 2 : N / 2 + 1;

    Real* restricted Y_1    = Y + Y_ld;
    Real* restricted Y_2    = Y + 2 * Y_ld;

    for (Integer i = 1, pos = 0; i < last; ++i, pos += 2)
    {
        Real ci             = cs[pos];
        Real si             = cs[pos+1];

        //This should be vectorized by compiler (checked)
        for (Integer j = 0; j < rows; ++j)
        {
            Real re         = Y_1[j];
            Real im         = Y_2[j];

            Y_1[j]          = re * si - im * ci;
            Y_2[j]          = re * ci + im * si;
        };

        Y_1                 += 2 * Y_ld;
        Y_2                 += 2 * Y_ld;
    };

    Integer k;
    if (even == true)
    {
        matcl::level1::ay<Real, Real, Rows>::eval(Y + (cols-1)*Y_ld, rows0, 0.5);

        Y_1                 = Y + (N-3) * Y_ld;
        Y_2                 = Y + (N-1) * Y_ld;
        k                   = N - 3;
    }
    else
    {
        Y_1                 = Y + (N-4) * Y_ld;
        Y_2                 = Y + (N-2) * Y_ld;
        k                   = N - 4;
    };

    for (; k >= 1; k -= 2)
    {
        matcl::level1::ypx<Real, Real, Rows>::eval(Y_1, Y_2, rows0);

        Y_1             -= 2*Y_ld;
        Y_2             -= 2*Y_ld;
    };
};

void recurrence_dct3_dim1::eval(Real* Y, Integer Y_ld, Integer rows, Integer cols, const Real* cs)
{
    cs                      += (rows/2)*3;
    Integer last_row        = rows/2;
    bool even               = (rows % 2 == 0);

    for (Integer i = 0; i < cols; ++i)
    {
        Integer k           = rows - 1;
        for (Integer j = 0; j < last_row; ++j, --k)
        {
            Real tmp1       = Y[j] + Y[k];
            Real tmp2       = Y[j] - Y[k];
            Real tmp3       = tmp2 * cs[j];

            Y[j]            = tmp1 + tmp3;
            Y[k]            = tmp1 - tmp3;
        };

        if (even == false)
        {
            Y[last_row]     = Y[last_row] * 2.0;
        };

        Y                   += Y_ld;
    };    
};
template<Integer Rows, Integer Y_LD>
void recurrence_dct3_dim2<Rows,Y_LD>::eval(Real* Y, Integer Y_ld0, Integer rows0, Integer cols, 
                                           const Real* cs)
{
    cs                      += (cols/2)*3;
    Integer last_row        = cols/2;
    bool even               = (cols % 2 == 0);
    Integer rows            = get_int<Rows>::eval(rows0);
    Integer Y_ld            = get_int<Y_LD>::eval(Y_ld0);

    Real* restricted Y_j    = Y;
    Real* restricted Y_k    = Y + (cols - 1) * Y_ld;

    for (Integer j = 0; j < last_row; ++j)
    {
        Real csj            = cs[j];

        //This should be vectorized by compiler (checked)
        for (Integer i = 0; i < rows; ++i)
        {
            Real tmp1       = Y_j[i] + Y_k[i];
            Real tmp2       = Y_j[i] - Y_k[i];
            Real tmp3       = tmp2 * csj;

            Y_j[i]          = tmp1 + tmp3;
            Y_k[i]          = tmp1 - tmp3;
        };

        Y_j                 += Y_ld;
        Y_k                 -= Y_ld;
    };

    if (even == false)
    {
        matcl::level1::ay<Real, Real, Rows>::eval(Y + last_row * Y_ld, rows0, 2.0);
    };
};



#define INSTANTIATE_DCT_3(class_name)\
template struct class_name<0,0,0>;  \
template struct class_name<1,0,0>;  \
template struct class_name<2,0,0>;  \
template struct class_name<3,0,0>;  \
template struct class_name<4,0,0>;  \
template struct class_name<1,1,0>;  \
template struct class_name<2,2,0>;  \
template struct class_name<3,3,0>;  \
template struct class_name<4,4,0>;  \
template struct class_name<1,0,1>;  \
template struct class_name<2,0,2>;  \
template struct class_name<3,0,3>;  \
template struct class_name<4,0,4>;  \
template struct class_name<1,1,1>;  \
template struct class_name<2,2,2>;  \
template struct class_name<3,3,3>;  \
template struct class_name<4,4,4>;

#define INSTANTIATE_DCT_2(class_name)\
template struct class_name<0,0>;    \
template struct class_name<1,0>;    \
template struct class_name<2,0>;    \
template struct class_name<3,0>;    \
template struct class_name<4,0>;    \
template struct class_name<1,1>;    \
template struct class_name<2,2>;    \
template struct class_name<3,3>;    \
template struct class_name<4,4>;

INSTANTIATE_DCT_3(make_permute_dct2_dim2)
INSTANTIATE_DCT_3(make_pack_dct3_dim2)
INSTANTIATE_DCT_3(mult_by_cos_dim2)
INSTANTIATE_DCT_3(prepare_dct1_dim2)
INSTANTIATE_DCT_3(prepare_dct2_dim2)
INSTANTIATE_DCT_3(prepare_dct3_dim2)

INSTANTIATE_DCT_2(extract_all_from_pack_with_mult_dim2)
INSTANTIATE_DCT_2(make_permute_dct3_dim2)
INSTANTIATE_DCT_2(make_permute_dct2_inpl_dim2)
INSTANTIATE_DCT_2(make_pack_dct3_inpl_dim2)
INSTANTIATE_DCT_2(recursion_dct4_dim2)
INSTANTIATE_DCT_2(mult_by_cos_inpl_dim2)
INSTANTIATE_DCT_2(recurrence_dct1_dim2)
INSTANTIATE_DCT_2(prepare_dct1_inpl_dim2)
INSTANTIATE_DCT_2(recurrence_dct2_dim2)
INSTANTIATE_DCT_2(prepare_dct2_inpl_dim2)
INSTANTIATE_DCT_2(prepare_dct3_inpl_dim2)
INSTANTIATE_DCT_2(recurrence_dct3_dim2)

}}
