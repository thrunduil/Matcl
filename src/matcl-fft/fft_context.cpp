#include "matcl-fft/fft_context.h"
#include "dct_helpers.h"
#include "mkl_wrapper/fft_wrapper.h"

namespace matcl { namespace fft
{

namespace details
{

struct problem_type
{
    Integer                 m_length;
    fft_context::fft_type   m_type;

    friend bool operator<(const problem_type& p1, const problem_type& p2)
    {
        if (p1.m_length < p2.m_length)
            return true;
        if (p1.m_length > p2.m_length)
            return true;

        return p1.m_type < p2.m_type;
    };
};

struct dct_data
{
    Real*   m_ptr;
    void    release();
};
struct fft_data
{
    void*   m_ptr;
    void    release();
};

class fft_context_impl
{
    private:
        using problem_map_dct   = std::map<problem_type, dct_data>;
        using problem_map_fft   = std::map<problem_type, fft_data>;

    private:
        problem_map_dct m_map_dct;
        problem_map_fft m_map_fft;

    public:
        fft_context_impl(){};
        ~fft_context_impl();

        void        initialize(Integer length, fft_context::fft_type type);
        const Real* get_table(Integer length, fft_context::fft_type type) const;
        void*       initialize_fft(Integer length, fft_context::fft_type type);

    private:
        fft_context_impl(const fft_context_impl&) = delete;
        fft_context_impl& operator=(const fft_context_impl&) = delete;
};

fft_context_impl::~fft_context_impl()
{
    for(auto pos: m_map_dct)
    {
        pos.second.release();
    };
    for(auto pos: m_map_fft)
    {
        pos.second.release();
    };
};

void fft_context_impl::initialize(Integer length, fft_context::fft_type type)
{
    problem_type pt{length,type};
    auto pos = m_map_dct.find(pt);

    if (pos != m_map_dct.end())
        return;

    switch(type)
    {
        case fft_context::dct_1:
        case fft_context::dct_23_1:
        case fft_context::dct_23_2:
            break;
        default:
            return;
    };

    initialize_fft(length, fft_context::fft_real);

    dct_data data;
    data.m_ptr      = create_dct_table(length, type);

    m_map_dct.insert(pos, problem_map_dct::value_type(pt,data));
    return;
};
const Real* fft_context_impl::get_table(Integer length, fft_context::fft_type type) const
{
    problem_type pt{length,type};
    auto pos = m_map_dct.find(pt);

    if (pos == m_map_dct.end())
    {
        throw std::runtime_error("dct problem not initialized");
    };

    return pos->second.m_ptr;
};
void* fft_context_impl::initialize_fft(Integer length, fft_context::fft_type type)
{
    problem_type pt{length,type};
    auto pos = m_map_fft.find(pt);

    if (pos != m_map_fft.end())
    {
        return pos->second.m_ptr;
    };

    details::mkl_fft::fft_wrapper wrapper;

    fft_data data;
    data.m_ptr = wrapper.make_mkl_context(length, type);

    m_map_fft.insert(pos, problem_map_fft::value_type(pt,data));

    return data.m_ptr;
};

void dct_data::release()
{
    destroy_dct_table(m_ptr);
};

void fft_data::release()
{
    details::mkl_fft::fft_wrapper wrapper;
    wrapper.destroy_mkl_context(m_ptr);
};


};

fft_context::fft_context()
    :m_impl(new details::fft_context_impl())
{};
fft_context::~fft_context()
{};

void fft_context::initialize(Integer length, fft_type type)
{
    m_impl->initialize(length,type);
};
const Real* fft_context::get_table(Integer length, fft_type type) const
{
    return m_impl->get_table(length,type);
};
void* fft_context::initialize_fft(Integer length, fft_type type)
{
    return m_impl->initialize_fft(length,type);
};

}}
