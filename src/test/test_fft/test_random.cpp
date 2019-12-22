#include "matcl-core/details/mpl.h"
#include "mkgen/details/utils/mpl.h"

#include <iostream>

namespace mk    = matcl::mkgen;
namespace mkl   = mk::list;
namespace mkd   = mk::details;

template<int cnt, class RS>
struct print_randoms
{
    static void print()
    {        
        std::cout << RS::value << " ";
        print_randoms<cnt - 1, RS::next >::print();
    }
};

template<class RS>
struct print_randoms<0, RS>
{
    static void print()
    {
        std::cout << RS::value << " " << std::endl;
    }
};

template<uint64_t ... V>
struct int_list;

void test_random()
{
    static const int seed   = mkd::get_time_seed();

    // Last template argument is the random seed. It is optional.
    using X     = mkd::initialize<mkd::LCE<>>::type;
    using RS0   = mkd::random_state<X>;

    // Prints the first random number    
    std::cout << RS0::value << std::endl;

    // next value
    using RS1   = RS0::next;
  
    // Prints the second random number
    std::cout << RS1::value << std::endl;
  
    using L = int_list<seed, RS0::value, RS1::value>;
    //L::print;

    print_randoms<15, RS1>::print();
}

