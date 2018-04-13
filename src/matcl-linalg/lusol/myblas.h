#pragma once

#include "matcl-matrep/general/config.h"
#include "commonlib.h"

namespace lusol
{

template<class T>
inline int amax(int n, const T* x)
{
    using TR = typename details::real_type<T>::type;

    if (n <= 1)
        return n;

    ++x;

    const int n4    = n/4;
    int i           = 0;

    TR val_max0     = 0;
    TR val_max1     = 0;

    int pos_max0    = 0;
    int pos_max1    = 0;

    TR tmp0;
    TR tmp1;
    TR tmp2;
    TR tmp3;

    for (i = 0; i < n4; i += 4)
    {
        tmp0 = abs(x[i+0]);
        tmp1 = abs(x[i+1]);
        tmp2 = abs(x[i+2]);
        tmp3 = abs(x[i+3]);

        if (tmp0 > val_max0)
        {
            val_max0 = tmp0;
            pos_max0 = i;
        };
        if (tmp1 > val_max0)
        {
            val_max0 = tmp1;
            pos_max0 = i+1;
        };

        if (tmp2 > val_max1)
        {
            val_max1 = tmp2;
            pos_max1 = i+2;
        };
        if (tmp3 > val_max1)
        {
            val_max1 = tmp3;
            pos_max1 = i+3;
        };
    };

    if (val_max1 > val_max0)
    {
        val_max0 = val_max1;
        pos_max0 = pos_max1;
    }

    for (; i < n; ++i)
    {
        tmp0 = abs(x[i]);
        if (tmp0 > val_max0)
        {
            val_max0 = tmp0;
            pos_max0 = i;
        };
    };

    return pos_max0 + 1;
};

};