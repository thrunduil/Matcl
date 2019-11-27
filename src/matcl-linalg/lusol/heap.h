#pragma once 

#include "matcl-linalg/lusol/commonlib.h"

namespace lusol
{

template<class Val>
struct heap
{
    private:
        Val*            HA;
        INT*            HJ;
        INT*            HK;

    private:
        void            HDOWN(INT N, INT K);
        void            HUP(INT K);
        void            HINSERT(INT N, Val V, INT JV);

    public:
        heap();

        void            change(INT N, INT K, Val V, INT JV);
        void            remove(INT *N, INT K);
        void            build(INT N);
        void            set(INT pos, INT j, Val val);

        bool            realloc(INT newsize,INT oldsize);
        void            clear(INT len);    

        const Val*      get_Ha() const          { return HA; };
        const INT*      get_Hj() const          { return HJ; };
        const INT*      get_Hk() const          { return HK; };
};

};
