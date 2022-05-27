

#pragma once

#include "mpi.h"

#include "pp.h"

namespace PP_NS
{
    class Ptrs
    {
    public:
        Ptrs(PP *ptr) :
            pp(ptr),
            mem(ptr->mem)
            {}

        virtual ~Ptrs() {}

    protected:
        PP *pp;
        Mem *&mem;
    };
}

