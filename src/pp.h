
/* Declaration of pointers used in the whole program. */

#pragma once
#include <string>
#include <vector>
#include "mpi.h"


namespace PP_NS
{
    class PP
    {
    public:


        class Mem *mem;
        class Powerspec *powerspec;
        //class PopTimer *poptimer;
        //class Config;
        PP(int, char **);
        ~PP();
        void create();
        void initialize();
        void run(int,char **);
        void finalize();
        int nprocs;
        int rank;

        int seed;

        double delta;
        double cutoff;
        int order;

        std::string task; // task to perform

    };
}

