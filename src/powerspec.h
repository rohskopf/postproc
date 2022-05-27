#pragma once

#include <vector>
#include <string>
#include "mpi.h"

#include <iostream>
#include <new>
#include <cstdlib>
#include "ptrs.h"

using namespace std;

namespace PP_NS
{
  class Powerspec: protected Ptrs
  {
  public:
    Powerspec(class PP *);
    ~Powerspec();

    // Constructor.
    FILE * fh_debug;
    int rank;

    // Functions.
    void initialize();
    void test(); // Test compare to the python ps script.
    
    // variables
    int ntimesteps; // number of timesteps DIVIDED by ndump, or how many timesteps per dump.
                    // e.g. if nsteps = 10000 in LAMMPS, and we dumped data every 5 timesteps, we put 2000 here.
    int natoms;
    double sampling_interval;
    double **data;
    double *time;
 
  };
}

