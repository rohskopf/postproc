#include <stdlib.h>
#include <iostream>
#include "pp.h"
#include "mpi.h"

using namespace PP_NS;

int main(int argc, char **argv)
{

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  /* Begin a MC instance */
  PP *pp = new PP(argc, argv);

  /* Delete the memory */
  delete pp;

  /* Close MPI */
  int MPI_Comm_free(MPI_Comm *comm);
  MPI_Finalize();

  return EXIT_SUCCESS;
}

