#include <iostream>
#include <iomanip>
#include <vector>

/* ModeCode (MC) include files */
#include "mem.h"
//#include "poptimer.h"
//#include "in.h"
//#include "verify.h"
//#include "ifc2mcc.h"
//#include "asr.h"
//#include "compute.h"
//#include "tic.h"
//#include "qhgk.h"
//#include "visualize.h"
#include "powerspec.h"
//#include "config.h"

/* MPI include file */
#include "mpi.h"

using namespace std;

using namespace PP_NS;

PP::PP(int narg, char **arg)
{

    /************************** Set up MPI settings **************************/

    int color,key,global,local;
    MPI_Comm comm;

    // Split the communicators so that multiple instances of LAMMPS can be run
    MPI_Comm_rank(MPI_COMM_WORLD, &global);
    color = global / 1; // Change "1" to 2 in order to use 2 procs per instance, etc..
    key = global; 
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
    MPI_Comm_rank(comm,&local);

    //  Get the number of processes.
    nprocs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group
    //  Get the individual process ID.
    rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

    /************************** Initial Screen Output **************************/
    if (rank == 0)
    {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                          Post Processor 0.0                     +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        //std::cout << " Running on " << nprocs << " procs" << std::endl;
    }
  
    //poptimer = new PopTimer(this);

    //if (rank == 0) std::cout << " Job started at " << poptimer->DateAndTime() << std::endl;

    /************************** Proceed **************************/

    // Extract desired task
    task = std::string(arg[1]);
    //std::cout << task << std::endl;

    // Dynamically allocated pointers
    create();

    // Grab user input to initialize settings
    initialize();
    
    // Run the code and desired task
    run(narg,arg);

    // Delete dynamically allocated pointers
    finalize();

    /*
    if (rank == 0) std::cout << std::endl << " Job finished at " 
        << poptimer->DateAndTime() << std::endl;
    */
}

void PP::create()
{
    mem = new Mem(this);

    if (task=="powerspec"){
        powerspec = new Powerspec(this);
    }
    //else{
    //    printf(" INVALID TASK.\n");
    //}
    

}

void PP::initialize()
{


}

void PP::run(int nargs, char **args)
{

  if (task=="powerspec"){
  
    if (rank==0) printf(" Calculating powerspectrum from data.\n");
    //if (rank==0) printf(" Reading EMAT.\n");
    
    //postproc->readEmat();
    //if (rank==0) printf(" Finished reading EMAT.\n");
    
    //int task_postproc = atoi(args[2]); // 1 - loop over ensembles, calculate FFT and |FFT|^2 (power spectrum) of all mode amplitudes and velocities, calculate DOS overlap for all pairs, then ensemble average.
    //postproc->task=task_postproc;
    
    if (nargs == 5){
        if (rank==0){
            printf(" ntimesteps, natoms, sampling_interval: %d %d %f\n", atoi(args[2]), atoi(args[3]), atof(args[4]));
        }
        
        //powerspec->test(); // Test powerspec calculation and compare with python example.
        powerspec->ntimesteps = atoi(args[2]); // ntimesteps DIVIDED by ndump, or how many timesteps per dump.
        powerspec->natoms = atoi(args[3]);
        powerspec->sampling_interval = atof(args[4]); // e.g. 5fs timestep with dump every 5 timesteps has sampling interval 0.0025 ps
                                                      //      1fs timestep with dump every 5 timesteps has sampling interval 0.0005 ps
        //powerspec->natoms = atoi(args[6]);
          
        powerspec->initialize();
    }
    else{
        printf("Need 5 args\n");
    }
    
  }
}

void PP::finalize()
{

    delete mem;

    if (task=="powerspec"){
        delete powerspec;
    }

}

PP::~PP()
{

};
