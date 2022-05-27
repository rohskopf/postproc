/*
This class post processes simulation data to analyze the modes.
Current functions include:

*/

#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <random>
#include "mpi.h"
#include <fftw3.h>

#include "powerspec.h"
#include "mem.h"
//#include "config.h"

#include <ctime>

using namespace std;

using namespace PP_NS;

Powerspec::Powerspec(PP *pp) : Ptrs(pp) {
    //fh_debug = fopen("DEBUG_PP","w");
    //fh_fc2 = fopen("FC2_ASR","w");
    
    rank = pp->rank;

}

Powerspec::~Powerspec() 
{

   //mem->deallocate(data);

};

/*
Initialize
*/

void Powerspec::initialize()
{

    // variables
    //int ntimesteps = 10000/5;
    //double sampling_interval = 0.0025;
    //int natoms = 54;
    
    //Post process the input settings to determine frequency axis of power spectrums.
    double sampling_frequency = 1.0/sampling_interval; // frequency of data collection
    double end_time = ntimesteps*sampling_interval;
    if (rank==0) printf(" End time: %e\n", end_time);
    // There will be (N/2)+1 output values from the FFT (real2complex does it symmetrically), where N = ntimesteps
    // So each index from 0 to (N/2)+1 represents a frequency. 
    // Need to divide these indices by the number of timesteps. 
    int nfreq = (ntimesteps/2)+1;
    if (rank==0) printf(" Number of freqency points: %d\n", nfreq);
    double *freq;
    mem->allocate(freq, (ntimesteps/2)+1);
    for (int f=0; f<(ntimesteps/2)+1; f++){
      freq[f] = f/(end_time);
    }
    // The maximum frequency will be ( (ntimesteps/2)/end_time )
    //printf("%d %f\n", ntimesteps, end_time);
    if (rank==0) printf(" Max frequency we can sample: %f\n", (ntimesteps/2)/end_time);
   


    // Read data
    
    char filename[1000];
    string data_dirname = "dump_spins.lammpstrj";
    sprintf(filename, "%s", data_dirname.c_str());
    if (rank==0) printf(" Data filename: %s\n", filename);
    
  // Split atoms across procs.
  int ndof = 3*natoms;
  int *nepp;
  mem->allocate(nepp, pp->nprocs); // number elements (MCC3s) per proc
  for (int p=0; p<pp->nprocs; p++){
      nepp[p] = natoms/pp->nprocs;
  }
  // divide up the remainder
  for (int p=0; p<(natoms % pp->nprocs); p++){
      nepp[p] += 1;
  }
  int start_indx = 0;
  for (int p=0; p<rank; p++){
      start_indx += nepp[p];
  }
  int end_indx = 0; //napp[0]-1;
  for (int p=0; p<rank+1; p++){
      end_indx += nepp[p];
  }
  end_indx=end_indx+1-1;
  //end_indx = end_indx+1-1;
  if (rank==0) printf("rank startindx endindx: %d %d %d\n", rank,start_indx, end_indx);
  if (rank==0){
      printf(" Splitting modes on procs like:\n");
      for (int p=0; p<pp->nprocs; p++){
          printf("  %d Atoms on proc %d.\n", nepp[p],p);
      }
  }
  
  
    // Allocate arrays for data
    mem->allocate(data, ntimesteps, nepp[rank]*3);
    mem->allocate(time, ntimesteps);
    
    
    
    if (rank==0) printf(" Opening %s\n", filename);
    //std::cout << filename_xm << std::endl;
    
    //int ndump = 5; // output dumped every this many timesteps
    int nrows_file = ((ntimesteps) + 1)*natoms;
    if (rank==0) printf(" nrows_file: %d\n", nrows_file);
    
    ifstream fh;
    fh.open(filename);
    
    double junk;
    if (fh.is_open()) {

      for (int t=0; t<(ntimesteps); t++){
        //printf("%d\n", t);
        //fh>>time[t];
        //for (int n=0; n<2; n++){
          //printf("  %d\n", n);
        //  fh>>data[t][n];
        //}
        int ii=0; //  atom index on a proc
        for (int i=0; i<natoms; i++){
        //for (int i=start_indx; i<end_indx; i++){
            //printf("i start end: %d %d %d\n", i, start_indx, end_indx);
            if ( (start_indx<=i) && (i<end_indx)){
                for (int a=0; a<3; a++){
                    fh>>data[t][3*ii+a];
                    //printf("%d %d %d\n", t,i,a);
                }
                ii++;
            }
            else{
                for (int a=0; a<3; a++){
                    fh>>junk;
                }
            }
            //printf("%d\n", ii);
        }
      }

      fh.close();
    }
    else {
      if (rank==0) printf("Unable to open %s.\n", filename);
    }
    
    // Check data
    /*
    char fname[1000];
    sprintf(fname, "data_%d.dat", rank);
    
    FILE * fh_dat;
    fh_dat = fopen(fname,"w");
    for (int t=0; t<ntimesteps; t++){
        for (int i=0; i<nepp[rank]; i++){
            for (int a=0; a<3; a++){
                fprintf(fh_dat, "%f ", data[t][3*i+a]);
            }
            fprintf(fh_dat, "\n");
        }
    }
    fclose(fh_dat);
    */
    
    // Calculate power spectrum of data.
    
    // Allocate arrays for storing PS
    double **ps;
    mem->allocate(ps, nfreq,nepp[rank]);
    
  for (int f=0; f<nfreq; f++){
    for (int n=0; n<nepp[rank]; n++){
      ps[f][n] = 0.0;
    }
  }
  
  // Declare the FFT output
  // out[i][0] and out[i][1] are the real and imaginary parts of a complex number.
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ntimesteps);
  // Declare the FFT input and plan
  double *in;
  mem->allocate(in, ntimesteps);
  // plan for Fourier transform, save the result in 'out'
  fftw_plan p; // plan for xn
  p = fftw_plan_dft_r2c_1d(ntimesteps, in, out, FFTW_ESTIMATE);  
  
  
  /*--------------------------------------------------------------
              Done declaring and allocating variables and arrays
  ---------------------------------------------------------------*/
  
    // Make first time to be zero.
    double start_time = time[0];
    for (int t=0; t<ntimesteps; t++){
      time[t] = time[t] - start_time;
    }
    
    // Subtract the mean from the data.
    for (int n=0; n<nepp[rank]; n++){
      double mean=0.0;
      for (int t=0; t<ntimesteps; t++){
        mean += data[t][n]/ntimesteps;
      }
      for (int t=0; t<ntimesteps; t++){
        data[t][n] = data[t][n]-mean;
      }
    }
    
    // Loop over mode amplitudes and velocities, calculate their power spectrum, and store it.
    //int modeindx = 11;
    //for (int n=modeindx; n<modeindx+1; n++){
    if (rank==0) printf(" Calculating power spectrums.\n");
    for (int n=0; n<nepp[rank]; n++){
      // Build the input array for FFT calculation.
      for (int t=0; t<ntimesteps; t++){
        in[t] = data[t][n];
        //printf("%f\n", in[t]);
      }
      
      // Calculate FFT
      fftw_execute(p); // amplitude
      
      // Calculate the power spectrum
      //printf("freq | powerspectrum\n");
      //fh = fopen("ps.dat","w");
      for (int f = 0; f < nfreq; f++){
        //ps1[f] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
        //ps1[i] = 1.0;
        //printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
        //fprintf(fh, "%e %e\n", freq[f],ps1[f]);
        ps[f][n] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
      }
      //fclose(fh);
    }
    
    // Now we can sum the power spectrums across atoms, and average.
    if (rank==0) printf(" Summing power spectrums.\n");
    double *psum;
    mem->allocate(psum, nfreq);
    for (int f = 0; f < nfreq; f++){
      psum[f]=0.0;
    }
    for (int n=0; n<nepp[rank]; n++){
      
      // Calculate the power spectrum
      //printf("freq | powerspectrum\n");
      //fh = fopen("ps.dat","w");
      for (int f = 0; f < nfreq; f++){
        //ps1[f] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
        //ps1[i] = 1.0;
        //printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
        //fprintf(fh, "%e %e\n", freq[f],ps1[f]);
        psum[f] += ps[f][n]/(3*nepp[rank]);
      }
      //fclose(fh);
    }
    
    // Print summed/average spectrums
    char fsum_name[1000];
    sprintf(fsum_name, "psum_%d.dat", rank);
    FILE * fh_sum;
    fh_sum = fopen(fsum_name,"w");
    for (int f=0; f<nfreq; f++){
        fprintf(fh_sum, "%e %e\n", freq[f], psum[f]);
    }
    fclose(fh_sum);
    
    // Deallocate
    fftw_free(out);
    fftw_cleanup();
    fftw_destroy_plan(p);
    mem->deallocate(in);
    mem->deallocate(data);
    mem->deallocate(time);
    mem->deallocate(freq);
    mem->deallocate(nepp);
    mem->deallocate(ps);
    mem->deallocate(psum);
    
}

void Powerspec::test()
{

    // variables
    int ntimesteps = 400000;
    double sampling_interval = 0.0025;
    int nind = 2;
    
    //Post process the input settings to determine frequency axis of power spectrums.
    double sampling_frequency = 1.0/sampling_interval; // frequency of data collection
    double end_time = ntimesteps*sampling_interval;
    if (rank==0) printf(" End time: %e\n", end_time);
    // There will be (N/2)+1 output values from the FFT (real2complex does it symmetrically), where N = ntimesteps
    // So each index from 0 to (N/2)+1 represents a frequency. 
    // Need to divide these indices by the number of timesteps. 
    int nfreq = (ntimesteps/2)+1;
    if (rank==0) printf(" Number of freqency points: %d\n", nfreq);
    double *freq;
    mem->allocate(freq, (ntimesteps/2)+1);
    for (int f=0; f<(ntimesteps/2)+1; f++){
      freq[f] = f/(end_time);
    }
    // The maximum frequency will be ( (ntimesteps/2)/end_time )
    if (rank==0) printf(" Max frequency we can sample: %f\n", (ntimesteps/2)/end_time);
    
    // Allocate arrays for data
    mem->allocate(data, ntimesteps, 2);
    mem->allocate(time, ntimesteps);
    
    // Read data
    
    char filename[1000];
    string data_dirname = "../100K/xm.dat";
    sprintf(filename, "%s", data_dirname.c_str());
    printf(" Data filename: %s\n", filename);
    
    if (rank==0) printf(" Opening %s\n", filename);
    //std::cout << filename_xm << std::endl;
    
    
    ifstream fh;
    fh.open(filename);
    
    double junk;
    if (fh.is_open()) {

      for (int t=0; t<ntimesteps; t++){
        //printf("%d\n", t);
        fh>>time[t];
        for (int n=0; n<2; n++){
          //printf("  %d\n", n);
          fh>>data[t][n];
        }
      }

      fh.close();
    }
    else {
      if (rank==0) printf("Unable to open %s.\n", filename);
    }
    
    // Allocate arrays for storing PS
    double **ps;
    mem->allocate(ps, nfreq,2);
    
  for (int f=0; f<nfreq; f++){
    for (int n=0; n<2; n++){
      ps[f][n] = 0.0;
    }
  }
  
  // Declare the FFT output
  // out[i][0] and out[i][1] are the real and imaginary parts of a complex number.
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ntimesteps);
  // Declare the FFT input and plan
  double *in;
  mem->allocate(in, ntimesteps);
  // plan for Fourier transform, save the result in 'out'
  fftw_plan p; // plan for xn
  p = fftw_plan_dft_r2c_1d(ntimesteps, in, out, FFTW_ESTIMATE);  
  
  
  /*--------------------------------------------------------------
              Done declaring and allocating variables and arrays
  ---------------------------------------------------------------*/
  
    // Make first time to be zero.
    double start_time = time[0];
    for (int t=0; t<ntimesteps; t++){
      time[t] = time[t] - start_time;
    }
    
    // Subtract the mean from the data.
    for (int n=0; n<nind; n++){
      double mean=0.0;
      for (int t=0; t<ntimesteps; t++){
        mean += data[t][n]/ntimesteps;
      }
      for (int t=0; t<ntimesteps; t++){
        data[t][n] = data[t][n]-mean;
      }
    }
    
    // Loop over mode amplitudes and velocities, calculate their power spectrum, and store it.
    //int modeindx = 11;
    //for (int n=modeindx; n<modeindx+1; n++){
    if (rank==0) printf(" Calculating power spectrums.\n");
    for (int n=0; n<nind; n++){
      // Build the input array for FFT calculation.
      for (int t=0; t<ntimesteps; t++){
        in[t] = data[t][n];
        //printf("%f\n", in[t]);
      }
      
      // Calculate FFT
      fftw_execute(p); // amplitude
      
      // Calculate the power spectrum
      //printf("freq | powerspectrum\n");
      //fh = fopen("ps.dat","w");
      for (int f = 0; f < nfreq; f++){
        //ps1[f] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
        //ps1[i] = 1.0;
        //printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
        //fprintf(fh, "%e %e\n", freq[f],ps1[f]);
        ps[f][n] = (out[f][0]*out[f][0] + out[f][1]*out[f][1])/(ntimesteps*ntimesteps);
      }
      //fclose(fh);
    }
    
    // Write the powerspectrums 
    FILE * fh_w;
    fh_w = fopen("ps.dat","w");
    //for (int n=0; n<nind; n++){
        for (int f=0; f<nfreq; f++){
            fprintf(fh_w, "%e %e\n", freq[f], ps[f][1]); //, ps[f][1]);
        }
    //}
    fclose(fh_w);
    
    // Deallocate
    fftw_free(out);
    fftw_cleanup();
    fftw_destroy_plan(p);
    mem->deallocate(in);
    mem->deallocate(data);
    mem->deallocate(time);
    mem->deallocate(freq);
    //mem->deallocate(nepp);
    mem->deallocate(ps);
}
