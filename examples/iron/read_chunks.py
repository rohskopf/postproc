import warnings
import numpy as np
from mpi4py import MPI

ntimesteps = 100 # number of timesteps per dump
#fname = "dump_00"

# MPI settings
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

fname = "dump_0" + str(rank)

print(f"rank: {rank}")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    N = int(np.genfromtxt(fname, skip_header=3, max_rows=1))
    box_data = np.genfromtxt(fname, skip_header=5, max_rows=3)
    data = np.genfromtxt(fname, skip_header=9, invalid_raise=False)

# Remove nans
#print(f"N: {N}")
data = data[~np.isnan(data).all(axis=1)]
data_split = np.split(data, ntimesteps)
data_split_flat = [x.flatten() for x in data_split]
print(data_split_flat)
stacked = np.vstack(data_split_flat)
#stacked = np.stack(data_split_flat) #, axis=1)
#print(stacked)

#print(data_split[1])

#print(data[0:100,:])
print(f"Nrows on rank {rank}:")
#print(np.shape(data)[0])
print(np.shape(stacked))
