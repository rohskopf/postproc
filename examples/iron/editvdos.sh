natoms=63 # natoms + 9 
dumps_per_worker=100
lines_per_worker=$(( $natoms*$dumps_per_worker ))
#awk 'NF==10 {print}' vdos_petn.traj > vdos5
time split -d -l $lines_per_worker dump_spins.lammpstrj dump_
