"""
Run Backward trajectory Code with Different Timesteps

Shawn Murdzek
sfm5282@psu.edu
Date Created: 17 June 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import os
import datetime as dt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Timesteps are in seconds

#steps = [60, 30, 15, 10, 5, 1] + [15]*3
#plist = [3]*6 + list(range(6, 13, 3))

steps = [15]*3
nprcl = list(range(6, 13, 3))


#---------------------------------------------------------------------------------------------------
# Create Trajectory Input Files and Run Program
#---------------------------------------------------------------------------------------------------

fptr = open('traject.input', 'r')
contents = fptr.readlines()
fptr.close()

for s, p in zip(steps, nprcl):

    contents[9]  = "  plist='plist%d.input', ! list of positions of parcels (m)\n" % p
    contents[11] = "  outfile='traj_test_%ds_p%d.nc', ! output filename\n" % (s, p)
    contents[17] = "  dt=%.1f, ! time step for parcel integration\n" % s

    fptr = open('traject.input', 'w')
    for l in contents:
        fptr.write(l)
    fptr.close()

    start = dt.datetime.now()
    os.system('../bin/traj < traject.input')
    print()
    print('dt = %d s, nprcl = %s' % (s, p))
    print('Elapsed time = ', str(dt.datetime.now() - start))


"""
End run_traj.py
""" 
