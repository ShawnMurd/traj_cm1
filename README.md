# traj_cm1

### Compute Backward Trajectories Using CM1 Output

This code is a stripped-down version of the trajectory code in Ryan Hastings' 
[cm1_pp](https://github.com/RyanHastings/cm1_pp) package. 
Unlike Ryan's code, only user-specified fields are interpolated 
to the parcel locations, forward parcel integrations are not performed, and 
only fourth-order Runge-Kutta (RK4) can be used to integrate parcels in time. 
This results in cleaner code that *should* run faster. 

Unlike Ryan's code, this program does perform time interpolation between CM1
file dumps. This means that the timestep used in the RK4 scheme can be less than
the time between CM1 file dumps, which results in more accurate trajectories
than those using the time between file dumps as the timestep. This is 
illustrated by the tests that used different timesteps (found in the `tests`
directory)

Although I have made several changes to this program, most of the code is still
Ryan's original code. I encourage visitors to also check out his 
[cm1_pp](https://github.com/RyanHastings/cm1_pp) program on GitHub.

### Running traj_cm1

Code can be compiled by navigating to the `src` directory and running

`make clean && make` 

This will create an executable in the `bin` directory. To run the program,
navigate to the `input` directory and run the command

`../bin/traj < traject.input`

### Input Files

traj_cm1 has four input files:
1. __cm1list.input:__ List of all CM1 3D output files needed for the backward trajectory 
integration. First line should contain the number of output files and each subsequent line
should contain a different output file in single quotes. Files should be listed in ascending 
order (i.e., earliest to latest).
2. __outlist.input:__ List of scalar variables to interpolate to the parcel locations
at each timestep. First line should contain the number of variables and each subsequent
line should contain a different variable in single quotes.
3. __plist.input:__ List of initial parcel locations. First line should contain the number
of parcels and each subsequent line should contain the starting coordinate for a different 
parcel (in m). Leave one space (and no commas) between the x, y, and z coordinate.
4. __traject.input:__ Specifies the various parameters for the program, including the parcel
initialization time, the time to integrate the parcels to, and the timestep.
