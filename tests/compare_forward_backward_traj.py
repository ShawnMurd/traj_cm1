"""
Compare Forward and Backward Parcel Trajectories

Shawn Murdzek
sfm5282@psu.edu
Date Created: 17 June 2022
"""

#%%-------------------------------------------------------------------------------------------------
# Input Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


#%%-------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# CM1 3D output at parcel initialization time

cm1_fname = ('/storage/home/sfm5282/scratch/mmp_supercell_ctrl_reruns/sims/prcl_sims/cLFC_lcl500/' +
             '3600s_restart/cm1out_000023.nc')

# NetCDF files holding backward parcel trajectory info and timesteps for these output files

dt = [60, 30, 15, 10, 5, 1]
back_fname = ['traj_test_%ds_p3.nc' % i for i in dt]

# Runtime for each backward trajectory calculation (one per output file, in s)

#runtime = [278, 298, 294, 312, 400, 1151]

# Forward parcel trajectory output and parcels used for comparisons with backward trajectories

for_fname = ('/storage/home/sfm5282/scratch/mmp_supercell_ctrl_reruns/sims/prcl_sims/cLFC_lcl500/' +
             '3600s_restart/cm1out_pdata.nc')
nprcl = [11979, 13648, 88850]
tinds = list(range(4, 41))

xlim = [-12, 0]
ylim = [-2, 10]
zind = 0
refzind = 20
clim = [0, 1000]

save_fname = 'forward_backward_compare.pdf'


#%%-------------------------------------------------------------------------------------------------
# Extract Data
#---------------------------------------------------------------------------------------------------

cm1 = xr.open_dataset(cm1_fname)
forward = xr.open_dataset(for_fname)

xind = np.where(np.logical_and(cm1['xh'] >= xlim[0], cm1['xh'] <= xlim[1]))[0]
yind = np.where(np.logical_and(cm1['yh'] >= ylim[0], cm1['yh'] <= ylim[1]))[0]

x2d, y2d = np.meshgrid(cm1['xh'][xind].values, cm1['yh'][yind].values)
z = cm1['zh'].values
u = cm1['uinterp'][0, zind, yind, xind].values
v = cm1['vinterp'][0, zind, yind, xind].values


#%%-------------------------------------------------------------------------------------------------
# Plot Trajectories
#---------------------------------------------------------------------------------------------------

pdf = PdfPages(save_fname)

for f in back_fname:

    back = xr.open_dataset(f)
    back = back.where(back['xp'] > -8e9)
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 7))
    
    # Plot CM1 fields
    
    cax = axes[0].contourf(x2d, y2d, cm1['dbz'][0, refzind, yind, xind], np.arange(5, 75, 5), 
                       cmap='Greys_r', extend='max')
    cbar = plt.colorbar(cax, ax=axes[0], orientation='horizontal', pad=0.06)
    cbar.set_label('z = %.3f km reflectivity (dBZ)' % z[zind], size=14)
    
    st = 3
    vect = axes[0].quiver(x2d[::st, ::st], y2d[::st, ::st], u[::st, ::st], v[::st, ::st])
    
    axes[0].set_title('CM1 Fields at %d s' % (np.float64(cm1['time'][0].values)*1e-9), size=16)
    
    # Plot parcel trajectories
    
    for j, (nb, nf) in enumerate(zip(back['parcel'], nprcl)):
        
        axes[0].plot(1e-3*back.xp[nb, -1], 1e-3*back.yp[nb, -1], 'bo', markersize=8)
        axes[0].plot(1e-3*back.xp[nb, :], 1e-3*back.yp[nb, :], 'b-')
        axes[0].plot(1e-3*forward.x[tinds, nf], 1e-3*forward.y[tinds, nf], 'r--')
        if j == 0:
            axes[1].plot(back.tp, 1e-3*back.zp[nb, :], 'b-', label='backward')
            axes[1].plot(forward.mtime[tinds, nf], 1e-3*forward.z[tinds, nf], 'r--', label='forward')
        else:
            axes[1].plot(back.tp, 1e-3*back.zp[nb, :], 'b-')
            axes[1].plot(forward.mtime[tinds, nf], 1e-3*forward.z[tinds, nf], 'r--')
    
    axes[1].set_ylim(bottom=0)
    axes[1].grid()
    axes[1].legend()
    axes[1].set_xlabel('time', size=14)
    axes[1].set_ylabel('height (km)', size=14)
    
    plt.suptitle(f, size=20)
    
    pdf.savefig(fig)
    #plt.show() 


#%%-------------------------------------------------------------------------------------------------
# Plot Runtimes
#---------------------------------------------------------------------------------------------------
'''
fig = plt.figure()

plt.plot(dt, runtime, 'b-')
plt.plot(dt, runtime, 'bo', markersize=8)

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('parcel integration timestep (s)', size=12)
plt.ylabel('runtime (s)', size=12)
plt.title('Runtime For 3 Parcels', size=16)
plt.grid()

pdf.savefig(fig)
plt.show()
'''
pdf.close()


"""
End plot_traj.py
"""
