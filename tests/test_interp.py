"""
Checking Interpolation of Scalars to Backward Parcel Trajectories

Shawn Murdzek
sfm5282@psu.edu
Date Created: 20 June 2022
"""

#%%-------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import xarray as xr
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)


#%%-------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

cm1_dir = ('/storage/home/sfm5282/scratch/mmp_supercell_ctrl_reruns/sims/prcl_sims/cLFC_lcl500/' +
           '3600s_restart/')

cm1_files = [cm1_dir + 'cm1out_000021.nc',
             cm1_dir + 'cm1out_000022.nc']

traj_file = 'traj_test_15s_p3.nc'
test_times = [4095, 4110, 4125]
fields = ['qv', 'qr', 'qc', 'qi', 'qi2']


#%%-------------------------------------------------------------------------------------------------
# Open Datasets
#---------------------------------------------------------------------------------------------------

cm1_ds = xr.open_mfdataset(cm1_files)
    
traj = xr.open_dataset(traj_file)
traj = traj.where(traj['xp'] > -8e9)


#%%-------------------------------------------------------------------------------------------------
# Test Spatial Interpolation
#---------------------------------------------------------------------------------------------------

for i, f in enumerate(cm1_files):
    print()
    print(f)
    time = np.float64(cm1_ds.time[i].values)*1e-9
    tind = np.where(np.isclose(traj.tp, time))[0][0]
    for field in fields:
        print('Field = %s' % field)
        for n in traj.parcel:
            if not np.isnan(traj[field][n, tind].values):
                truth = cm1_ds[field][i, :, :, :].interp(zh=(1e-3*traj.zp[n, tind].values), 
                                                         yh=(1e-3*traj.yp[n, tind].values),
                                                         xh=(1e-3*traj.xp[n, tind].values))
                print('percent diff = %.5f (truth = %.5f)' % 
                      (100*((traj[field][n, tind].values - truth) / truth), truth))
    

#%%-------------------------------------------------------------------------------------------------
# Test Temporal and Spatial Interpolation
#---------------------------------------------------------------------------------------------------

print('------------------------------------------------------')

for t in test_times:
    print()
    print('Time = %.1f s' % t)
    tind = np.where(np.isclose(traj.tp, t))[0][0]
    time = np.timedelta64(t, 's')
    for field in fields:
        print('Field = %s' % field)
        for n in traj.parcel:
            if not np.isnan(traj[field][n, tind].values):
                truth = cm1_ds[field].interp(zh=(1e-3*traj.zp[n, tind].values), 
                                             yh=(1e-3*traj.yp[n, tind].values),
                                             xh=(1e-3*traj.xp[n, tind].values),
                                             time=time)
                print('percent diff = %.5f (truth = %.5f)' % 
                      (100*((traj[field][n, tind].values - truth) / truth), truth))


#%%-------------------------------------------------------------------------------------------------
# Digging Into a Specific Value
#---------------------------------------------------------------------------------------------------

tprcl = 4140.
nprcl = 1
cm1_tind = 1
prcl_tind = np.where(np.isclose(tprcl, traj.tp.values))[0][0]
field = 'qr'

xp = 1e-3*traj.xp[nprcl, prcl_tind].values
yp = 1e-3*traj.yp[nprcl, prcl_tind].values
zp = 1e-3*traj.zp[nprcl, prcl_tind].values

truth = cm1_ds[field].interp(zh=zp, yh=yp, xh=xp, time=np.timedelta64(int(tprcl), 's'))

# Determine lower left corner grid indices

dx = np.round_(cm1_ds.xh[1].values - cm1_ds.xh[0].values, decimals=3)
dy = np.round_(cm1_ds.yh[1].values - cm1_ds.yh[0].values, decimals=3)
dz = np.round_(cm1_ds.zh[1].values - cm1_ds.zh[0].values, decimals=3)

xind = int(np.floor((xp - cm1_ds.xh[0].values) / dx))
yind = int(np.floor((yp - cm1_ds.yh[0].values) / dy))
zind = int(np.floor((zp - cm1_ds.zh[0].values) / dz))

print('------------------------------------------------------')
print()
print('Digging Into a Specific Parce')
print('parcel location = (%.3f, %.3f, %.3f) km' % (xp, yp, zp))
print('parcel value for %s = %.5f' % (field, traj[field][nprcl, prcl_tind]))
print('CM1 lower-left corner = (%.3f, %.3f, %.3f) km' % (cm1_ds.xh[xind], cm1_ds.yh[yind],
                                                         cm1_ds.zh[zind]))
print('CM1 interpolated value for %s = %.5f' % (field, truth))

print()
print('CM1 values at z = %.3f km' % cm1_ds.zh[zind])
print('%.5f  %.5f' % (cm1_ds[field][cm1_tind, zind, yind+1, xind], 
                      cm1_ds[field][cm1_tind, zind, yind+1, xind+1]))
print('%.5f  %.5f' % (cm1_ds[field][cm1_tind, zind, yind, xind], 
                      cm1_ds[field][cm1_tind, zind, yind, xind+1]))

print()
print('CM1 values at z = %.3f km' % cm1_ds.zh[zind+1])
print('%.5f  %.5f' % (cm1_ds[field][cm1_tind, zind+1, yind+1, xind], 
                      cm1_ds[field][cm1_tind, zind+1, yind+1, xind+1]))
print('%.5f  %.5f' % (cm1_ds[field][cm1_tind, zind+1, yind, xind], 
                      cm1_ds[field][cm1_tind, zind+1, yind, xind+1]))


"""
End test_interp.py
"""
