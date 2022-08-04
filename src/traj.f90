program traj
  !#############################################################################!
  !#############################################################################!
  !####                                                                     ####!
  !####                         PROGRAM TRAJ                                ####!
  !####                                                                     ####!
  !####  This program computes trajectories of parcels from specified       ####!
  !####  starting positions. The CM1 output must be netCDF with each his-   ####!
  !####  tory dump in a different file. The timesteo for the parcel         ####!
  !####  integration does not need to equal the time between CM1 file       ####!
  !####  dumps.                                                             ####!
  !####                                                                     ####!
  !####=====================================================================####!
  !#### Shawn Murdzek, 21 Jun 2022                                          ####!
  !#### Ryan Hastings, 8 Sep 2011                                           ####!
  !####---------------------------------------------------------------------####!
  !#### v1.0 completed on 10 December 2016                                  ####!
  !#############################################################################!
  !#############################################################################!
  use globals
  use nc_inout
  use traj_utils
  implicit none
  !=============================================================================!
  !#############################################################################!
  !################################# VARIABLES #################################!
  !#############################################################################!
  !=============================================================================!
  !-----------------------------------------------------------------------------!
  ! configuration variables
  real                                        :: begintime ! time to integrate back to
  real                                        :: parcelt ! initialization time

  character(len=200),dimension(:),allocatable :: infiles ! array to hold list of CM1 files
  character(len=10),dimension(:),allocatable  :: outvars ! array to hold list of output variables
  character(len=200)                          :: cm1list ! name of file containing list of
                                                         ! CM1 output files
  character(len=200)                          :: plist ! file containing list of parcel
                                                       ! initialization positions
  character(len=200)                          :: varlist ! name of file containing list of
                                                         ! scalar variables to interpolate
  character(len=200)                          :: outfile ! output filename

  real :: vcm1 ! CM1 version

  !-----------------------------------------------------------------------------!
  ! dimensional variables
  integer                       :: npstep ! number of integration steps
  integer                       :: nfiles ! number of CM1 3D output files
  integer                       :: nparcels ! number of parcels
  integer                       :: nvars ! number of variables in varlist
  integer                       :: np ! looping variable

  real                          :: tt ! time step
  real,dimension(:),allocatable :: zzh ! scalar height array (see below for explanation)

  real,dimension(:),allocatable :: times ! times of files in infiles
  real,dimension(:),allocatable :: tp ! integration times
  integer                       :: ti, tstep ! looping variables

  !-----------------------------------------------------------------------------!
  ! netcdf variables
  integer                          :: rcode, ncid_out, varid, ncid
  integer                          :: tpid, xpid, ypid, zpid, upid, vpid, wpid
  integer,dimension(:),allocatable :: varids

  !-----------------------------------------------------------------------------!
  ! field variables
  real,dimension(:,:,:),allocatable   :: u, v, w ! wind (m/s)
  real,dimension(:,:,:),allocatable   :: u1, v1, w1 ! wind (m/s)
  real,dimension(:,:,:),allocatable   :: u2, v2, w2 ! wind (m/s)
  real,dimension(:,:,:),allocatable   :: utem, vtem, wtem ! wind (m/s)

  real,dimension(:,:,:,:),allocatable :: varcm1, varcm11, varcm12 ! scalar variables to interpolate 

  !-----------------------------------------------------------------------------!
  ! parcel variables
  real,dimension(:,:),allocatable   :: xpc, ypc, zpc ! parcel position (m)
  real,dimension(:,:),allocatable   :: up, vp, wp ! parcel winds

  real,dimension(:,:,:),allocatable :: varp ! parcel scalar variables

  !-----------------------------------------------------------------------------!
  ! looping and dummy variables
  integer i, t
  integer,dimension(1) :: tmp


  !-----------------------------------------------------------------------------!
  ! namelist

  namelist /filenames/  cm1list, plist, varlist, outfile
  namelist /parameters/ begintime, parcelt, dt, vcm1

  !=============================================================================!
  !#############################################################################!
  !################################# MAIN BODY #################################!
  !#############################################################################!
  !=============================================================================!

  write(*,*)
  write(*,*) '                  TRAJECT_CM1'
  write(*,*)

  !##############################################################################!
  !#                              INITIALIZATION                                #!
  !##############################################################################!

  !------------------------------------------------------------------------------!
  ! read namelist.input                                                          !
  !------------------------------------------------------------------------------!

  write(*,*) 'opening traject.input'
  open(8,file='traject.input')

  read(8,nml=filenames)
  write(*,*) 'cm1list    =',cm1list     ! list of cm1 files
  write(*,*) 'plist      =',plist       ! file with list of parcel starting positions
  write(*,*) 'varlist    =',varlist     ! list of scalar variables to interpolate
  write(*,*) 'outfile    =',outfile     ! outfile

  read(8,nml=parameters)
  write(*,*) 'begintime  =',begintime   ! beginning time for trajectory run
  write(*,*) 'parcelt    =',parcelt      ! parcel beginning time
  write(*,*) 'dt         =',dt          ! time step
  write(*,*) 'vcm1       =',vcm1          ! CM1 version

  close(8)

  !------------------------------------------------------------------------------!
  ! set up arrays for times                                                      !
  !------------------------------------------------------------------------------!
  ! set up array for times of parcels
  write(*,*)
  write(*,*) 'setting up arrays for parcel times'
  npstep     = int( 1 + (parcelt-begintime)/dt ) ! total number of times

  allocate( tp(npstep) )          ! parcel times
  do t=1,npstep
    tp(t) = begintime + (t-1)*dt
    write(*,*) 't,tp(t)=',t,tp(t)
  enddo

  !------------------------------------------------------------------------------!
  ! set up arrays for times of files

  ! read infilelist
  write(*,*) 'reading list of cm1 files:'
  open(8,file=cm1list)

  read(8,*) nfiles
  allocate( infiles(nfiles) )
  do i=1,nfiles
    read(8,*) infiles(i)
    write(*,*) infiles(i) 
  enddo
  close(8)

  write(*,*)
  ! read times
  allocate( times(nfiles) )
  write(*,*) 'reading times:'
  do i=1,nfiles
    call cm1_nc_gettime_traj( infiles(i), times(i) )
    write(*,*) infiles(i),times(i)
  enddo

  write(*,*)
  ! read varlist
  write(*,*) 'reading list of variables to interpolate:'
  open(8,file=varlist)
  read(8,*) nvars

  allocate( outvars(nvars) )
  do i=1,nvars
    read(8,*) outvars(i)
    write(*,*) outvars(i)
  enddo
  close(8)

  !------------------------------------------------------------------------------!
  ! get dimensions and set up arrays                                             !
  !------------------------------------------------------------------------------!
  ! get dimensions

  write(*,*) 'reading dimensions..'
  ncid = ncopn( infiles(1), NCNOWRIT, rcode )

  call cm1_nc_getdims( ncid,vcm1 )

  !------------------------------------------------------------------------------!
  ! set up spacial arrays

  write(*,*) 'reading space..'
  allocate( xh(ni) )
  allocate( xf(nip1) )
  allocate( yh(nj) )
  allocate( yf(njp1) )
  allocate( zh(nk) )
  allocate( zf(nkp1) )

  ! we will be assuming free-slip conditions, so the horizontal winds will remain
  ! constant below the lowest grid point.  zzh(0) represents the surface for
  ! mass grid points, then.  
  allocate( zzh(0:nk) )

  allocate( mzh(1:nk) )
  allocate( mzf(1:nk+1) )

  write(*,*) 'reading space-files..'
  call cm1_nc_getspace( ncid,vcm1 )

  dz = xh(2)-xh(1)
  rdz = 1.0/dz

  write(*,*) 'setting up stretching arrays'
  do i=2,nk
    mzh(i) = dz/(zf(i+1)-zf(i))
  enddo
  mzh(1)=1.0
  do i=2,nk
    mzf(i) = dz/(zh(i)-zh(i-1))
  enddo
  mzf(1)=1.0
  mzf(nk+1) = 1.0

  zzh(1:nk) = zh(1:nk)
  
  zzh(0)    = -zzh(1)

!  stop

  !---------------------------------------------------------------------------!
  ! set up other arrays, now that we have dimensions                          !
  !---------------------------------------------------------------------------!
  write(*,*) 'setting up some wind arrays'
  allocate( u(1:ni+1,1:nj,0:nk) )
  allocate( v(1:ni,1:nj+1,0:nk) )
  allocate( w(1:ni,1:nj,1:nk+1) )
  allocate( utem(1:ni+1,1:nj,0:nk) )
  allocate( vtem(1:ni,1:nj+1,0:nk) )
  allocate( wtem(1:ni,1:nj,1:nk+1) )
  allocate( u1(1:ni+1,1:nj,0:nk) )
  allocate( v1(1:ni,1:nj+1,0:nk) )
  allocate( w1(1:ni,1:nj,1:nk+1) )
  allocate( u2(1:ni+1,1:nj,0:nk) )
  allocate( v2(1:ni,1:nj+1,0:nk) )
  allocate( w2(1:ni,1:nj,1:nk+1) )


  !------------------------------------------!
  ! scalar variables to be interpolated

  allocate( varcm1(ni,nj,0:nk,nvars) )
  allocate( varcm11(ni,nj,0:nk,nvars) )
  allocate( varcm12(ni,nj,0:nk,nvars) )
  allocate( varids(nvars) )

  !----------------------------------------------------------!
  ! read in parcels
  write(*,*)
  write(*,*) 'reading parcels'
  open(8,file=plist)

  read(8,*) nparcels
  write(*,*) nparcels,' parcels'

  !----------------------------------------------------------!
  ! create netcdf file
  call nc_traject_vardef( outfile, ncid_out, npstep, nparcels, nvars, &
      outvars, tpid, xpid, ypid, zpid, upid, vpid, wpid, varids )

  !-----------------------------------------------------------!
  ! set up trajectory arrays
  write(*,*) 'setting up trajectory arrays'
  allocate( xpc(npstep,nparcels) )
  allocate( ypc(npstep,nparcels) )
  allocate( zpc(npstep,nparcels) )

  xpc=missing_value
  ypc=missing_value
  zpc=missing_value

  allocate( up(npstep,nparcels) )
  allocate( vp(npstep,nparcels) )
  allocate( wp(npstep,nparcels) )

  allocate( varp(npstep,nparcels,nvars) )

  write(*,*) 'number of backward integration steps = ',npstep

  !----------------------------------------------------------------------------!
  ! read in parcels
  do t=1,nparcels
    read(8,*) xpc(npstep,t), ypc(npstep,t), zpc(npstep,t)
    write(*,*) 'parcel ',t,' starts at x,y,z=',xpc(npstep,t),'m,',ypc(npstep,t),'m,',zpc(npstep,t),'m'
  enddo

  write(*,*)
  close(8)

  !##############################################################################!
  !#                           TRAJECTORY CALCULATIONS                          #!
  !##############################################################################!

  !------------------------------------------------------------------------------!
  ! do backward trajectory calculations

  write(*,*) '################################################################'
  write(*,*) 'starting backward trajectory calculations'
  write(*,*) 

  dt=-dt
  ti=nfiles-1

  !----------------- set up initial times ------------------------!

  write(*,*) 'opening ',infiles(ti+1),'...time ',times(ti+1)

  call read_cm1_winds( infiles(ti+1),  u2, v2, w2 )
  call read_cm1_scalars( infiles(ti+1), nvars, outvars, varcm12 )

  write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

  call read_cm1_winds( infiles(ti),  u1, v1, w1 )
  call read_cm1_scalars( infiles(ti), nvars, outvars, varcm11 )


  !----------------- loop through backwards times ----------------!
  do tstep=npstep,2,-1

    tt=tp(tstep)

    ! Check that parcel timestep lies between the two CM1 output times
    if(tt.lt.times(ti))then
      ti=ti-1

      u2 = u1
      v2 = v1
      w2 = w1
      varcm12 = varcm11
  
      write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

      call read_cm1_winds( infiles(ti),  u1, v1, w1 )
      call read_cm1_scalars( infiles(ti), nvars, outvars, varcm11 )

    endif
    
    write(*,*) '--------------------------------------------------------------'
    write(*,*) 'time step ',tstep,' (',tt,'s)'

    ! Interpolate winds and scalars to parcel timestep

    if (tstep.eq.npstep) then
      ! Only need to do this during the first timestep b/c u, v, and w are
      ! saved from the previous loop iteration
      call interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, times(ti), times(ti+1), tt )
    endif

    call interpintime_winds( u1, u2, utem, v1, v2, vtem, w1, w2, wtem, times(ti), &
      times(ti+1), tp(tstep-1) )

    call interpintime_scalars( nvars, varcm11, varcm12, varcm1, times(ti), times(ti+1), tt )

    ! Interpolate winds and scalars to each parcel location, then advance parcel back in time
    do np=1,nparcels

      if (xpc(tstep,np).ne.missing_value) then
        ! Parcel has not exited domain yet
        call interp_all( xpc(tstep,np), ypc(tstep,np), zpc(tstep,np), zzh, u, up(tstep,np), v, vp(tstep,np), &
          w, wp(tstep,np), nvars, varcm1, varp(tstep,np,:) )
      endif

      if (xpc(tstep,np).ne.missing_value) then
        ! Parcel has not exited domain yet
        call rk4( xpc(tstep-1:tstep,np), ypc(tstep-1:tstep,np), zpc(tstep-1:tstep,np), &
          up(tstep,np), vp(tstep,np), wp(tstep,np), u, v, w, utem, vtem, wtem, dt )
      endif

    enddo

    u = utem
    v = vtem
    w = wtem

  enddo ! do tstep=npstep,2,-1
  tt=tp(1)
  !--------------- do first time ---------------------------!

    call interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, times(1), times(2), tt )
    call interpintime_scalars( nvars, varcm11, varcm12, varcm1, times(1), times(2), tt )

  do np=1,nparcels

    call interp_all( xpc(1,np), ypc(1,np), zpc(1,np), zzh, u, up(1,np), v, vp(1,np), &
      w, wp(1,np), nvars, varcm1, varp(1,np,:) )

  enddo


  !##############################################################################!
  !#                                WRITE OUTPUT                                #!
  !##############################################################################!

  write(*,*) '###################################################################'
  write(*,*) '# write out netcdf file'
  write(*,*)


  call nc_traject_writeout( npstep, nparcels, ncid_out, tpid, tp, xpid, xpc, ypid, ypc, zpid, zpc, &
      upid, up, vpid, vp, wpid, wp, nvars, varids, varp )


  write(*,*) 'done writing out netcdf file'
  write(*,*) 'program terminated normally'

end program
