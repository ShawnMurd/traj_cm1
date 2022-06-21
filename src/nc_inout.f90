module nc_inout
  !#############################################################################!
  !############################################################################!
  !####                                                                    ####!
  !####                          MODULE NC_INOUT                           ####!
  !####                                                                    ####!
  !####  Module containing subroutines for input and output from netCDF    ####!
  !####  files useful for CM1_PP.                                          ####!
  !####====================================================================####!
  !####  Ryan Hastings, 7 September 2011                                   ####!
  !####                                                                    ####!
  !####--------------------------------------------------------------------####!
  !#### Modified by Ryan Hastings on 26 September 2012                     ####!
  !####      added nc_error and modified existing subroutines to use       ####!
  !####      F77 netcdf libraries (i.e., nf_inq_dimlen rather than ncdinq, ####!
  !####      which doesn't seem to be working)                             ####!
  !####                                                                    ####!
  !#### Modified by Shawn Murdzek on 15 June 2022                          ####!
  !####--------------------------------------------------------------------####!
  !#### v0.8 28 Mar 2013                                                   ####!
  !####--------------------------------------------------------------------####!
  !#### v1.0 13 Dec 2016                                                   ####!
  !############################################################################!
  !############################################################################!
  use globals
  implicit none
  include 'netcdf.inc'

  contains

    !###########################################################################!
    !#                    SUBROUTINE CM1_NC_GETTIME_TRAJ                       #!
    !#                                                                         #!
    !#  Gets the time from a cm1 output file.                                  #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_gettime_traj( filename, tt )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ------------------------------!
      character(len=200) :: filename
      !------------------------- output variables ------------------------------!
      real tt

      !-------------------------- local variables ------------------------------!
      integer ncid, rcode, timedid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      ncid = ncopn( filename, NCNOWRIT, rcode )
      timedid = ncvid( ncid, 'time', rcode )
      call ncvgt( ncid, timedid, 1, 1, tt, rcode )
      call ncclos( ncid, rcode )

    end subroutine

    !#########################################################################!
    !#                        SUBROUTINE NC_HANDLE_ERROR                     #!
    !#                                                                       #!
    !#  Handle error.                                                        #!
    !#-----------------------------------------------------------------------#!
    !# Ryan Hastings, 26 September, 2012                                     #!
    !#########################################################################!
      subroutine nc_handle_error( retval )
      implicit none

      integer retval

      if (retval.ne.nf_noerr) then
        write(*,*) 'error: ',nf_strerror(retval)
        stop 2
      endif

    end subroutine



    !###########################################################################!
    !#                        SUBROUTINE CM1_NC_GETDIMS                        #!
    !#                                                                         #!
    !#  Get dimensions from CM1 output file.                                   #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_getdims( ncid,vcm1 )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!

      !--------------------------- input variables -----------------------------!
      integer ncid
      real vcm1

      !--------------------------- local variables -----------------------------!
      integer rcode, nidid, njdid, nkdid, nip1did, njp1did, nkp1did
      character(len=2) char2
      character(len=4) char4

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*)
      write(*,*) 'reading dimensions'
      write(*,*)

      if (vcm1.lt.20) then
        nidid   = ncdid(ncid,'ni',rcode)
        nip1did = ncdid(ncid,'nip1',rcode)
        njdid   = ncdid(ncid,'nj',rcode)
        njp1did = ncdid(ncid,'njp1',rcode)
        nkdid   = ncdid(ncid,'nk',rcode)
        nkp1did = ncdid(ncid,'nkp1',rcode)
      else
        nidid   = ncdid(ncid,'xh',rcode)
        nip1did = ncdid(ncid,'xf',rcode)
        njdid   = ncdid(ncid,'yh',rcode)
        njp1did = ncdid(ncid,'yf',rcode)
        nkdid   = ncdid(ncid,'zh',rcode)
        nkp1did = ncdid(ncid,'zf',rcode)
      endif

      call nc_handle_error( nf_inq_dimlen(ncid,nidid,ni) )
      call nc_handle_error( nf_inq_dimlen(ncid,nip1did,nip1) )
      call nc_handle_error( nf_inq_dimlen(ncid,njdid,nj) )
      call nc_handle_error( nf_inq_dimlen(ncid,njp1did,njp1) )
      call nc_handle_error( nf_inq_dimlen(ncid,nkdid,nk) )
      call nc_handle_error( nf_inq_dimlen(ncid,nkp1did,nkp1) )

      write(*,*)
      write(*,*) 'ni,nj,nk=',ni,nj,nk
      write(*,*) 'nip1,njp1,nkp1=',nip1,njp1,nkp1
      write(*,*)

      return

    end subroutine
    !===========================================================================!
    !###########################################################################!
    !#                        SUBROUTINE CM1_NC_GETTIME                        #!
    !#                                                                         #!
    !#  Gets the time from a cm1 output file.                                  #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_gettime( ncid )
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ------------------------------!
      integer ncid

      !-------------------------- local variables ------------------------------!
      integer rcode, timedid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      timedid = ncvid( ncid, 'time', rcode )
      call ncvgt( ncid, timedid, 1, 1, time, rcode )

    end subroutine
   !===========================================================================!
    !===========================================================================!
    !###########################################################################!
    !#                        SUBROUTINE CM1_NC_GETSPACE                       #!
    !#                                                                         #!
    !#  Reads in spacial variables from cm1 file:  xh, yh, zh, zf.  Calculates #!
    !# mzh and mfh, too.  These are the Jacobian functions. 'h' indicates pos- #!
    !# ition on mass grid, 'f' on momentum grid.  mzh(k) is the Jacobian for   #!
    !# the kth point on zh, mzf for zf.                                        #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_getspace( ncid,vcm1 )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ------------------------------!
      integer ncid
      real vcm1

      !-------------------------- local variables ------------------------------! 
      integer i, j, k
      integer rcode
      integer xhid, xfid, yhid, yfid, zhid, zfid, timeid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*)
      write(*,*) 'reading in spatial variables'
      write(*,*)

      !-------------------------------------------------------------------------!
      ! get xh, yh, zh, zf

            write(*,*) '....xh....'
      if (vcm1.gt.20) then
        xhid = ncvid(ncid, 'xh', rcode)
      else
        xhid = ncvid(ncid, 'ni', rcode)
      endif
      call ncvgt( ncid, xhid, 1,   ni, xh(1:ni), rcode )

      write(*,*) '....xf....'
      if (vcm1.gt.20) then
        xfid = ncvid(ncid, 'xf', rcode)
      else
        xfid = ncvid(ncid, 'nip1', rcode)
      endif
      call ncvgt( ncid, xfid, 1, nip1, xf(1:nip1), rcode )

      if(nj.eq.1)then
        yh(1)=0.5
        yh(2)=1.5
        yf(1)=0.0
        yf(2)=1.0
        yf(3)=2.0
      else

        if (vcm1.gt.20) then

          write(*,*) '....yh....'
          yhid = ncvid(ncid, 'yh', rcode)
          call ncvgt( ncid, yhid, 1, nj,   yh(1:nj),   rcode )

          write(*,*) '....yf....'
          yfid = ncvid(ncid, 'yf', rcode)
          call ncvgt( ncid, yfid, 1, njp1, yf(1:njp1),   rcode )

          write(*,*) '....zh....'
          zhid = ncvid(ncid, 'zh', rcode)
          call ncvgt( ncid, zhid, 1, nk,   zh(1:nk),   rcode )

          write(*,*) '....zf....'
          zfid = ncvid(ncid, 'zf', rcode)
          call ncvgt( ncid, zfid, 1, nk+1, zf(1:nkp1), rcode )

        else

          write(*,*) '....yh....'
          yhid = ncvid(ncid, 'nj', rcode)
          call ncvgt( ncid, yhid, 1, nj,   yh(1:nj),   rcode )

          write(*,*) '....yf....'
          yfid = ncvid(ncid, 'njp1', rcode)
          call ncvgt( ncid, yfid, 1, njp1, yf(1:njp1),   rcode )

          write(*,*) '....zh....'
          zhid = ncvid(ncid, 'nk', rcode)
          call ncvgt( ncid, zhid, 1, nk,   zh(1:nk),   rcode )

          write(*,*) '....zf....'
          zfid = ncvid(ncid, 'nkp1', rcode)
          call ncvgt( ncid, zfid, 1, nk+1, zf(1:nkp1), rcode )

        endif

        ! convert km to m
        xh=xh*1000
        xf=xf*1000
        yh=yh*1000
        yf=yf*1000
        zh=zh*1000
        zf=zf*1000

      endif

      ! calculate mzh and mzf
      write(*,*) '....mzh....'
      do k=1,nk
        mzh(k) = dz/(zf(k+1)-zf(k))
      enddo

      write(*,*) '....mzf....'
      do k=2,nk
        mzf(k) = dz/(zh(k)-zh(k-1))
      enddo
      mzf(1)=1.0
      mzf(nk+1)=1.0

      ! write out variables
      !write(*,*)
      !do i=1,ni
      !  write(*,*) 'i,xh(i),xf(i)=',i,xh(i),xf(i)
      !enddo
      !write(*,*) 'i,xf(i)=',nip1,xf(nip1)
      !write(*,*)
      !do j=1,nj
      !  write(*,*) 'j,yh(j),yf(j)=',j,yh(j),yf(j)
      !enddo
      !write(*,*) 'j,yf(j)=',njp1,yf(njp1)
      !write(*,*)
      !do k=1,nk
      !  write(*,*) 'k,zh(k),zf(k)=',k,zh(k),zf(k)
      !enddo
      !write(*,*) 'k,zh(k),zf(k)=',nk+1,zf(nk+1)
      !write(*,*) '---------------------------------------------------------'

      write(*,*)
      write(*,*) '.......done reading in space'
      write(*,*)

    end subroutine

    !===========================================================================!
    !###########################################################################!
    !#                           SUBROUTINE NC_GETVAR3                         #!
    !#                                                                         #!
    !#  Get a generic 3-dimensional variable from an netCDF file and store it  #!
    !#  in a 3-dimensional variable.                                           #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !###########################################################################!
    subroutine nc_getvar3( ncid, varname, strt, cnt, var )
      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      integer                 ncid
      character(len=10)        varname
      integer,dimension(3) :: strt, cnt

      !------------------------- output variables ------------------------------!
      real,dimension(ni,nj,nk) :: var

      !-------------------------- local variables ------------------------------!
      integer rcode, varid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      varid = ncvid( ncid, varname, rcode )
      call ncvgt( ncid, varid, strt, cnt, var, rcode )

    end subroutine

    !===========================================================================!
    !###########################################################################!
    !#                          SUBROUTINE NC_GETVAR4                          #!
    !#                                                                         #!
    !#  Get a generic 4-dimensional variable and store it in a 3-dimensional   #!
    !# variable.                                                               #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !###########################################################################!
    subroutine nc_getvar4( ncid, varname, strt, cnt, var )
      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      integer                  ncid
      character(len=10)         varname
      integer,dimension(4)  :: strt, cnt

      !------------------------- output variables ------------------------------!
      real,dimension(ni,nj,nk) :: var

      !-------------------------- local variables ------------------------------!
      real,dimension(ni,nj,nk,1) :: dummyvar
      integer                       rcode, varid, i, j, k

      !#########################################################################!
      !############################## MAIN BODY ################################!

      varid = ncvid( ncid, varname, rcode )
      call ncvgt( ncid, varid, strt, cnt, dummyvar, rcode )
      do i=1,ni
      do j=1,nj
      do k=1,nk
        var(i,j,k) = dummyvar(i,j,k,1)
      enddo
      enddo
      enddo

    end subroutine


    !##########################################################################!
    !#                      SUBROUTINE NC_TRAJECT_VARDEF                      #!
    !#                                                                        #!
    !#  Sets up netCDF output file for traj.  Defines all variables.          #!
    !#  Passed variable:                                                      #!
    !#    character outfile : name of output file                             #!
    !#  Returned variables:                                                   #!
    !#    integers that are netCDF variable identifiers                       #!
    !#========================================================================#!
    !# v1.0, Ryan Hastings, 14 Dec 2016                                       #!
    !##########################################################################!
    subroutine nc_traject_vardef( outfile, ncid, ntimes, nparcels, nvars, & 
      outvars, tpid, xpid, ypid, zpid, upid, vpid, wpid, varids )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ----------------------------!
      character(len=200)                 :: outfile
      integer                            :: nparcels, ntimes, nvars
      character(len=10),dimension(nvars) :: outvars
  
      !---------------------------- output variables ---------------------------!
      integer                            :: ncid
      integer                            :: tpid, xpid, ypid, zpid, upid, vpid, wpid
      integer,dimension(nvars)           :: varids

      !---------------------------- local variables -------------------------!
      integer rcode, i
      integer tid, npid
      integer,dimension(2) :: dimid2


     !#########################################################################!
     !############################## MAIN BODY ################################!
     !----------------------------------------------------------!
     ! create netcdf file
     write(*,*) 'creating ',outfile
     ncid = nccre( outfile, NCNOCLOB, rcode )

     !----------------------------------------------------------!
     ! define dimensions
       write(*,*) 'ntimes=',ntimes
       tid  = ncddef( ncid, 'time',    ntimes,   rcode )
       write(*,*) 'nparcels=',nparcels
       npid = ncddef( ncid, 'parcel', ncunlim, rcode )

    dimid2 = (/tid,npid/)

    !----------------------------------------------------------!
    ! define variables
    tpid = ncvdef( ncid, 'tp', nf_real, 1,    tid, rcode )

    ! spacial variables
    xpid = ncvdef( ncid, 'xp', nf_real, 2, dimid2, rcode )
    ypid = ncvdef( ncid, 'yp', nf_real, 2, dimid2, rcode )
    zpid = ncvdef( ncid, 'zp', nf_real, 2, dimid2, rcode )

    ! wind
    upid = ncvdef( ncid, 'up', nf_real, 2, dimid2, rcode )
    vpid = ncvdef( ncid, 'vp', nf_real, 2, dimid2, rcode )
    wpid = ncvdef( ncid, 'wp', nf_real, 2, dimid2, rcode )
 
    ! scalar variables to interpolate
    do i=1,nvars
      varids(i) = ncvdef( ncid, outvars(i), nf_real, 2, dimid2, rcode )
    enddo

    ! end definitions
    call ncendf( ncid, rcode )

  end subroutine

  !############################################################################!
  !#                       SUBROUTINE READ_CM1_WINDS                          #!
  !#                                                                          #!
  !#  Reads winds from CM1 file for traj.                                     #!
  !############################################################################!
  subroutine read_cm1_winds( infile, u, v, w )

    implicit none

    character(len=200),intent(in) :: infile

    real,intent(inout),dimension(1:ni+1,1:nj,0:nk) :: u
    real,intent(inout),dimension(1:ni,1:nj+1,0:nk) :: v
    real,intent(inout),dimension(1:ni,1:nj,1:nk+1) :: w

    integer              :: ncid, rcode, varid
    integer,dimension(4) :: st

    ncid = ncopn( infile, NCNOWRIT, rcode )

    st=(/1,1,1,1/)

    varid = ncvid( ncid, 'u', rcode )
    call ncvgt( ncid, varid, st, (/nip1,nj,nk,1/), u(:,:,1:nk), rcode )
    varid = ncvid( ncid, 'v', rcode )
    call ncvgt( ncid, varid, st, (/ni,njp1,nk,1/), v(:,:,1:nk), rcode )
    varid = ncvid( ncid, 'w', rcode )
    call ncvgt( ncid, varid, st, (/ni,nj,nkp1,1/), w(:,:,:), rcode )

    u(:,:,0)=u(:,:,1)
    v(:,:,0)=v(:,:,1)

  end subroutine
  

  !############################################################################!
  !#                     SUBROUTINE READ_CM1_SCALARS                          #!
  !#                                                                          #!
  !#  Reads scalars from CM1 file for traj.                                   #!
  !############################################################################!
  subroutine read_cm1_scalars( infile, nvars, outvars, varcm1 )

    implicit none

    character(len=200),intent(in) :: infile
    integer,intent(in) :: nvars
    character(len=10),dimension(nvars),intent(in) :: outvars

    real,intent(inout),dimension(1:ni,1:nj,0:nk,nvars) :: varcm1

    integer i, ncid, rcode

    ncid = ncopn( infile, NCNOWRIT, rcode )

    do i=1,nvars
      call nc_readvar4( ncid, outvars(i), varcm1(:,:,:,i) )
    enddo

  end subroutine


  !############################################################################!
  !#                   SUBROUTINE NC_READVAR4                                 #!
  !#                                                                          #!
  !#  Reads 4-D netcdf variable on scalar grid.  This differs from            #!
  !#  CM1_NC_GETVAR4 in that it reads to an extended scalar vertical grid.    #!
  !#  I.e., the lowest grid point is below the ground.  Assumes values remain #!
  !#  constant below lowest level.                                            #!
  !############################################################################!
  subroutine nc_readvar4( ncid, varname, var )

    implicit none

    integer,intent(in) :: ncid
    character(len=10),intent(in) :: varname
    real,intent(out),dimension(ni,nj,0:nk) :: var

    integer rcode, varid
    integer i, j

    varid = ncvid( ncid, varname, rcode )
    call ncvgt( ncid, varid, (/1,1,1,1/), (/ni,nj,nk,1/), var(:,:,1:nk), rcode )

    do i=1,ni
    do j=1,nj
      var(i,j,0) = var(i,j,1)
    enddo
    enddo
 
  end subroutine

  !############################################################################!
  !#                     SUBROUTINE NC_TRAJECT_WRITEOUT                       #!
  !#                                                                          #!
  !#  Writes out variables for traj.                                          #!
  !############################################################################!
  subroutine nc_traject_writeout( ntimes, nparcels, ncid, tpid, tp, &
    xpid, xpc, ypid, ypc, zpid, zpc, upid, up, vpid, vp, wpid, wp, &
    nvars, varids, varp )
    
    implicit none

    integer,intent(in) :: ntimes, nparcels
    integer,intent(in) :: ncid
    integer,intent(in) :: tpid, xpid, ypid, zpid, upid, vpid, wpid
    integer,intent(in) :: nvars
    integer,intent(in),dimension(nvars) :: varids

    real,intent(in),dimension(ntimes) :: tp
    real,intent(in),dimension(ntimes,nparcels) :: xpc, ypc, zpc, up, vp, wp
    real,intent(in),dimension(ntimes,nparcels,nvars) :: varp

    integer,dimension(2) :: s2, c2
    integer rcode, i


    s2 = (/1,1/)
    c2 = (/ntimes,nparcels/)

    call ncvpt( ncid, tpid, 1, ntimes, tp, rcode )
    call ncvpt( ncid, xpid, s2, c2, xpc, rcode )
    call ncvpt( ncid, ypid, s2, c2, ypc, rcode )
    call ncvpt( ncid, zpid, s2, c2, zpc, rcode )

    call ncvpt( ncid, upid, s2, c2, up, rcode )
    call ncvpt( ncid, vpid, s2, c2, vp, rcode )
    call ncvpt( ncid, wpid, s2, c2, wp, rcode )

    do i=1,nvars
      call ncvpt( ncid, varids(i), s2, c2, varp(:,:,i), rcode )
    enddo

    call ncclos( ncid, rcode )

  end subroutine

end module
