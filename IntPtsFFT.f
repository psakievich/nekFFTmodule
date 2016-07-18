C***********************************************************************
C     ROUTINES FOR PERFORMING FFT's INSIDE NEK5000
C     REQUIRES FFTW (www.fftw.org), BUILT AND TESTED WITH v3.3.4
C
C     CODE WRITEN BY PHIL SAKIEVICH (psakievi@asu.edu)
C
C     UTILIZES DEFINITIONS IN fftw3.f by Matteo Frigo and Steven Johnson
C      (http://people.sc.fsu.edu/~jburkardt/f77_src/fftw3/fftw3.html)
C
C     INCLUDE FILE:  fftw3.f
C     LINK LIBS:    -lfftw3 -lw
C
C     IMPORTANT DOCUMENTATION:
C     http://www.fftw.org/fftw3_doc/Calling-FFTW-from-Legacy-Fortran.html#Calling-FFTW-from-Legacy-Fortran
C***********************************************************************

C     OVERVIEW:
C     These subroutines are designed to allow one to set up a set of
C     points on a given processor get their values from somewhere in the
C     parallel envirnoment using intpts and then perform FFT's on them
C     locally using routines from fftw3.3.4.
C
C     It is the users responsibility to ensure that the:
C
C     1) The way the points are defined are compatible with the FFT they
C        intend to perform.

C     SUBROUTINE MYFFT
      subroutine MyFFT()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'


      integer nFFTSetup  !variable to determine if setup has been called
      data nFFTSetup /0/ !initialize value to zero
      save nFFTSetup     !save value between subsequent calls

      ! 1) Make sure a valid number of processors are present
      if(nFFTProcs<np) then
         if(nid.eq.0)write(6,*),"ERROR FFT: FFT Procs< Total processors"
         call exitt()
      end if

      ! 2) Perform setup procedures
      if(nFFTSetup.eq.0) then
           call FFT_Define_Points()
           call FFT_Find_Points()
           call FFT_Create_Plan()
           nFFTSetup=1
      end if

      ! 3) Perform Interpolation
           call FFT_Interp_Points()
      ! 4) If desired write to file

      return
      end

C     SUBROUTINE DEFINE POINTS
C     Users should use this to define the sampling points they want for
C     their FFT's. This example will be for points in a cylinder with
C     FFT in the theta direction
      subroutine FFT_Define_Points()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      parameter,real::PI=4.0*atan(1.0)
      integer i,j,k,ii

      real dR=1.0/nFFTlx1, dTheta=2.0*PI/nFFTly1, dZ=1.0/nFFTlz1

      ! Good idea to zero out any thing that won't be using an FFT
      if(nid.gt.nFFTprocs) then
       do i=1,nFFTtotal
          rFFTpts(1,i)=0.0
          rFFTpts(2,i)=0.0
          if(if3d) rFFTpts(3,1)=0.0
       end do
      else !Initialize fields for processors of interest
       do i=1,nFFTlx1
          do j=1,nFFTly1
             do k=1,nFFTlz1
                ii=k+j*nFFTlz1+i*nFFTlz1*nFFTly1
                rFFTpts(1,ii)=(i-1)*dR
                rFFTpts(2,ii)=(j-1)*dTheta
                if(if3d) rFFTpts(3,ii)=(k-1)*dZ
             end do
          end do
       end do
      end if

      if(nio.eq.0) write(6,*) 'done::FFT points declared'

      return
      end

C     SUBROUTINE FIND POINTS
C     Find the points that will be used for the FFT and interpolate them
C     to the array rFFTvals.  Note that rFFTvals is type real, and when
C     the FFT is performed the complex values will be held in cFFTvals.
C     DO NOT confuse rFFTvals and cFFTvals
      subroutine FFT_FIND_POINTS()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      !integer nFFTsetup !Integer for if initialization been performed
      !integer inth_hpts !Handle for intpts routines
      integer nxyz
      integer ntot
      !save    inth_hpts
      data nxyz,ntot /nx1*ny1*nz1,nx1*ny1*nz1*nelt/

      nFFTitp_handle=0

      !Perform setup on first call
      call inpts_setup(-1.0,nFFTitp_handle)

      ! pack working array
      if(ifvo) then
        call copy(rFFTwrk(1,ndim),vx,ntot)
        call copy(rFFTwrk(1,ndim),vy,ntot)
        if(if3d) call copy(rFFTwrk(1,ndim),vz,ntot)
        nflds = ndim
      endif
      if(ifpo) then
        nflds = nflds + 1
        call copy(rFFTwrk(1,nflds),pm1,ntot)
      endif
      if(ifto) then
        nflds = nflds + 1
        call copy(rFFTwrk(1,nflds),t,ntot)
      endif
      do i = 1,ldimt
         if(ifpsco(i)) then
           nflds = nflds + 1
           call copy(rFFTwrk(1,nflds),T(1,1,1,1,i+1),ntot)
         endif
      enddo

      !find points
      call findpts(nFFTitp_handle,nFFTrcode,1,
     &                 nFFTproc,1,
     &                 nFFTelid,1,
     &                 rFFTrst,ndim,
     &                 rFFTdist,1,
     &                 rFFTpts(1,1),ndim,
     &                 rFFTpts(2,1),ndim,
     &                 rFFTpts(3,1),ndim,nFFTtotal)
      !check return codes
      do i=1,nFFTtotal
           ! check return code
           if(nFFTrcode(i).eq.1) then
             if (rFFTdist(i).gt.1e-12) then
                nfail = nfail + 1
                IF (NFAIL.LE.5) WRITE(6,'(a,1p4e15.7)')
     &     ' WARNING: point on boundary or outside the mesh xy[z]d^2:'
     &     ,(rFFTpts(k,i),k=1,ndim),rFFTdist(i)
             endif
           elseif(nFFTrcode(i).eq.2) then
             nfail = nfail + 1
             if (nfail.le.5) write(6,'(a,1p3e15.7)')
     &        ' WARNING: point not within mesh xy[z]: !',
     &        (rFFTpts(k,i),k=1,ndim)
           endif
        enddo

      if(nio.eq.0) write(6,*) 'done::FFT points found'

      return
      end

C     SUBROUTINE FFT INTERP POINTS
C     This is the routine where the actual interpolation takes place
      subroutine FFT_INTERP_POINTS()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure

      !Map pressure to grid1
      call prepost_map(0)

      !evaluate field at given points
      do ifld = 1,nFFTflds
         call findpts_eval(nFFTitp_handle,rFFTvals(ifld,1),nFFTflds,
     &                     nFFTrcode,1,
     &                     nFFTproc,1,
     &                     nFFTelid,1,
     &                     rFFTrst,ndim,nFFTtotal,
     &                     rFFTwrk(1,ifld))
      enddo

      !Map pressure back to grid2
      call prepost_map(1)
      if(nio.eq.0) write(6,*) 'done::FFT points interpolation'

      return
      end
C     SUBROUTINE FFT CREATE PLAN
C     FFTW requires plans for performing FFT's to be generated. These
C     must also be destroyed later.  Documenation can be found in the
C     fftw resources online.
      subroutine FFT_CREATE_PLAN()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      ! I will use one of the most general plans fftw_plan_many_dft_r2c
      ! it requires the following parameters:


      return
      end
C     SUBROUTINE PERFORM FFT
C     SUBROUTINE PRINT FFT DATA TO TEXT
