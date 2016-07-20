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
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE MYFFT
      subroutine MyFFT()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'


      integer nFFTSetup  !variable to determine if setup has been called
      data nFFTSetup /0/ !initialize value to zero
      save nFFTSetup     !save value between subsequent calls


      ! 1) Make sure a valid number of processors are present
      if(nFFTp2c>np) then
         if(nid.eq.0)write(6,*),"ERROR FFT: FFTp2c> Total processors"
         call exitt()
      end if

      ! 2) Perform setup procedures
      if(nFFTSetup.eq.0) then
           call FFT_Define_Points()
           call FFT_Find_Points()
      end if
      if(nFFTdauto.ne.1.or.nFFTsetup.eq.0)then
           call FFT_Create_Plan()
      endif
      nFFTSetup=1

      ! 3) Perform Interpolation
           call FFT_Interp_Points()
      ! 4) Perform Transform
           call FFT_Transform()
      ! 5) Destroy Plan
          if(nFFTdeauto.ne.1)then
           call FFT_Destroy_Plan()
          end if
      ! 6) If desired write to file

      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE DEFINE POINTS
C     Users should use this to define the sampling points they want for
C     their FFT's. This example will be for points in a cylinder with
C     FFT in the theta direction
      subroutine FFT_Define_Points()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      real PI

      integer i,j,k,ii

      real dX, dY

      PI=4.0*atan(1.0)
      dX=2.0*PI/nFFTlx1
      dY=2.0*PI/nFFTly1

      ! Good idea to zero out any thing that won't be using an FFT
      if(nid.gt.nFFTp2c) then
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
                rFFTpts(1,ii)=(i-1)*dX
                rFFTpts(2,ii)=(j-1)*dY
                if(if3d) rFFTpts(3,ii)=(k-1)*dZ
             end do
          end do
       end do
      end if

      if(nio.eq.0) write(6,*) 'done::FFT points declared'

      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE FIND POINTS
C     Find the points that will be used for the FFT and interpolate them
C     to the array rFFTvals.  Note that rFFTvals is type real, and when
C     the FFT is performed the complex values will be held in cFFTvals.
C     DO NOT confuse rFFTvals and cFFTvals
      subroutine FFT_FIND_POINTS()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'
      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure
      integer nxyz, nflds
      integer ntot

      nxyz=nx1*ny1*nz1
      ntot=nxyz*nelt

      nFFTitp_handle=0
      !Map pressure to grid1
      call prepost_map(0)

      !Perform setup on first call
      call intpts_setup(-1.0,nFFTitp_handle)

      ! pack working array
      if(ifvo) then
        call copy(rFFTwrk(1,ndim),vx,ntot)
        call copy(rFFTwrk(1,ndim),vy,ntot)
        if(if3d) call copy(rFFTwrk(1,ndim),vz,ntot)
        nflds = ndim
      endif

      if(ifpo) then
        nflds = nflds + 1
        call copy(rFFTwrk(1,nFFTflds),pm1,ntot)
      endif

      if(ifto) then
        nflds = nflds + 1
        call copy(rFFTwrk(1,nFFTflds),t,ntot)
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

      !Map pressure back to grid2
      call prepost_map(1)

      if(nio.eq.0) write(6,*) 'done::FFT points found'

      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

      integer rank,n(nFFTorder),howmany, idist, odist, istride, ostride,
     $  inembed(nFFTorder), enembed(nFFTorder), swap

      !This shares the memory for n and the embedding.  My routine
      !does not function with embedded FFT for padding etc. If you need
      !embedding you will need to modify this routine
      equivalence (n,inembed)
      equivalence (n,enembed)

      !Set up FFT plan parameters
      rank=nFFTorder !order of FFT 1d, 2d, 3d

      !Sort nFFTd2t to make sure it is ordered lowest to highest
      if (rank.gt.1) then
        do i=1,rank-1
            if(n(i).gt.nFFTd2t(i+1))then
               swap=nFFTd2t(i+1)
               nFFTd2t(i+1)=nFFTd2t(i)
               nFFTd2t(i)=swap
            end if
        end do
      end if


      ! 1-D FFT parameters----------------------------------------
      if(rank.eq.1)then
        if(nFFTd2t(1).eq.1)then
           n(1)=nFFTlx1 !size of FFT
           istride=1  !distance between points in memory
           idist=nFFTlx1 !dist between first elements of different FFTs
           howmany=nFFTtotal/nFFTlx1*nFFTflds !total num FFTs
         else if(nFFTd2t(1).eq.2) then
           n(1)=nFFTly1
           istride=nFFTlx1
           idist=1
           howmany=nFFTtotal/nFFTly1*nFFTflds
         else if(if3d.and.nFFTd2t(1).eq.3) then
           n(1)=nFFTlz1
           istride=nFFTlx1*nFFTly1
           idist=1
           howmany=nFFTtotal/nFFTlz1*nFFTflds
         else
           go to 100
         endif
      ! 2-D FFT parameters----------------------------------------
      else if(rank.eq.2)then
        if(nFFTd2t(1).eq.1) then
           n(1)=nFFTlx1
           istride=1
           if(nFFTdt2(2).eq.2)then
           ! FFT in dimensions 1 and 2 of pts array
              n(2)=nFFTly1
              idist=nFFTlx1*nFFTly1
              howmany=nFFTlz1*nFFTflds
           else
           ! FFT in dimensions 1 and 3 of pts array
           ! currently not sure how this would work since
           ! it would require 2 different strides. Therefore,
           ! I choose not to support it
              go to 100
           endif
        else
        ! FFT in dimensions 2 and 3 of pts array
          n(1)=nFFTly1
          n(2)=nFFTlz1
          idist=1
          istride=nFFTlx1
          howmany=nFFTlx1*nFFTflds
        endif
      else if (rank.eq.3)then
      ! 3-D FFT parameters----------------------------------------
        n(1)=nFFTlx1
        n(2)=nFFTly1
        n(3)=nFFTlz1
        idist=1
        istride=1
        howmany=nFFTflds
      else
        go to 100
      end if

      ostride=istride
      odist=idist

      call dfftw_plan_many_dft_r2c(nFFTplan,rank,n,howmany,rFFTvals,
     $                              inembed,istride,idist,cFFTvals,
     $                              onembed,ostride,odist,FFTW_ESTIMATE)

      if(nio.eq.0) write(6,*) 'done::FFT plan creation'
      return
 100  if(nid.eq.0) write(6,*) 'ERROR:: unsupported nFFTd2t entry'
      call exitt()
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE PERFORM FFT
      subroutine FFT_TRANSFORM()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      call dfftw_execute_(nFFTplan,rFFTvals,cFFTvals)

      if(nio.eq.0) write(6,*) 'done::FFT transform'

      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE DESTROY FFT PLAN
      subroutine FFT_DESTROY_PLAN()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      call dfftw_destroy_plan_(nFFTplan)
      if(nio.eq.0) write(6,*) 'done::FFT plan destroyed'
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE PRINT FFT DATA TO TEXT
