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

      character*32 chFilename
      integer nFFTSetup  !variable to determine if setup has been called
      integer nFFToutstep
      data nFFTSetup,nFFToutstep /0,0/ !initialize value to zero
      save nFFTSetup,nFFToutstep     !save value between subsequent calls


      ! 1) Make sure a valid number of processors are present
      if(nFFTp2c.gt.np) then
         if(nid.eq.0)write(6,*),"ERROR FFT: FFTp2c> Total processors"
         call exitt()
      end if

      ! 2) Perform setup procedures
      if(nFFTSetup.eq.0) then
           call FFT_Define_Points()
          ! call FFT_Find_Points()
      end if
      !if(nFFTdmanual.ne.1.or.nFFTsetup.eq.0)then
      !     call FFT_Create_Plan()
      !endif
      nFFTSetup=1
      ! 3) Perform Interpolation
           call FFT_Find_Points()
          ! call FFT_Interp_Points()
      if(nFFTdmanual.ne.1.or.nFFTsetup.eq.0)then
           call FFT_Create_Plan()
      endif
      if(nid.lt.nFFTp2c)then
       write(chFilename,"(A4)")"test"
       call dwritevts(nid,nFFTdims,nFFTflds,rFFTpts,rFFTvals,chFilename)
      endif
      ! 3a) Convert velocity to cylindrical coordinates
           call FFT_Cart2Cyl_Vel()
      ! 4) Perform Transform
           call FFT_Transform()
      ! 5) Destroy Plan
          if(nFFTdmanual.ne.1)then
           call FFT_Destroy_Plan()
          end if
      ! 6) If desired write to file
      !     call FFT_ASCII_PRINT()
           call FFT_OUTPUT_WAVENUMBERS(nFFToutstep)
           call FFT_ENERGY_REPORT(nFFToutstep)
      nFFToutstep=nFFToutstep+1
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
      common /myDomainRange/rMax,zMax,zMin,
     $ rRPnt(nFFTGy),rRWgt(nFFTGy),rZPnt(nFFTGz),rZWgt(nFFTGz)
      real PI

      integer i,j,k,ii,ntot
      integer iG,jG,kG
      real dR,dTheta,dZ,Rval,Tval,Zval,Xval,Yval
      real rMax,zMax,zMin

      ntot=lx1*ly1*lz1*lelv
      zMax=glmax(zm1,ntot)
      zMin=glmin(zm1,ntot)
      rMax=glmax(xm1,ntot)
      PI=4.0*atan(1.0)

      !Determine R locations
      call ZWGJD(rRPnt,rRWgt,nFFTGy,0.,0.)
      !Determin Z locations
      call ZWGJD(rZPnt,rZWgt,nFFTGz,0.,0.)
C      dR=xMax/(nFFTly1*nFFTbly-1)
      dTheta=2.0*PI/(nFFTlx1*nFFTblx)
C      dZ=1.0/(nFFTlz1*nFFTblZ-1)

      ! Good idea to zero out any thing that won't be using an FFT
      if(nid.ge.nFFTp2c) then
       do i=1,nFFTtotal
          rFFTpts(1,i)=0.0
          rFFTpts(2,i)=0.0
          if(if3d) rFFTpts(3,i)=0.0
       end do
      else !Initialize fields for processors of interest
       do i=1,nFFTlz1
          do j=1,nFFTly1
             do k=1,nFFTlx1
                ii=(k-1)+(j-1)*nFFTlx1+(i-1)*nFFTlx1*nFFTly1+1
                call FFT_L2G(i,j,k,iG,jG,kG,nid)
                
                Rval=rMax*(0.5*rRPnt(jG)+0.5)
                Tval=dTheta*(k-1)
                Zval=zMin+(zMax-zMin)*(0.5*rZPnt(iG)+0.5)

                rFFTpts(1,ii)=Rval*cos(Tval)
                rFFTpts(2,ii)=Rval*sin(Tval)
                rFFTpts(3,ii)=ZVal

             end do
          end do
       end do
      end if

      if(nid.eq.0) write(6,*) 'done::FFT points declared'

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
      integer nxyz, nflds, iCalled
      integer ntot
      save iCalled
      data iCalled/0/

      nxyz=nx1*ny1*nz1
      ntot=nxyz*nelt

      nFFTitp_handle=0
      nFlds=0 
      !Map pressure to grid1
      call prepost_map(0)

      !Perform setup on first call
      if(iCalled.eq.0)call intpts_setup(-1.0,nFFTitp_handle)

      ! pack working array
      if(ifvo) then
        call copy(rFFTwrk(1,1),vx,ntot)
        call copy(rFFTwrk(1,2),vy,ntot)
        if(if3d) call copy(rFFTwrk(1,3),vz,ntot)
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
      if(icalled.eq.0)then
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
       icalled=1
      endif
      !Map pressure back to grid2
      !call prepost_map(1)

      if(nid.eq.0) write(6,*) 'done::FFT points found'

      !return
      !end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE FFT INTERP POINTS
C     This is the routine where the actual interpolation takes place
      !subroutine FFT_INTERP_POINTS()

      !include 'SIZE'
      !include 'TOTAL'
      !include 'MYFFT'

      !common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure

      !Map pressure to grid1
      !call prepost_map(0)

      !evaluate field at given points
      do ifld = 1,nFFTflds
         call findpts_eval(nFFTitp_handle,rFFTvals(ifld,1),nFFTflds,
     &                     nFFTrcode,1,
     &                     nFFTproc,1,
     &                     nFFTelid,1,
     &                     rFFTrst,ndim,nFFTtotal,
     &                     rFFTwrk(1,ifld))
         do j=1,nFFTtotal
            cFFTvals(j,ifld)=DCMPLX(rFFTvals(ifld,j),0.d0)
         end do
      enddo

      !Map pressure back to grid2
      call prepost_map(1)
      if(nid.eq.0) write(6,*) 'done::FFT points interpolation'

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
     $  inembed(nFFTorder), onembed(nFFTorder), swap

      !This shares the memory for n and the embedding.  My routine
      !does not function with embedded FFT for padding etc. If you need
      !embedding you will need to modify this routine
      equivalence (n,inembed)
      equivalence (n,onembed)

      !Set up FFT plan parameters
      rank=nFFTorder !order of FFT 1d, 2d, 3d

      ! 1-D FFT parameters----------------------------------------
      if(rank.eq.1)then
        if(bFFTd2t(1))then
           n(1)=nFFTlx1 !size of FFT
           istride=1  !distance between points in memory
           idist=nFFTlx1 !dist between first elements of different FFTs
           howmany=nFFTtotal/nFFTlx1*nFFTflds !total num FFTs
         else if(bFFTd2t(2)) then
           n(1)=nFFTly1
           istride=nFFTlx1
           idist=1
           howmany=nFFTtotal/nFFTly1*nFFTflds
         else if(bFFTd2t(3)) then
           n(1)=nFFTlz1
           istride=nFFTlx1*nFFTly1
           idist=nFFTtotal
           howmany=nFFTtotal/nFFTlz1*nFFTflds
         else
           go to 100
         endif
      else
        go to 100
      end if
      
      ostride=istride
      odist=idist
      if(nid.eq.0)then
        write(6,*)"FFT rank:",rank 
        write(6,*)"FFT n:",n 
        write(6,*)"FFT howmany:",howmany 
        write(6,*)"FFT inembed:",inembed 
        write(6,*)"FFT istride:",istride 
        write(6,*)"FFT idist:",idist 
        write(6,*)"FFT onembed:",onembed 
        write(6,*)"FFT ostride:",ostride 
        write(6,*)"FFT odist:",odist 
       ! write(6,*)"FFT rVals:",rFFTvals
      endif
      call dfftw_plan_many_dft(nFFTplan,rank,n,howmany,cFFTvals,
     $                              inembed,istride,idist,cFFTvals,
     $                              onembed,ostride,odist,
     $                              FFTW_FORWARD,FFTW_ESTIMATE)

      if(nid.eq.0) write(6,*) 'done::FFT plan creation'
      return
 100  if(nid.eq.0) write(6,*) 'ERROR:: unsupported bFFTd2t entry'
      call exitt()
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE PERFORM FFT
      subroutine FFT_TRANSFORM()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      call dfftw_execute_dft(nFFTplan,cFFTvals,cFFTvals)
      if(nid.eq.0) write(6,*) 'done::FFT transform'

      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE DESTROY FFT PLAN
      subroutine FFT_DESTROY_PLAN()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      call dfftw_destroy_plan(nFFTplan)
      if(nid.eq.0) write(6,*) 'done::FFT plan destroyed'
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE PRINT FFT DATA TO TEXT FILE
      subroutine FFT_ASCII_PRINT()

      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      integer nFileNum,nFileErr
      character*32 filename
      data nFileNum /1/
      save nFileNum

      !Each processor will write to the same file one at a time
      do i=1,nFFTp2c
       !if i=my rank then write data
       if(nid.eq.i-1) then
         write(filename,"('myFFT',I0,'.dat')")nFileNum

         !if rank=0 create new file, else open old file
         if(i.eq.1)then
           open(unit=10,file=filename,iostat=nFileErr,status='REPLACE')
         else
           open(unit=10,file=filename,iostat=nFileErr,status='OLD',
     $      access='APPEND')
         end if

         do j=1,nFFTtotal
           write(10,*) nid,GetAngle(rFFTpts(1,j),rFFTpts(2,j)),
     $                        (real(cFFTvals(j,ii)),ii=5,5)
         end do

         close(unit=10)
       endif
       call nekgsync()
      end do
      if(nid.eq.0) write(6,*) 'done::FFT results printed to file'
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     subroutine FFT_OFFSET
C     Determines local to global array offset based on lexigraphical
C     ordering i.e. x+y*nX+z*nX*nY
C     Must use even number of processors divisible for 1d
      subroutine FFT_OFFSET(nMyR,nXoff,nYoff,nZoff)
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'
      !INPUT
      integer nMyR
      !OUTPUT
      integer nXoff,nYoff,nZoff

      nXoff=0
      nYoff=0
      nZoff=0

      if(nFFTblX.gt.1) nXoff=mod(nMyR,nFFTbLX)
      if(nFFTblY.gt.1.and.nFFTblX.gt.1) then
          nYoff=nMyR/nFFTblX
      else
          if(nFFTblY.gt.1) nYoff=mod(nMyR,nFFTblY)
      endif
      if(nFFTblz.gt.1) nZoff=nMyR/(nFFTblX*nFFTblY)

      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE FFT_OUTPUT_WAVENUMBERS()
C     This subroutine is to output the data for each wave number in a
C     seperate file.  Data is collected onto rank0 and written there
      subroutine FFT_OUTPUT_WAVENUMBERS(nFFToutstep)
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      integer,parameter::nInFFT=nFFTlx1
      integer,parameter::nInSlice=nFFTly1*nFFTbly*nFFTlz1*nFFTblZ

      !Define Size of array for wave number
      real dataOutReal(nFFTflds,nInSlice)
      real dataOutComp(nFFTflds,nInSlice)
      real dataWork(nFFTflds,nInSlice)
      !TODO determine parameters for size
      real dataOutPts(3,nInSlice)
      real dataWrkPts(3,nInSlice)
      real theta
      !Define number of wave numbers to output
      integer nMyWave,nFFToutstep
      integer,dimension(3):: nDimension
      character*32 chFileName

      data nDimension /nFFTGx,1,nFFTGz/
      !Loop over the wave numbers
      do i=1,16
        !Zero out workign arrays
        call rzero(dataOutReal,nInSlice*nFFTflds)
        call rzero(dataOutComp,nInSlice*nFFTflds)
        call rzero(dataWork,nInSlice*nFFTflds)
        call rzero(dataOutPts,nInSlice*3)
        call rzero(dataWrkPts,nInSlice*3)
        !populate local spot in the array
        if(nid.lt.nFFTp2c)then
        do j=1,nFFTlz1
          do k=1,nFFTly1
            kg=(k)+mod(nid,nFFTblY)*nFFTly1 
            jg=(j)+(nid/nFFTblY)*nFFTlz1
            iii=(i)+(k-1)*nFFTlx1+(j-1)*nFFTlx1*nFFTly1
            ii=(kg)+(jg-1)*nFFTly1*nFFTbly 
         !   if(nid.eq.1) write(6,*) i,j,k,ii,iii,"INDEX"
            if(i.lt.nInFFT/2)then
               nMyWave=i-1
            else
               nMyWave=nInFFT-i
            endif
            theta=GetAngle(rFFTpts(1,iii),rFFTpts(2,iii))
            !Convert to cylindrical coordinates
            dataOutPts(1,ii)=sqrt(rFFTpts(1,iii)**2+rFFTpts(2,iii)**2)
            dataOutPts(2,ii)=0.d0
            dataOutPts(3,ii)=rFFTpts(3,iii)

            do jj=1,nFFTflds
               dataOutReal(jj,ii)=real(cFFTvals(iii,jj))/dble(nInFFT)
               dataOutComp(jj,ii)=dimag(cFFTvals(iii,jj))/dble(nInFFT)
            end do

          enddo
        enddo
        endif
        !gather procedure
        call gop(dataOutPts,dataWrkPts,'+  ',3*nInSlice)
        call gop(dataOutReal,dataWork,'+  ',nInSlice*nFFTflds)
        call gop(dataOutComp,dataWork,'+  ',nInSlice*nFFTflds)
        !write data to file on rank0
         write(chFileName,"('WAVE_',I0,'_TSTEP')")
     $                         nMyWave
        if(nid.eq.0) call dwritevtsc(nFFToutstep,nDimension,nFFTflds,
     $     dataOutPts,
     $                    dataOutReal,dataOutComp,chFileName)
      end do
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine FFT_Cart2Cyl_Vel()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'
      
      real locTheta,velX,velY
      
      do i=1,nFFTtotal
         velX=cFFTvals(i,1)
         velY=cFFTvals(i,2)
             locTheta=getangle(rFFTpts(1,i),rFFTpts(2,i))
         cFFTvals(i,1)=velX*cos(locTheta)+velY*sin(locTheta)
         cFFTvals(i,2)=-velX*sin(locTheta)+velY*cos(locTheta)
      end do
  
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function GetAngle(x,y)
      real x,y
         if(x.gt.0.0)then
            GetAngle=atan(y/x)
         else
            if(x.lt.0.0)then
              GetAngle=atan(y/x)+4.0*atan(1.0)
            else
               if(y.ge.0.0)then
                  GetAngle=2.0*atan(1.0)
               else
                  GetAngle=-2.0*atan(1.0)
               end if
            endif
         endif
      return 
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine FFT_ENERGY_REPORT(nFFToutstep)
C     OUTPUT THE ENERGY FOR EACH WAVE NUMBER INEGRATED OVER THE R-Z
C     PLANE USEING NUMERICAL INTEGRATION
C      ---> \int_0^R\int_0^H A(r,z)*A(r,z) r dr dz (*is complex conj)
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'
      common /myDomainRange/rMax,zMax,zMin,
     $ rRPnt(nFFTGy),rRWgt(nFFTGy),rZPnt(nFFTGz),rZWgt(nFFTGz)
      real EnergyWave(nFFTflds,nFFTlx1),EnergyWork(nFFTflds,nFFTlx1)
      real radius
      integer nFFToutstep,ii,iii,iG,jG,kG
      character*80 fileName
      call rzero(EnergyWave,nFFTflds*nFFTlx1)
      do i=1,nFFTlz1
      do ii=1,nFFTly1
        do j=1,nFFTlx1
           do k=1,nFFTflds
              call FFT_L2G(i,ii,j,iG,jG,kG,nid)
              iii=(ii-1)+(i-1)*nFFTly1
              radius=sqrt(rFFTpts(1,j+iii*nFFTlx1)**2+
     $             rFFTpts(2,j+iii*nFFTlx1)**2)

C             Ek=sum_z sum_r phi(r,z)*CC[phi(r,z)] *wr*R/2*wz*H/2
C             Divide by nFFTlx1 to normalize the Fourier Coefficient
              EnergyWave(k,j)=EnergyWave(k,j)+
     $         real(dconjg(cFFTvals(j+iii*nFFTlx1,k)/dble(nFFTlx1))
     $               *(cFFTvals(j+iii*nFFTlx1,k)/dble(nFFTlx1)))
     $               *radius*rRWgt(jG)*rZWgt(iG)
     $               *0.25*rMax*(zMax-zMin)
           end do
        end do
      end do
      end do
      call gop(EnergyWave,EnergyWork,'+  ',nFFTflds*nFFTlx1)

      write(fileName,"('EnergyReport_',I0,'.dat')")nFFToutstep
      if(nid.eq.0)then
        open(unit=10,file=fileName,status='REPLACE')
        do i =1,nFFTlx1/2
          write(10,*)i-1,(EnergyWave(k,i),k=1,nFFTflds)
        end do
        do i =nFFTlx1/2+1,nFFTlx1
          write(10,*)i-(nFFTlx1+1),(EnergyWave(k,i),k=1,nFFTflds)
        end do
        close(unit=10) 
      end if
      return
      end 
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine FFT_L2G(iL,jL,kL,iG,jG,kG,me)
C     MAP LOCAL INDEX IN GRID TO GLOBAL INDEX
C     K-Innermost loop corresponding to nFFTlx1
C     J-Middle loop corresponding to nFFTly1
C     I-Outermost loop corresponding to nFFTlz1
C     me-mpi rank
      include 'SIZE'
      include 'MYFFT'
      integer iL,jL,kL,iG,jG,kG,me
c      write(6,*)'INSIDE FFT_L2G',nFFTbly,me
c      write(6,*) iL,jL,kL,iG,jG,kG,me
      iG=iL+(me/nFFTblY)*nFFTlz1
      jG=jL+mod(me,nFFTblY)*nFFTly1
      kG=kL
c      write(6,*) iL,jL,kL,iG,jG,kG,me
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function CGL_POINT(N,I)
C     Chebyshev-Gauss-Lobatto quadrature point
      !N is total number of points (N-1 polynomial)
      !I is current point in series from 1:N
      integer N,I
      CGL_POINT=-cos(4.0*atan(1.0)*dble(I-1)/dble(N-1))
      return      
      end 
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function CGL_WEIGHT(N,I)
C     Chebyshev-Gauss-Lobatto quadrature weight
C     i from 1:N
      integer N,I
      if(I.eq.1.or.I.eq.N)then
        CGL_WEIGHT=4.0*atan(1.0)*0.5/dble(N-1)
      else
        CGL_WEIGHT=4.0*atan(1.0)/dble(N-1)
      end if
      CGL_WEIGHT=CGL_WEIGHT*sqrt(1.0-CGL_POINT(N,I)**2)
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function CG_POINT(N,I)
C     Chebyshev-Gauss quadrature point
C     i from 1:N
      integer N,I
      PI=4.0*atan(1.0)
      CG_POINT=-cos((2.0*i-1)*PI/(2.0*dble(N)))
      return
      end 
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function CG_WEIGHT(N,I)
      integer N,I
      CG_WEIGHT=4.0*atan(1.0)/dble(N)*sqrt(1.0-CG_POINT(N,I)**2)
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function QuadPnt(N,I)
      integer N,I
      QuadPnt=CG_POINT(N,I)
      !QuadPnt=CGL_POINT(N,I)
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function QuadWgt(N,I)
      integer N,I
      QuadWgt=CG_Weight(N)
      !QuadWgt=CGL_Weight(N,I)
      return
      end
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
