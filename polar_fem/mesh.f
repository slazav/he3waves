
C =====================================================================
      subroutine draw_mesh(filename, desc)
        include 'th.fh'
        character*(*) filename, desc
        Write(*,*) '   Writing ', desc,' mesh into ', filename
        Write(*,*) nt, ' triangles and ', nv, 'vertices.'
        Write(*,*) nc, ' curves ', nb, 'boundary edges'
        call graph(nv,vrt, nt,tri, filename)
        return
      end

C ====== make initial mesh from cell dimensions
      subroutine create_mesh(Lx,Ly,R)
        include 'th.fh'
        Integer   aft2dfront
        EXTERNAL  aft2dfront
        double precision vbr(2,nbmax), tmp(2)
        real*8  Lx,Ly,R
        integer nbr, dummy

        integer  ipIRE,ipWork, MaxWiWork, iERR
        DATA nvfix/0/, fixedV/0/,  nbfix/0/, ntfix/0/, nc/0/

        Write (*,*) 'Creating mesh.'

        nbr=1

        vbr(1,nbr) = 0D0
        vbr(2,nbr) = 0D0

        nbr=nbr+1
        vbr(1,nbr) = 0D0
        vbr(2,nbr) = Ly

        nbr=nbr+1
        vbr(1,nbr) = Lx
        vbr(2,nbr) = Ly

        nbr=nbr+1
        vbr(1,nbr) = Lx
        vbr(2,nbr) = 0D0

        nbr=nbr+1
        vbr(1,nbr) = R
        vbr(2,nbr) = 0D0

        nbr=nbr+1
        vbr(1,nbr) = 0D0
        vbr(2,nbr) = 0D0

C Generate a mesh
        Write(*,5101) Nbr
        iERR=aft2dfront(
     &           0, dummy, nbr, vbr,   ! segment data
     &           nv, vrt,              ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
        If (iERR.ne.0) stop ' error in function aft2dfront'

        labelB(1) = 5     ! center line, a'=0
        labelB(2) = 4     ! upper edge, a=0
        labelB(3) = 4     ! right edge, a=0
        labelB(4) = 4     ! cut, a=0
        labelB(5) = 3     ! a=pi/2

        nc = 0

        return

c =======
 5101   format(
     &  '  The initial front has ', I4, ' edges.')

      end

C =====================================================================
      subroutine refine_mesh()
        include 'th.fh'
        Integer  nEStar, iERR
        Integer  control(6)
        Real*8   Quality
        Integer  metric_func
        External metric_func

c === generate adaptive mesh
        nEStar  = 200      !  desired number of triangles (not used?)
        control(1) = nEStar/10   !  MaxSkipE
        control(2) = nEStar*10   !  MaxQItr
        control(3) = 1       !  status
        control(4) = 1       !  flagAuto
        control(5) = 0       !  iPrint:   average level of output information
        control(6) = 0       !  iErrMesgt: only critical termination allowed

        Quality = 1D0      !  request shape-regular triangles in metric

      call mbaAnalytic(
     &      nv, nvfix, nvmax, vrt, labelV, fixedV,
     &      nb, nbfix, nbmax, bnd, labelB, fixedB,
     &      nc,               crv, labelC, 0,
     &      nt, ntfix, ntmax, tri, labelT, fixedT,
     &      nEStar, Quality, control, metric_func,
     &      MaxWr, MaxWi, rW, iW, iERR)
      write (*,*) 'quality after refining: ', Quality

      end
C =====================================================================
      Integer Function metric_func(x, y, M)
C =====================================================================
C  This routine creates a metric at the given point (x,y). The
C  metric is a 2x2 positive definite symmetric tensor.
C  Only the upper triangular part of 2x2 array Metric may be defined.
C =====================================================================
      include 'th.fh'
      Real*8  x, y, M(2, 2)

c      Metric(1,1) = (D/2.0-x)*(L/2.0-y)/D/L*4
        M(1,1) = 1
        M(1,2) = 0
        M(2,1) = 0
        M(2,2) = M(1,1)
        metric_func = 0
      End

C =====================================================================
      subroutine adapt_mesh(Quality, F)
        include 'th.fh'

        Real*8   F(nfmax)

        Real*8   Lp
        Real*8   Metric(3,nvmax)

        Integer  control(6), nEStar, iERR
        Real*8   Quality

c  ===  generate metric (from SOL) optimal for the L_p norm
c       Lp = 0             ! maximum norm
        Lp = 1             ! L_1 norm
        Call Nodal2MetricVAR(F,
     &          vrt, nv, tri, nt, bnd, nb, Metric,
     &          MaxWr, rW, MaxWi, iW)

        If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


        Write(*,*) 'generate the adaptive mesh...'
c === generate the adaptive mesh to u
         nEStar = TRI_NUM
         control(1) = nEStar/10  ! MaxSkipE
         control(2) = nEStar*10  ! MaxQItr
         control(3) = 32+1    ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 0       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         Quality = 0.6

         Call mbaNodal(
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, 0,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
     &        nEStar, Quality, control, Metric,
     &        MaxWr, MaxWi, rW, iW, iERR)

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error if mbaNodal')
      write (*,*) 'mesh quality: ', Quality

        return
      end
