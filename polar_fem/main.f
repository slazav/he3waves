c ======================================================================
      Program Th
c ======================================================================
      include 'th.fh'

      Integer  iLoop, nLOOPs, i
      Real*8   Quality

c ======================================================================
c number of adaptive loops
      nLOOPs = 5

      call read_cfg
      call create_mesh(4D0,4D0,1.8D0)
      Call draw_mesh('ps/mesh0.ps', 'initial')
      call refine_mesh
      Call draw_mesh('ps/mesh1.ps', 'refined')

c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> U LOOP: ', iLoop
         call solve_u
         If(iLoop.ne.nLOOPs) call adapt_mesh(Quality, SOL_U)
      End do

      Call draw_mesh('ps/mesh2.ps', 'final')
      Call draw_u('ps/sol_u.ps')
      call write_zlines(SOL_U, 'u')


      Stop
      End


c ======================================================================


      subroutine write_line(x1,x2,y1,y2,np,u, fd)
        include 'th.fh'
        real*8 x1,x2,y1,y2,dx,dy
        real*8 xy(2,np), out(np)
        real*8 u(*)
        integer np,info(3)
        integer i,fd

        do i=1,np
          xy(1,i)=x1 + (x2-x1)*(i-1)/(np-1)
          xy(2,i)=y1 + (y2-y1)*(i-1)/(np-1)
        end do
        info(1)=1
        info(2)=1
        call LINTRP(nt,tri,nv,vrt,1,u, np,xy, out, iW,rW, info)
        do i=1,np
          write (fd,*) xy(1,i), xy(2,i), out(i)
        enddo
      end




      subroutine write_zlines(u,par)
        include 'th.fh'
        real*8 x(5),y1,y2
        double precision u(*)
        character*(64) filename
        character*(*) par

        integer i, np

        x(1)=0D0
        x(2)=DIM_Dd*0.25
        x(3)=DIM_Dd*0.40
        x(4)=DIM_Dd*0.48
        x(5)=DIM_Dd*0.50
        y1=0D0
        y2=0.005D0
        np=500

        do i=1,5
          write(filename,'(A,A,I1,A)') 'dat/solz_', par, i,'.dat'
          open (57, FILE=filename)
          call write_line(x(i),x(i),y1,y2, np, u, 57)
          close(57)
        enddo
      end

      subroutine write_rlines(u,par)
        include 'th.fh'
        real*8 y(5),x1,x2
        double precision u(*)
        character*(64) filename
        character*(*) par

        integer i, np

        y(1)=0D0
        y(2)=0.0005D0
        y(3)=0.001D0
        y(4)=0.002D0
        y(5)=0.004D0
        x1=0D0
        x2=0.002D0
        np=500

        do i=1,5
          write(filename,'(A,A,I1,A)') 'dat/solr_', par, i,'.dat'
          open (57, FILE=filename)
          call write_line(x1,x2,y(i),y(i), np, u, 57)
          close(57)
        enddo
      end

      subroutine write_res(nx,ny)
        include 'th.fh'
        integer i,j, nx,ny
        real*8 x,y

        open (57, FILE='dat/sol.dat')
        write (57,*) '## r z T'
        do i=0,nx
          x=(DIM_Dd/2.0D0*i)/nx
          call write_line(x,x,0,DIM_L/2.0D0-DIM_Lw, ny, SOL_N, 57)
          write (57,*) ''
        enddo
        close(57)

      end
