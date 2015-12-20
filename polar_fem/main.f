c ======================================================================
      Program Th
c ======================================================================
      include 'th.fh'

      Integer  i, nLoops
      Real*8   Quality

c ======================================================================
c number of adaptive loops

        call create_mesh(8D0,8D0,4D0, 500)
        call draw_mesh('ps/mesh0.ps', 'initial')

        ! begin adaptive iterative loop
        nLoops = 10
        do i = 1, nLoops
          write(*,'(/,A,I2)') '===> U LOOP: ', i
          call solve_u
          if(i.ne.nLoops) call adapt_mesh(Quality, SOL_U, 1000)
        end do

        call draw_mesh('ps/mesh2.ps', 'final')
        call draw_u('ps/sol_u.ps')

      End


c ======================================================================


