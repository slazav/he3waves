      Subroutine draw_u(filename)
        include 'th.fh'
        character*(*) filename
        Call isolines(SOL_U, nv,vrt, nt,tri, nb,bnd, filename, 50)
      end

C ======================================================================
      Subroutine solve_u
        include 'th.fh'
        include 'fem2Dtri.fd'
        include 'assemble.fd'

        Integer  IA(nfmax), JA(namax)
        Real*8    A(namax), RHS(nfmax), RES(nfmax)

        Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3)
        Real*8   dDATAFEM(1)

        Integer  Dbc_u
        External FEM2Dext_u, Dbc_u

        Integer  symbolic(2), numeric(2), sys
        Real*8   lucontrol(20), luinfo(90)

        Integer  i, nRow, nCol
        Real*8   rmax, h, x, y, eBC(2)

        Write(*,*) 'finding U'

c === no data is provided for the user subroutine Ddiff
        dDATAFEM(1) = 0D0
        iDATAFEM(1) = 0

c mark the Dirichlet points with the maximal edge color
        Call markDIR(nv, vrt, labelV, nb, bnd, labelB,
     &                Dbc_u, dDATAFEM, iDATAFEM, iSYS)

c === general sparse matrix in a 0-based CSC format used in UMFPACK
        controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
        controlFEM(2) = 1 ! verbosity
 
        Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext_u, dDATAFEM, iDATAFEM, controlFEM,
     &        nfmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c  ===   call the driver for LU factorization and solution
        Call CSC2CSC0(nCol, IA, JA)

c set up default control parameters & print only error messages
        Call umf4def(lucontrol)
        lucontrol(1) = 1

c pre-order and symbolic analysis
        Call umf4sym(nCol, nCol, IA,JA,A, symbolic,lucontrol,luinfo)
        If(luinfo(1).LT.0) Goto 5001

c numeric factorization
        Call umf4num(IA,JA,A, symbolic,numeric,lucontrol,luinfo)
        If(luinfo(1).LT.0) Goto 5002

c free the symbolic analysis data
        Call umf4fsym(symbolic)

c solve Ax=b, without iterative refinement
        sys = 0
        Call umf4sol(sys, SOL_U, RHS, numeric, lucontrol,luinfo)
        If(luinfo(1).LT.0) Goto 5003

c free the numeric factorization data
        Call umf4fnum (numeric)

c check the residual
        Call mulAcsc0(nRow, IA, JA, A, SOL_U, RES)
        rmax = 0
        Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
        End do
        Write(*,'(A,E12.6)') 'LU:  maximal norm of residual: ', rmax

         return
c error messages
 5001   Continue
        Write(*,*) 'Error occurred in umf4sym: ', luinfo(1)
        Stop 911

 5002   Continue
        Write(*,*) 'Error occurred in umf4num: ', luinfo(1)
        Stop 911

 5003   Continue
        Write(*,*) 'Error occurred in umf4sol: ', luinfo(1)
        Stop 911


      end


C ======================================================================
      Subroutine FEM2Dext_u(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Implicit none
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*), LDA, nRow, nCol

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff_u, Drhs_u, Dbc_u
      External Ddiff_u, Drhs_u, Dbc_u

      Integer  i,k, ir,ic, label, ibc
      Real*8   x, y, eBC(1)
      Logical  ifXbc

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates 
      Do i = 1, 3
         templateR(i) = Vdof
         templateC(i) = Vdof
      End do

c ... compute the stiffness matrix A
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff_u, dDATA, iDATA, iSYS, 1,
     &              LDA, A, ir, ic)

c ... compute right hand side F
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs_u, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)

c ... impose boundary conditions (assume nRow = nCol)
      Do k = 1, 3
         If(lbP(k).GT.0) Then
            If(k.EQ.1) Then
               x = XY1(1)
               y = XY1(2)
            ElseIf(k.EQ.2) Then
               x = XY2(1)
               y = XY2(2) 
            ElseIf(k.EQ.3) Then
               x = XY3(1)
               y = XY3(2)
            End if

            ibc = Dbc_u(x, y, lbP(k), dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC)
            End if
         End if
      End do

      Return
      End



C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff_u(x, y, label, dDATA, iDATA, iSYS, Coef)
      Include 'th.fh'
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, 4)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      if (labelT(iSYS(3)).eq.1) then
        Coef(1,1) = PHYS_SIGMA * x
      else
        Coef(1,1) = PHYS_SIGMAd * x
      endif

      Ddiff_u = TENSOR_SCALAR

      Return
      End



C ======================================================================
C Boundary condition
C ======================================================================
      Integer Function Dbc_u(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'th.fh'
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If (label.EQ.2) Then
         Dbc_u = BC_DIRICHLET
         eBC(1,1) = 0D0
      Else If (label.EQ.3) Then
         Dbc_u = BC_DIRICHLET
         eBC(1,1) = PHYS_U0/2.0D0
      Else If (label.EQ.1) Then
         Dbc_u = BC_NULL
         eBC(1,1) = 0D0
      Else
         Dbc_u = BC_NEUMANN
         eBC(1,1) = 0D0
      End if

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs_u(x, y, label, dDATA, iDATA, iSYS, F)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, F(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      F(1, 1) = 0D0
      Drhs_u = TENSOR_SCALAR

      Return
      End
