      subroutine qp ( xn, n, zsmall, lenz, inform ) 

      implicit           double precision (a-h,o-z)
      integer            n, lenz, inform
      double precision   xn(n+1), zsmall(lenz)
*     ------------------------------------------------------------------
      integer            nwcore
      parameter          (nwcore = 10000000)
*     ------------------------------------------------------------------
      integer            m, nb, ne, nname,
     &                   nncon, nnobj, nnjac, iobj,  
     &                   mincor, ns, ninf, 
     &                   ka(n+1), name1, name2,
     &                   iprint, isumm, ispecs, i
      integer*4          ha(n), hs(n+1)
      double precision   objadd, sinf, obj,
     &                   a(n), bl(n+1), bu(n+1),
     &                   pi(1), rc(n+1), z(nwcore)
      character*8        names(5)

      double precision   zero,             one
      parameter         (zero   = 0.0d+0,  one    = 1.0d+0)

      iprint = 0   ! The MINOS PRINT   file.
      isumm  = 0   ! The MINOS SUMMARY file.
      ispecs = 0   ! The MINOS SPECS   file.
      
      call mistart( iprint, isumm, ispecs )  ! Initialize MINOS and open
      call miopti( 'Workspace (user) ', lenz, 0, 0, inform )
      call miopti( 'LOG FREQUENCY ', 0, 0, 0, inform )
      call miopti( 'PRINT LEVEL ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FILE ', 0, 0, 0, inform )
      call miopti( 'SUMMARY FREQUENCY ', 0, 0, 0, inform )
*     ------------------------------------------------------------------
*     Now set parameters for moniss
*     ------------------------------------------------------------------
      do i = 1, lenz
         z(i) = zsmall(i)
      end do
      m = 1
      nb = n+1
      ne = n
      nname = 1
      nncon = 0
      nnobj = n
      nnjac = 0
      iobj = 0
      objadd = zero
      do i = 1, ne
         ha(i) = 1
         ka(i) = i
         a(i) = one
      end do
      ka(n+1) = ne+1
      do i = 1, n
         bl(i) = zero
         bu(i) = one
      end do
      bl(nb) = -z(lenz)
      bu(nb) = -z(lenz)
      do i = 1, nb
         hs(i) = 0
      end do
      pi(1) = zero

      call minoss( 'Cold', m, n, nb, ne, nname,
     $             nncon, nnobj, nnjac,
     $             iobj, objadd, names,
     $             a, ha, ka, bl, bu, name1, name2,
     $             hs, xn, pi, rc, 
     $             inform, mincor, ns, ninf, sinf, obj,
     $             z, nwcore )

      zsmall(lenz) = z(lenz)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine funobj( mode, n, x, f, g, nstate, nprob, z, nwcore )

      implicit           double precision (a-h,o-z)
      integer            mode, n, nstate, nprob, nwcore
      double precision   x(n), f, g(n), z(nwcore)

      integer            i, j, ii
      double precision   Q(n,n), cvec(n), grad, ddot  
      double precision   zero,          one,          two
      parameter         (zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0)

      mode=mode
      nstate=nstate
      nprob=nprob

      ii = 0
      do i = 1, n
         do j = 1, n
            ii = ii+1
            Q(j,i) = z(ii)
         end do
      end do
      do i = 1, n
         ii = ii+1
         cvec(i) = z(ii)
      end do
      f = zero
      do i = 1, n
         grad = zero
         do j = 1, n
            grad = grad + Q(i,j)*x(j)
         end do
         g(i) = grad
      end do
      f = ddot(n,x,1,g,1)/two + ddot(n,cvec,1,x,1)
      call daxpy(n,one,cvec,1,g,1)
      z(ii+1) = f
      return
      end

******************************************************************************
*     inform  says what happened; see Chapter 6.3 of the User's Guide.       
*             A summary of possible values follows:                          
*                                                                            
*             inform   Meaning                                               
*
*                0     Optimal solution found.
*                1     The problem is infeasible.
*                2     The problem is unbounded (or badly scaled).
*                3     Too many iterations.
*                4     Apparent stall.  The solution has not changed
*                      for a large number of iterations (e.g. 1000).
*                5     The Superbasics limit is too small.
*                6     Subroutine funobj or funcon requested termination
*                      by returning mode < 0.
*                7     Subroutine funobj seems to be giving incorrect
*                      gradients.
*                8     Subroutine funcon seems to be giving incorrect
*                      gradients.
*                9     The current point cannot be improved.
*               10     Numerical error in trying to satisfy the linear
*                      constraints (or the linearized nonlinear
*                      constraints).  The basis is very ill-conditioned.
*               11     Cannot find a superbasic to replace a basic
*                      variable.
*               12     Basis factorization requested twice in a row.
*                      Should probably be treated as inform = 9.
*               13     Near-optimal solution found.
*                      Should probably be treated as inform = 9.
*
*               20     Not enough storage for the basis factorization.
*               21     Error in basis package.
*               22     The basis is singular after several attempts to
*                      factorize it (and add slacks where necessary).
*
*               30     An OLD BASIS file had dimensions that did not
*                      match the current problem.
*               32     System error.  Wrong number of basic variables.
*
*               40     Fatal errors in the MPS file.
*               41     Not enough storage to read the MPS file.
*               42     Not enough storage to solve the problem.
*
******************************************************************************
