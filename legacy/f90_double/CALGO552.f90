!--**--CH751--552--C:ID--1:8:1999
module l1_calgo552
   !! L1 solution of (constrained) linear systems: ACM TOMS algorithm 552 (Barrodale & Roberts).
   use, intrinsic :: iso_c_binding, only: c_int
   use l1_precision, only: wp, dp
   implicit none (type, external)
   private

   public :: CL1

contains

      ! The legacy 18-argument public interface is kept unchanged for now, and the routine is still
      ! one monolithic simplex loop; splitting it into procedures is planned (structure layer).
      ! allow(too-many-arguments, too-complex)
      subroutine CL1(K,L,M,N,KLMD,KLM2D,NKLMD,N2D,Q,KODE,TOLER,iter,X, &
                     RES,ERROR,CU,IU,S) bind(C, name="cl1_")
! THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX
! METHOD OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION
! TO A K BY N SYSTEM OF LINEAR EQUATIONS
!             AX=B
! SUBJECT TO L LINEAR EQUALITY CONSTRAINTS
!             CX=D
! AND M LINEAR INEQUALITY CONSTRAINTS
!             EX.LE.F.
! DESCRIPTION OF PARAMETERS
! K      NUMBER OF ROWS OF THE MATRIX A (K.GE.1).
! L      NUMBER OF ROWS OF THE MATRIX C (L.GE.0).
! M      NUMBER OF ROWS OF THE MATRIX E (M.GE.0).
! N      NUMBER OF COLUMNS OF THE MATRICES A,C,E (N.GE.1).
! KLMD   SET TO AT LEAST K+L+M FOR ADJUSTABLE DIMENSIONS.
! KLM2D  SET TO AT LEAST K+L+M+2 FOR ADJUSTABLE DIMENSIONS.
! NKLMD  SET TO AT LEAST N+K+L+M FOR ADJUSTABLE DIMENSIONS.
! N2D    SET TO AT LEAST N+2 FOR ADJUSTABLE DIMENSIONS
! Q      TWO DIMENSIONAL REAL ARRAY WITH KLM2D ROWS AND
!        AT LEAST N2D COLUMNS.
!        ON ENTRY THE MATRICES A,C AND E, AND THE VECTORS
!        B,D AND F MUST BE STORED IN THE FIRST K+L+M ROWS
!        AND N+1 COLUMNS OF Q AS FOLLOWS
!             A B
!         Q = C D
!             E F
!        THESE VALUES ARE DESTROYED BY THE SUBROUTINE.
! KODE   A CODE USED ON ENTRY TO, AND EXIT
!        FROM, THE SUBROUTINE.
!        ON ENTRY, THIS SHOULD NORMALLY BE SET TO 0.
!        HOWEVER, IF CERTAIN NONNEGATIVITY CONSTRAINTS
!        ARE TO BE INCLUDED IMPLICITLY, RATHER THAN
!        EXPLICITLY IN THE CONSTRAINTS EX.LE.F, THEN KODE
!        SHOULD BE SET TO 1, AND THE NONNEGATIVITY
!        CONSTRAINTS INCLUDED IN THE ARRAYS X AND
!        RES (SEE BELOW).
!        ON EXIT, KODE HAS ONE OF THE
!        FOLLOWING VALUES
!             0- OPTIMAL SOLUTION FOUND,
!             1- NO FEASIBLE SOLUTION TO THE
!                CONSTRAINTS,
!             2- CALCULATIONS TERMINATED
!                PREMATURELY DUE TO ROUNDING ERRORS,
!             3- MAXIMUM NUMBER OF ITERATIONS REACHED.
! TOLER  A SMALL POSITIVE TOLERANCE. EMPIRICAL
!        EVIDENCE SUGGESTS TOLER = 10**(-D*2/3),
!        WHERE D REPRESENTS THE NUMBER OF DECIMAL
!        DIGITS OF ACCURACY AVAILABLE. ESSENTIALLY,
!        THE SUBROUTINE CANNOT DISTINGUISH BETWEEN ZERO
!        AND ANY QUANTITY WHOSE MAGNITUDE DOES NOT EXCEED
!        TOLER. IN PARTICULAR, IT WILL NOT PIVOT ON ANY
!        NUMBER WHOSE MAGNITUDE DOES NOT EXCEED TOLER.
! ITER   ON ENTRY ITER MUST CONTAIN AN UPPER BOUND ON
!        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
!        A SUGGESTED VALUE IS 10*(K+L+M). ON EXIT ITER
!        GIVES THE NUMBER OF SIMPLEX ITERATIONS.
! X      ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST N2D.
!        ON EXIT THIS ARRAY CONTAINS A
!        SOLUTION TO THE L1 PROBLEM. IF KODE=1
!        ON ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE
!        SIMPLE NONNEGATIVITY CONSTRAINTS ON THE
!        VARIABLES. THE VALUES -1, 0, OR 1
!        FOR X(J) INDICATE THAT THE J-TH VARIABLE
!        IS RESTRICTED TO BE .LE.0, UNRESTRICTED,
!        OR .GE.0 RESPECTIVELY.
! RES    ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST KLMD.
!        ON EXIT THIS CONTAINS THE RESIDUALS B-AX
!        IN THE FIRST K COMPONENTS, D-CX IN THE
!        NEXT L COMPONENTS (THESE WILL BE =0),AND
!        F-EX IN THE NEXT M COMPONENTS. IF KODE=1 ON
!        ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE SIMPLE
!        NONNEGATIVITY CONSTRAINTS ON THE RESIDUALS
!        B-AX. THE VALUES -1, 0, OR 1 FOR RES(I)
!        INDICATE THAT THE I-TH RESIDUAL (1.LE.I.LE.K) IS
!        RESTRICTED TO BE .LE.0, UNRESTRICTED, OR .GE.0
!        RESPECTIVELY.
! ERROR  ON EXIT, THIS GIVES THE MINIMUM SUM OF
!        ABSOLUTE VALUES OF THE RESIDUALS.
! CU     A TWO DIMENSIONAL REAL ARRAY WITH TWO ROWS AND
!        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.
! IU     A TWO DIMENSIONAL INTEGER ARRAY WITH TWO ROWS AND
!        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.
! S      INTEGER ARRAY OF SIZE AT LEAST KLMD, USED FOR
!        WORKSPACE.
! IF YOUR FORTRAN COMPILER PERMITS A SINGLE COLUMN OF A TWO
! DIMENSIONAL ARRAY TO BE PASSED TO A ONE DIMENSIONAL ARRAY
! THROUGH A SUBROUTINE CALL, CONSIDERABLE SAVINGS IN
! EXECUTION TIME MAY BE ACHIEVED THROUGH THE USE OF THE
! FOLLOWING SUBROUTINE, WHICH OPERATES ON COLUMN VECTORS.
!     SUBROUTINE COL(V1, V2, XMLT, NOTROW, K)
! THIS SUBROUTINE ADDS TO THE VECTOR V1 A MULTIPLE OF THE
! VECTOR V2 (ELEMENTS 1 THROUGH K EXCLUDING NOTROW).
!     DIMENSION V1(K), V2(K)
!     KEND = NOTROW - 1
!     KSTART = NOTROW + 1
!     IF (KEND .LT. 1) GO TO 20
!     DO 10 I=1,KEND
!        V1(I) = V1(I) + XMLT*V2(I)
!  10 CONTINUE
!     IF(KSTART .GT. K) GO TO 40
!  20 DO 30 I=KSTART,K
!       V1(I) = V1(I) + XMLT*V2(I)
!  30 CONTINUE
!  40 RETURN
!     END
! SEE COMMENTS FOLLOWING STATEMENT LABELLED 440 FOR
! INSTRUCTIONS ON THE IMPLEMENTATION OF THIS MODIFICATION.
!
! INITIALIZATION.
!
!     .. Scalar Arguments ..
      integer(c_int), intent(in) :: K,L,M,N,KLMD,KLM2D,NKLMD,N2D
      integer(c_int), intent(inout) :: KODE
      integer(c_int), intent(inout) :: iter
      real(wp), intent(in) :: TOLER
      real(wp), intent(out) :: ERROR
!     ..
!     .. Array Arguments ..
      real(wp), intent(inout) :: Q(KLM2D,N2D)
      real(wp), intent(inout) :: X(N2D)
      real(wp), intent(inout) :: RES(KLMD)
      real(wp), intent(inout) :: CU(2,NKLMD)
      integer(c_int), intent(inout) :: IU(2,NKLMD)
      integer(c_int), intent(inout) :: S(KLMD)
!     ..
!     .. Local Scalars ..
      real(dp) :: SUM
      real(wp) :: CUV,PIVOT,SN,TPIVOT,XMAX,XMIN,tmp1,ZU,ZV
      integer :: I,iq,II,IIMN,IINEG,IN,IOUT,IPHASE,J,JMN,JPN,JS,KFORCE,KK, &
              KLM,KLM1,KLM2,max_iterter,N1,N2,NK,NK1,NKL,NKL1,NKLM
      logical :: to_phase2,at_optimum,pivot_found
!     ..
      max_iterter = iter
      N1 = N + 1
      N2 = N + 2
      NK = N + K
      NK1 = NK + 1
      NKL = NK + L
      NKL1 = NKL + 1
      KLM = K + L + M
      KLM1 = KLM + 1
      KLM2 = KLM + 2
      NKLM = N + KLM
      KFORCE = 1
      iter = 0
      JS = 1
      iq = 0
! SET UP LABELS IN Q.
      do J = 1,N
          Q(KLM2,J) = real(J, wp)
      end do
      label_rows: do I = 1,KLM
          Q(I,N2) = real(N + I, wp)
          if (Q(I,N1) >= 0.0_wp) cycle label_rows
          do J = 1,N2
              Q(I,J) = -Q(I,J)
          end do
      end do label_rows
! SET UP PHASE 1 COSTS.
      IPHASE = 2
      do J = 1,NKLM
          CU(1,J) = 0.0_wp
          CU(2,J) = 0.0_wp
          IU(1,J) = 0
          IU(2,J) = 0
      end do
      if (L /= 0) then
          do J = NK1,NKL
              CU(1,J) = 1.0_wp
              CU(2,J) = 1.0_wp
              IU(1,J) = 1
              IU(2,J) = 1
          end do
          IPHASE = 1
      end if
      if (M /= 0) then
          do J = NKL1,NKLM
              CU(2,J) = 1.0_wp
              IU(2,J) = 1
              JMN = J - N
              if (Q(JMN,N2) < 0.0_wp) IPHASE = 1
          end do
      end if
      if (KODE /= 0) then
! NONNEGATIVITY RESTRICTIONS ON X (SIGN IN X(J)) AND ON THE RESIDUALS (SIGN IN RES(J)).
          do J = 1,N
              if (X(J) < 0.0_wp) then
                  CU(1,J) = 1.0_wp
                  IU(1,J) = 1
              else if (X(J) /= 0.0_wp) then
                  CU(2,J) = 1.0_wp
                  IU(2,J) = 1
              end if
          end do
          do J = 1,K
              JPN = J + N
              if (RES(J) < 0.0_wp) then
                  CU(1,JPN) = 1.0_wp
                  IU(1,JPN) = 1
                  if (Q(J,N2) > 0.0_wp) IPHASE = 1
              else if (RES(J) /= 0.0_wp) then
                  CU(2,JPN) = 1.0_wp
                  IU(2,JPN) = 1
                  if (Q(J,N2) < 0.0_wp) IPHASE = 1
              end if
          end do
      end if
!
! MAIN LOOP.  EACH PASS (RE)STARTS WITH THE PHASE 2 COSTS (ONLY WHEN TO_PHASE2 IS SET) AND
! THE MARGINAL COSTS, THEN PIVOTS UNTIL THE CURRENT PHASE IS OPTIMAL.
!
      to_phase2 = (IPHASE == 2)
      simplex: do
          if (to_phase2) then
! SET UP PHASE 2 COSTS.
              IPHASE = 2
              do J = 1,NKLM
                  CU(1,J) = 0.0_wp
                  CU(2,J) = 0.0_wp
              end do
              do J = N1,NK
                  CU(1,J) = 1.0_wp
                  CU(2,J) = 1.0_wp
              end do
              phase2_costs: do I = 1,KLM
                  II = int(Q(I,N2))
                  if (II > 0) then
                      if (IU(1,II) == 0) cycle phase2_costs
                      CU(1,II) = 0.0_wp
                  else
                      II = -II
                      if (IU(2,II) == 0) cycle phase2_costs
                      CU(2,II) = 0.0_wp
                  end if
                  iq = iq + 1
                  do J = 1,N2
                      tmp1 = Q(iq,J)
                      Q(iq,J) = Q(I,J)
                      Q(I,J) = tmp1
                  end do
              end do phase2_costs
          end if
! COMPUTE THE MARGINAL COSTS.
          do J = JS,N1
              SUM = 0.0_dp
              do I = 1,KLM
                  II = int(Q(I,N2))
                  if (II < 0) then
                      IINEG = -II
                      tmp1 = CU(2,IINEG)
                  else
                      tmp1 = CU(1,II)
                  end if
                  SUM = SUM + real(Q(I,J), dp)*real(tmp1, dp)
              end do
              Q(KLM1,J) = real(SUM, wp)
          end do
          do J = JS,N
              II = int(Q(KLM2,J))
              if (II < 0) then
                  IINEG = -II
                  tmp1 = CU(2,IINEG)
              else
                  tmp1 = CU(1,II)
              end if
              Q(KLM1,J) = Q(KLM1,J) - tmp1
          end do

          iterate: do
! DETERMINE THE VECTOR TO ENTER THE BASIS.
              XMAX = 0.0_wp
              at_optimum = (JS > N)
              if (.not. at_optimum) then
                  entering_candidates: do J = JS,N
                      ZU = Q(KLM1,J)
                      II = int(Q(KLM2,J))
                      if (II > 0) then
                          ZV = -ZU - CU(1,II) - CU(2,II)
                      else
                          II = -II
                          ZV = ZU
                          ZU = -ZU - CU(1,II) - CU(2,II)
                      end if
                      if (KFORCE == 1 .and. II > N) cycle entering_candidates
                      candidate_u: block
                          if (IU(1,II) == 1) exit candidate_u
                          if (ZU <= XMAX) exit candidate_u
                          XMAX = ZU
                          IN = J
                      end block candidate_u
                      if (IU(2,II) == 1) cycle entering_candidates
                      if (ZV <= XMAX) cycle entering_candidates
                      XMAX = ZV
                      IN = J
                  end do entering_candidates
                  at_optimum = (XMAX <= TOLER)
              end if
!
! TEST FOR OPTIMALITY.
!
              if (at_optimum) then
                  if (KFORCE == 0) then
                      if (IPHASE == 1) then
                          if (Q(KLM1,N1) <= TOLER) then
                              to_phase2 = .true.
                              cycle simplex
                          end if
                          KODE = 1
                      else
                          KODE = 0
                      end if
                      exit simplex
                  end if
                  if (IPHASE == 1 .and. Q(KLM1,N1) <= TOLER) then
                      to_phase2 = .true.
                      cycle simplex
                  end if
                  KFORCE = 0
                  cycle iterate
              end if
              if (Q(KLM1,IN) /= XMAX) then
                  do I = 1,KLM2
                      Q(I,IN) = -Q(I,IN)
                  end do
                  Q(KLM1,IN) = XMAX
              end if
!
! DETERMINE THE VECTOR TO LEAVE THE BASIS.
!
              pivot_found = .false.
              artificial_pivot: block
                  if (IPHASE == 1 .or. iq == 0) exit artificial_pivot
                  XMAX = 0.0_wp
                  find_artificial_pivot: do I = 1,iq
                      tmp1 = abs(Q(I,IN))
                      if (tmp1 <= XMAX) cycle find_artificial_pivot
                      XMAX = tmp1
                      IOUT = I
                  end do find_artificial_pivot
                  if (XMAX <= TOLER) exit artificial_pivot
                  do J = 1,N2
                      tmp1 = Q(iq,J)
                      Q(iq,J) = Q(IOUT,J)
                      Q(IOUT,J) = tmp1
                  end do
                  IOUT = iq
                  iq = iq - 1
                  PIVOT = Q(IOUT,IN)
                  pivot_found = .true.
              end block artificial_pivot

              if (.not. pivot_found) then
                  KK = 0
                  ratio_candidates: do I = 1,KLM
                      tmp1 = Q(I,IN)
                      if (tmp1 <= TOLER) cycle ratio_candidates
                      KK = KK + 1
                      RES(KK) = Q(I,N1)/tmp1
                      S(KK) = I
                  end do ratio_candidates
! RATIO TEST: PICK THE SMALLEST RATIO, OR BYPASS THE VERTEX AND PICK AGAIN.
                  bypass: do
                      if (KK <= 0) then
                          KODE = 2
                          exit simplex
                      end if
                      XMIN = RES(1)
                      IOUT = S(1)
                      J = 1
                      if (KK /= 1) then
                          smallest_ratio: do I = 2,KK
                              if (RES(I) >= XMIN) cycle smallest_ratio
                              J = I
                              XMIN = RES(I)
                              IOUT = S(I)
                          end do smallest_ratio
                          RES(J) = RES(KK)
                          S(J) = S(KK)
                      end if
                      KK = KK - 1
                      PIVOT = Q(IOUT,IN)
                      II = int(Q(IOUT,N2))
                      if (IPHASE /= 1) then
                          if (II < 0) then
                              IINEG = -II
                              if (IU(1,IINEG) == 1) exit bypass
                          else
                              if (IU(2,II) == 1) exit bypass
                          end if
                      end if
                      II = abs(II)
                      CUV = CU(1,II) + CU(2,II)
                      if (Q(KLM1,IN)-PIVOT*CUV <= TOLER) exit bypass
! BYPASS INTERMEDIATE VERTICES.
                      do J = JS,N1
                          tmp1 = Q(IOUT,J)
                          Q(KLM1,J) = Q(KLM1,J) - tmp1*CUV
                          Q(IOUT,J) = -tmp1
                      end do
                      Q(IOUT,N2) = -Q(IOUT,N2)
                  end do bypass
              end if
!
! GAUSS-JORDAN ELIMINATION.
!
              if (iter >= max_iterter) then
                  KODE = 3
                  exit simplex
              end if
              iter = iter + 1
              do J = JS,N1
                  if (J /= IN) Q(IOUT,J) = Q(IOUT,J)/PIVOT
              end do
! IF PERMITTED, USE SUBROUTINE COL OF THE DESCRIPTION
! SECTION AND REPLACE THE FOLLOWING SEVEN STATEMENTS DOWN
! TO AND INCLUDING STATEMENT NUMBER 460 BY..
!     DO 460 J=JS,N1
!        IF(J .EQ. IN) GO TO 460
!        Z = -Q(IOUT,J)
!        CALL COL(Q(1,J), Q(1,IN), Z, IOUT, KLM1)
! 460 CONTINUE
              eliminate_columns: do J = JS,N1
                  if (J == IN) cycle eliminate_columns
                  tmp1 = -Q(IOUT,J)
                  do I = 1,KLM1
                      if (I /= IOUT) Q(I,J) = Q(I,J) + tmp1*Q(I,IN)
                  end do
              end do eliminate_columns
              TPIVOT = -PIVOT
              do I = 1,KLM1
                  if (I /= IOUT) Q(I,IN) = Q(I,IN)/TPIVOT
              end do
              Q(IOUT,IN) = 1.0_wp/PIVOT
              tmp1 = Q(IOUT,N2)
              Q(IOUT,N2) = Q(KLM2,IN)
              Q(KLM2,IN) = tmp1
              II = int(abs(tmp1))
              if (IU(1,II) == 0 .or. IU(2,II) == 0) cycle iterate
              do I = 1,KLM2
                  tmp1 = Q(I,IN)
                  Q(I,IN) = Q(I,JS)
                  Q(I,JS) = tmp1
              end do
              JS = JS + 1
          end do iterate
      end do simplex
!
! PREPARE OUTPUT.
!
      SUM = 0.0_dp
      do J = 1,N
          X(J) = 0.0_wp
      end do
      do I = 1,KLM
          RES(I) = 0.0_wp
      end do
      do I = 1,KLM
          II = int(Q(I,N2))
          if (II > 0) then
              SN = 1.0_wp
          else
              II = -II
              SN = -1.0_wp
          end if
          if (II <= N) then
              X(II) = SN*Q(I,N1)
          else
              IIMN = II - N
              RES(IIMN) = SN*Q(I,N1)
              if (II >= N1 .and. II <= NK) SUM = SUM + real(Q(I,N1), dp)
          end if
      end do
      ERROR = real(SUM, wp)

      end subroutine CL1

end module l1_calgo552
