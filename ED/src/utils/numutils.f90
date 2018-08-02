!==========================================================================================!
!==========================================================================================!
!  Change Log                                                                              !
!  2.0.0                                                                                   !
!                                                                                          !
!------------------------------------------------------------------------------------------!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!

!******************************************************************
subroutine atob(n,a,b)
   implicit none
   integer, intent(in)                :: n
   real   , intent(in) , dimension(n) :: a
   real   , intent(out), dimension(n) :: b
   integer :: i
   do i=1,n
     b(i)=a(i)
   end do
   return
end subroutine atob





!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from smallest to largest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_up(nmax,variable,ranking)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                    :: nmax
   real    , intent(in)  , dimension(nmax) :: variable
   integer , intent(out) , dimension(nmax) :: ranking
   !----- Local variables. ----------------------------------------------------------------!
   logical ,               dimension(nmax) :: unlocked
   integer                                 :: n
   integer                                 :: locmin
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.
   ranking (:) = 0
   do n=1,nmax
     locmin           = minloc(variable,1,unlocked)
     unlocked(locmin) = .false.
     ranking (locmin) = n
   end do

   return
end subroutine rank_up
!==========================================================================================!
!==========================================================================================!


!==========================================================================================!
!==========================================================================================!
!     This sub-routine sorts the elements of vector a from smallest to largest.            !
!------------------------------------------------------------------------------------------!
subroutine sort_up(a,n)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                  :: n
   integer , intent(inout), dimension(n) :: a
   !----- Local variables. ----------------------------------------------------------------!
   logical ,                dimension(n) :: unlocked
   integer                               :: atmp
   integer                               :: imin
   integer                               :: k
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.

   do k=1,n
      imin        = minloc(a,1,unlocked)
      atmp        = a(imin)
      a(imin)     = a(k)
      a(k)        = atmp
      unlocked(k) = .false.
   end do
   return
end subroutine sort_up
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This sub-routine ranks the elements of vector a from largest to smallest.            !
!------------------------------------------------------------------------------------------!
subroutine rank_down(nmax,variable,ranking)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer , intent(in)                    :: nmax
   real    , intent(in)  , dimension(nmax) :: variable
   integer , intent(out) , dimension(nmax) :: ranking
   !----- Local variables. ----------------------------------------------------------------!
   logical ,               dimension(nmax) :: unlocked
   integer                                 :: n
   integer                                 :: locmax
   !---------------------------------------------------------------------------------------!

   unlocked(:) = .true.
   ranking (:) = 0
   do n=1,nmax
     locmax           = maxloc(variable,1,unlocked)
     unlocked(locmax) = .false.
     ranking (locmax) = n
   end do

   return
end subroutine rank_down
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function finds the element of the rank array that has a given rank.             !
!------------------------------------------------------------------------------------------!
integer function find_rank(ranking,nmax,rankarray)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)                  :: ranking
   integer, intent(in)                  :: nmax
   integer, intent(in), dimension(nmax) :: rankarray
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: n
   !---------------------------------------------------------------------------------------!
   find_rank=-1
   do n=1,nmax
      if (rankarray(n) == ranking) then
         find_rank=n
         return
      end if
   end do
   if (find_rank < 0) call fatal_error('Index not found','find_rank','numutils.f90')
   return
end function find_rank
!==========================================================================================!
!==========================================================================================!







!==========================================================================================!
!==========================================================================================!
!   This function simply computes the cubic root of all numbers, including the negative    !
! ones.                                                                                    !
!------------------------------------------------------------------------------------------!
real function cbrt(x)
   use consts_coms, only: onethird
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in) :: x
   !---------------------------------------------------------------------------------------!

   if (x > 0.0) then
     cbrt=x**onethird
   else
     cbrt=-((-x)**onethird)
   end if

   return
end function cbrt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function simply computes the cubic root of all numbers, including the negative    !
! ones, for a double precision number.                                                     !
!------------------------------------------------------------------------------------------!
real(kind=8) function cbrt8(x)
   use consts_coms, only: onethird8
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   !---------------------------------------------------------------------------------------!
   if (x > 0.d0) then
     cbrt8 = x**onethird8
   else
     cbrt8 = -((-x)**onethird8)
   end if 

   return
end function cbrt8
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!    This function converts the double precision variable into single, in a way to prevent !
! floating point exception when they are tiny.  In case the number is too small, less than !
! off, then the output value is flushed to 0.                                              !
!------------------------------------------------------------------------------------------!
real function sngloff(x,off)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8), intent(in) :: x
   real(kind=8), intent(in) :: off
   !---------------------------------------------------------------------------------------!
   
   if (abs(x) < off) then
      sngloff = 0.
   else
      sngloff = sngl(x)
   end if
   return
end function sngloff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine returns the accumulated sum of a given vector.                        !
!------------------------------------------------------------------------------------------!
subroutine cumsum(nsiz,vec)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                 , intent(in)    :: nsiz
   real   , dimension(nsiz), intent(inout) :: vec
   !----- Local variables. ----------------------------------------------------------------!
   integer                :: n
   !---------------------------------------------------------------------------------------!
   do n=2,nsiz
      vec(n) = vec(n) + vec(n-1)
   end do

   return
end subroutine cumsum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine is the double precision version of the linear system solver above.   !
! It will solve the linear system AA . X = Y for given AA and Y, using the Gaussian        !
! elimination method with partial pivoting and back-substitution.  This subroutine is      !
! based on:                                                                                !
!                                                                                          !
! Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery: 1992. Numerical recipes !
!    in Fortran 77.  Cambridge University Press.                                           !
!------------------------------------------------------------------------------------------!
subroutine lisys_solver8(nsiz,AA,Y,X,sing)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)  :: nsiz  ! matrix and vector size
   real(kind=8), dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
   real(kind=8), dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
   real(kind=8), dimension(nsiz)     , intent(out) :: X     ! unknown vector
   logical                           , intent(out) :: sing  ! The matrix was singular [T|F]
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
   real(kind=8), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
   real(kind=8), dimension(nsiz)                   :: dumvec ! Dummy vector (row swapping)
   real(kind=8)                                    :: pivot  ! The pivot
   real(kind=8)                                    :: multip ! Multiplier
   integer                                         :: r      ! Row index
   integer                                         :: b      ! Row below index
   integer                                         :: p      ! Pivot index
   real(kind=8)                                    :: dumsca ! Dummy scalar (row swapping)
   !----- Local parameters. ---------------------------------------------------------------!
   real(kind=8)                      , parameter   :: tinyoff=1.d-20
   !---------------------------------------------------------------------------------------!
   
   !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
   EE(:,:) = AA(:,:)
   Z (:)   = Y (:)
   dumvec  = 0.d0
   dumsca  = 0.d0
   !---------------------------------------------------------------------------------------!
   !     We initialise X with a huge, non-sense value, which will become the answer when   !
   ! the matrix is singular.                                                               !
   !---------------------------------------------------------------------------------------!
   X (:)   = -huge(1.d0)
   !----- We first assume that everything will be fine. -----------------------------------!
   sing    = .false.

   !---------------------------------------------------------------------------------------!
   ! 1. Main elimination loop, done row by row.                                            !
   !---------------------------------------------------------------------------------------!
   elimloop: do r = 1, nsiz-1
      !------ 1a. Finding the largest element, which will become our pivot ----------------!
      p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)
      
      pivot = maxval(abs(EE(r:nsiz,r)))
      !------------------------------------------------------------------------------------!
      ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
      !     singular or almost singular, and we cannot solve it, so we switch the flag and !
      !     return.                                                                        !
      !------------------------------------------------------------------------------------!
      if (pivot < tinyoff) then
         sing = .true.
         return
      end if
      
      !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
      if (p /= r) then
         dumvec(r:nsiz) = EE(r,r:nsiz)
         dumsca         = Z(r)
         EE(r,r:nsiz)   = EE(p,r:nsiz)
         Z(r)           = Z(p)
         EE(p,r:nsiz)   = dumvec(r:nsiz)
         Z(p)           = dumsca
      end if

      !------------------------------------------------------------------------------------!
      ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
      !      zero (we won't compute that, but they will be.).                              !
      !------------------------------------------------------------------------------------!
      belowloop: do b=r+1,nsiz
         multip = EE(b,r)/EE(r,r)
         EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
         Z(b)         = Z(b)         - multip * Z(r)
      end do belowloop
   end do elimloop

   !---------------------------------------------------------------------------------------!
   ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
   !    check the last pivot too.                                                          ! 
   !---------------------------------------------------------------------------------------!
   if (abs(EE(nsiz,nsiz)) < tinyoff) then
      sing = .true.
      return
   end if

   !---------------------------------------------------------------------------------------!
   ! 3. We now perform the back-substitution, to find the solution.                        !
   !---------------------------------------------------------------------------------------!
   X(nsiz) = Z(nsiz) / EE(nsiz,nsiz)
   backsubloop: do r=nsiz-1,1,-1
      b    = r+1
      X(r) = (Z(r) - sum(EE(r,b:nsiz)*x(b:nsiz))) / EE(r,r)
   end do backsubloop

   return
end subroutine lisys_solver8
!==========================================================================================!
!==========================================================================================!

!==========================================================================================!
!==========================================================================================!
!   Function to compute the great circle distance between two points: the s suffix denotes !
! source point, and f denotes the destination - "forepoint"). The results are given in     !
! metres. The formula is intended to be accurate for both small and large distances and    !
! uses double precision to avoid ill-conditioned behaviour of sin and cos for numbers      !
! close to the n*pi/2.                                                                     !
!------------------------------------------------------------------------------------------!
real function dist_gc(slons,slonf,slats,slatf)
   use consts_coms, only : erad    & ! intent(in)
                         , pio1808 ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real, intent(in) :: slons
   real, intent(in) :: slonf
   real, intent(in) :: slats
   real, intent(in) :: slatf
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)     :: lons
   real(kind=8)     :: lonf
   real(kind=8)     :: lats
   real(kind=8)     :: latf
   real(kind=8)     :: dlon
   real(kind=8)     :: dlat
   real(kind=8)     :: x
   real(kind=8)     :: y
   !---------------------------------------------------------------------------------------!

   !----- Convert the co-ordinates to double precision and to radians. --------------------!
   lons = dble(slons) * pio1808
   lonf = dble(slonf) * pio1808
   lats = dble(slats) * pio1808
   latf = dble(slatf) * pio1808
   dlon = lonf - lons
   dlat = latf - lats

   !----- Find the arcs. ------------------------------------------------------------------!
   x    = dsin(lats) * dsin(latf) + dcos(lats) * dcos(latf) * dcos(dlon)
   y    = dsqrt( (dcos(latf)*dsin(dlon)) * (dcos(latf)*dsin(dlon))                         &
               + (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon))                  &
               * (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon)) )

   !----- Convert the arcs to actual distance. --------------------------------------------!
   dist_gc = erad*sngl(datan2(y,x))

   return
end function dist_gc
!==========================================================================================!
!==========================================================================================!

