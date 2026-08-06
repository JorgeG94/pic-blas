module test_pic_lapack_interfaces_getrf
   !! Tests for pic_getrf / pic_getrs / pic_getri / pic_gecon.
   !!
   !! These exist to replace LINPACK's DGEFA/DGESL/DGEDI/DGECO, which are
   !! pasted into more than one legacy code. The tests therefore pin the two
   !! things that differ between the two libraries and would otherwise be
   !! found the hard way: the sign convention of the stored L, and the fact
   !! that the determinant has to be recovered by hand because LAPACK has no
   !! routine for it.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs, pic_getri, pic_gecon
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_lapack_getrf_tests

   real(dp), parameter :: tol = 1.0e-11_dp

contains

   subroutine collect_pic_lapack_getrf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("getrf_getrs_solves", test_solve), &
                  new_unittest("getrs_transpose_solves", test_solve_trans), &
                  new_unittest("getrs_multiple_rhs", test_multi_rhs), &
                  new_unittest("getri_inverts", test_inverse), &
                  new_unittest("determinant_from_u_and_pivots", test_determinant), &
                  new_unittest("stored_L_is_unnegated", test_l_sign), &
                  new_unittest("singular_matrix_reports_info", test_singular), &
                  new_unittest("numerically_singular_needs_rcond", test_numerically_singular), &
                  new_unittest("gecon_identity_is_well_conditioned", test_gecon), &
                  new_unittest("getrf_single_precision", test_single) &
                  ]

   end subroutine collect_pic_lapack_getrf_tests

   !  A well conditioned 3x3 with an exactly known inverse and determinant.
   !      A = [ 2 1 1 ; 4 -6 0 ; -2 7 2 ],  det(A) = -16
   subroutine make_a(A)
      real(dp), intent(out) :: A(3, 3)
      A(1, :) = [2.0_dp, 1.0_dp, 1.0_dp]
      A(2, :) = [4.0_dp, -6.0_dp, 0.0_dp]
      A(3, :) = [-2.0_dp, 7.0_dp, 2.0_dp]
   end subroutine make_a

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), B(3, 1)
      integer(default_int) ::  info
      integer(default_int) :: ipiv(3)

      call make_a(A)
      !  A * [1,2,3]' = [7, -8, 18]'
      B(:, 1) = [7.0_dp, -8.0_dp, 18.0_dp]
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return
      call pic_getrs(A, ipiv, B, info=info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, maxval(abs(B(:, 1) - [1.0_dp, 2.0_dp, 3.0_dp])) < tol)
   end subroutine test_solve

   subroutine test_solve_trans(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), Aorig(3, 3), B(3, 1), rhs(3)
      integer(default_int) ::  info, i, j
      integer(default_int) :: ipiv(3)

      call make_a(A)
      Aorig = A
      !  build A' * [1,2,3]' explicitly so the expected answer is independent
      do i = 1, 3
         rhs(i) = 0.0_dp
         do j = 1, 3
            rhs(i) = rhs(i) + Aorig(j, i)*real(j, dp)
         end do
      end do
      B(:, 1) = rhs
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return
      call pic_getrs(A, ipiv, B, "T", info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, maxval(abs(B(:, 1) - [1.0_dp, 2.0_dp, 3.0_dp])) < tol)
   end subroutine test_solve_trans

   subroutine test_multi_rhs(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), Aorig(3, 3), B(3, 2), X(3, 2), r(3)
      integer(default_int) ::  info, i, j, k
      integer(default_int) :: ipiv(3)

      call make_a(A)
      Aorig = A
      X(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
      X(:, 2) = [-1.0_dp, 0.5_dp, 4.0_dp]
      do k = 1, 2
         do i = 1, 3
            r(i) = 0.0_dp
            do j = 1, 3
               r(i) = r(i) + Aorig(i, j)*X(j, k)
            end do
         end do
         B(:, k) = r
      end do
      call pic_getrf(A, ipiv, info)
      call pic_getrs(A, ipiv, B, info=info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, maxval(abs(B - X)) < tol)
   end subroutine test_multi_rhs

   subroutine test_inverse(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), Aorig(3, 3), prod(3, 3)
      integer(default_int) ::  info, i, j, k
      integer(default_int) :: ipiv(3)

      call make_a(A)
      Aorig = A
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return
      call pic_getri(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return

      do i = 1, 3
         do j = 1, 3
            prod(i, j) = 0.0_dp
            do k = 1, 3
               prod(i, j) = prod(i, j) + Aorig(i, k)*A(k, j)
            end do
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            if (i == j) then
               call check(error, abs(prod(i, j) - 1.0_dp) < tol)
            else
               call check(error, abs(prod(i, j)) < tol)
            end if
            if (allocated(error)) return
         end do
      end do
   end subroutine test_inverse

   !  LAPACK has no determinant routine, so callers recover it from the
   !  factors: product of U's diagonal, negated once per pivot that moved a
   !  row. U is the same whichever way L was signed, so this works on either
   !  library's output -- which is what lets a LINPACK DGEDI caller keep its
   !  determinant code unchanged while the factorisation underneath changes.
   subroutine test_determinant(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) ::  det
      real(dp) :: A(3, 3)
      integer(default_int) ::  info, i
      integer(default_int) :: ipiv(3)

      call make_a(A)
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return

      det = 1.0_dp
      do i = 1, 3
         if (ipiv(i) /= i) det = -det
         det = det*A(i, i)
      end do
      call check(error, abs(det - (-16.0_dp)) < 1.0e-10_dp)
   end subroutine test_determinant

   !  The difference that silently breaks a mixed LINPACK/LAPACK code.
   !  For [2 1; 6 5] the pivot is 6, so after the row swap the multiplier is
   !  2/6 = 1/3. LAPACK stores +1/3; LINPACK's DGEFA would store -1/3.
   subroutine test_l_sign(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2)
      integer(default_int) ::  info
      integer(default_int) :: ipiv(2)

      A(1, :) = [2.0_dp, 1.0_dp]
      A(2, :) = [6.0_dp, 5.0_dp]
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, ipiv(1) == 2)
      if (allocated(error)) return
      call check(error, abs(A(2, 1) - (1.0_dp/3.0_dp)) < tol)
      if (allocated(error)) return
      !  and positive, which is the whole point
      call check(error, A(2, 1) > 0.0_dp)
   end subroutine test_l_sign

   !  info is only set for an *exactly* zero pivot -- LAPACK compares against
   !  zero, and so does LINPACK's DGEFA. The multipliers here are halves, so
   !  the elimination cancels exactly and the zero is real rather than 1e-16.
   !  Row 2 is twice row 1.
   subroutine test_singular(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3)
      integer(default_int) ::  info
      integer(default_int) :: ipiv(3)

      A(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]
      A(2, :) = [2.0_dp, 4.0_dp, 6.0_dp]
      A(3, :) = [1.0_dp, 1.0_dp, 1.0_dp]
      call pic_getrf(A, ipiv, info)
      call check(error, info > 0)
   end subroutine test_singular

   !  A matrix that is singular in exact arithmetic but whose elimination
   !  does not cancel exactly reports info == 0 and produces a factorisation
   !  that looks fine. This is not a defect and it is not something a caller
   !  can detect from info -- it is the reason a condition estimate exists,
   !  and the reason LINPACK shipped DGECO alongside DGEFA. rcond must come
   !  back tiny where info came back clean.
   subroutine test_numerically_singular(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) ::  rcond, anorm
      real(dp) :: A(3, 3)
      integer(default_int) ::  info
      integer(default_int) :: ipiv(3)

      !  row 3 = row 1 + row 2, with multipliers that are not representable
      A(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]
      A(2, :) = [4.0_dp, 5.0_dp, 6.0_dp]
      A(3, :) = [5.0_dp, 7.0_dp, 9.0_dp]
      anorm = 18.0_dp   ! 1-norm: largest column sum, |3|+|6|+|9|
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)          ! not flagged, as expected
      if (allocated(error)) return
      call pic_gecon(A, anorm, rcond, "1", info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, rcond < 1.0e-14_dp)
   end subroutine test_numerically_singular

   subroutine test_gecon(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) ::  rcond, anorm
      real(dp) :: A(4, 4)
      integer(default_int) ::  info, i
      integer(default_int) :: ipiv(4)

      A = 0.0_dp
      do i = 1, 4
         A(i, i) = 1.0_dp
      end do
      anorm = 1.0_dp        ! 1-norm of the identity
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return
      call pic_gecon(A, anorm, rcond, "1", info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, abs(rcond - 1.0_dp) < 1.0e-12_dp)
   end subroutine test_gecon

   subroutine test_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: A(3, 3), B(3, 1)
      integer(default_int) ::  info
      integer(default_int) :: ipiv(3)

      A(1, :) = [2.0_sp, 1.0_sp, 1.0_sp]
      A(2, :) = [4.0_sp, -6.0_sp, 0.0_sp]
      A(3, :) = [-2.0_sp, 7.0_sp, 2.0_sp]
      B(:, 1) = [7.0_sp, -8.0_sp, 18.0_sp]
      call pic_getrf(A, ipiv, info)
      call check(error, info == 0)
      if (allocated(error)) return
      call pic_getrs(A, ipiv, B, info=info)
      call check(error, maxval(abs(B(:, 1) - [1.0_sp, 2.0_sp, 3.0_sp])) < 1.0e-4_sp)
   end subroutine test_single

end module test_pic_lapack_interfaces_getrf
