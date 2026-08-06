module test_pic_blas_interfaces_ger
   !! Tests for pic_ger / pic_ger_x / pic_gemv_x.
   !!
   !! GER is the partner GEMV needs: together they turn a Gram-Schmidt written
   !! as a loop of DOTs and AXPYs over one column at a time into two calls over
   !! the whole trailing block. The explicit-dimension forms exist because the
   !! callers that want them hold a block inside a larger array, so every test
   !! here works on a padded array and checks the padding is left alone.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_blas_interfaces, only: pic_ger, pic_ger_x, pic_gemv_x, pic_gemv, &
                                 pic_dgemv_x, pic_dger_x
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_ger_tests

   real(dp), parameter :: tol = 1.0e-12_dp

contains

   subroutine collect_pic_ger_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ger_rank_one_update", test_ger), &
                  new_unittest("ger_alpha", test_ger_alpha), &
                  new_unittest("ger_x_on_a_subblock", test_ger_x), &
                  new_unittest("gemv_x_matches_gemv", test_gemv_x), &
                  new_unittest("gemv_x_transpose", test_gemv_x_trans), &
                  new_unittest("gemv_x_and_ger_x_do_one_mgs_step", test_mgs_step), &
                  new_unittest("ger_single_precision", test_ger_single) &
                  ]

   end subroutine collect_pic_ger_tests

   subroutine test_ger(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 2), x(3), y(2)
      integer(default_int) :: i, j

      A = 1.0_dp
      x = [1.0_dp, 2.0_dp, 3.0_dp]
      y = [10.0_dp, 20.0_dp]
      call pic_ger(A, x, y)
      do j = 1, 2
         do i = 1, 3
            call check(error, abs(A(i, j) - (1.0_dp + x(i)*y(j))) < tol)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_ger

   subroutine test_ger_alpha(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2), x(2), y(2)

      A = 0.0_dp
      x = [1.0_dp, 2.0_dp]
      y = [3.0_dp, 4.0_dp]
      call pic_ger(A, x, y, -2.0_dp)
      call check(error, abs(A(1, 1) - (-6.0_dp)) < tol)
      if (allocated(error)) return
      call check(error, abs(A(2, 2) - (-16.0_dp)) < tol)
   end subroutine test_ger_alpha

   !  the whole point of the _x form: a block inside a padded array
   subroutine test_ger_x(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 7, m = 3, n = 2
      real(dp) :: A(ld, n), x(ld), y(n)
      integer(default_int) :: i, j

      A = -5.0_dp
      do i = 1, ld
         x(i) = real(i, dp)
      end do
      y = [2.0_dp, 3.0_dp]
      call pic_ger_x(m, n, 1.0_dp, x, 1_default_int, y, 1_default_int, A, ld)
      do j = 1, n
         do i = 1, m
            call check(error, abs(A(i, j) - (-5.0_dp + x(i)*y(j))) < tol)
            if (allocated(error)) return
         end do
      end do
      !  rows past m are outside the update and must be untouched
      do j = 1, n
         do i = m + 1, ld
            call check(error, A(i, j) == -5.0_dp)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_ger_x

   subroutine test_gemv_x(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 9, m = 4, n = 3
      real(dp) :: A(ld, n), x(n), y(ld)
      real(dp) :: As(m, n), xs(n), ys(m)
      integer(default_int) :: i, j

      do j = 1, n
         do i = 1, ld
            A(i, j) = real(i*2 + j, dp)
         end do
         x(j) = real(j, dp)
      end do
      y = -7.0_dp
      do j = 1, n
         do i = 1, m
            As(i, j) = A(i, j)
         end do
      end do
      xs = x
      ys = 0.0_dp

      call pic_gemv(As, xs, ys)
      call pic_gemv_x("N", m, n, 1.0_dp, A, ld, x, 1_default_int, 0.0_dp, y, 1_default_int)

      do i = 1, m
         call check(error, abs(y(i) - ys(i)) < tol)
         if (allocated(error)) return
      end do
      do i = m + 1, ld
         call check(error, y(i) == -7.0_dp)
         if (allocated(error)) return
      end do
   end subroutine test_gemv_x

   subroutine test_gemv_x_trans(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 8, m = 5, n = 2
      real(dp) :: A(ld, n), x(ld), y(n), expect(n)
      integer(default_int) :: i, j

      do j = 1, n
         do i = 1, ld
            A(i, j) = real(i + 3*j, dp)
         end do
      end do
      do i = 1, ld
         x(i) = real(i, dp)*0.5_dp
      end do
      !  y := A**T x over the leading m rows only
      do j = 1, n
         expect(j) = 0.0_dp
         do i = 1, m
            expect(j) = expect(j) + A(i, j)*x(i)
         end do
      end do
      y = 0.0_dp
      call pic_gemv_x("T", m, n, 1.0_dp, A, ld, x, 1_default_int, 0.0_dp, y, 1_default_int)
      do j = 1, n
         call check(error, abs(y(j) - expect(j)) < tol)
         if (allocated(error)) return
      end do
   end subroutine test_gemv_x_trans

   !  One step of the blocked modified Gram-Schmidt these two exist for:
   !  project the trailing columns onto the current one and subtract. Compared
   !  against the column-at-a-time form it replaces.
   subroutine test_mgs_step(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 12, n = 6, ncol = 4
      real(dp) :: V(ld, ncol), W(ld, ncol), x(ncol)
      integer(default_int) :: i, j, k
      real(dp) :: d

      do j = 1, ncol
         do i = 1, ld
            V(i, j) = real(mod(i*5 + j*3, 17), dp)
         end do
      end do
      W = V

      !  blocked: one GEMV then one GER over columns 2..ncol
      !  Array elements as the matrix arguments -- the legacy calling style,
      !  and the reason these entries exist. That means the *specific* names:
      !  sequence association lets V(1,2) stand for the submatrix starting at
      !  column 2, but generic resolution matches on rank and will not accept
      !  a scalar actual for an array dummy.
      call pic_dgemv_x("T", n, ncol - 1, 1.0_dp, V(1, 2), ld, V(1, 1), 1_default_int, &
                       0.0_dp, x, 1_default_int)
      call pic_dger_x(n, ncol - 1, -1.0_dp, V(1, 1), 1_default_int, x, 1_default_int, &
                      V(1, 2), ld)

      !  column at a time, which is what it replaces
      do j = 2, ncol
         d = 0.0_dp
         do k = 1, n
            d = d + W(k, 1)*W(k, j)
         end do
         do k = 1, n
            W(k, j) = W(k, j) - d*W(k, 1)
         end do
      end do

      do j = 2, ncol
         do i = 1, n
            call check(error, abs(V(i, j) - W(i, j)) < 1.0e-11_dp)
            if (allocated(error)) return
         end do
      end do
      !  rows past n were not part of the step
      do j = 1, ncol
         do i = n + 1, ld
            call check(error, V(i, j) == W(i, j))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_mgs_step

   subroutine test_ger_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: A(2, 2), x(2), y(2)

      A = 0.0_sp
      x = [1.0_sp, 2.0_sp]
      y = [3.0_sp, 4.0_sp]
      call pic_ger(A, x, y)
      call check(error, abs(A(2, 2) - 8.0_sp) < 1.0e-5_sp)
   end subroutine test_ger_single

end module test_pic_blas_interfaces_ger
