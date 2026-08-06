module test_pic_blas_interfaces_level1_x
   !! Tests for the explicit-length Level-1 entries.
   !!
   !! Two things they add over the assumed-shape forms. They cost nothing --
   !! the assumed-shape ones build an array descriptor per call, which is
   !! about 4 ns and therefore 1.26x on a DOT of length 10 -- and they take
   !! strides, so a row of a matrix goes in directly instead of through a
   !! section that would be packed. The stride tests below are the ones that
   !! matter; a wrapper that ignored incx would pass every contiguous test.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_blas_interfaces, only: pic_dot_x, pic_axpy_x, pic_copy_x, pic_scal_x, &
                                  pic_ddot_x, pic_daxpy_x, pic_dcopy_x, pic_dscal_x
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_level1_x_tests

   real(dp), parameter :: tol = 1.0e-12_dp
   integer(default_int), parameter :: one = 1_default_int

contains

   subroutine collect_pic_level1_x_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dot_x_contiguous", test_dot), &
                  new_unittest("dot_x_strided_row", test_dot_stride), &
                  new_unittest("axpy_x_contiguous", test_axpy), &
                  new_unittest("axpy_x_strided_row", test_axpy_stride), &
                  new_unittest("copy_x_contiguous_and_strided", test_copy), &
                  new_unittest("scal_x_contiguous_and_strided", test_scal), &
                  new_unittest("zero_length_is_a_no_op", test_zero), &
                  new_unittest("level1_x_single_precision", test_single) &
                  ]

   end subroutine collect_pic_level1_x_tests

   subroutine test_dot(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: x(4), y(4)
      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      y = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
      call check(error, abs(pic_dot_x(4_default_int, x, one, y, one) - 300.0_dp) < tol)
   end subroutine test_dot

   !  A row of a column-major matrix has stride equal to the leading
   !  dimension. This is the case the assumed-shape form cannot express
   !  without a packed temporary.
   subroutine test_dot_stride(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 5, n = 3
      real(dp) :: A(ld, n), v(n), expect
      integer(default_int) :: i, j
      do j = 1, n
         do i = 1, ld
            A(i, j) = real(i*10 + j, dp)
         end do
         v(j) = real(j, dp)
      end do
      !  dot of row 2 of A with v
      expect = 0.0_dp
      do j = 1, n
         expect = expect + A(2, j)*v(j)
      end do
      call check(error, abs(pic_ddot_x(n, A(2, 1), ld, v, one) - expect) < tol)
   end subroutine test_dot_stride

   subroutine test_axpy(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: x(3), y(3)
      x = [1.0_dp, 2.0_dp, 3.0_dp]
      y = [10.0_dp, 20.0_dp, 30.0_dp]
      call pic_axpy_x(3_default_int, 2.0_dp, x, one, y, one)
      call check(error, maxval(abs(y - [12.0_dp, 24.0_dp, 36.0_dp])) < tol)
   end subroutine test_axpy

   subroutine test_axpy_stride(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 4, n = 3
      real(dp) :: A(ld, n), v(n)
      integer(default_int) :: i, j
      do j = 1, n
         do i = 1, ld
            A(i, j) = real(i, dp)
         end do
         v(j) = real(j*100, dp)
      end do
      !  add v into row 3 of A, leaving every other row alone
      call pic_daxpy_x(n, 1.0_dp, v, one, A(3, 1), ld)
      do j = 1, n
         call check(error, abs(A(3, j) - (3.0_dp + real(j*100, dp))) < tol)
         if (allocated(error)) return
      end do
      do j = 1, n
         do i = 1, ld
            if (i /= 3) then
               call check(error, A(i, j) == real(i, dp))
               if (allocated(error)) return
            end if
         end do
      end do
   end subroutine test_axpy_stride

   subroutine test_copy(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 4, n = 3
      real(dp) :: x(3), y(3), A(ld, n)
      integer(default_int) :: j
      x = [7.0_dp, 8.0_dp, 9.0_dp]
      y = 0.0_dp
      call pic_copy_x(3_default_int, x, one, y, one)
      call check(error, maxval(abs(y - x)) < tol)
      if (allocated(error)) return
      A = -1.0_dp
      call pic_dcopy_x(n, x, one, A(2, 1), ld)
      do j = 1, n
         call check(error, A(2, j) == x(j))
         if (allocated(error)) return
         call check(error, A(1, j) == -1.0_dp)
         if (allocated(error)) return
      end do
   end subroutine test_copy

   subroutine test_scal(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: ld = 4, n = 3
      real(dp) :: x(3), A(ld, n)
      integer(default_int) :: i, j
      x = [1.0_dp, 2.0_dp, 3.0_dp]
      call pic_scal_x(3_default_int, 3.0_dp, x, one)
      call check(error, maxval(abs(x - [3.0_dp, 6.0_dp, 9.0_dp])) < tol)
      if (allocated(error)) return
      do j = 1, n
         do i = 1, ld
            A(i, j) = 2.0_dp
         end do
      end do
      call pic_dscal_x(n, 5.0_dp, A(4, 1), ld)
      do j = 1, n
         call check(error, A(4, j) == 10.0_dp)
         if (allocated(error)) return
         call check(error, A(1, j) == 2.0_dp)
         if (allocated(error)) return
      end do
   end subroutine test_scal

   subroutine test_zero(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: x(2), y(2)
      x = [1.0_dp, 2.0_dp]
      y = [3.0_dp, 4.0_dp]
      call check(error, abs(pic_dot_x(0_default_int, x, one, y, one)) < tol)
      if (allocated(error)) return
      call pic_axpy_x(0_default_int, 5.0_dp, x, one, y, one)
      call check(error, y(1) == 3.0_dp .and. y(2) == 4.0_dp)
   end subroutine test_zero

   subroutine test_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: x(3), y(3)
      x = [1.0_sp, 2.0_sp, 3.0_sp]
      y = [1.0_sp, 1.0_sp, 1.0_sp]
      call check(error, abs(pic_dot_x(3_default_int, x, one, y, one) - 6.0_sp) < 1.0e-5_sp)
      if (allocated(error)) return
      call pic_scal_x(3_default_int, 2.0_sp, x, one)
      call check(error, abs(x(3) - 6.0_sp) < 1.0e-5_sp)
   end subroutine test_single

end module test_pic_blas_interfaces_level1_x
