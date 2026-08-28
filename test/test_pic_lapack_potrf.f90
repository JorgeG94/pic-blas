module test_pic_lapack_interfaces_potrf
   !! Tests for pic_potrf.
   !!
   !! The reason this wrapper exists is the `info` path, not the happy one: a
   !! Cholesky is chosen over an eigendecomposition for speed and has to be
   !! able to say when it was the wrong choice. So the tests pin what `info`
   !! reports on an indefinite matrix and on a positive semi-definite one, as
   !! well as that the factor multiplies back to the original.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_lapack_interfaces, only: pic_potrf
   use pic_blas_interfaces, only: pic_gemm
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_lapack_potrf_tests

   real(dp), parameter :: tol = 1.0e-12_dp

contains

   subroutine collect_pic_lapack_potrf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("upper_factor_reproduces_matrix", test_upper), &
                  new_unittest("lower_factor_reproduces_matrix", test_lower), &
                  new_unittest("untouched_triangle_is_left_alone", test_other_triangle), &
                  new_unittest("indefinite_reports_info", test_indefinite), &
                  new_unittest("semidefinite_reports_info", test_semidefinite), &
                  new_unittest("potrf_single_precision", test_single) &
                  ]
   end subroutine collect_pic_lapack_potrf_tests

   !  A symmetric positive definite 3x3.
   subroutine make_spd(A)
      real(dp), intent(out) :: A(3, 3)
      A(1, :) = [4.0_dp, 2.0_dp, -2.0_dp]
      A(2, :) = [2.0_dp, 10.0_dp, 2.0_dp]
      A(3, :) = [-2.0_dp, 2.0_dp, 6.0_dp]
   end subroutine make_spd

   subroutine test_upper(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), u(3, 3), back(3, 3), original(3, 3)
      integer(default_int) :: info, i, j

      call make_spd(A)
      original = A
      call pic_potrf(A, uplo="U", info=info)
      call check(error, info == 0)
      if (allocated(error)) return

      ! A = U**T U, with U the upper triangle of what came back.
      u = 0.0_dp
      do j = 1, 3
         do i = 1, j
            u(i, j) = A(i, j)
         end do
      end do
      call pic_gemm(u, u, back, transa="T")
      call check(error, maxval(abs(back - original)) < tol)
   end subroutine test_upper

   subroutine test_lower(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), l(3, 3), back(3, 3), original(3, 3)
      integer(default_int) :: info, i, j

      call make_spd(A)
      original = A
      call pic_potrf(A, uplo="L", info=info)
      call check(error, info == 0)
      if (allocated(error)) return

      l = 0.0_dp
      do j = 1, 3
         do i = j, 3
            l(i, j) = A(i, j)
         end do
      end do
      call pic_gemm(l, l, back, transb="T")
      call check(error, maxval(abs(back - original)) < tol)
   end subroutine test_lower

   subroutine test_other_triangle(error)
      !! LAPACK neither reads nor writes the triangle it was not asked for.
      !! Worth pinning because a caller that forgets is handed the *original*
      !! matrix's entries there, not zeros, and they look like factor entries.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), original(3, 3)
      integer(default_int) :: info, i, j

      call make_spd(A)
      original = A
      call pic_potrf(A, uplo="U", info=info)
      call check(error, info == 0)
      if (allocated(error)) return

      do j = 1, 3
         do i = j + 1, 3
            call check(error, abs(A(i, j) - original(i, j)) < tol)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_other_triangle

   subroutine test_indefinite(error)
      !! Symmetric but with a negative eigenvalue: info names the order of the
      !! leading minor that failed.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2)
      integer(default_int) :: info

      A(1, :) = [1.0_dp, 2.0_dp]
      A(2, :) = [2.0_dp, 1.0_dp]
      call pic_potrf(A, uplo="U", info=info)
      call check(error, info > 0)
   end subroutine test_indefinite

   subroutine test_semidefinite(error)
      !! Exactly singular, which is the case a fitting metric degrades towards.
      !! It must be reported rather than producing a factor full of infinities.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2)
      integer(default_int) :: info

      A(1, :) = [1.0_dp, 1.0_dp]
      A(2, :) = [1.0_dp, 1.0_dp]
      call pic_potrf(A, uplo="U", info=info)
      call check(error, info > 0)
   end subroutine test_semidefinite

   subroutine test_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: A(2, 2)
      integer(default_int) :: info

      A(1, :) = [4.0_sp, 2.0_sp]
      A(2, :) = [2.0_sp, 5.0_sp]
      call pic_potrf(A, uplo="U", info=info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, abs(A(1, 1) - 2.0_sp) < 1.0e-5_sp)
   end subroutine test_single

end module test_pic_lapack_interfaces_potrf
