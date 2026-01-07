module test_pic_lapack_interfaces_syevd
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_lapack_interfaces, only: pic_syevd
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_lapack_syevd_tests

contains

   subroutine collect_pic_lapack_syevd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_ssyevd_basic", test_ssyevd_basic), &
                  new_unittest("test_dsyevd_basic", test_dsyevd_basic), &
                  new_unittest("test_dsyevd_eigenvalues_only", test_dsyevd_eigenvalues_only), &
                  new_unittest("test_dsyevd_diagonal", test_dsyevd_diagonal), &
                  new_unittest("test_dsyevd_identity", test_dsyevd_identity) &
                  ]

   end subroutine collect_pic_lapack_syevd_tests

   subroutine test_ssyevd_basic(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: A(3, 3), W(3)
      real(sp), parameter :: tol = 1.0e-5_sp
      integer(default_int) :: info

      ! Symmetric 3x3 matrix (stored in upper triangle)
      A(1, 1) = 2.0_sp; A(1, 2) = 1.0_sp; A(1, 3) = 0.0_sp
      A(2, 1) = 1.0_sp; A(2, 2) = 2.0_sp; A(2, 3) = 1.0_sp
      A(3, 1) = 0.0_sp; A(3, 2) = 1.0_sp; A(3, 3) = 2.0_sp

      call pic_syevd(A, W, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Eigenvalues: 2 - sqrt(2), 2, 2 + sqrt(2)
      call check(error, abs(W(1) - (2.0_sp - sqrt(2.0_sp))) < tol)
      if (allocated(error)) return
      call check(error, abs(W(2) - 2.0_sp) < tol)
      if (allocated(error)) return
      call check(error, abs(W(3) - (2.0_sp + sqrt(2.0_sp))) < tol)
      if (allocated(error)) return
   end subroutine test_ssyevd_basic

   subroutine test_dsyevd_basic(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), W(3)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Symmetric 3x3 matrix
      A(1, 1) = 2.0_dp; A(1, 2) = 1.0_dp; A(1, 3) = 0.0_dp
      A(2, 1) = 1.0_dp; A(2, 2) = 2.0_dp; A(2, 3) = 1.0_dp
      A(3, 1) = 0.0_dp; A(3, 2) = 1.0_dp; A(3, 3) = 2.0_dp

      call pic_syevd(A, W, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Eigenvalues: 2 - sqrt(2), 2, 2 + sqrt(2)
      call check(error, abs(W(1) - (2.0_dp - sqrt(2.0_dp))) < tol)
      if (allocated(error)) return
      call check(error, abs(W(2) - 2.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(W(3) - (2.0_dp + sqrt(2.0_dp))) < tol)
      if (allocated(error)) return
   end subroutine test_dsyevd_basic

   subroutine test_dsyevd_eigenvalues_only(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2), W(2)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Simple 2x2 symmetric matrix
      A(1, 1) = 3.0_dp; A(1, 2) = 1.0_dp
      A(2, 1) = 1.0_dp; A(2, 2) = 3.0_dp

      call pic_syevd(A, W, jobz='N', info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Eigenvalues: 2, 4
      call check(error, abs(W(1) - 2.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(W(2) - 4.0_dp) < tol)
      if (allocated(error)) return
   end subroutine test_dsyevd_eigenvalues_only

   subroutine test_dsyevd_diagonal(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), W(3)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Diagonal matrix - eigenvalues are the diagonal elements
      A = 0.0_dp
      A(1, 1) = 5.0_dp
      A(2, 2) = 2.0_dp
      A(3, 3) = 8.0_dp

      call pic_syevd(A, W, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Eigenvalues sorted ascending: 2, 5, 8
      call check(error, abs(W(1) - 2.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(W(2) - 5.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(W(3) - 8.0_dp) < tol)
      if (allocated(error)) return
   end subroutine test_dsyevd_diagonal

   subroutine test_dsyevd_identity(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(4, 4), W(4)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info, i

      ! Identity matrix - all eigenvalues should be 1
      A = 0.0_dp
      do i = 1, 4
         A(i, i) = 1.0_dp
      end do

      call pic_syevd(A, W, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! All eigenvalues should be 1
      do i = 1, 4
         call check(error, abs(W(i) - 1.0_dp) < tol)
         if (allocated(error)) return
      end do
   end subroutine test_dsyevd_identity

end module test_pic_lapack_interfaces_syevd
