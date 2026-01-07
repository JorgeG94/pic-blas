module test_pic_lapack_interfaces_gesvd
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_lapack_interfaces, only: pic_gesvd
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_lapack_gesvd_tests

contains

   subroutine collect_pic_lapack_gesvd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_sgesvd_basic", test_sgesvd_basic), &
                  new_unittest("test_dgesvd_basic", test_dgesvd_basic), &
                  new_unittest("test_dgesvd_singular_values_only", test_dgesvd_singular_values_only), &
                  new_unittest("test_dgesvd_rectangular", test_dgesvd_rectangular), &
                  new_unittest("test_dgesvd_identity", test_dgesvd_identity), &
                  new_unittest("test_dgesvd_full_decomposition", test_dgesvd_full_decomposition) &
                  ]

   end subroutine collect_pic_lapack_gesvd_tests

   subroutine test_sgesvd_basic(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: A(2, 2), S(2)
      real(sp), parameter :: tol = 1.0e-5_sp
      integer(default_int) :: info

      ! Simple 2x2 matrix
      A(1, 1) = 3.0_sp; A(1, 2) = 0.0_sp
      A(2, 1) = 0.0_sp; A(2, 2) = 4.0_sp

      call pic_gesvd(A, S, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Singular values are sorted descending: 4, 3
      call check(error, abs(S(1) - 4.0_sp) < tol)
      if (allocated(error)) return
      call check(error, abs(S(2) - 3.0_sp) < tol)
      if (allocated(error)) return
   end subroutine test_sgesvd_basic

   subroutine test_dgesvd_basic(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2), S(2)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Simple 2x2 diagonal matrix
      A(1, 1) = 3.0_dp; A(1, 2) = 0.0_dp
      A(2, 1) = 0.0_dp; A(2, 2) = 4.0_dp

      call pic_gesvd(A, S, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Singular values sorted descending: 4, 3
      call check(error, abs(S(1) - 4.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(S(2) - 3.0_dp) < tol)
      if (allocated(error)) return
   end subroutine test_dgesvd_basic

   subroutine test_dgesvd_singular_values_only(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), S(3)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Diagonal matrix
      A = 0.0_dp
      A(1, 1) = 5.0_dp
      A(2, 2) = 2.0_dp
      A(3, 3) = 8.0_dp

      call pic_gesvd(A, S, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Singular values sorted descending: 8, 5, 2
      call check(error, abs(S(1) - 8.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(S(2) - 5.0_dp) < tol)
      if (allocated(error)) return
      call check(error, abs(S(3) - 2.0_dp) < tol)
      if (allocated(error)) return
   end subroutine test_dgesvd_singular_values_only

   subroutine test_dgesvd_rectangular(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 2), S(2)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Rectangular 3x2 matrix (column of 1s and column of 0s)
      A(1, 1) = 1.0_dp; A(1, 2) = 0.0_dp
      A(2, 1) = 1.0_dp; A(2, 2) = 0.0_dp
      A(3, 1) = 1.0_dp; A(3, 2) = 0.0_dp

      call pic_gesvd(A, S, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Singular values: sqrt(3), 0
      call check(error, abs(S(1) - sqrt(3.0_dp)) < tol)
      if (allocated(error)) return
      call check(error, abs(S(2)) < tol)
      if (allocated(error)) return
   end subroutine test_dgesvd_rectangular

   subroutine test_dgesvd_identity(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(3, 3), S(3)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info, i

      ! Identity matrix - all singular values should be 1
      A = 0.0_dp
      do i = 1, 3
         A(i, i) = 1.0_dp
      end do

      call pic_gesvd(A, S, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! All singular values should be 1
      do i = 1, 3
         call check(error, abs(S(i) - 1.0_dp) < tol)
         if (allocated(error)) return
      end do
   end subroutine test_dgesvd_identity

   subroutine test_dgesvd_full_decomposition(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: A(2, 2), A_orig(2, 2), S(2), U(2, 2), VT(2, 2)
      real(dp) :: reconstructed(2, 2)
      real(dp), parameter :: tol = 1.0e-10_dp
      integer(default_int) :: info

      ! Test matrix
      A(1, 1) = 1.0_dp; A(1, 2) = 2.0_dp
      A(2, 1) = 3.0_dp; A(2, 2) = 4.0_dp
      A_orig = A

      call pic_gesvd(A, S, U, VT, info=info)

      call check(error, info == 0)
      if (allocated(error)) return

      ! Reconstruct A from SVD: A = U * S * VT
      ! For 2x2: reconstructed = U(:,1)*S(1)*VT(1,:) + U(:,2)*S(2)*VT(2,:)
      reconstructed(1, 1) = U(1, 1)*S(1)*VT(1, 1) + U(1, 2)*S(2)*VT(2, 1)
      reconstructed(1, 2) = U(1, 1)*S(1)*VT(1, 2) + U(1, 2)*S(2)*VT(2, 2)
      reconstructed(2, 1) = U(2, 1)*S(1)*VT(1, 1) + U(2, 2)*S(2)*VT(2, 1)
      reconstructed(2, 2) = U(2, 1)*S(1)*VT(1, 2) + U(2, 2)*S(2)*VT(2, 2)

      call check(error, all(abs(reconstructed - A_orig) < tol))
      if (allocated(error)) return
   end subroutine test_dgesvd_full_decomposition

end module test_pic_lapack_interfaces_gesvd
