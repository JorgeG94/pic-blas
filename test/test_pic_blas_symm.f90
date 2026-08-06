module test_pic_blas_interfaces_symm
   !! Tests for pic_symm, pic_syrk and pic_syr2k.
   !!
   !! The point of all three is that only one triangle of the symmetric
   !! operand is read. Every test here deliberately fills the *other*
   !! triangle with values that would give a different answer if it were
   !! touched, which is the property that lets a caller skip mirroring.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_blas_interfaces, only: pic_symm, pic_syrk, pic_syr2k, pic_gemm
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_symm_tests

   real(dp), parameter :: tol = 1.0e-12_dp

contains

   subroutine collect_pic_symm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("symm_matches_gemm_on_mirrored", test_symm_vs_gemm), &
                  new_unittest("symm_ignores_unread_triangle", test_symm_ignores), &
                  new_unittest("symm_right_side", test_symm_right), &
                  new_unittest("symm_alpha_beta", test_symm_alpha_beta), &
                  new_unittest("syrk_matches_gemm", test_syrk_vs_gemm), &
                  new_unittest("syrk_transpose", test_syrk_trans), &
                  new_unittest("syr2k_half_alpha_gives_atb", test_syr2k_identity), &
                  new_unittest("symm_single_precision", test_symm_single) &
                  ]

   end subroutine collect_pic_symm_tests

   !  A symmetric matrix and its mirror, plus a version whose lower triangle
   !  is poisoned. Anything reading the lower triangle gets a wrong answer.
   subroutine make_sym(S, poisoned, n)
      integer(default_int), intent(in) :: n
      real(dp), intent(out) :: S(n, n), poisoned(n, n)
      integer(default_int) :: i, j
      do j = 1, n
         do i = 1, j
            S(i, j) = real(i + 2*j, dp)*0.5_dp
            S(j, i) = S(i, j)
         end do
      end do
      poisoned = S
      do j = 1, n
         do i = j + 1, n
            poisoned(i, j) = -1000.0_dp
         end do
      end do
   end subroutine make_sym

   subroutine test_symm_vs_gemm(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 6, k = 4
      real(dp) :: S(n, n), P(n, n), B(n, k), C1(n, k), C2(n, k)
      integer(default_int) :: i, j

      call make_sym(S, P, n)
      do j = 1, k
         do i = 1, n
            B(i, j) = real(i, dp) - 0.25_dp*real(j, dp)
         end do
      end do
      C1 = 0.0_dp
      C2 = 0.0_dp
      call pic_symm(S, B, C1, "L", "U")
      call pic_gemm(S, B, C2)
      call check(error, maxval(abs(C1 - C2)) < tol)
   end subroutine test_symm_vs_gemm

   subroutine test_symm_ignores(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 6, k = 4
      real(dp) :: S(n, n), P(n, n), B(n, k), C1(n, k), C2(n, k)
      integer(default_int) :: i, j

      call make_sym(S, P, n)
      do j = 1, k
         do i = 1, n
            B(i, j) = real(i, dp) - 0.25_dp*real(j, dp)
         end do
      end do
      C1 = 0.0_dp
      C2 = 0.0_dp
      call pic_symm(S, B, C1, "L", "U")
      call pic_symm(P, B, C2, "L", "U")
      call check(error, maxval(abs(C1 - C2)) < tol)
   end subroutine test_symm_ignores

   subroutine test_symm_right(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: m = 5, n = 4
      real(dp) :: S(n, n), P(n, n), B(m, n), C1(m, n), C2(m, n)
      integer(default_int) :: i, j

      call make_sym(S, P, n)
      do j = 1, n
         do i = 1, m
            B(i, j) = 0.75_dp*real(i, dp) + real(j, dp)
         end do
      end do
      C1 = 0.0_dp
      C2 = 0.0_dp
      call pic_symm(P, B, C1, "R", "U")
      call pic_gemm(B, S, C2)
      call check(error, maxval(abs(C1 - C2)) < tol)
   end subroutine test_symm_right

   subroutine test_symm_alpha_beta(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 5, k = 3
      real(dp) :: S(n, n), P(n, n), B(n, k), C(n, k), C0(n, k), ref(n, k)
      integer(default_int) :: i, j

      call make_sym(S, P, n)
      do j = 1, k
         do i = 1, n
            B(i, j) = real(i*j, dp)*0.1_dp
            C0(i, j) = real(i - j, dp)
         end do
      end do
      ref = 0.0_dp
      call pic_gemm(S, B, ref)
      ref = 2.5_dp*ref + 3.0_dp*C0
      C = C0
      call pic_symm(P, B, C, "L", "U", 2.5_dp, 3.0_dp)
      call check(error, maxval(abs(C - ref)) < tol)
   end subroutine test_symm_alpha_beta

   subroutine test_syrk_vs_gemm(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 5, k = 3
      real(dp) :: A(n, k), C(n, n), ref(n, n)
      integer(default_int) :: i, j

      do j = 1, k
         do i = 1, n
            A(i, j) = real(i, dp) - 0.5_dp*real(j, dp)
         end do
      end do
      ref = 0.0_dp
      call pic_gemm(A, A, ref, "N", "T")
      C = 0.0_dp
      call pic_syrk(A, C, "U", "N")
      !  only the named triangle is written, so compare only that
      do j = 1, n
         do i = 1, j
            call check(error, abs(C(i, j) - ref(i, j)) < tol)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_syrk_vs_gemm

   subroutine test_syrk_trans(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 4, k = 6
      real(dp) :: A(k, n), C(n, n), ref(n, n)
      integer(default_int) :: i, j

      do j = 1, n
         do i = 1, k
            A(i, j) = 0.3_dp*real(i, dp) + real(j, dp)
         end do
      end do
      ref = 0.0_dp
      call pic_gemm(A, A, ref, "T", "N")
      C = 0.0_dp
      call pic_syrk(A, C, "L", "T")
      do j = 1, n
         do i = j, n
            call check(error, abs(C(i, j) - ref(i, j)) < tol)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_syrk_trans

   !  The identity the congruence-transform rewrites depend on: when A'B is
   !  symmetric, SYR2K with alpha=1/2 produces its triangle, because it forms
   !  alpha*(A'B + B'A) = alpha*2*A'B.
   subroutine test_syr2k_identity(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 5, k = 5
      real(dp) :: S(n, n), P(n, n), X(n, n), AB(n, n), Y(n, n), C(n, n)
      integer(default_int) :: i, j

      !  X'SX is symmetric for symmetric S, which is the shape that arises
      call make_sym(S, P, n)
      do j = 1, n
         do i = 1, n
            X(i, j) = real(i, dp)*0.2_dp - real(j, dp)*0.1_dp
         end do
      end do
      Y = 0.0_dp
      call pic_gemm(S, X, Y)          ! Y = S X
      AB = 0.0_dp
      call pic_gemm(X, Y, AB, "T", "N")   ! AB = X' S X, symmetric

      C = 0.0_dp
      call pic_syr2k(X, Y, C, "U", "T", 0.5_dp, 0.0_dp)
      do j = 1, n
         do i = 1, j
            call check(error, abs(C(i, j) - AB(i, j)) < 1.0e-10_dp)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_syr2k_identity

   subroutine test_symm_single(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 4, k = 2
      real(sp) :: S(n, n), B(n, k), C1(n, k), C2(n, k)
      integer(default_int) :: i, j

      do j = 1, n
         do i = 1, j
            S(i, j) = real(i + j, sp)
            S(j, i) = S(i, j)
         end do
      end do
      do j = 1, k
         do i = 1, n
            B(i, j) = real(i - j, sp)
         end do
      end do
      C1 = 0.0_sp
      C2 = 0.0_sp
      call pic_symm(S, B, C1, "L", "U")
      call pic_gemm(S, B, C2)
      call check(error, maxval(abs(C1 - C2)) < 1.0e-4_sp)
   end subroutine test_symm_single

end module test_pic_blas_interfaces_symm
