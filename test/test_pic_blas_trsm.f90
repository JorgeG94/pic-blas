module test_pic_blas_interfaces_trsm
   !! Tests for pic_trsm.
   !!
   !! The property that matters to callers is that the solve and an explicit
   !! multiplication by the inverse agree, because the whole point of reaching
   !! for this is to avoid forming that inverse. Each side/transpose
   !! combination is checked against a matrix product built independently.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_blas_interfaces, only: pic_trsm, pic_gemm
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_trsm_tests

   real(dp), parameter :: tol = 1.0e-11_dp

contains

   subroutine collect_pic_trsm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("left_upper_solves", test_left_upper), &
                  new_unittest("right_upper_solves", test_right_upper), &
                  new_unittest("transpose_solves", test_transpose), &
                  new_unittest("lower_triangle_solves", test_lower), &
                  new_unittest("unit_diagonal_ignores_stored_diagonal", test_unit_diag), &
                  new_unittest("alpha_scales_the_solution", test_alpha), &
                  new_unittest("untouched_triangle_is_not_read", test_other_triangle), &
                  new_unittest("trsm_single_precision", test_single) &
                  ]
   end subroutine collect_pic_trsm_tests

   subroutine make_upper(u)
      real(dp), intent(out) :: u(3, 3)
      u = 0.0_dp
      u(1, :) = [2.0_dp, 1.0_dp, -1.0_dp]
      u(2, 2:) = [3.0_dp, 2.0_dp]
      u(3, 3) = 4.0_dp
   end subroutine make_upper

   subroutine make_b(b)
      real(dp), intent(out) :: b(3, 2)
      b(1, :) = [1.0_dp, 4.0_dp]
      b(2, :) = [2.0_dp, -1.0_dp]
      b(3, :) = [3.0_dp, 0.5_dp]
   end subroutine make_b

   subroutine test_left_upper(error)
      !! X = U**-1 B, checked by multiplying back: U X must equal B.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: u(3, 3), b(3, 2), original(3, 2), back(3, 2)

      call make_upper(u)
      call make_b(b)
      original = b
      call pic_trsm(u, b, side="L", uplo="U")
      call pic_gemm(u, b, back)
      call check(error, maxval(abs(back - original)) < tol)
   end subroutine test_left_upper

   subroutine test_right_upper(error)
      !! X = B U**-1, so X U must equal B. This is the shape density fitting
      !! uses -- a tall B against a small triangle.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: u(3, 3), b(4, 3), original(4, 3), back(4, 3)
      integer :: i

      call make_upper(u)
      do i = 1, 4
         b(i, :) = [real(i, dp), 2.0_dp*i - 1.0_dp, 0.5_dp*i]
      end do
      original = b
      call pic_trsm(u, b, side="R", uplo="U")
      call pic_gemm(b, u, back)
      call check(error, maxval(abs(back - original)) < tol)
   end subroutine test_right_upper

   subroutine test_transpose(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: u(3, 3), b(3, 2), original(3, 2), back(3, 2)

      call make_upper(u)
      call make_b(b)
      original = b
      call pic_trsm(u, b, side="L", uplo="U", transa="T")
      call pic_gemm(u, b, back, transa="T")
      call check(error, maxval(abs(back - original)) < tol)
   end subroutine test_transpose

   subroutine test_lower(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: l(3, 3), b(3, 2), original(3, 2), back(3, 2)

      l = 0.0_dp
      l(1, 1) = 2.0_dp
      l(2, :2) = [1.0_dp, 3.0_dp]
      l(3, :) = [-1.0_dp, 2.0_dp, 4.0_dp]
      call make_b(b)
      original = b
      call pic_trsm(l, b, side="L", uplo="L")
      call pic_gemm(l, b, back)
      call check(error, maxval(abs(back - original)) < tol)
   end subroutine test_lower

   subroutine test_unit_diag(error)
      !! diag="U" must not read the stored diagonal. Poisoning it with values
      !! that would change the answer is the only way to prove that.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: u(3, 3), poisoned(3, 3), b(3, 2), first(3, 2), second(3, 2)

      call make_upper(u)
      poisoned = u
      poisoned(1, 1) = 99.0_dp
      poisoned(2, 2) = -7.0_dp
      poisoned(3, 3) = 0.0_dp

      call make_b(b)
      first = b
      call pic_trsm(u, first, side="L", uplo="U", diag="U")
      call make_b(b)
      second = b
      call pic_trsm(poisoned, second, side="L", uplo="U", diag="U")
      call check(error, maxval(abs(first - second)) < tol)
   end subroutine test_unit_diag

   subroutine test_alpha(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: u(3, 3), b(3, 2), plain(3, 2), scaled(3, 2)

      call make_upper(u)
      call make_b(b)
      plain = b
      call pic_trsm(u, plain, side="L", uplo="U")
      call make_b(b)
      scaled = b
      call pic_trsm(u, scaled, side="L", uplo="U", alpha=2.5_dp)
      call check(error, maxval(abs(scaled - 2.5_dp*plain)) < tol)
   end subroutine test_alpha

   subroutine test_other_triangle(error)
      !! Junk below the diagonal must not reach the answer. A caller handing
      !! this a full symmetric matrix and naming one triangle is the normal
      !! case, so reading the other one would be silently wrong.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: u(3, 3), junked(3, 3), b(3, 2), first(3, 2), second(3, 2)

      call make_upper(u)
      junked = u
      junked(2, 1) = 13.0_dp
      junked(3, 1) = -5.0_dp
      junked(3, 2) = 8.0_dp

      call make_b(b)
      first = b
      call pic_trsm(u, first, side="L", uplo="U")
      call make_b(b)
      second = b
      call pic_trsm(junked, second, side="L", uplo="U")
      call check(error, maxval(abs(first - second)) < tol)
   end subroutine test_other_triangle

   subroutine test_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: u(2, 2), b(2, 1), back(2, 1)

      u = 0.0_sp
      u(1, :) = [2.0_sp, 1.0_sp]
      u(2, 2) = 4.0_sp
      b(1, 1) = 3.0_sp
      b(2, 1) = 8.0_sp
      call pic_trsm(u, b, side="L", uplo="U")
      call pic_gemm(u, b, back)
      call check(error, abs(back(1, 1) - 3.0_sp) < 1.0e-5_sp)
      if (allocated(error)) return
      call check(error, abs(back(2, 1) - 8.0_sp) < 1.0e-5_sp)
   end subroutine test_single

end module test_pic_blas_interfaces_trsm
