module test_pic_lapack_interfaces_unpack
   !! Tests for pic_unpack / pic_unpack_x.
   !!
   !! The routine dispatches on n across three regimes, so a test on a 4x4
   !! proves only that one of the three loop nests is right. Every case here
   !! that can be run at more than one size is run at one n per regime, taken
   !! from either side of the two boundaries.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_lapack_interfaces, only: pic_unpack, pic_unpack_x
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_lapack_unpack_tests

   !  One n in each regime, and both sides of each boundary. The dispatch
   !  points are 56 and 2560; getting these wrong is how a rewrite of one
   !  branch passes while another stays broken.
   integer(default_int), parameter :: regime_sizes(6) = &
      [7_default_int, 55_default_int, 56_default_int, &
       200_default_int, 2559_default_int, 2560_default_int]

contains

   subroutine collect_pic_lapack_unpack_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("unpack_symmetric_upper_known", test_sym_upper_known), &
                  new_unittest("unpack_antisymmetric_upper_known", test_anti_upper_known), &
                  new_unittest("unpack_symmetric_lower_known", test_sym_lower_known), &
                  new_unittest("unpack_antisymmetric_lower_known", test_anti_lower_known), &
                  new_unittest("unpack_all_regimes_vs_reference", test_all_regimes), &
                  new_unittest("unpack_x_matches_unpack", test_x_matches), &
                  new_unittest("unpack_defaults_to_symmetric_upper", test_defaults), &
                  new_unittest("unpack_antisymmetric_diagonal_is_zero", test_anti_diagonal), &
                  new_unittest("unpack_block_size_does_not_change_result", test_nb_invariant), &
                  new_unittest("unpack_single_precision", test_single) &
                  ]

   end subroutine collect_pic_lapack_unpack_tests

   !  Reference expansion, written from the storage definition rather than
   !  from the implementation, so a shared misreading of the packing cannot
   !  make both agree.
   subroutine reference(AP, R, n, mode, uplo)
      real(dp), intent(in) :: AP(:)
      real(dp), intent(out) :: R(:, :)
      integer(default_int), intent(in) :: n
      character(len=1), intent(in) :: mode, uplo
      integer(default_int) :: i, j
      real(dp) :: s
      s = 1.0_dp
      if (mode == "A") s = -1.0_dp
      if (uplo == "U") then
         do j = 1, n
            do i = 1, j
               R(i, j) = AP((j*(j - 1))/2 + i)
               R(j, i) = s*AP((j*(j - 1))/2 + i)
            end do
         end do
      else
         do j = 1, n
            do i = j, n
               R(i, j) = AP(((j - 1)*(2*n - j))/2 + i)
               R(j, i) = s*AP(((j - 1)*(2*n - j))/2 + i)
            end do
         end do
      end if
      if (mode == "A") then
         do j = 1, n
            R(j, j) = 0.0_dp
         end do
      end if
   end subroutine reference

   subroutine fill(AP)
      real(dp), intent(out) :: AP(:)
      integer(default_int) :: i
      do i = 1, size(AP, kind=default_int)
         AP(i) = real(i, dp)*0.25_dp - 7.0_dp
      end do
   end subroutine fill

   !  A 3x3 worked by hand, so the tests do not rest entirely on a reference
   !  that could share a mistake with the code. Packed upper is
   !  [a11, a12, a22, a13, a23, a33].
   subroutine test_sym_upper_known(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: AP(6), A(3, 3)
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = -99.0_dp
      call pic_unpack(AP, A, "S", "U")
      call check(error, all(A(1, :) == [1.0_dp, 2.0_dp, 4.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(2, :) == [2.0_dp, 3.0_dp, 5.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(3, :) == [4.0_dp, 5.0_dp, 6.0_dp]))
   end subroutine test_sym_upper_known

   subroutine test_anti_upper_known(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: AP(6), A(3, 3)
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = -99.0_dp
      call pic_unpack(AP, A, "A", "U")
      call check(error, all(A(1, :) == [0.0_dp, 2.0_dp, 4.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(2, :) == [-2.0_dp, 0.0_dp, 5.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(3, :) == [-4.0_dp, -5.0_dp, 0.0_dp]))
   end subroutine test_anti_upper_known

   !  Packed lower is [a11, a21, a31, a22, a32, a33].
   subroutine test_sym_lower_known(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: AP(6), A(3, 3)
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = -99.0_dp
      call pic_unpack(AP, A, "S", "L")
      call check(error, all(A(:, 1) == [1.0_dp, 2.0_dp, 3.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(:, 2) == [2.0_dp, 4.0_dp, 5.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(:, 3) == [3.0_dp, 5.0_dp, 6.0_dp]))
   end subroutine test_sym_lower_known

   subroutine test_anti_lower_known(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: AP(6), A(3, 3)
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = -99.0_dp
      call pic_unpack(AP, A, "A", "L")
      call check(error, all(A(:, 1) == [0.0_dp, 2.0_dp, 3.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(:, 2) == [-2.0_dp, 0.0_dp, 5.0_dp]))
      if (allocated(error)) return
      call check(error, all(A(:, 3) == [-3.0_dp, -5.0_dp, 0.0_dp]))
   end subroutine test_anti_lower_known

   !  Both modes and both packings at a size in each dispatch regime. The
   !  expansion is a copy, so these must agree exactly, not to a tolerance.
   subroutine test_all_regimes(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: AP(:), A(:, :), R(:, :)
      character(len=1) :: modes(2) = ["S", "A"], uplos(2) = ["U", "L"]
      integer(default_int) :: n
      integer :: t, im, iu

      do t = 1, size(regime_sizes)
         n = regime_sizes(t)
         allocate (AP(n*(n + 1)/2), A(n, n), R(n, n))
         call fill(AP)
         do im = 1, 2
            do iu = 1, 2
               A = -huge(1.0_dp)
               call reference(AP, R, n, modes(im), uplos(iu))
               call pic_unpack(AP, A, modes(im), uplos(iu))
               call check(error, all(A == R))
               if (allocated(error)) then
                  deallocate (AP, A, R)
                  return
               end if
            end do
         end do
         deallocate (AP, A, R)
      end do
   end subroutine test_all_regimes

   !  The explicit-dimension entry runs the same kernel and must agree bit
   !  for bit with the assumed-shape one, in every regime.
   subroutine test_x_matches(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: AP(:), A(:, :), B(:, :)
      character(len=1) :: modes(2) = ["S", "A"], uplos(2) = ["U", "L"]
      integer(default_int) :: n
      integer :: t, im, iu

      do t = 1, size(regime_sizes)
         n = regime_sizes(t)
         allocate (AP(n*(n + 1)/2), A(n, n), B(n, n))
         call fill(AP)
         do im = 1, 2
            do iu = 1, 2
               A = -huge(1.0_dp)
               B = huge(1.0_dp)
               call pic_unpack(AP, A, modes(im), uplos(iu))
               call pic_unpack_x(n, AP, B, n, modes(im), uplos(iu), 8_default_int)
               call check(error, all(A == B))
               if (allocated(error)) then
                  deallocate (AP, A, B)
                  return
               end if
            end do
         end do
         deallocate (AP, A, B)
      end do
   end subroutine test_x_matches

   subroutine test_defaults(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: AP(6), A(3, 3), B(3, 3)
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      call pic_unpack(AP, A)
      call pic_unpack(AP, B, "S", "U")
      call check(error, all(A == B))
   end subroutine test_defaults

   !  The antisymmetric diagonal is zero regardless of what the packed
   !  diagonal held, which is the part a mirror-and-negate gets wrong.
   subroutine test_anti_diagonal(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: AP(:), A(:, :)
      integer(default_int) :: n, j
      integer :: t

      do t = 1, size(regime_sizes)
         n = regime_sizes(t)
         allocate (AP(n*(n + 1)/2), A(n, n))
         AP = 3.5_dp
         call pic_unpack(AP, A, "A", "U")
         do j = 1, n
            call check(error, A(j, j) == 0.0_dp)
            if (allocated(error)) then
               deallocate (AP, A)
               return
            end if
         end do
         deallocate (AP, A)
      end do
   end subroutine test_anti_diagonal

   !  Block size is a performance knob and must not be able to change an
   !  answer. Only the largest regime uses it, so this runs there.
   subroutine test_nb_invariant(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: AP(:), A(:, :), B(:, :)
      integer(default_int), parameter :: n = 2600_default_int

      allocate (AP(n*(n + 1)/2), A(n, n), B(n, n))
      call fill(AP)
      call pic_unpack(AP, A, "A", "L", 8_default_int)
      call pic_unpack(AP, B, "A", "L", 97_default_int)
      call check(error, all(A == B))
      deallocate (AP, A, B)
   end subroutine test_nb_invariant

   subroutine test_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: AP(6), A(3, 3)
      AP = [1.0_sp, 2.0_sp, 3.0_sp, 4.0_sp, 5.0_sp, 6.0_sp]
      call pic_unpack(AP, A, "A", "U")
      call check(error, all(A(1, :) == [0.0_sp, 2.0_sp, 4.0_sp]))
      if (allocated(error)) return
      call check(error, all(A(3, :) == [-4.0_sp, -5.0_sp, 0.0_sp]))
   end subroutine test_single

end module test_pic_lapack_interfaces_unpack
