module test_pic_lapack_interfaces_tpttr
   !! Tests for pic_tpttr / pic_trttp, the packed<->full conversions.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_lapack_interfaces, only: pic_tpttr, pic_trttp
   use pic_types, only: sp, dp, default_int
   implicit none
   private
   public :: collect_pic_lapack_tpttr_tests

contains

   subroutine collect_pic_lapack_tpttr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("tpttr_upper_known", test_tpttr_upper), &
                  new_unittest("tpttr_lower_known", test_tpttr_lower), &
                  new_unittest("tpttr_leaves_other_triangle_alone", test_tpttr_leaves), &
                  new_unittest("trttp_is_inverse_of_tpttr", test_roundtrip), &
                  new_unittest("tpttr_single_precision", test_single) &
                  ]

   end subroutine collect_pic_lapack_tpttr_tests

   subroutine test_tpttr_upper(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 3_default_int
      real(dp) :: AP(n*(n + 1)/2), A(n, n)
      integer(default_int) :: info
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = 0.0_dp
      call pic_tpttr(AP, A, "U", info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, A(1, 1) == 1.0_dp .and. A(1, 2) == 2.0_dp .and. A(2, 2) == 3.0_dp)
      if (allocated(error)) return
      call check(error, A(1, 3) == 4.0_dp .and. A(2, 3) == 5.0_dp .and. A(3, 3) == 6.0_dp)
   end subroutine test_tpttr_upper

   subroutine test_tpttr_lower(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 3_default_int
      real(dp) :: AP(n*(n + 1)/2), A(n, n)
      integer(default_int) :: info
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = 0.0_dp
      call pic_tpttr(AP, A, "L", info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, A(1, 1) == 1.0_dp .and. A(2, 1) == 2.0_dp .and. A(3, 1) == 3.0_dp)
      if (allocated(error)) return
      call check(error, A(2, 2) == 4.0_dp .and. A(3, 2) == 5.0_dp .and. A(3, 3) == 6.0_dp)
   end subroutine test_tpttr_lower

   !  The documented difference from pic_unpack: TPTTR writes one triangle
   !  and does not touch the other. Callers rely on that when the result
   !  feeds SYMM or SYRK, so it is worth pinning down.
   subroutine test_tpttr_leaves(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 3_default_int
      real(dp) :: AP(n*(n + 1)/2), A(n, n)
      integer(default_int) :: info
      AP = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      A = -42.0_dp
      call pic_tpttr(AP, A, "U", info)
      call check(error, A(2, 1) == -42.0_dp)
      if (allocated(error)) return
      call check(error, A(3, 1) == -42.0_dp)
      if (allocated(error)) return
      call check(error, A(3, 2) == -42.0_dp)
   end subroutine test_tpttr_leaves

   subroutine test_roundtrip(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), parameter :: n = 17_default_int
      real(dp) :: AP(n*(n + 1)/2), BP(n*(n + 1)/2), A(n, n)
      integer(default_int) :: i, info
      character(len=1), parameter :: uplos(2) = ["U", "L"]
      integer :: k

      do k = 1, 2
         do i = 1, size(AP, kind=default_int)
            AP(i) = real(i, dp)*1.5_dp - 4.0_dp
         end do
         A = 0.0_dp
         BP = 0.0_dp
         call pic_tpttr(AP, A, uplos(k), info)
         call check(error, info == 0)
         if (allocated(error)) return
         call pic_trttp(A, BP, uplos(k), info)
         call check(error, info == 0)
         if (allocated(error)) return
         call check(error, all(AP == BP))
         if (allocated(error)) return
      end do
   end subroutine test_roundtrip

   subroutine test_single(error)
      type(error_type), allocatable, intent(out) :: error
      real(sp) :: AP(3), A(2, 2)
      integer(default_int) :: info
      AP = [1.0_sp, 2.0_sp, 3.0_sp]
      A = 0.0_sp
      call pic_tpttr(AP, A, "U", info)
      call check(error, info == 0)
      if (allocated(error)) return
      call check(error, A(1, 1) == 1.0_sp .and. A(1, 2) == 2.0_sp .and. A(2, 2) == 3.0_sp)
   end subroutine test_single

end module test_pic_lapack_interfaces_tpttr
