! SPDX-License-Identifier: MIT
! Copyright (c) 2025 Jorge Luis Galvez Vallejo
!! This file contains the interfaces for LAPACK eigenvalue and SVD routines
!! Following the same pattern as pic_blas_interfaces.F90

module pic_lapack_interfaces
   !! pic_lapack_interfaces.F90 provides the interfaces for LAPACK routines
   !! the idea is to have a two level interface, first pic_lapack_xyz which
   !! is the way programmers will use LAPACK, it'll do some checks and then
   !! call the "overloaded" LAPACK interfaces to call the correct routine
   use pic_types, only: sp, dp, default_int
   implicit none
   private

   ! Public overloaded interfaces
   public :: pic_syev, pic_syevd, pic_gesvd

   interface pic_syev
      !! General interface for LAPACK SYEV routines (symmetric eigenvalue problem)
      !!
      !! Usage: call pic_syev(A, W, [optional] jobz, [optional] uplo, [optional] info)
      !!
      !! where A is a symmetric matrix, W is the output eigenvalue array.
      !! If jobz='V' (default), A is overwritten with eigenvectors on output.
      !! If jobz='N', only eigenvalues are computed.
      !! uplo='U' (default) means upper triangle of A is stored.
      !!
      module procedure :: pic_ssyev
      module procedure :: pic_dsyev
   end interface pic_syev

   interface pic_syevd
      !! General interface for LAPACK SYEVD routines (divide-and-conquer eigenvalue)
      !!
      !! Usage: call pic_syevd(A, W, [optional] jobz, [optional] uplo, [optional] info)
      !!
      !! Same as pic_syev but uses divide-and-conquer algorithm.
      !! Faster for large matrices but uses more memory.
      !!
      module procedure :: pic_ssyevd
      module procedure :: pic_dsyevd
   end interface pic_syevd

   interface pic_gesvd
      !! General interface for LAPACK GESVD routines (singular value decomposition)
      !!
      !! Usage: call pic_gesvd(A, S, [optional] U, [optional] VT, [optional] jobu, [optional] jobvt, [optional] info)
      !!
      !! Computes the SVD of A: A = U * SIGMA * VT
      !! where S contains the singular values (diagonal of SIGMA).
      !! A is destroyed on output. U and VT are optional outputs.
      !!
      module procedure :: pic_sgesvd
      module procedure :: pic_dgesvd
   end interface pic_gesvd

   ! Low-level LAPACK interfaces (not public)
   interface lapack_syev
      !! Explicit interface for LAPACK SYEV routines
      subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         real(sp), intent(inout) :: a(lda, *)
         integer(default_int), intent(in) :: lda
         real(sp), intent(out) :: w(*)
         real(sp), intent(out) :: work(*)
         integer(default_int), intent(in) :: lwork
         integer(default_int), intent(out) :: info
      end subroutine ssyev
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         real(dp), intent(inout) :: a(lda, *)
         integer(default_int), intent(in) :: lda
         real(dp), intent(out) :: w(*)
         real(dp), intent(out) :: work(*)
         integer(default_int), intent(in) :: lwork
         integer(default_int), intent(out) :: info
      end subroutine dsyev
   end interface lapack_syev

   interface lapack_syevd
      !! Explicit interface for LAPACK SYEVD routines
      subroutine ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         real(sp), intent(inout) :: a(lda, *)
         integer(default_int), intent(in) :: lda
         real(sp), intent(out) :: w(*)
         real(sp), intent(out) :: work(*)
         integer(default_int), intent(in) :: lwork
         integer(default_int), intent(out) :: iwork(*)
         integer(default_int), intent(in) :: liwork
         integer(default_int), intent(out) :: info
      end subroutine ssyevd
      subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         real(dp), intent(inout) :: a(lda, *)
         integer(default_int), intent(in) :: lda
         real(dp), intent(out) :: w(*)
         real(dp), intent(out) :: work(*)
         integer(default_int), intent(in) :: lwork
         integer(default_int), intent(out) :: iwork(*)
         integer(default_int), intent(in) :: liwork
         integer(default_int), intent(out) :: info
      end subroutine dsyevd
   end interface lapack_syevd

   interface lapack_gesvd
      !! Explicit interface for LAPACK GESVD routines
      subroutine sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: jobu
         character(len=1), intent(in) :: jobvt
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         real(sp), intent(inout) :: a(lda, *)
         integer(default_int), intent(in) :: lda
         real(sp), intent(out) :: s(*)
         real(sp), intent(out) :: u(ldu, *)
         integer(default_int), intent(in) :: ldu
         real(sp), intent(out) :: vt(ldvt, *)
         integer(default_int), intent(in) :: ldvt
         real(sp), intent(out) :: work(*)
         integer(default_int), intent(in) :: lwork
         integer(default_int), intent(out) :: info
      end subroutine sgesvd
      subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: jobu
         character(len=1), intent(in) :: jobvt
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         real(dp), intent(inout) :: a(lda, *)
         integer(default_int), intent(in) :: lda
         real(dp), intent(out) :: s(*)
         real(dp), intent(out) :: u(ldu, *)
         integer(default_int), intent(in) :: ldu
         real(dp), intent(out) :: vt(ldvt, *)
         integer(default_int), intent(in) :: ldvt
         real(dp), intent(out) :: work(*)
         integer(default_int), intent(in) :: lwork
         integer(default_int), intent(out) :: info
      end subroutine dgesvd
   end interface lapack_gesvd

contains

   ! ============================================================================
   ! SYEV wrappers (symmetric eigenvalue problem)
   ! ============================================================================

   subroutine pic_ssyev(A, W, jobz, uplo, info)
      !! Single precision symmetric eigenvalue problem
      !! Computes eigenvalues (and optionally eigenvectors) of symmetric matrix A
      real(sp), intent(inout) :: A(:, :)
      real(sp), intent(out) :: W(:)
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_jobz, l_uplo
      integer(default_int) :: n, lda, lwork, l_info
      real(sp), allocatable :: work(:)
      real(sp) :: work_query(1)

      ! Set defaults
      l_jobz = 'V'
      if (present(jobz)) l_jobz = jobz
      l_uplo = 'U'
      if (present(uplo)) l_uplo = uplo

      ! Get dimensions
      n = size(A, 1)
      lda = max(1, n)

      ! Query optimal workspace size
      lwork = -1
      call lapack_syev(l_jobz, l_uplo, n, A, lda, W, work_query, lwork, l_info)
      lwork = int(work_query(1))

      ! Allocate workspace and compute
      allocate (work(lwork))
      call lapack_syev(l_jobz, l_uplo, n, A, lda, W, work, lwork, l_info)

      if (present(info)) info = l_info
   end subroutine pic_ssyev

   subroutine pic_dsyev(A, W, jobz, uplo, info)
      !! Double precision symmetric eigenvalue problem
      !! Computes eigenvalues (and optionally eigenvectors) of symmetric matrix A
      real(dp), intent(inout) :: A(:, :)
      real(dp), intent(out) :: W(:)
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_jobz, l_uplo
      integer(default_int) :: n, lda, lwork, l_info
      real(dp), allocatable :: work(:)
      real(dp) :: work_query(1)

      ! Set defaults
      l_jobz = 'V'
      if (present(jobz)) l_jobz = jobz
      l_uplo = 'U'
      if (present(uplo)) l_uplo = uplo

      ! Get dimensions
      n = size(A, 1)
      lda = max(1, n)

      ! Query optimal workspace size
      lwork = -1
      call lapack_syev(l_jobz, l_uplo, n, A, lda, W, work_query, lwork, l_info)
      lwork = int(work_query(1))

      ! Allocate workspace and compute
      allocate (work(lwork))
      call lapack_syev(l_jobz, l_uplo, n, A, lda, W, work, lwork, l_info)

      if (present(info)) info = l_info
   end subroutine pic_dsyev

   ! ============================================================================
   ! SYEVD wrappers (divide-and-conquer eigenvalue problem)
   ! ============================================================================

   subroutine pic_ssyevd(A, W, jobz, uplo, info)
      !! Single precision divide-and-conquer eigenvalue problem
      !! Computes eigenvalues (and optionally eigenvectors) of symmetric matrix A
      real(sp), intent(inout) :: A(:, :)
      real(sp), intent(out) :: W(:)
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_jobz, l_uplo
      integer(default_int) :: n, lda, lwork, liwork, l_info
      real(sp), allocatable :: work(:)
      integer(default_int), allocatable :: iwork(:)
      real(sp) :: work_query(1)
      integer(default_int) :: iwork_query(1)

      ! Set defaults
      l_jobz = 'V'
      if (present(jobz)) l_jobz = jobz
      l_uplo = 'U'
      if (present(uplo)) l_uplo = uplo

      ! Get dimensions
      n = size(A, 1)
      lda = max(1, n)

      ! Query optimal workspace size
      lwork = -1
      liwork = -1
      call lapack_syevd(l_jobz, l_uplo, n, A, lda, W, work_query, lwork, iwork_query, liwork, l_info)
      lwork = int(work_query(1))
      liwork = iwork_query(1)

      ! Allocate workspace and compute
      allocate (work(lwork), iwork(liwork))
      call lapack_syevd(l_jobz, l_uplo, n, A, lda, W, work, lwork, iwork, liwork, l_info)

      if (present(info)) info = l_info
   end subroutine pic_ssyevd

   subroutine pic_dsyevd(A, W, jobz, uplo, info)
      !! Double precision divide-and-conquer eigenvalue problem
      !! Computes eigenvalues (and optionally eigenvectors) of symmetric matrix A
      real(dp), intent(inout) :: A(:, :)
      real(dp), intent(out) :: W(:)
      character(len=1), intent(in), optional :: jobz
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_jobz, l_uplo
      integer(default_int) :: n, lda, lwork, liwork, l_info
      real(dp), allocatable :: work(:)
      integer(default_int), allocatable :: iwork(:)
      real(dp) :: work_query(1)
      integer(default_int) :: iwork_query(1)

      ! Set defaults
      l_jobz = 'V'
      if (present(jobz)) l_jobz = jobz
      l_uplo = 'U'
      if (present(uplo)) l_uplo = uplo

      ! Get dimensions
      n = size(A, 1)
      lda = max(1, n)

      ! Query optimal workspace size
      lwork = -1
      liwork = -1
      call lapack_syevd(l_jobz, l_uplo, n, A, lda, W, work_query, lwork, iwork_query, liwork, l_info)
      lwork = int(work_query(1))
      liwork = iwork_query(1)

      ! Allocate workspace and compute
      allocate (work(lwork), iwork(liwork))
      call lapack_syevd(l_jobz, l_uplo, n, A, lda, W, work, lwork, iwork, liwork, l_info)

      if (present(info)) info = l_info
   end subroutine pic_dsyevd

   ! ============================================================================
   ! GESVD wrappers (singular value decomposition)
   ! ============================================================================

   subroutine pic_sgesvd(A, S, U, VT, jobu, jobvt, info)
      !! Single precision SVD: A = U * SIGMA * VT
      !! A is destroyed on output
      real(sp), intent(inout) :: A(:, :)
      real(sp), intent(out) :: S(:)
      real(sp), intent(out), optional :: U(:, :)
      real(sp), intent(out), optional :: VT(:, :)
      character(len=1), intent(in), optional :: jobu
      character(len=1), intent(in), optional :: jobvt
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_jobu, l_jobvt
      integer(default_int) :: m, n, lda, ldu, ldvt, lwork, l_info
      real(sp), allocatable :: work(:), u_local(:, :), vt_local(:, :)
      real(sp) :: work_query(1)

      ! Get dimensions
      m = size(A, 1)
      n = size(A, 2)
      lda = max(1, m)

      ! Set defaults for job options based on presence of output arrays
      if (present(jobu)) then
         l_jobu = jobu
      else if (present(U)) then
         l_jobu = 'S'  ! Compute first min(m,n) columns of U
      else
         l_jobu = 'N'  ! Don't compute U
      end if

      if (present(jobvt)) then
         l_jobvt = jobvt
      else if (present(VT)) then
         l_jobvt = 'S'  ! Compute first min(m,n) rows of VT
      else
         l_jobvt = 'N'  ! Don't compute VT
      end if

      ! Set up U workspace
      if (present(U)) then
         ldu = max(1, size(U, 1))
      else
         ldu = 1
         allocate (u_local(1, 1))
      end if

      ! Set up VT workspace
      if (present(VT)) then
         ldvt = max(1, size(VT, 1))
      else
         ldvt = 1
         allocate (vt_local(1, 1))
      end if

      ! Query optimal workspace size
      lwork = -1
      if (present(U) .and. present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, &
              & work_query, lwork, l_info)
      else if (present(U)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, vt_local, ldvt, &
              & work_query, lwork, l_info)
      else if (present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, VT, ldvt, &
              & work_query, lwork, l_info)
      else
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, vt_local, ldvt, &
              & work_query, lwork, l_info)
      end if
      lwork = int(work_query(1))

      ! Allocate workspace and compute
      allocate (work(lwork))
      if (present(U) .and. present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, &
              & work, lwork, l_info)
      else if (present(U)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, vt_local, ldvt, &
              & work, lwork, l_info)
      else if (present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, VT, ldvt, &
              & work, lwork, l_info)
      else
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, vt_local, ldvt, &
              & work, lwork, l_info)
      end if

      if (present(info)) info = l_info
   end subroutine pic_sgesvd

   subroutine pic_dgesvd(A, S, U, VT, jobu, jobvt, info)
      !! Double precision SVD: A = U * SIGMA * VT
      !! A is destroyed on output
      real(dp), intent(inout) :: A(:, :)
      real(dp), intent(out) :: S(:)
      real(dp), intent(out), optional :: U(:, :)
      real(dp), intent(out), optional :: VT(:, :)
      character(len=1), intent(in), optional :: jobu
      character(len=1), intent(in), optional :: jobvt
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_jobu, l_jobvt
      integer(default_int) :: m, n, lda, ldu, ldvt, lwork, l_info
      real(dp), allocatable :: work(:), u_local(:, :), vt_local(:, :)
      real(dp) :: work_query(1)

      ! Get dimensions
      m = size(A, 1)
      n = size(A, 2)
      lda = max(1, m)

      ! Set defaults for job options based on presence of output arrays
      if (present(jobu)) then
         l_jobu = jobu
      else if (present(U)) then
         l_jobu = 'S'  ! Compute first min(m,n) columns of U
      else
         l_jobu = 'N'  ! Don't compute U
      end if

      if (present(jobvt)) then
         l_jobvt = jobvt
      else if (present(VT)) then
         l_jobvt = 'S'  ! Compute first min(m,n) rows of VT
      else
         l_jobvt = 'N'  ! Don't compute VT
      end if

      ! Set up U workspace
      if (present(U)) then
         ldu = max(1, size(U, 1))
      else
         ldu = 1
         allocate (u_local(1, 1))
      end if

      ! Set up VT workspace
      if (present(VT)) then
         ldvt = max(1, size(VT, 1))
      else
         ldvt = 1
         allocate (vt_local(1, 1))
      end if

      ! Query optimal workspace size
      lwork = -1
      if (present(U) .and. present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, &
              & work_query, lwork, l_info)
      else if (present(U)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, vt_local, ldvt, &
              & work_query, lwork, l_info)
      else if (present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, VT, ldvt, &
              & work_query, lwork, l_info)
      else
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, vt_local, ldvt, &
              & work_query, lwork, l_info)
      end if
      lwork = int(work_query(1))

      ! Allocate workspace and compute
      allocate (work(lwork))
      if (present(U) .and. present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, VT, ldvt, &
              & work, lwork, l_info)
      else if (present(U)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, U, ldu, vt_local, ldvt, &
              & work, lwork, l_info)
      else if (present(VT)) then
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, VT, ldvt, &
              & work, lwork, l_info)
      else
         call lapack_gesvd(l_jobu, l_jobvt, m, n, A, lda, S, u_local, ldu, vt_local, ldvt, &
              & work, lwork, l_info)
      end if

      if (present(info)) info = l_info
   end subroutine pic_dgesvd

end module pic_lapack_interfaces
