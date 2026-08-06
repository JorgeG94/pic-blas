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

   !> Size dispatch for UNPACK. Three regimes, each measured against the
   !> obvious interleaved loop -- the one LAPACK does not provide and every
   !> legacy code writes by hand, storing A(i,j) and A(j,i) in the same
   !> iteration. All figures are single-threaded gfortran -O2, which is how
   !> this library is normally built.
   !>
   !>   n < 56     nothing beats the loop itself. It is two stores per
   !>              element running at about one store per cycle, which is the
   !>              store port's limit, so there is no headroom to trade.
   !>              Eight variants were tried: unrolling two columns at a time
   !>              (0.93-1.01x), filling rows before mirroring (0.86-0.99x),
   !>              walking the destination column (0.44-0.61x) and tiling the
   !>              mirror (0.40-0.60x), the last two because a loop nest costs
   !>              more than a 16x16 matrix is worth. So this range runs the
   !>              same loop, and the only thing that matters here is not
   !>              spending anything around it -- see the note on the diagonal
   !>              in the small branch below.
   !>
   !>   n < 2560   walk the destination column instead of the source triangle.
   !>              Every store is then contiguous and only the load is
   !>              strided, which the prefetcher handles far better than a
   !>              scattered store handles anything. 1.1x at n=56, 1.6x at
   !>              64-96, 4.5x at 128, 8.5x at 256, 11.0x at 512, 10.4x at
   !>              1024, 7.0x at 2048.
   !>
   !>   otherwise  fill the stored triangle contiguously, then mirror it
   !>              in tiles. Past about 2560 the column walk's strided load
   !>              spans more than cache can hold and decays -- 3.1x at 4096,
   !>              2.6x at 6144 -- while this holds 4.1-4.9x.
   !>
   !> Block size for the mirror. 8 doubles is one 64-byte cache line, so an
   !> 8x8 tile fills each line of the transposed stream completely before it
   !> can be evicted. Measured against 16/24/32/48/64/128/256 at n from 256 to
   !> 4096: 8 wins at every size, and the penalty for getting this wrong is
   !> large -- 128 gives only 1.03x where 8 gives 2.6x.
   integer, parameter, private :: unpack_small_max = 56
   integer, parameter, private :: unpack_colwise_max = 2560
   integer, parameter, private :: unpack_nb = 8

   !> Threading pays only well above the serial crossovers, and pic is
   !> usually built without OpenMP, so the serial paths above are what
   !> normally runs. Applied through the if() clause so there is one code
   !> path rather than two, and so it degrades to a comment without -fopenmp.
   integer, parameter, private :: unpack_par_min = 256

   ! Public overloaded interfaces
   public :: pic_syev, pic_syevd, pic_gesvd
   public :: pic_tpttr, pic_trttp, pic_unpack, pic_unpack_x

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

   interface pic_tpttr
      !! packed triangular storage -> full storage
      !!
      !! Usage: call pic_tpttr(AP, A, [uplo], [info])
      !!
      !! Note what this does *not* do: it writes only the triangle named by
      !! uplo, and leaves the rest of A alone. It does not mirror into a full
      !! symmetric matrix.
      !!
      !! That is usually what you want. Code that unpacks a symmetric matrix
      !! and mirrors it is generally doing so because the next call is a GEMM,
      !! which needs every element. Reach for pic_symm, pic_syrk or pic_syr2k
      !! instead and only one triangle is ever read, so the mirror is wasted
      !! work on top of a multiply that was already doing twice what it needed.
      module procedure :: pic_stpttr
      module procedure :: pic_dtpttr
   end interface pic_tpttr

   interface pic_trttp
      !! full storage -> packed triangular storage
      !!
      !! Usage: call pic_trttp(A, AP, [uplo], [info])
      !!
      !! The inverse of pic_tpttr; reads the uplo triangle of A.
      module procedure :: pic_strttp
      module procedure :: pic_dtrttp
   end interface pic_trttp

   interface pic_unpack
      !! packed triangular storage -> a *full* matrix, mirrored
      !!
      !! Usage: call pic_unpack(AP, A, [mode], [uplo], [nb])
      !!
      !!   mode "S"  symmetric:     A(j,i) =  A(i,j)
      !!   mode "A"  antisymmetric: A(j,i) = -A(i,j), diagonal zeroed
      !!
      !! LAPACK has no equivalent -- TPTTR fills one triangle and stops, which
      !! is the right thing when the result feeds SYMM/SYRK/SYR2K. This is for
      !! the cases that genuinely need every element, and for the
      !! antisymmetric case, which has no packed BLAS operation at all.
      !!
      !! Dispatches on n across three regimes; see unpack_small_max for what
      !! each does and what it was measured at. Against the loop it replaces,
      !! serial: parity to n=48, then 1.1x at 56, 1.6x at 64-96, 4.5x at 128,
      !! 8.5x at 256, 11.0x at 512, 10.4x at 1024, 4.1x at 4096.
      !!
      !! At n below about 32 the array descriptors this interface builds cost
      !! more than the kernel does, and it runs at 0.7-0.97x. pic_unpack_x
      !! takes explicit dimensions instead and holds 0.93x at n=4 and 1.0x
      !! from n=8; call that one from a loop over small blocks.
      module procedure :: pic_sunpack
      module procedure :: pic_dunpack
   end interface pic_unpack

   interface pic_unpack_x
      !! pic_unpack with explicit dimensions, in the shape the BLAS uses
      !!
      !! Usage: call pic_unpack_x(n, AP, A, lda, mode, uplo, nb)
      !!
      !! Same operation and same size dispatch, but order and leading
      !! dimension are passed and nothing is optional, so no array descriptor
      !! is built at the call. Worth about 5 ns a call, which matters only
      !! when n is small enough for that to be a large share of the work: at
      !! n=4 it is the difference between 0.69x and 0.93x, and by n=64 the
      !! two are indistinguishable.
      module procedure :: pic_sunpack_x
      module procedure :: pic_dunpack_x
   end interface pic_unpack_x

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

   interface lapack_tpttr
      !! not a public interface, used internally by pic_tpttr
      pure subroutine stpttr(uplo, n, ap, a, lda, info)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         real(sp), intent(in) :: ap(*)
         real(sp), intent(inout) :: a(lda, *)
         integer(default_int), intent(out) :: info
      end subroutine stpttr
      pure subroutine dtpttr(uplo, n, ap, a, lda, info)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: a(lda, *)
         integer(default_int), intent(out) :: info
      end subroutine dtpttr
   end interface lapack_tpttr

   interface lapack_trttp
      !! not a public interface, used internally by pic_trttp
      pure subroutine strttp(uplo, n, a, lda, ap, info)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(out) :: ap(*)
         integer(default_int), intent(out) :: info
      end subroutine strttp
      pure subroutine dtrttp(uplo, n, a, lda, ap, info)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(out) :: ap(*)
         integer(default_int), intent(out) :: info
      end subroutine dtrttp
   end interface lapack_trttp

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
      l_jobz = "V"
      if (present(jobz)) l_jobz = jobz
      l_uplo = "U"
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
      l_jobz = "V"
      if (present(jobz)) l_jobz = jobz
      l_uplo = "U"
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
      l_jobz = "V"
      if (present(jobz)) l_jobz = jobz
      l_uplo = "U"
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
      l_jobz = "V"
      if (present(jobz)) l_jobz = jobz
      l_uplo = "U"
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
         l_jobu = "S"  ! Compute first min(m,n) columns of U
      else
         l_jobu = "N"  ! Don't compute U
      end if

      if (present(jobvt)) then
         l_jobvt = jobvt
      else if (present(VT)) then
         l_jobvt = "S"  ! Compute first min(m,n) rows of VT
      else
         l_jobvt = "N"  ! Don't compute VT
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
         l_jobu = "S"  ! Compute first min(m,n) columns of U
      else
         l_jobu = "N"  ! Don't compute U
      end if

      if (present(jobvt)) then
         l_jobvt = jobvt
      else if (present(VT)) then
         l_jobvt = "S"  ! Compute first min(m,n) rows of VT
      else
         l_jobvt = "N"  ! Don't compute VT
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

   subroutine pic_stpttr(AP, A, uplo, info)
      !! packed triangular -> full storage; writes only the uplo triangle
      real(sp), intent(in) :: AP(:)
      real(sp), intent(inout) :: A(:, :)
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_uplo
      integer(default_int) :: n, lda, l_info

      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if

      n = size(A, 2)
      lda = max(1, size(A, 1))

      call lapack_tpttr(l_uplo, n, AP, A, lda, l_info)

      if (present(info)) info = l_info

   end subroutine pic_stpttr

   subroutine pic_dtpttr(AP, A, uplo, info)
      !! packed triangular -> full storage; writes only the uplo triangle
      real(dp), intent(in) :: AP(:)
      real(dp), intent(inout) :: A(:, :)
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_uplo
      integer(default_int) :: n, lda, l_info

      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if

      n = size(A, 2)
      lda = max(1, size(A, 1))

      call lapack_tpttr(l_uplo, n, AP, A, lda, l_info)

      if (present(info)) info = l_info

   end subroutine pic_dtpttr

   subroutine pic_strttp(A, AP, uplo, info)
      !! full storage -> packed triangular; reads only the uplo triangle
      real(sp), intent(in) :: A(:, :)
      real(sp), intent(out) :: AP(:)
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_uplo
      integer(default_int) :: n, lda, l_info

      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if

      n = size(A, 2)
      lda = max(1, size(A, 1))

      call lapack_trttp(l_uplo, n, A, lda, AP, l_info)

      if (present(info)) info = l_info

   end subroutine pic_strttp

   subroutine pic_dtrttp(A, AP, uplo, info)
      !! full storage -> packed triangular; reads only the uplo triangle
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(out) :: AP(:)
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(out), optional :: info
      character(len=1) :: l_uplo
      integer(default_int) :: n, lda, l_info

      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if

      n = size(A, 2)
      lda = max(1, size(A, 1))

      call lapack_trttp(l_uplo, n, A, lda, AP, l_info)

      if (present(info)) info = l_info

   end subroutine pic_dtrttp

   subroutine pic_sunpack_x(n, AP, A, lda, mode, uplo, nb)
      !! expand packed triangular storage into a full mirrored matrix,
      !! explicit dimensions
      !!
      !! The same operation as pic_unpack, in the shape the BLAS uses: order
      !! and leading dimension passed, arrays assumed-size, nothing optional.
      !! That costs the caller a little clarity and saves an array descriptor
      !! at every call, which is the whole runtime at small n -- the friendly
      !! interface is about 9 ns slower per call, against a kernel that is
      !! only 13 ns at n=4. Call this one from inside a loop over small
      !! blocks; call pic_unpack everywhere else.
      integer(default_int), intent(in) :: n, lda
      real(sp), intent(in) :: AP(*)
      real(sp), intent(out) :: A(lda, *)
      character(len=1), intent(in) :: mode, uplo
      integer(default_int), intent(in) :: nb
      logical :: upper, anti
      integer(default_int) :: ii, jj, ihi, jhi, i, j, ij, base
      real(sp) :: v, sgn

      if (n <= 0) return

      upper = (uplo == "U" .or. uplo == "u")
      anti = (mode == "A" .or. mode == "a")

      !  The sign is carried as a multiply rather than by duplicating every
      !  loop. Measured: identical to two decimal places at every size and
      !  every block size, because these loops are bound by stores and not by
      !  arithmetic. Two code paths here would buy nothing.
      sgn = 1.0_sp
      if (anti) sgn = -1.0_sp

      if (n < unpack_small_max) then
         !
         !  Small. The interleaved loop, which is what the hand-written
         !  version does. Kept deliberately: at this size it already runs at
         !  the store port's limit and nothing measured beat it.
         !
         !  Four loops rather than two with a branch inside. At n of a few
         !  the branch per column and the separately written diagonal were
         !  enough to lose 30% against the hand-written loop; hoisted out,
         !  the symmetric case below is that loop exactly, including writing
         !  the diagonal twice rather than special-casing it.
         ij = 0
         if (upper .and. .not. anti) then
            do j = 1, n
               do i = 1, j - 1
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = v
               end do
               ij = ij + 1
               A(j, j) = AP(ij)
            end do
         else if (upper) then
            do j = 1, n
               do i = 1, j - 1
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = -v
               end do
               ij = ij + 1
               A(j, j) = 0.0_sp
            end do
         else if (.not. anti) then
            do j = 1, n
               ij = ij + 1
               A(j, j) = AP(ij)
               do i = j + 1, n
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = v
               end do
            end do
         else
            do j = 1, n
               ij = ij + 1
               A(j, j) = 0.0_sp
               do i = j + 1, n
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = -v
               end do
            end do
         end if

      else if (n < unpack_colwise_max) then
         !
         !  Medium. Walk the destination column rather than the source
         !  triangle. Part of each column comes from a contiguous run of AP
         !  and the rest from a strided one, but every store is contiguous
         !  either way -- which is what matters, since a scattered store
         !  cannot be combined in the write buffer the way a strided load can
         !  be prefetched.
         !
         if (upper) then
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do i = 1, n
               base = (i*(i - 1))/2
               do j = 1, i - 1
                  A(j, i) = AP(base + j)
               end do
               if (anti) then
                  A(i, i) = 0.0_sp
               else
                  A(i, i) = AP(base + i)
               end if
               do j = i + 1, n
                  A(j, i) = sgn*AP((j*(j - 1))/2 + i)
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do i = 1, n
               base = ((i - 1)*(2*n - i))/2
               do j = 1, i - 1
                  A(j, i) = sgn*AP(((j - 1)*(2*n - j))/2 + i)
               end do
               if (anti) then
                  A(i, i) = 0.0_sp
               else
                  A(i, i) = AP(base + i)
               end if
               do j = i + 1, n
                  A(j, i) = AP(base + j)
               end do
            end do
            !$omp end parallel do
         end if

      else
         !
         !  Large. Past the point where the column walk's strided load
         !  outruns cache, fill the stored triangle contiguously and then
         !  mirror it in tiles, so both streams of the mirror stay resident.
         !
         if (upper) then
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do j = 1, n
               base = (j*(j - 1))/2
               do i = 1, j
                  A(i, j) = AP(base + i)
               end do
            end do
            !$omp end parallel do
            !$omp parallel do private(ii, jhi, ihi, i, j) &
            !$omp             schedule(dynamic) if(n >= unpack_par_min)
            do jj = 1, n, nb
               jhi = min(jj + nb - 1, n)
               do ii = 1, jhi, nb
                  ihi = min(ii + nb - 1, n)
                  do j = jj, jhi
                     do i = ii, min(ihi, j - 1)
                        A(j, i) = sgn*A(i, j)
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do j = 1, n
               base = ((j - 1)*(2*n - j))/2
               do i = j, n
                  A(i, j) = AP(base + i)
               end do
            end do
            !$omp end parallel do
            !$omp parallel do private(ii, jhi, ihi, i, j) &
            !$omp             schedule(dynamic) if(n >= unpack_par_min)
            do jj = 1, n, nb
               jhi = min(jj + nb - 1, n)
               do ii = jj, n, nb
                  ihi = min(ii + nb - 1, n)
                  do j = jj, jhi
                     do i = max(ii, j + 1), ihi
                        A(j, i) = sgn*A(i, j)
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end if

         if (anti) then
            do j = 1, n
               A(j, j) = 0.0_sp
            end do
         end if
      end if

   end subroutine pic_sunpack_x

   subroutine pic_sunpack(AP, A, mode, uplo, nb)
      !! expand packed triangular storage into a full mirrored matrix
      real(sp), intent(in), contiguous :: AP(:)
      real(sp), intent(out), contiguous :: A(:, :)
      character(len=1), intent(in), optional :: mode
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(in), optional :: nb
      character(len=1) :: l_mode, l_uplo
      integer(default_int) :: l_nb

      l_mode = "S"
      if (present(mode)) l_mode = mode
      l_uplo = "U"
      if (present(uplo)) l_uplo = uplo
      l_nb = unpack_nb
      if (present(nb)) l_nb = nb

      call pic_sunpack_x(size(A, 2, kind=default_int), AP, A, &
                         size(A, 1, kind=default_int), l_mode, l_uplo, l_nb)

   end subroutine pic_sunpack

   subroutine pic_dunpack_x(n, AP, A, lda, mode, uplo, nb)
      !! expand packed triangular storage into a full mirrored matrix,
      !! explicit dimensions
      !!
      !! The same operation as pic_unpack, in the shape the BLAS uses: order
      !! and leading dimension passed, arrays assumed-size, nothing optional.
      !! That costs the caller a little clarity and saves an array descriptor
      !! at every call, which is the whole runtime at small n -- the friendly
      !! interface is about 9 ns slower per call, against a kernel that is
      !! only 13 ns at n=4. Call this one from inside a loop over small
      !! blocks; call pic_unpack everywhere else.
      integer(default_int), intent(in) :: n, lda
      real(dp), intent(in) :: AP(*)
      real(dp), intent(out) :: A(lda, *)
      character(len=1), intent(in) :: mode, uplo
      integer(default_int), intent(in) :: nb
      logical :: upper, anti
      integer(default_int) :: ii, jj, ihi, jhi, i, j, ij, base
      real(dp) :: v, sgn

      if (n <= 0) return

      upper = (uplo == "U" .or. uplo == "u")
      anti = (mode == "A" .or. mode == "a")

      !  The sign is carried as a multiply rather than by duplicating every
      !  loop. Measured: identical to two decimal places at every size and
      !  every block size, because these loops are bound by stores and not by
      !  arithmetic. Two code paths here would buy nothing.
      sgn = 1.0_dp
      if (anti) sgn = -1.0_dp

      if (n < unpack_small_max) then
         !
         !  Small. The interleaved loop, which is what the hand-written
         !  version does. Kept deliberately: at this size it already runs at
         !  the store port's limit and nothing measured beat it.
         !
         !  Four loops rather than two with a branch inside. At n of a few
         !  the branch per column and the separately written diagonal were
         !  enough to lose 30% against the hand-written loop; hoisted out,
         !  the symmetric case below is that loop exactly, including writing
         !  the diagonal twice rather than special-casing it.
         ij = 0
         if (upper .and. .not. anti) then
            do j = 1, n
               do i = 1, j - 1
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = v
               end do
               ij = ij + 1
               A(j, j) = AP(ij)
            end do
         else if (upper) then
            do j = 1, n
               do i = 1, j - 1
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = -v
               end do
               ij = ij + 1
               A(j, j) = 0.0_dp
            end do
         else if (.not. anti) then
            do j = 1, n
               ij = ij + 1
               A(j, j) = AP(ij)
               do i = j + 1, n
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = v
               end do
            end do
         else
            do j = 1, n
               ij = ij + 1
               A(j, j) = 0.0_dp
               do i = j + 1, n
                  ij = ij + 1
                  v = AP(ij)
                  A(i, j) = v
                  A(j, i) = -v
               end do
            end do
         end if

      else if (n < unpack_colwise_max) then
         !
         !  Medium. Walk the destination column rather than the source
         !  triangle. Part of each column comes from a contiguous run of AP
         !  and the rest from a strided one, but every store is contiguous
         !  either way -- which is what matters, since a scattered store
         !  cannot be combined in the write buffer the way a strided load can
         !  be prefetched.
         !
         if (upper) then
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do i = 1, n
               base = (i*(i - 1))/2
               do j = 1, i - 1
                  A(j, i) = AP(base + j)
               end do
               if (anti) then
                  A(i, i) = 0.0_dp
               else
                  A(i, i) = AP(base + i)
               end if
               do j = i + 1, n
                  A(j, i) = sgn*AP((j*(j - 1))/2 + i)
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do i = 1, n
               base = ((i - 1)*(2*n - i))/2
               do j = 1, i - 1
                  A(j, i) = sgn*AP(((j - 1)*(2*n - j))/2 + i)
               end do
               if (anti) then
                  A(i, i) = 0.0_dp
               else
                  A(i, i) = AP(base + i)
               end if
               do j = i + 1, n
                  A(j, i) = AP(base + j)
               end do
            end do
            !$omp end parallel do
         end if

      else
         !
         !  Large. Past the point where the column walk's strided load
         !  outruns cache, fill the stored triangle contiguously and then
         !  mirror it in tiles, so both streams of the mirror stay resident.
         !
         if (upper) then
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do j = 1, n
               base = (j*(j - 1))/2
               do i = 1, j
                  A(i, j) = AP(base + i)
               end do
            end do
            !$omp end parallel do
            !$omp parallel do private(ii, jhi, ihi, i, j) &
            !$omp             schedule(dynamic) if(n >= unpack_par_min)
            do jj = 1, n, nb
               jhi = min(jj + nb - 1, n)
               do ii = 1, jhi, nb
                  ihi = min(ii + nb - 1, n)
                  do j = jj, jhi
                     do i = ii, min(ihi, j - 1)
                        A(j, i) = sgn*A(i, j)
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do private(i, j, base) schedule(static) &
            !$omp             if(n >= unpack_par_min)
            do j = 1, n
               base = ((j - 1)*(2*n - j))/2
               do i = j, n
                  A(i, j) = AP(base + i)
               end do
            end do
            !$omp end parallel do
            !$omp parallel do private(ii, jhi, ihi, i, j) &
            !$omp             schedule(dynamic) if(n >= unpack_par_min)
            do jj = 1, n, nb
               jhi = min(jj + nb - 1, n)
               do ii = jj, n, nb
                  ihi = min(ii + nb - 1, n)
                  do j = jj, jhi
                     do i = max(ii, j + 1), ihi
                        A(j, i) = sgn*A(i, j)
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end if

         if (anti) then
            do j = 1, n
               A(j, j) = 0.0_dp
            end do
         end if
      end if

   end subroutine pic_dunpack_x

   subroutine pic_dunpack(AP, A, mode, uplo, nb)
      !! expand packed triangular storage into a full mirrored matrix
      real(dp), intent(in), contiguous :: AP(:)
      real(dp), intent(out), contiguous :: A(:, :)
      character(len=1), intent(in), optional :: mode
      character(len=1), intent(in), optional :: uplo
      integer(default_int), intent(in), optional :: nb
      character(len=1) :: l_mode, l_uplo
      integer(default_int) :: l_nb

      l_mode = "S"
      if (present(mode)) l_mode = mode
      l_uplo = "U"
      if (present(uplo)) l_uplo = uplo
      l_nb = unpack_nb
      if (present(nb)) l_nb = nb

      call pic_dunpack_x(size(A, 2, kind=default_int), AP, A, &
                         size(A, 1, kind=default_int), l_mode, l_uplo, l_nb)

   end subroutine pic_dunpack

end module pic_lapack_interfaces
