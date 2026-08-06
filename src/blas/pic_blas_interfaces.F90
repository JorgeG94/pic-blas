! SPDX-License-Identifier: MIT
! Copyright (c) 2025 Jorge Luis Galvez Vallejo
!! this file contains the interfaces for the BLAS routines of all levels
!! I might consider splitting them up later but alas, I don't have the time now
!! the idea of this file is to provide something akin to
!! interface blas_gemm
!!  subroutine sgemm()
!!  subroutine dgemm()
!!   ... etc,
!! end interface blas_gemm
!! so that I can use the same interface for all BLAS routines

module pic_blas_interfaces
  !! pic_blas_interfaces.F90 provides the interfaces for the BLAS routines
  !! the idea is to have a two level interface, first pic_blas_xyz which
  !! is the way programmers will use BLAS, it'll do some checks and then
  !! call the "overloaded" BLAS interfaces to call the correct BLAS routine
   use pic_types, only: sp, dp, default_int
   implicit none
   private

   ! these are the cool overloaded interfaces, the pic_xyz function
   ! has the procedures pic_(type)xyz which will call the correct BLAS routine
   ! depending on the data type of the arguments
   ! this _needs_ allocatable arrays since we deduce shapes from the arrays themselves
   public :: pic_gemm, pic_gemv, pic_asum, pic_axpy, pic_copy, pic_dot, pic_scal, pic_iamax
   public :: pic_symm, pic_syrk, pic_syr2k
   public :: pic_gemm_x, pic_gemv_x, pic_ger, pic_ger_x
   public :: pic_sgemv_x, pic_dgemv_x, pic_sger_x, pic_dger_x
   public :: pic_sgemm_x, pic_dgemm_x

   ! tested
   interface pic_gemm
      !! general interface of the BLAS GEMM routines, will call SGEMM, DGEMM, CGEMM
      !!
      !! Usage: call pic_gemm(A, B, C, [optional] transa, [optional] transb, [optional] alpha, [optional] beta)
      !!
      !! where A, B, C are matrices, transa and transb are optional transpose options,
      !! alpha and beta are optional scaling factors
      !!
      !! By default, if not specified transA and transB are "N" (no transpose),
      !! and alpha and beta are 1.0 and 0.0 respectively.
      !!
      !! The matrices A, B, C must be allocatable arrays, we deduce the shapes from them.
      module procedure :: pic_sgemm
      module procedure :: pic_dgemm
   end interface pic_gemm

   interface pic_gemv
      !! general interface of the BLAS GEMV routines, will call SGEMV, DGEMV, CGEMV, ZGEMV
      !!
      !! Usage: call pic_gemv(A, x, y, [optional] transa, [optional] alpha, [optional] beta)
      !!
      !! where A is a matrix, x and y are vectors, transa is an optional transpose option,
      !! alpha and beta are optional scaling factors.
      !!
      !! The matrix A must be an allocatable array, we deduce the shapes from it.
      !! TransA is "N" (no transpose) by default. And alpha and beta are 1.0 and 0.0 respectively.
      !!
      module procedure :: pic_sgemv
      module procedure :: pic_dgemv
   end interface pic_gemv
   ! tested
   interface pic_asum
      !! general interface of the BLAS ASUM routines, will call SASUM, DASUM, SCASUM, DZASUM
      !!
      !! Usage: result = pic_asum(x, incx)
      !!
      !! where x is a vector and incx is the increment, this will return the sum of the absolute values
      !! of the elements of x.
      !!
      !! The vector x must be an allocatable array, we deduce the shape from it.
      !! The increment incx is 1 by default.
      !!
      module procedure :: pic_sasum
      module procedure :: pic_dasum
      module procedure :: pic_scasum
      module procedure :: pic_dzasum
   end interface pic_asum

   interface pic_axpy
      !! general interface of the BLAS AXPY routines, will call SAXPY, DAXPY, CAXPY, ZAXPY
      !!
      !! Usage: call pic_axpy(n, alpha, x, incx, y, incy)
      !!
      !! where n is the number of elements, alpha is the scaling factor,
      !! x is the input vector, incx is the increment for x, y is the output vector,
      !! and incy is the increment for y.
      !!
      !! The vectors x and y must be allocatable arrays, we deduce the shapes from them.
      !! The increments incx and incy are 1 by default.
      !!
      module procedure :: pic_saxpy
      module procedure :: pic_daxpy
   end interface pic_axpy

   interface pic_copy
      !! general interface of the BLAS COPY routines, will call SCOPY, DCOPY, CCOPY, ZCOPY
      !!
      !! Usage: call pic_copy(x, y)
      !!
      !! where x is the input vector, y is the output vector.
      !! The vectors x and y must be allocatable arrays, we deduce the shapes from them.
      !!
      module procedure :: pic_scopy
      module procedure :: pic_dcopy
   end interface pic_copy

   interface pic_dot
      !! general interface of the BLAS DOT routines, will call SDOT, DDOT, CDOTC, ZDOTC
      !!
      !! Usage: result = pic_dot(x, y)
      !!
      !! where x is the input vector, y is the output vector.
      !! The vectors x and y must be allocatable arrays, we deduce the shapes from them.
      !!
      module procedure :: pic_sdot
      module procedure :: pic_ddot
   end interface pic_dot

   interface pic_scal
      !! general interface of the BLAS SCAL routines, will call SSCAL, DSCAL, CSCAL, ZSCAL
      !!
      !! Usage: call pic_scal(x, [optional] alpha)
      !!
      !! where x is the input vector, alpha is the scaling factor.
      !! The vector x must be an allocatable array, we deduce the shape from it.
      !! The scaling factor alpha is 1.0 by default.
      !!
      module procedure :: pic_sscal
      module procedure :: pic_dscal
   end interface pic_scal

   interface pic_iamax
      !! general interface of the BLAS IAMAX routines, will call ISAMAX, IDAMAX, ICAMAX, IZAMAX
      !!
      !! Usage: idx = pic_iamax(x, incx)
      !!
      !! where x is the input vector, incx is the increment.
      !! The vector x must be an allocatable array, we deduce the shape from it.
      !! The increment incx is 1 by default.
      !!
      module procedure :: pic_isamax
      module procedure :: pic_idamax
   end interface pic_iamax

   interface pic_symm
      !! general interface of the BLAS SYMM routines, will call SSYMM, DSYMM
      !!
      !! Usage: call pic_symm(A, B, C, [side], [uplo], [alpha], [beta])
      !!
      !! C := alpha*A*B + beta*C with A symmetric (side "L", the default), or
      !! C := alpha*B*A + beta*C (side "R"). Only the triangle of A named by
      !! uplo ("U" by default) is read.
      !!
      !! This is the one to reach for when a symmetric matrix has been
      !! expanded to a full square and would otherwise go through GEMM: it
      !! does half the work for the same result.
      module procedure :: pic_ssymm
      module procedure :: pic_dsymm
   end interface pic_symm

   interface pic_syrk
      !! general interface of the BLAS SYRK routines, will call SSYRK, DSYRK
      !!
      !! Usage: call pic_syrk(A, C, [uplo], [trans], [alpha], [beta])
      !!
      !! C := alpha*A*A**T + beta*C  (trans "N", the default)
      !! C := alpha*A**T*A + beta*C  (trans "T")
      !!
      !! C is symmetric and only the uplo triangle is written. Half the flops
      !! of the equivalent GEMM, and the result is symmetric by construction
      !! rather than by assumption.
      module procedure :: pic_ssyrk
      module procedure :: pic_dsyrk
   end interface pic_syrk

   interface pic_syr2k
      !! general interface of the BLAS SYR2K routines, will call SSYR2K, DSYR2K
      !!
      !! Usage: call pic_syr2k(A, B, C, [uplo], [trans], [alpha], [beta])
      !!
      !! C := alpha*A*B**T + alpha*B*A**T + beta*C  (trans "N", the default)
      !! C := alpha*A**T*B + alpha*B**T*A + beta*C  (trans "T")
      !!
      !! Worth knowing: when A**T*B is itself symmetric, B**T*A is the same
      !! matrix, so trans="T" with alpha=0.5 yields exactly the triangle of
      !! A**T*B for half a GEMM's work. That is the shape a congruence
      !! transformation C**T*Y*C takes once Y*C has been formed.
      module procedure :: pic_ssyr2k
      module procedure :: pic_dsyr2k
   end interface pic_syr2k

   interface blas_asum
      !! this is the interface for the BLAS ASUM routines, it will call SASUM, DASUM, SCASUM, DZASUM
      !! Usage: result = blas_asum(x, incx)
      !! where x is the input vector, incx is the increment.
      !!
      !! This is not a public interface, it is used internally by pic_asum
      pure function sasum(n, x, incx) result(res_sasum)
         import :: sp, default_int
         implicit none
         real(sp) :: res_sasum
         real(sp), intent(in) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end function sasum
      pure function dasum(n, x, incx) result(res_dasum)
         import :: dp, default_int
         implicit none
         real(dp) :: res_dasum
         real(dp), intent(in) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end function dasum
      pure function scasum(n, x, incx) result(res_scasum)
         import :: sp, default_int
         implicit none
         real(sp) :: res_scasum
         complex(sp), intent(in) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end function scasum
      pure function dzasum(n, x, incx) result(res_dzasum)
         import :: dp, default_int
         implicit none
         real(dp) :: res_dzasum
         complex(dp), intent(in) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end function dzasum
   end interface blas_asum

   interface blas_axpy
      !! explicit interface for BLAS AXPY routines
      !!
      !! Usage: call blas_axpy(n, alpha, x, incx, y, incy)
      !!
      !! This is not a public interface, it is used internally by pic_axpy
      pure subroutine saxpy(n, alpha, x, incx, y, incy)
         import :: sp, default_int
         implicit none
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end subroutine saxpy
      pure subroutine daxpy(n, alpha, x, incx, y, incy)
         import :: dp, default_int
         implicit none
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end subroutine daxpy
      pure subroutine caxpy(n, alpha, x, incx, y, incy)
         import :: sp, default_int
         implicit none
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end subroutine caxpy
      pure subroutine zaxpy(n, alpha, x, incx, y, incy)
         import :: dp, default_int
         implicit none
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end subroutine zaxpy
   end interface blas_axpy

   interface blas_copy
      !! explicit interface for BLAS COPY routines
      !!
      !! Usage: call blas_copy(x, y)
      !!
      !! This is not a public interface, it is used internally by pic_copy
      pure subroutine scopy(n, x, incx, y, incy)
         import :: sp, default_int
         implicit none
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end subroutine scopy
      pure subroutine dcopy(n, x, incx, y, incy)
         import :: dp, default_int
         implicit none
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end subroutine dcopy
   end interface blas_copy

   interface blas_dot
      !! explicit interface for BLAS DOT routines
      !!
      !! Usage: result = blas_dot(x, y, incx, incy, n)
      !! This is not a public interface, it is used internally by pic_dot
      !!
      pure function sdot(n, x, incx, y, incy) result(res)
         import :: sp, default_int
         implicit none
         real(sp) :: res
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end function sdot
      pure function ddot(n, x, incx, y, incy) result(res)
         import :: dp, default_int
         implicit none
         real(dp) :: res
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end function ddot
      pure function cdotc(n, x, incx, y, incy) result(res)
         import :: sp, default_int
         implicit none
         complex(sp) :: res
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end function cdotc
      pure function zdotc(n, x, incx, y, incy) result(res)
         import :: dp, default_int
         implicit none
         complex(dp) :: res
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         integer(default_int), intent(in) :: n
      end function zdotc
   end interface blas_dot

   interface blas_scal
      !! explicit interface for BLAS SCAL routines
      !!
      !! Usage: call blas_scal(n, alpha, x, incx)
      !!
      !! This is not a public interface, it is used internally by pic_scal
      pure subroutine sscal(n, alpha, x, incx)
         import :: sp, default_int
         implicit none
         real(sp), intent(in) :: alpha
         real(sp), intent(inout) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end subroutine sscal
      pure subroutine dscal(n, alpha, x, incx)
         import :: dp, default_int
         implicit none
         real(dp), intent(in) :: alpha
         real(dp), intent(inout) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end subroutine dscal
   end interface blas_scal

   interface blas_iamax
      !! explicit interface for BLAS IAMAX routines
      !!
      !! Usage: idx = blas_iamax(x, incx)
      !!
      !! This is not a public interface, it is used internally by pic_iamax
      pure function isamax(n, x, incx) result(idx)
         import :: sp, default_int
         implicit none
         integer(default_int) :: idx
         real(sp), intent(in) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end function isamax
      pure function idamax(n, x, incx) result(idx)
         import :: dp, default_int
         implicit none
         integer(default_int) :: idx
         real(dp), intent(in) :: x(*)
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: n
      end function idamax
   end interface blas_iamax

   interface blas_gemv
      !! explicit interface for BLAS GEMV routines
      !!
      !! Usage: call blas_gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      !!
      !! This is not a public interface, it is used internally by pic_gemv
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: trans
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: trans
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
      end subroutine dgemv
      pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: trans
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
      end subroutine cgemv
      pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: incx
         integer(default_int), intent(in) :: incy
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: trans
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
      end subroutine zgemv
   end interface blas_gemv

   interface pic_gemm_x
      !! GEMM with explicit dimensions, in the shape the BLAS itself uses
      !!
      !! Usage: call pic_gemm_x(transa, transb, m, n, k, alpha, A, lda, &
      !!                        B, ldb, beta, C, ldc)
      !!
      !! pic_gemm derives m, n, k and the leading dimensions from the shapes
      !! of its arguments, which is the pleasant way to call it and is right
      !! whenever the arrays are exactly the matrices being multiplied. It is
      !! the wrong way when they are not.
      !!
      !! A caller holding a block inside a larger array has to pass a section.
      !! If the leading dimension exceeds the number of rows that section is
      !! strided, and handing it to an explicit-shape dummy means the compiler
      !! packs a temporary -- for C, which is read and written, on the way in
      !! and again on the way out. Measured against calling DGEMM directly,
      !! with the leading dimension larger than the block, that costs:
      !!
      !!     16x16x16    1.59x        256x256x256   1.72x
      !!     32x32x32    1.40x        512x512x512   1.48x
      !!     64x64x64    1.20x        100x20x400    1.35x
      !!
      !! It does not wash out with size, because the packing grows with the
      !! matrices too. This entry forwards straight to the BLAS and adds
      !! nothing. Reach for it when the arrays are bigger than the multiply;
      !! reach for pic_gemm otherwise.
      module procedure :: pic_sgemm_x
      module procedure :: pic_dgemm_x
   end interface pic_gemm_x

   interface pic_gemv_x
      !! GEMV with explicit dimensions, in the shape the BLAS itself uses
      !!
      !! Usage: call pic_gemv_x(trans, m, n, alpha, A, lda, x, incx, &
      !!                        beta, y, incy)
      !!
      !! Same reasoning as pic_gemm_x: pic_gemv takes the matrix as
      !! assumed-shape and so needs a section when the caller holds a block
      !! inside a larger array, and a section whose leading dimension exceeds
      !! its row count is packed into a temporary before the call. Legacy
      !! callers pass a leading dimension separately as a matter of course --
      !! in GAMESS's mthlib, all seven GEMV calls do.
      !!
      !! This also takes incx and incy, so a row of a matrix can be passed
      !! directly rather than through a strided section.
      !!
      !! Calling convention, and it is a trap: legacy code names a submatrix by
      !! its first element -- A(1,j) meaning "from column j on, with leading
      !! dimension lda". Sequence association allows that, but only for a
      !! specific procedure. Generic resolution matches on rank and rejects a
      !! scalar actual argument for an array dummy, so a call through this
      !! generic fails to compile with "no specific subroutine". Call
      !! pic_dgemv_x or pic_sgemv_x by name from that style of caller; the
      !! generic is for whole arrays and rank-2 sections. Both specific names
      !! are public for exactly this reason.
      module procedure :: pic_sgemv_x
      module procedure :: pic_dgemv_x
   end interface pic_gemv_x

   interface pic_ger
      !! rank-one update: A := alpha*x*y**T + A
      !!
      !! Usage: call pic_ger(A, x, y, [alpha])
      !!
      !! The natural partner to GEMV, and what turns a Gram-Schmidt written
      !! as a loop of AXPYs into two calls over the whole trailing block.
      module procedure :: pic_sger
      module procedure :: pic_dger
   end interface pic_ger

   interface pic_ger_x
      !! rank-one update with explicit dimensions
      !!
      !! Usage: call pic_ger_x(m, n, alpha, x, incx, y, incy, A, lda)
      module procedure :: pic_sger_x
      module procedure :: pic_dger_x
   end interface pic_ger_x

   interface blas_ger
      !! not a public interface, used internally by pic_ger
      pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
         import :: sp, default_int
         implicit none
         integer(default_int), intent(in) :: m, n, incx, incy, lda
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         real(sp), intent(inout) :: a(lda, *)
      end subroutine sger
      pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
         import :: dp, default_int
         implicit none
         integer(default_int), intent(in) :: m, n, incx, incy, lda
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         real(dp), intent(inout) :: a(lda, *)
      end subroutine dger
   end interface blas_ger

   interface blas_gemm
      !! explicit interface for BLAS GEMM routines
      !!
      !! Usage: call blas_gemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      !!
      !! This is not a public interface, it is used internally by pic_gemm
      pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
           & beta, c, ldc)
         import :: sp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
      end subroutine sgemm
      pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
           & beta, c, ldc)
         import :: dp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
      end subroutine dgemm
      pure subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
           & beta, c, ldc)
         import :: sp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
      end subroutine cgemm
      pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
           & beta, c, ldc)
         import :: dp, default_int
         implicit none
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
      end subroutine zgemm
   end interface blas_gemm

   interface blas_symm
      !! not a public interface, used internally by pic_symm
      pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssymm
      pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         integer(default_int), intent(in) :: m
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsymm
   end interface blas_symm

   interface blas_syrk
      !! not a public interface, used internally by pic_syrk
      pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldc
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssyrk
      pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldc
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsyrk
   end interface blas_syrk

   interface blas_syr2k
      !! not a public interface, used internally by pic_syr2k
      pure subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: sp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssyr2k
      pure subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: dp, default_int
         implicit none
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         integer(default_int), intent(in) :: n
         integer(default_int), intent(in) :: k
         integer(default_int), intent(in) :: lda
         integer(default_int), intent(in) :: ldb
         integer(default_int), intent(in) :: ldc
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsyr2k
   end interface blas_syr2k

contains

   pure subroutine pic_sgemm(A, B, C, transa, transb, alpha, beta)
      !! interface for single precision matrix multiplication
      real(sp), intent(in) :: A(:, :)
      real(sp), intent(in) :: B(:, :)
      real(sp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: transa
      character(len=1), intent(in), optional :: transb
      real(sp), intent(in), optional :: alpha
      real(sp), intent(in), optional :: beta
      character(len=1) :: OP_A, OP_B
      real(sp) :: l_alpha, l_beta
      integer(default_int) :: m, n, k, lda, ldb, ldc

      ! first check for the constants
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      ! check the OP options, maybe this should not be optional
      if (present(transa)) then
         OP_A = transa
      else
         OP_A = "N"
      end if
      if (present(transb)) then
         OP_B = transb
      else
         OP_B = "N"
      end if

      ! check for the dimensions now
      if ((OP_A == "N" .or. OP_A == "n")) then
         k = size(A, 2)
      else
         k = size(A, 1)
      end if

      ! get LDA, LDB, and LDC
      lda = max(1, size(A, 1))
      ldb = max(1, size(B, 1))
      ldc = max(1, size(C, 1))
      m = size(C, 1)
      n = size(C, 2)

      call blas_gemm(OP_A, OP_B, m, n, k, l_alpha, A, lda, B, ldb, l_beta, C, ldc)

   end subroutine pic_sgemm

   pure subroutine pic_dgemm(A, B, C, transa, transb, alpha, beta)
      !! interface for single precision matrix multiplication
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(in) :: B(:, :)
      real(dp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: transa
      character(len=1), intent(in), optional :: transb
      real(dp), intent(in), optional :: alpha
      real(dp), intent(in), optional :: beta
      character(len=1) :: OP_A, OP_B
      real(dp) :: l_alpha, l_beta
      integer(default_int) :: m, n, k, lda, ldb, ldc

      ! first check for the constants
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      ! check the OP options, maybe this should not be optional
      if (present(transa)) then
         OP_A = transa
      else
         OP_A = "N"
      end if
      if (present(transb)) then
         OP_B = transb
      else
         OP_B = "N"
      end if

      ! check for the dimensions now
      if ((OP_A == "N" .or. OP_A == "n")) then
         k = size(A, 2)
      else
         k = size(A, 1)
      end if

      ! get LDA, LDB, and LDC
      lda = max(1, size(A, 1))
      ldb = max(1, size(B, 1))
      ldc = max(1, size(C, 1))
      m = size(C, 1)
      n = size(C, 2)

      call blas_gemm(OP_A, OP_B, m, n, k, l_alpha, A, lda, B, ldb, l_beta, C, ldc)

   end subroutine pic_dgemm

   pure subroutine pic_ssymm(A, B, C, side, uplo, alpha, beta)
      !! symmetric matrix times general matrix
      real(sp), intent(in) :: A(:, :)
      real(sp), intent(in) :: B(:, :)
      real(sp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: side
      character(len=1), intent(in), optional :: uplo
      real(sp), intent(in), optional :: alpha
      real(sp), intent(in), optional :: beta
      character(len=1) :: l_side, l_uplo
      real(sp) :: l_alpha, l_beta
      integer(default_int) :: m, n, lda, ldb, ldc

      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      if (present(side)) then
         l_side = side
      else
         l_side = "L"
      end if
      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if

      ! the shape of C fixes the problem: A is m x m for side "L" and n x n
      ! for side "R", which the caller has already had to get right for the
      ! shapes to be conformable at all
      m = size(C, 1)
      n = size(C, 2)
      lda = max(1, size(A, 1))
      ldb = max(1, size(B, 1))
      ldc = max(1, size(C, 1))

      call blas_symm(l_side, l_uplo, m, n, l_alpha, A, lda, B, ldb, l_beta, C, ldc)

   end subroutine pic_ssymm

   pure subroutine pic_dsymm(A, B, C, side, uplo, alpha, beta)
      !! symmetric matrix times general matrix
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(in) :: B(:, :)
      real(dp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: side
      character(len=1), intent(in), optional :: uplo
      real(dp), intent(in), optional :: alpha
      real(dp), intent(in), optional :: beta
      character(len=1) :: l_side, l_uplo
      real(dp) :: l_alpha, l_beta
      integer(default_int) :: m, n, lda, ldb, ldc

      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_dp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_dp
      end if
      if (present(side)) then
         l_side = side
      else
         l_side = "L"
      end if
      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if

      ! the shape of C fixes the problem: A is m x m for side "L" and n x n
      ! for side "R", which the caller has already had to get right for the
      ! shapes to be conformable at all
      m = size(C, 1)
      n = size(C, 2)
      lda = max(1, size(A, 1))
      ldb = max(1, size(B, 1))
      ldc = max(1, size(C, 1))

      call blas_symm(l_side, l_uplo, m, n, l_alpha, A, lda, B, ldb, l_beta, C, ldc)

   end subroutine pic_dsymm

   pure subroutine pic_ssyrk(A, C, uplo, trans, alpha, beta)
      !! symmetric rank-k update
      real(sp), intent(in) :: A(:, :)
      real(sp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: uplo
      character(len=1), intent(in), optional :: trans
      real(sp), intent(in), optional :: alpha
      real(sp), intent(in), optional :: beta
      character(len=1) :: l_uplo, l_trans
      real(sp) :: l_alpha, l_beta
      integer(default_int) :: n, k, lda, ldc

      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if
      if (present(trans)) then
         l_trans = trans
      else
         l_trans = "N"
      end if

      n = size(C, 1)
      if (l_trans == "N" .or. l_trans == "n") then
         k = size(A, 2)
      else
         k = size(A, 1)
      end if
      lda = max(1, size(A, 1))
      ldc = max(1, size(C, 1))

      call blas_syrk(l_uplo, l_trans, n, k, l_alpha, A, lda, l_beta, C, ldc)

   end subroutine pic_ssyrk

   pure subroutine pic_dsyrk(A, C, uplo, trans, alpha, beta)
      !! symmetric rank-k update
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: uplo
      character(len=1), intent(in), optional :: trans
      real(dp), intent(in), optional :: alpha
      real(dp), intent(in), optional :: beta
      character(len=1) :: l_uplo, l_trans
      real(dp) :: l_alpha, l_beta
      integer(default_int) :: n, k, lda, ldc

      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_dp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_dp
      end if
      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if
      if (present(trans)) then
         l_trans = trans
      else
         l_trans = "N"
      end if

      n = size(C, 1)
      if (l_trans == "N" .or. l_trans == "n") then
         k = size(A, 2)
      else
         k = size(A, 1)
      end if
      lda = max(1, size(A, 1))
      ldc = max(1, size(C, 1))

      call blas_syrk(l_uplo, l_trans, n, k, l_alpha, A, lda, l_beta, C, ldc)

   end subroutine pic_dsyrk

   pure subroutine pic_ssyr2k(A, B, C, uplo, trans, alpha, beta)
      !! symmetric rank-2k update
      real(sp), intent(in) :: A(:, :)
      real(sp), intent(in) :: B(:, :)
      real(sp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: uplo
      character(len=1), intent(in), optional :: trans
      real(sp), intent(in), optional :: alpha
      real(sp), intent(in), optional :: beta
      character(len=1) :: l_uplo, l_trans
      real(sp) :: l_alpha, l_beta
      integer(default_int) :: n, k, lda, ldb, ldc

      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if
      if (present(trans)) then
         l_trans = trans
      else
         l_trans = "N"
      end if

      n = size(C, 1)
      if (l_trans == "N" .or. l_trans == "n") then
         k = size(A, 2)
      else
         k = size(A, 1)
      end if
      lda = max(1, size(A, 1))
      ldb = max(1, size(B, 1))
      ldc = max(1, size(C, 1))

      call blas_syr2k(l_uplo, l_trans, n, k, l_alpha, A, lda, B, ldb, l_beta, C, ldc)

   end subroutine pic_ssyr2k

   pure subroutine pic_dsyr2k(A, B, C, uplo, trans, alpha, beta)
      !! symmetric rank-2k update
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(in) :: B(:, :)
      real(dp), intent(inout) :: C(:, :)
      character(len=1), intent(in), optional :: uplo
      character(len=1), intent(in), optional :: trans
      real(dp), intent(in), optional :: alpha
      real(dp), intent(in), optional :: beta
      character(len=1) :: l_uplo, l_trans
      real(dp) :: l_alpha, l_beta
      integer(default_int) :: n, k, lda, ldb, ldc

      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_dp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_dp
      end if
      if (present(uplo)) then
         l_uplo = uplo
      else
         l_uplo = "U"
      end if
      if (present(trans)) then
         l_trans = trans
      else
         l_trans = "N"
      end if

      n = size(C, 1)
      if (l_trans == "N" .or. l_trans == "n") then
         k = size(A, 2)
      else
         k = size(A, 1)
      end if
      lda = max(1, size(A, 1))
      ldb = max(1, size(B, 1))
      ldc = max(1, size(C, 1))

      call blas_syr2k(l_uplo, l_trans, n, k, l_alpha, A, lda, B, ldb, l_beta, C, ldc)

   end subroutine pic_dsyr2k

   pure subroutine pic_sgemv(A, x, y, trans_a, alpha, beta)
      !! interface for single precision matrix-vector multiplication
      real(sp), intent(in) :: A(:, :)
      real(sp), intent(in) :: x(:)
      real(sp), intent(inout) :: y(:)
      character(len=1), intent(in), optional :: trans_a
      real(sp), intent(in), optional :: alpha
      real(sp), intent(in), optional :: beta
      real(sp) :: l_alpha, l_beta
      character(len=1) :: l_trans_a
      integer(default_int) :: incx, incy, m, n, lda
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      if (present(trans_a)) then
         l_trans_a = trans_a
      else
         l_trans_a = "n"
      end if
      incx = 1
      incy = 1
      lda = max(1, size(A, 1))
      m = size(A, 1)
      n = size(A, 2)
      call blas_gemv(l_trans_a, m, n, l_alpha, A, lda, x, incx, l_beta, y, incy)
   end subroutine pic_sgemv

   pure subroutine pic_dgemv(A, x, y, trans_a, alpha, beta)
      !! interface for double precision matrix-vector multiplication
      real(dp), intent(in) :: A(:, :)
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      character(len=1), intent(in), optional :: trans_a
      real(dp), intent(in), optional :: alpha
      real(dp), intent(in), optional :: beta
      real(dp) :: l_alpha, l_beta
      character(len=1) :: l_trans_a
      integer(default_int) :: incx, incy, m, n, lda
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      if (present(beta)) then
         l_beta = beta
      else
         l_beta = 0.0_sp
      end if
      if (present(trans_a)) then
         l_trans_a = trans_a
      else
         l_trans_a = "n"
      end if
      incx = 1
      incy = 1
      lda = max(1, size(A, 1))
      m = size(A, 1)
      n = size(A, 2)
      call blas_gemv(l_trans_a, m, n, l_alpha, A, lda, x, incx, l_beta, y, incy)
   end subroutine pic_dgemv

   function pic_sasum(x) result(res)
      !! interface for single precision absolute sum
      real(sp), intent(in) :: x(:)
      real(sp) :: res
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      res = blas_asum(n, x, incx)
   end function pic_sasum

   function pic_dasum(x) result(res)
      !! interface for double precision absolute sum
      real(dp), intent(in) :: x(:)
      real(dp) :: res
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      res = blas_asum(n, x, incx)
   end function pic_dasum

   function pic_scasum(x) result(res)
      !! interface for single precision complex absolute sum
      complex(sp), intent(in) :: x(:)
      real(sp) :: res
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      res = blas_asum(n, x, incx)
   end function pic_scasum

   function pic_dzasum(x) result(res)
      !! interface for double precision complex absolute sum
      complex(dp), intent(in) :: x(:)
      real(dp) :: res
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      res = blas_asum(n, x, incx)
   end function pic_dzasum

   subroutine pic_saxpy(x, y, alpha)
      !! interface for single precision AXPY
      real(sp), intent(in) :: x(:)
      real(sp), intent(inout) :: y(:)
      real(sp), intent(in), optional :: alpha
      real(sp) :: l_alpha
      integer(default_int) :: n, incx, incy
      n = size(x)
      incx = 1
      incy = 1
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      call blas_axpy(n, l_alpha, x, incx, y, incy)
   end subroutine pic_saxpy

   subroutine pic_daxpy(x, y, alpha)
      !! interface for double precision AXPY
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      real(dp), intent(in), optional :: alpha
      real(dp) :: l_alpha
      integer(default_int) :: n, incx, incy
      n = size(x)
      incx = 1
      incy = 1
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_dp
      end if
      call blas_axpy(n, l_alpha, x, incx, y, incy)
   end subroutine pic_daxpy

   subroutine pic_scopy(x, y)
      !! interface for single precision copy
      real(sp), intent(in) :: x(:)
      real(sp), intent(inout) :: y(:)
      integer(default_int) :: n, incx, incy
      n = size(x)
      incx = 1
      incy = 1
      call blas_copy(n, x, incx, y, incy)
   end subroutine pic_scopy

   subroutine pic_dcopy(x, y)
      !! interface for double precision copy
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      integer(default_int) :: n, incx, incy
      n = size(x)
      incx = 1
      incy = 1
      call blas_copy(n, x, incx, y, incy)
   end subroutine pic_dcopy

   function pic_sdot(x, y) result(res)
      !! interface for single precision dot product
      real(sp), intent(in) :: x(:)
      real(sp), intent(in) :: y(:)
      real(sp) :: res
      integer(default_int) :: n, incx, incy
      n = size(x)
      incx = 1
      incy = 1
      res = blas_dot(n, x, incx, y, incy)
   end function pic_sdot

   function pic_ddot(x, y) result(res)
      !! interface for double precision dot product
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp) :: res
      integer(default_int) :: n, incx, incy
      n = size(x)
      incx = 1
      incy = 1
      res = blas_dot(n, x, incx, y, incy)
   end function pic_ddot

   subroutine pic_sscal(x, alpha)
      !! interface for single precision scaling
      real(sp), intent(inout) :: x(:)
      real(sp), intent(in), optional :: alpha
      real(sp) :: l_alpha
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_sp
      end if
      call blas_scal(n, l_alpha, x, incx)
   end subroutine pic_sscal

   subroutine pic_dscal(x, alpha)
      !! interface for double precision scaling
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in), optional :: alpha
      real(dp) :: l_alpha
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      if (present(alpha)) then
         l_alpha = alpha
      else
         l_alpha = 1.0_dp
      end if
      call blas_scal(n, l_alpha, x, incx)
   end subroutine pic_dscal

   function pic_isamax(x) result(idx)
      !! interface for single precision index of maximum absolute value
      real(sp), intent(in) :: x(:)
      integer(default_int) :: idx
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      idx = blas_iamax(n, x, incx)
   end function pic_isamax
   function pic_idamax(x) result(idx)
      !! interface for double precision index of maximum absolute value
      real(dp), intent(in) :: x(:)
      integer(default_int) :: idx
      integer(default_int) :: n, incx
      n = size(x)
      incx = 1
      idx = blas_iamax(n, x, incx)
   end function pic_idamax
   pure subroutine pic_sgemm_x(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      !! GEMM with explicit dimensions; forwards to the BLAS unchanged
      character(len=1), intent(in) :: transa, transb
      integer(default_int), intent(in) :: m, n, k, lda, ldb, ldc
      real(sp), intent(in) :: alpha, beta
      real(sp), intent(in) :: A(lda, *)
      real(sp), intent(in) :: B(ldb, *)
      real(sp), intent(inout) :: C(ldc, *)

      call blas_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
   end subroutine pic_sgemm_x

   pure subroutine pic_dgemm_x(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
      !! GEMM with explicit dimensions; forwards to the BLAS unchanged
      character(len=1), intent(in) :: transa, transb
      integer(default_int), intent(in) :: m, n, k, lda, ldb, ldc
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: A(lda, *)
      real(dp), intent(in) :: B(ldb, *)
      real(dp), intent(inout) :: C(ldc, *)

      call blas_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
   end subroutine pic_dgemm_x

   pure subroutine pic_sgemv_x(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
      !! GEMV with explicit dimensions; forwards to the BLAS unchanged
      character(len=1), intent(in) :: trans
      integer(default_int), intent(in) :: m, n, lda, incx, incy
      real(sp), intent(in) :: alpha, beta
      real(sp), intent(in) :: A(lda, *)
      real(sp), intent(in) :: x(*)
      real(sp), intent(inout) :: y(*)

      call blas_gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
   end subroutine pic_sgemv_x

   pure subroutine pic_sger_x(m, n, alpha, x, incx, y, incy, A, lda)
      !! rank-one update with explicit dimensions; forwards unchanged
      integer(default_int), intent(in) :: m, n, incx, incy, lda
      real(sp), intent(in) :: alpha
      real(sp), intent(in) :: x(*)
      real(sp), intent(in) :: y(*)
      real(sp), intent(inout) :: A(lda, *)

      call blas_ger(m, n, alpha, x, incx, y, incy, A, lda)
   end subroutine pic_sger_x

   pure subroutine pic_sger(A, x, y, alpha)
      !! rank-one update, dimensions taken from the arguments
      real(sp), intent(inout) :: A(:, :)
      real(sp), intent(in) :: x(:)
      real(sp), intent(in) :: y(:)
      real(sp), intent(in), optional :: alpha
      real(sp) :: l_alpha
      integer(default_int) :: m, n, lda

      l_alpha = 1.0_sp
      if (present(alpha)) l_alpha = alpha
      m = size(A, 1, kind=default_int)
      n = size(A, 2, kind=default_int)
      lda = max(1_default_int, m)
      call blas_ger(m, n, l_alpha, x, 1_default_int, y, 1_default_int, A, lda)
   end subroutine pic_sger

   pure subroutine pic_dgemv_x(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
      !! GEMV with explicit dimensions; forwards to the BLAS unchanged
      character(len=1), intent(in) :: trans
      integer(default_int), intent(in) :: m, n, lda, incx, incy
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: A(lda, *)
      real(dp), intent(in) :: x(*)
      real(dp), intent(inout) :: y(*)

      call blas_gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
   end subroutine pic_dgemv_x

   pure subroutine pic_dger_x(m, n, alpha, x, incx, y, incy, A, lda)
      !! rank-one update with explicit dimensions; forwards unchanged
      integer(default_int), intent(in) :: m, n, incx, incy, lda
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      real(dp), intent(in) :: y(*)
      real(dp), intent(inout) :: A(lda, *)

      call blas_ger(m, n, alpha, x, incx, y, incy, A, lda)
   end subroutine pic_dger_x

   pure subroutine pic_dger(A, x, y, alpha)
      !! rank-one update, dimensions taken from the arguments
      real(dp), intent(inout) :: A(:, :)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in), optional :: alpha
      real(dp) :: l_alpha
      integer(default_int) :: m, n, lda

      l_alpha = 1.0_dp
      if (present(alpha)) l_alpha = alpha
      m = size(A, 1, kind=default_int)
      n = size(A, 2, kind=default_int)
      lda = max(1_default_int, m)
      call blas_ger(m, n, l_alpha, x, 1_default_int, y, 1_default_int, A, lda)
   end subroutine pic_dger

end module pic_blas_interfaces
