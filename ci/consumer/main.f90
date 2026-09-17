program consumer
!! Minimal downstream consumer of pic-blas, exercised by CI only.
   use pic_blas_interfaces, only: pic_gemm
   use pic_types, only: dp, default_int
   implicit none

   integer(default_int), parameter :: n = 2
   real(dp) :: a(n, n), b(n, n), c(n, n)

   a = 1.0_dp
   b = 2.0_dp
   c = 0.0_dp

   call pic_gemm(a, b, c)

   if (abs(c(1, 1) - 4.0_dp) > 1.0e-10_dp) then
      error stop "pic_gemm through pic-blas::pic-blas gave an unexpected result"
   end if

   print *, "consumer linked pic-blas and pic modules resolved"

end program consumer
