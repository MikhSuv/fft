! Module providing precision parameters for floating-point arithmetic.
module precision_mod

   implicit none
   private

   public :: dp

   ! Double precision real kind: 15 significant digits, exponent range [-307,307]
   integer, parameter :: dp = selected_real_kind(15, 307)

end module precision_mod
