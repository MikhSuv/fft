! Fast Fourier transform of complex data, decimation in frequency, in place.
!
! The transform is the unitary DFT.
!
!   y_k = (1/sqrt(n)) * sum_{j=0}^{n-1} x_j * exp(-+ 2*pi*i*j*k/n),
!
! where the "-" sign corresponds to the forward transform,
! "+" to the inverse one. 
module fft
  use precision_mod, only: dp
  implicit none
  private
  public :: make_twiddles, fft_transform, bit_reverse_permute

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  ! Imaginary unit
  complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

contains

  ! True when n is a positive power of two
  pure logical function is_pow2(n)
    integer, intent(in) :: n
    is_pow2 = n >= 1 .and. iand(n, n - 1) == 0
  end function is_pow2

  ! Builds the table of twiddle factors w(k+1) = exp(sgn * 2*pi*i * k/n)
  ! for k = 0, ..., n/2 - 1, i.e. W_n^k of the assignment notation.
  !
  !   inverse = .false. -> sgn = -1: forward transform
  !   inverse = .true.  -> sgn = +1: inverse transform
  function make_twiddles(n, inverse) result(w)
    integer, intent(in) :: n
    logical, intent(in) :: inverse
    complex(dp), allocatable :: w(:)
    real(dp) :: sgn, coef
    integer :: k

    if (n < 2 .or. .not. is_pow2(n)) then
      error stop 'make_twiddles: n must be a power of two, n >= 2'
    end if

    ! Sign of the exponent: +1 for the inverse transform, -1 for the forward one
    sgn = merge(1.0_dp, -1.0_dp, inverse)
    coef = sgn * 2.0_dp * pi / real(n, dp)

    allocate(w(n / 2))
    !$omp parallel do schedule(static)
    do k = 0, n / 2 - 1
      w(k + 1) = exp(iu * coef * real(k, dp))
    end do
    !$omp end parallel do
  end function make_twiddles

  ! Decimation-in-frequency butterflies, in place.
  !
  ! Block of length m with half h = m/2, for j = 0..h-1:
  !   x(j)   <- x(j) + x(j+h)
  !   x(j+h) <- (x(j) - x(j+h)) * W_m^j,   W_m^j = W_n^(j*n/m)
  ! The stages run for m = n, n/2, ..., 2.
  ! The output is in bit-reversed order.
  ! Table layout: w(k+1) = W_n^k.
  subroutine dif_stages(x, w)
    complex(dp), intent(inout) :: x(:)
    complex(dp), intent(in) :: w(:)
    integer :: n, m, h, stride
    integer :: s, j
    complex(dp) :: a, b

    n = size(x)
    m = n

    do while (m >= 2)
      h = m / 2
      stride = n / m

      !$omp parallel do collapse(2) private(j, a, b) schedule(static)
      do s = 0, n - 1, m
        do j = 0, h - 1
          a = x(s + j + 1)
          b = x(s + j + h + 1)
          x(s + j + 1) = a + b
          x(s + j + h + 1) = (a - b) * w(j * stride + 1)
        end do
      end do
      !$omp end parallel do
      m = h
    end do
  end subroutine dif_stages

  ! Bit-reversal permutation of x, in place.
  ! Precondition: size(x) is a power of two (checked).
  subroutine bit_reverse_permute(x)
    complex(dp), intent(inout) :: x(:)
    integer :: n, i, j
    complex(dp) :: t

    n = size(x)
    if (.not. is_pow2(n)) then
      error stop 'bit_reverse_permute: size(x) must be a power of two'
    end if

    !$omp parallel do private(j, t) schedule(static)
    do i = 0, n - 1
      j = reversed_bits(i, n)
      if (i < j) then
        t = x(i + 1)
        x(i + 1) = x(j + 1)
        x(j + 1) = t
      end if
    end do
    !$omp end parallel do
  end subroutine bit_reverse_permute

  ! Reverses the p significant bit positions of i, where n = 2^p is the
  ! transform length and i = 0..n-1, i.e. D_p(i) of the bit-reversal
  ! permutation (the leading zeros of i land in the high positions)
  pure integer function reversed_bits(i, n) result(r)
    integer, intent(in) :: i, n
    integer :: k, p

    ! p = log2(n)
    p = 0
    k = n
    do while (k > 1)
      p = p + 1
      k = shiftr(k, 1)
    end do

    ! Reversal over all p bits, low-order zeros included
    r = 0
    k = i
    do while (p > 0)
      r = shiftl(r, 1) + iand(k, 1)
      k = shiftr(k, 1)
      p = p - 1
    end do
  end function reversed_bits

  ! fft_transform(x, inverse):
  !   inverse = .false. -> forward transform (sign "-")
  !   inverse = .true.  -> inverse transform (sign "+")
  ! Both directions are scaled by the unitary factor 1/sqrt(n).
  subroutine fft_transform(x, inverse)
    complex(dp), intent(inout) :: x(:)
    logical, intent(in) :: inverse
    complex(dp), allocatable :: w(:)
    integer :: n

    n = size(x)
    if (n == 1) return ! Degenerate one-point transform is the identity
    if (.not. is_pow2(n)) then
      error stop 'fft_transform: size(x) must be a power of two'
    end if

    w = make_twiddles(n, inverse)
    call dif_stages(x, w)
    call bit_reverse_permute(x)

    ! Unitary scaling
    x = x / sqrt(real(n, dp))
  end subroutine fft_transform

end module fft
