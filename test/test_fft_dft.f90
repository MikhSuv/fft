! Compares fft_transform against the direct O(n^2) evaluation of the DFT
! from the assignment,
!
!   y_k = (1/sqrt(n)) * sum_j x_j * exp(-+ 2*pi*i*j*k/n),
!
! for pseudo-random complex data in both directions. A check fails either
! when the error exceeds the tolerance or when the routine crashes
! (bounds/FPE are caught by the Debug build flags). The program exits with
! a non-zero status if at least one check failed, so it works with both
! `ctest` and `fpm test`.
program test_fft_dft
  use precision_mod, only: dp
  use fft, only: fft_transform
  implicit none

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

  integer :: nfail = 0
  integer :: npass = 0

  call random_init(repeatable=.true., image_distinct=.true.)
  call test_sizes(inverse=.false.)
  call test_sizes(inverse=.true.)

  print '(a,i0,a,i0)', 'passed: ', npass, ', failed: ', nfail
  if (nfail > 0) then
    error stop 'fft vs naive DFT tests FAILED'
  end if
  print '(a)', 'All fft-vs-DFT tests passed'

contains

  subroutine test_sizes(inverse)
    integer, parameter :: nsizes(*) = [2, 4, 8, 16, 32, 64, 128, 256]
    logical, intent(in) :: inverse
    integer :: j

    do j = 1, size(nsizes)
      call test_one(nsizes(j), inverse)
    end do
  end subroutine test_sizes

  subroutine test_one(n, inverse)
    integer, intent(in) :: n
    logical, intent(in) :: inverse
    complex(dp), allocatable :: x(:), y(:), ref(:)
    real(dp) :: scale, tol, err
    character(len=9) :: dir

    dir = merge('inverse ', 'forward ', inverse)
    allocate(x(n), y(n))

    call fill_random(x)
    y = x
    call fft_transform(y, inverse)
    ref = naive_dft(x, inverse)

    ! A few hundred double-precision ulps scaled by n and by the value range
    scale = max(maxval(abs(x)), maxval(abs(ref)), 1.0_dp)
    tol = 100.0_dp * epsilon(1.0_dp) * real(n, dp) * scale
    err = maxval(abs(y - ref))

    if (err <= tol) then
      npass = npass + 1
      print '(a,a,a,i0)', 'PASS: matches naive DFT, ', trim(dir), ', n = ', n
    else
      nfail = nfail + 1
      print '(a,a,a,i0)', 'FAIL: matches naive DFT, ', trim(dir), ', n = ', n
      print '(a,es14.6)', '  max error = ', err
      print '(a,es14.6)', '  tolerance = ', tol
    end if

    deallocate(x, y, ref)
  end subroutine test_one

  ! Random values uniform on the square [-1, 1) x [-1, 1)
  subroutine fill_random(x)
    complex(dp), intent(out) :: x(:)
    real(dp) :: u
    integer :: k

    do k = 1, size(x)
      call random_number(u)
      x(k) = cmplx(2.0_dp * u - 1.0_dp, 0.0_dp, kind=dp)
      call random_number(u)
      x(k) = x(k) + cmplx(0.0_dp, 2.0_dp * u - 1.0_dp, kind=dp)
    end do
  end subroutine fill_random

  ! Direct sum of the assignment formula; sgn picks the exponent sign
  function naive_dft(x, inverse) result(y)
    complex(dp), intent(in) :: x(:)
    logical, intent(in) :: inverse
    complex(dp) :: y(size(x))
    real(dp) :: sgn, coef
    integer :: n, j, k

    sgn = merge(1.0_dp, -1.0_dp, inverse)
    n = size(x)
    coef = sgn * 2.0_dp * pi / real(n, dp)

    do k = 1, n
      y(k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do j = 1, n
        y(k) = y(k) + x(j) * exp(iu * coef * real((j - 1) * (k - 1), dp))
      end do
      y(k) = y(k) / sqrt(real(n, dp))
    end do
  end function naive_dft

end program test_fft_dft
