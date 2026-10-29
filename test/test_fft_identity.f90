! Self-consistency checks of fft_transform that do not require the
! O(n^2) reference:
!
!   round trip: ifft(fft(x)) == x (both directions scaled by 1/sqrt(n))
!   impulse:    x = e0 -> every |y_k| = 1/sqrt(n)
!   constant:   x = 1/sqrt(n) -> y_0 = 1, the rest of the spectrum is zero
!   Parseval:   sum |y_k|^2 == sum |x_j|^2 (unitary transform)
!   linearity:  fft(a*x + b*y) == a*fft(x) + b*fft(y)
!
! The program exits with a non-zero status if at least one check failed,
! so it works with both `ctest` and `fpm test`.
program test_fft_identity
  use precision_mod, only: dp
  use fft, only: fft_transform
  implicit none

  integer :: nfail = 0
  integer :: npass = 0

  call random_init(repeatable=.true., image_distinct=.true.)
  call test_round_trip([2, 16, 65536, 1048576])
  call test_impulse([2, 8, 64, 1024])
  call test_constant([4, 8, 64, 1024])
  call test_parseval([16, 1024, 262144])
  call test_linearity([8, 64, 256])

  print '(a,i0,a,i0)', 'passed: ', npass, ', failed: ', nfail
  if (nfail > 0) then
    error stop 'fft identity tests FAILED'
  end if
  print '(a)', 'All fft-identity tests passed'

contains

  ! Fills x with pseudo-random values uniform on the square [-1,1)^2
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

  ! The same tolerance as in test_fft_dft: a few hundred ulps scaled by
  ! n and by the value range, which is very generous for real FFT bugs
  function tolerance(n, scale) result(tol)
    integer, intent(in) :: n
    real(dp), intent(in) :: scale
    real(dp) :: tol

    tol = 100.0_dp * epsilon(1.0_dp) * real(n, dp) * scale
  end function tolerance

  ! Applies fft_transform twice with opposite signs and compares with the
  ! original data. Includes n = 2^20 as a large-size smoke test.
  subroutine test_round_trip(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x(:), y(:)
    real(dp) :: scale, err, tol
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x(n), y(n))
      call fill_random(x)
      y = x
      call fft_transform(y, .false.)
      call fft_transform(y, .true.)

      scale = max(maxval(abs(x)), 1.0e-2_dp)
      tol = tolerance(n, scale)
      err = maxval(abs(y - x))

      if (err <= tol) then
        npass = npass + 1
        print '(a,i0)', 'PASS: round trip, n = ', n
      else
        nfail = nfail + 1
        print '(a,i0)', 'FAIL: round trip, n = ', n
        print '(a,es14.6)', '  max error = ', err
        print '(a,es14.6)', '  tolerance = ', tol
      end if
      deallocate(x, y)
    end do
  end subroutine test_round_trip

  ! Forward transform of the unit vector e0: y_k = exp(-2*pi*i*0*k/n)/sqrt(n)
  subroutine test_impulse(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x(:), y(:)
    real(dp) :: target, err, tol
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x(n), y(n))
      x = cmplx(0.0_dp, 0.0_dp, kind=dp)
      x(1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      call fft_transform(x, .false.)
      y = x

      target = 1.0_dp / sqrt(real(n, dp))
      tol = tolerance(n, 1.0_dp)
      ! Every |y_k| must equal the target
      err = maxval(abs(abs(y) - target))

      if (err <= tol) then
        npass = npass + 1
        print '(a,i0)', 'PASS: impulse spectrum, n = ', n
      else
        nfail = nfail + 1
        print '(a,i0)', 'FAIL: impulse spectrum, n = ', n
        print '(a,es14.6)', '  max error = ', err
        print '(a,es14.6)', '  tolerance = ', tol
      end if
      deallocate(x, y)
    end do
  end subroutine test_impulse

  subroutine test_constant(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x(:), y(:)
    real(dp) :: err, tol
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x(n), y(n))
      x = cmplx(1.0_dp, 0.0_dp, kind=dp) / sqrt(real(n, dp))
      call fft_transform(x, .false.)
      y = x

      tol = tolerance(n, 1.0_dp)
      err = max(abs(y(1) - cmplx(1.0_dp, 0.0_dp, kind=dp)), &
                maxval(abs(y(2:))))

      if (err <= tol) then
        npass = npass + 1
        print '(a,i0)', 'PASS: constant spectrum, n = ', n
      else
        nfail = nfail + 1
        print '(a,i0)', 'FAIL: constant spectrum, n = ', n
        print '(a,es14.6)', '  max error = ', err
        print '(a,es14.6)', '  tolerance = ', tol
      end if
      deallocate(x, y)
    end do
  end subroutine test_constant

  ! Parseval identity for the unitary transform
  subroutine test_parseval(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x(:)
    real(dp) :: ex, ey, scale, tol, err
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x(n))
      call fill_random(x)
      ex = sum(abs(x)**2)
      call fft_transform(x, .false.)
      ey = sum(abs(x)**2)

      scale = max(ex, ey)
      tol = tolerance(n, 1.0_dp) * ex
      err = abs(ey - ex)

      if (err <= tol) then
        npass = npass + 1
        print '(a,i0)', 'PASS: Parseval identity, n = ', n
      else
        nfail = nfail + 1
        print '(a,i0)', 'FAIL: Parseval identity, n = ', n
        print '(a,es14.6)', '  energy error = ', err
        print '(a,es14.6)', '  tolerance    = ', tol
      end if
      deallocate(x)
    end do
  end subroutine test_parseval

  ! Superposition: fft(a*x1 + b*x2) == a*fft(x1) + b*fft(x2)
  subroutine test_linearity(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x1(:), x2(:), z(:), y1(:), y2(:), t(:)
    complex(dp), parameter :: ca = (0.5_dp, 0.0_dp)
    complex(dp), parameter :: cb = (-0.3_dp, 0.1_dp)
    real(dp) :: scale, err, tol
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x1(n), x2(n), z(n), y1(n), y2(n), t(n))
      call fill_random(x1)
      call fill_random(x2)

      y1 = x1
      y2 = x2
      call fft_transform(y1, .false.)
      call fft_transform(y2, .false.)

      z = ca * x1 + cb * x2
      call fft_transform(z, .false.)
      t = ca * y1 + cb * y2

      scale = max(maxval(abs(z)), maxval(abs(t)), 1.0_dp)
      tol = tolerance(n, scale)
      err = maxval(abs(z - t))

      if (err <= tol) then
        npass = npass + 1
        print '(a,i0)', 'PASS: linearity, n = ', n
      else
        nfail = nfail + 1
        print '(a,i0)', 'FAIL: linearity, n = ', n
        print '(a,es14.6)', '  max error = ', err
        print '(a,es14.6)', '  tolerance = ', tol
      end if
      deallocate(x1, x2, z, y1, y2, t)
    end do
  end subroutine test_linearity

end program test_fft_identity
