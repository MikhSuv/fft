! Unit tests for fft::make_twiddles().
!
! A check fails either when a wrong value is produced or when the routine
! crashes (bounds/FPE are caught by the Debug build flags). The program
! exits with a non-zero status if at least one check failed, so it works
! with both `ctest` and `fpm test`.
program test_make_twiddles
  use precision_mod
  use fft, only: make_twiddles
  implicit none

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  ! Absolute tolerance: a few hundred double-precision ulps
  real(dp), parameter :: tol = 100.0_dp * epsilon(1.0_dp)

  integer :: nfail = 0
  integer :: npass = 0

  call test_sizes()
  call test_known_values()
  call test_unit_modulus()
  call test_against_cos_sin()
  call test_conjugate_pair()
  call test_recurrence()
  call test_last_element()
  call test_odd_size()

  print '(a,i0,a,i0)', 'passed: ', npass, ', failed: ', nfail
  if (nfail > 0) then
    error stop 'make_twiddles tests FAILED'
  end if
  print '(a)', 'All make_twiddles tests passed'

contains

  ! The table length must be n/2 for both transform directions
  subroutine test_sizes()
    integer, parameter :: nsizes(*) = [2, 4, 8, 16, 1024]
    integer :: j, n
    complex(dp), allocatable :: w(:)

    do j = 1, size(nsizes)
      n = nsizes(j)

      w = make_twiddles(n, .false.)
      call check_ok(size(w) == n/2, 'size(w) == n/2, forward, n = ', n)
      deallocate(w)

      w = make_twiddles(n, .true.)
      call check_ok(size(w) == n/2, 'size(w) == n/2, inverse, n = ', n)
      deallocate(w)
    end do
  end subroutine test_sizes

  ! Closed-form values that do not depend on the implementation of exp()
  subroutine test_known_values()
    complex(dp), allocatable :: w(:)

    ! n = 4, forward: exp(-i*pi/2) = -i, exp(-i*pi) = -1
    w = make_twiddles(4, .false.)
    call check_close(w(1), cmplx(0.0_dp, -1.0_dp, kind=dp), 'n=4 fwd w(1) = -i')
    call check_close(w(2), cmplx(-1.0_dp, 0.0_dp, kind=dp), 'n=4 fwd w(2) = -1')
    deallocate(w)

    ! n = 4, inverse: exp(+i*pi/2) = +i, exp(+i*pi) = -1
    w = make_twiddles(4, .true.)
    call check_close(w(1), cmplx(0.0_dp, 1.0_dp, kind=dp), 'n=4 inv w(1) = +i')
    call check_close(w(2), cmplx(-1.0_dp, 0.0_dp, kind=dp), 'n=4 inv w(2) = -1')
    deallocate(w)

    ! n = 2, forward: exp(-i*pi) = -1
    w = make_twiddles(2, .false.)
    call check_close(w(1), cmplx(-1.0_dp, 0.0_dp, kind=dp), 'n=2 fwd w(1) = -1')
    deallocate(w)

    ! n = 8, forward: exp(-i*pi/4) = (sqrt(2)/2) * (1 - i)
    w = make_twiddles(8, .false.)
    call check_close(w(1), &
      cmplx(sqrt(0.5_dp), -sqrt(0.5_dp), kind=dp), 'n=8 fwd w(1) = (1-i)/sqrt(2)')
    deallocate(w)
  end subroutine test_known_values

  ! Every twiddle factor lies on the unit circle
  subroutine test_unit_modulus()
    integer, parameter :: n = 64
    complex(dp), allocatable :: w(:)
    logical :: ok
    integer :: k

    w = make_twiddles(n, .false.)
    ok = .true.
    do k = 1, n/2
      if (abs(abs(w(k)) - 1.0_dp) > tol) ok = .false.
    end do
    call check_true(ok, '|w(k)| == 1')
    deallocate(w)
  end subroutine test_unit_modulus

  ! Compare against an independent cos/sin evaluation of the same exponent
  subroutine test_against_cos_sin()
    integer, parameter :: nsizes(*) = [16, 64]
    complex(dp), allocatable :: w(:), expected(:)
    logical :: ok
    integer :: j, k, n
    real(dp) :: th

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(expected(n/2))

      do k = 1, n/2
        th = 2.0_dp * pi * real(k, dp) / real(n, dp)
        expected(k) = cmplx(cos(th), -sin(th), kind=dp)
      end do
      w = make_twiddles(n, .false.)
      ok = size(w) == size(expected)
      if (ok) ok = all(abs(w - expected) <= tol)
      call check_ok(ok, 'matches cos/sin, forward, n = ', n)
      deallocate(w)

      ! Inverse factors are the mirror image: exp(+i*theta)
      do k = 1, n/2
        th = 2.0_dp * pi * real(k, dp) / real(n, dp)
        expected(k) = cmplx(cos(th), sin(th), kind=dp)
      end do
      w = make_twiddles(n, .true.)
      ok = size(w) == size(expected)
      if (ok) ok = all(abs(w - expected) <= tol)
      call check_ok(ok, 'matches cos/sin, inverse, n = ', n)
      deallocate(w)

      deallocate(expected)
    end do
  end subroutine test_against_cos_sin

  ! Inverse twiddle factors are the complex conjugates of the forward ones
  subroutine test_conjugate_pair()
    integer, parameter :: n = 32
    complex(dp), allocatable :: wf(:), wi(:)
    logical :: ok
    integer :: k

    wf = make_twiddles(n, .false.)
    wi = make_twiddles(n, .true.)

    ok = size(wf) == size(wi)
    if (ok) then
      do k = 1, size(wf)
        if (abs(conjg(wf(k)) - wi(k)) > tol) ok = .false.
      end do
    end if
    call check_true(ok, 'conj(w_fwd) == w_inv')
    deallocate(wf, wi)
  end subroutine test_conjugate_pair

  ! Twiddle factors form a geometric progression: w(k+1) = w(k) * w(1)
  subroutine test_recurrence()
    integer, parameter :: n = 128
    complex(dp), allocatable :: w(:)
    logical :: ok
    integer :: k

    w = make_twiddles(n, .false.)
    ok = .true.
    do k = 1, n/2 - 1
      if (abs(w(k + 1) - w(k) * w(1)) > tol) ok = .false.
    end do
    call check_true(ok, 'w(k+1) == w(k)*w(1)')
    deallocate(w)
  end subroutine test_recurrence

  ! For even n the half-turn factor w(n/2) = exp(+/-i*pi) equals -1
  subroutine test_last_element()
    integer, parameter :: nsizes(*) = [16, 1024]
    complex(dp), allocatable :: w(:)
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)

      w = make_twiddles(n, .false.)
      call check_close(w(n/2), cmplx(-1.0_dp, 0.0_dp, kind=dp), &
        'w(n/2) == -1, forward, n = ' // trim(itoa(n)))
      deallocate(w)

      w = make_twiddles(n, .true.)
      call check_close(w(n/2), cmplx(-1.0_dp, 0.0_dp, kind=dp), &
        'w(n/2) == -1, inverse, n = ' // trim(itoa(n)))
      deallocate(w)
    end do
  end subroutine test_last_element

  ! Documented truncation: for an odd n only floor(n/2) factors are produced
  subroutine test_odd_size()
    complex(dp), allocatable :: w(:)

    w = make_twiddles(6, .false.)
    call check_true(size(w) == 3, 'odd n=6 gives floor(n/2) = 3 factors')
    deallocate(w)
  end subroutine test_odd_size

  subroutine check_true(cond, name)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: name

    if (cond) then
      npass = npass + 1
      print '(a,a)', 'PASS: ', name
    else
      nfail = nfail + 1
      print '(a,a)', 'FAIL: ', name
    end if
  end subroutine check_true

  ! Same as check_true, but appends an integer label (usually n)
  subroutine check_ok(cond, name, ival)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: name
    integer, intent(in) :: ival

    if (cond) then
      npass = npass + 1
      print '(a,a,i0)', 'PASS: ', name, ival
    else
      nfail = nfail + 1
      print '(a,a,i0)', 'FAIL: ', name, ival
    end if
  end subroutine check_ok

  subroutine check_close(got, expected, name)
    complex(dp), intent(in) :: got, expected
    character(len=*), intent(in) :: name

    if (abs(got - expected) <= tol) then
      npass = npass + 1
      print '(a,a)', 'PASS: ', name
    else
      nfail = nfail + 1
      print '(a,a)', 'FAIL: ', name
      print '(a,2(es23.15,1x))', '  expected = ', expected%re, expected%im
      print '(a,2(es23.15,1x))', '  got      = ', got%re, got%im
    end if
  end subroutine check_close

  ! Minimal integer to string helper for composing test names
  function itoa(n) result(s)
    integer, intent(in) :: n
    character(len=16) :: s

    write(s, '(i0)') n
  end function itoa

end program test_make_twiddles
