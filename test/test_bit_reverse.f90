! Unit tests of fft::bit_reverse_permute()
program test_bit_reverse
  use precision_mod, only: dp
  use fft, only: bit_reverse_permute
  implicit none

  integer :: nfail = 0
  integer :: npass = 0

  call random_init(repeatable=.true., image_distinct=.true.)
  call test_matches_reference([2, 8, 64, 4096])
  call test_involution([2, 8, 1024])

  print '(a,i0,a,i0)', 'passed: ', npass, ', failed: ', nfail
  if (nfail > 0) then
    error stop 'bit_reverse tests FAILED'
  end if
  print '(a)', 'All bit_reverse tests passed'

contains

  ! Reference: the highest input bit lands at the lowest output position
  ! and vice versa, weights descending against weights ascending
  pure integer function ref_bitrev(i, n) result(r)
    integer, intent(in) :: i, n
    integer :: k, w

    r = 0
    w = 1
    k = n / 2
    do while (k > 0)
      if (iand(i, k) /= 0) r = r + w
      w = shiftl(w, 1)
      k = shiftr(k, 1)
    end do
  end function ref_bitrev

  ! After the permutation each input value must sit exactly at the
  ! bit-reversed position of its index  subroutine test_matches_reference(nsizes)
  subroutine test_matches_reference(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x(:), y(:)
    logical :: ok
    integer :: i, j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x(n), y(n))

      call fill_random(x)
      y = x
      call bit_reverse_permute(y)

      ok = .true.
      do i = 0, n - 1
        if (y(ref_bitrev(i, n) + 1) /= x(i + 1)) ok = .false.
      end do
      call check_true(ok, 'matches reference bitrev, n = ' // trim(itoa(n)))

      deallocate(x, y)
    end do
  end subroutine test_matches_reference

  ! Applying the permutation twice must return the original order
  subroutine test_involution(nsizes)
    integer, intent(in) :: nsizes(:)
    complex(dp), allocatable :: x(:), y(:)
    integer :: j, n

    do j = 1, size(nsizes)
      n = nsizes(j)
      allocate(x(n), y(n))

      call fill_random(x)
      y = x
      call bit_reverse_permute(y)
      call bit_reverse_permute(y)

      call check_ok(all(y == x), 'involution (twice = identity), n = ', n)
      deallocate(x, y)
    end do
  end subroutine test_involution

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

  ! Minimal integer to string helper for composing test names
  function itoa(n) result(s)
    integer, intent(in) :: n
    character(len=16) :: s

    write(s, '(i0)') n
  end function itoa

end program test_bit_reverse
